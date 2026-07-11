"""Displayed-formula generator.
Auto-layout from SMILES: horizontal backbone, C=O up, all H shown, Arial, monochrome.
Acyclic molecules (chains + one functional group + simple branches). Aromatic handled separately."""
import math, collections, re, statistics
from rdkit import Chem
from rdkit.Chem import rdDepictor
from PIL import Image, ImageDraw, ImageFont

ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
CALIBRI = "/Applications/Microsoft Word.app/Contents/Resources/DFonts/Calibri.ttf"
FONT_PATH = ARIAL              # label font: Arial for clean, even structure labels
FONT_NAME = "Arial"            # font-family used in SVG output
U = 46; STROKE = 3; GAP = 15; DBL = 6; TRP = 7
HEAVY_LEN = 1.95; H_LEN = 1.0
FONT = ImageFont.truetype(FONT_PATH, 34)
DIRS = {'up': (0, 1), 'down': (0, -1), 'left': (-1, 0), 'right': (1, 0)}
OPP = {'up': 'down', 'down': 'up', 'left': 'right', 'right': 'left'}


def _is_carbonyl_O(a):
    return (a.GetAtomicNum() == 8 and a.GetDegree() == 1
            and any(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds()))


def _longest_path(adj, weight=None, edgewt=None):
    if not adj:
        return []
    w = weight or (lambda n: 0)
    ew = edgewt or (lambda a, b: 0)
    def pw(path):
        return sum(w(n) for n in path)
    def ep(path):                                        # unsaturation on the path: keep C=C / C#C ON the main chain, not as a substituent
        return sum(ew(path[i], path[i + 1]) for i in range(len(path) - 1))
    def bfs(s):                                          # priority: most multiple bonds, then longest, then highest node-weight (=O / heteroatom)
        seen = {s: [s]}; q = collections.deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in seen:
                    seen[v] = seen[u] + [v]; q.append(v)
        far = max(seen, key=lambda k: (ep(seen[k]), len(seen[k]), pw(seen[k])))
        return far, seen[far]
    a, _ = bfs(next(iter(adj)))
    _, path = bfs(a)
    return path


def _elem(a):
    """Element symbol, mapping a SMILES wildcard/dummy atom to the house R-group glyph."""
    return _RGROUP.get(a.GetAtomMapNum(), "R") if a.GetAtomicNum() == 0 else a.GetSymbol()


def _condense(molH, idx, exclude, visited=None):
    """Full condensed-formula string for the branch rooted at `idx`, preserving EVERY atom (no lossy
    collapse). `exclude` = atoms outside the branch (the ring / the stereocentre)."""
    if visited is None:
        visited = set(exclude)
    visited = visited | {idx}
    a = molH.GetAtomWithIdx(idx)
    if a.GetAtomicNum() == 0:
        return _elem(a)
    nbrs = [n for n in a.GetNeighbors() if n.GetIdx() not in visited]
    nH = sum(1 for n in nbrs if n.GetAtomicNum() == 1)
    heavy = [n for n in nbrs if n.GetAtomicNum() > 1]
    core = a.GetSymbol() + ("H" if nH >= 1 else "") + (str(nH) if nH >= 2 else "") + _charge_suffix(a.GetFormalCharge())
    if not heavy:
        return core
    subs = []
    for n in heavy:
        bt = molH.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondType()
        bsym = {Chem.BondType.DOUBLE: "=", Chem.BondType.TRIPLE: "#"}.get(bt, "")
        subs.append((bsym, _condense(molH, n.GetIdx(), exclude, visited)))
    if len(subs) == 1:
        return core + subs[0][0] + subs[0][1]
    if len({s for _, s in subs}) == 1 and all(b == "" for b, _ in subs):
        return core + "(" + subs[0][1] + ")" + str(len(subs))    # e.g. CH(CH3)2
    subs.sort(key=lambda bs: len(bs[1]))                          # shortest branches parenthesised, longest inline
    *paren, inline = subs
    out = core
    for (b, s), k in collections.Counter(paren).items():
        out += "(" + b + s + ")" + (str(k) if k > 1 else "")
    return out + inline[0] + inline[1]


def _sub_label(molH, ringset, root):
    """Condensed label for a substituent rooted at `root` (attached across `ringset`). Canonical short
    forms for the common groups; otherwise a full condensed formula that preserves every atom."""
    a = molH.GetAtomWithIdx(root)
    Z = a.GetAtomicNum()
    nbrs = [n for n in a.GetNeighbors() if n.GetIdx() not in ringset]
    heavy = [n for n in nbrs if n.GetAtomicNum() > 1]
    nH = sum(1 for n in nbrs if n.GetAtomicNum() == 1)
    ch = _charge_suffix(a.GetFormalCharge())
    if Z == 17: return "Cl"
    if Z == 35: return "Br"
    if Z == 9:  return "F"
    if Z == 53: return "I"
    if Z == 8:  return _condense(molH, root, ringset) if heavy else ("OH" if nH else "O") + ch
    if Z == 7:
        if sum(1 for n in heavy if n.GetAtomicNum() == 8) >= 2:
            return "NO2"
        if not heavy:
            return ("NH2" if nH == 2 else "NH" if nH == 1 else "N") + ch
        return _condense(molH, root, ringset)                    # e.g. NHCOCH3, preserves atoms
    if Z == 6:
        if any(molH.GetBondBetweenAtoms(root, n.GetIdx()).GetBondType() == Chem.BondType.TRIPLE
               and n.GetAtomicNum() == 7 for n in heavy):
            return "CN"
        dO = next((n for n in heavy if n.GetAtomicNum() == 8 and
                   molH.GetBondBetweenAtoms(root, n.GetIdx()).GetBondType() == Chem.BondType.DOUBLE), None)
        if dO is not None:
            others = [n for n in heavy if n.GetIdx() != dO.GetIdx()]
            if not others:
                return "CHO"                                     # -CHO aldehyde
            if len(others) == 1:                                 # acid / ester / amide / ketone
                o = others[0]
                if o.GetAtomicNum() == 8 and any(nb.GetAtomicNum() == 1 for nb in o.GetNeighbors()):
                    return "COOH"
                return "CO" + _condense(molH, o.GetIdx(), ringset | {root})   # COOR / CONH2 / COCH3, preserves atoms
            return _condense(molH, root, ringset)
        return _condense(molH, root, ringset)                    # general C substituent, preserves every atom
    if Z == 0:
        return _elem(a)
    return a.GetSymbol() + ch


def _alkyl_branch(mol, ringset, root):
    """True if the substituent rooted at `root` is a pure alkyl branch (only C/H, no double bonds)."""
    seen = set(ringset); stack = [root]
    while stack:
        u = stack.pop()
        if u in seen:
            continue
        seen.add(u); a = mol.GetAtomWithIdx(u)
        if a.GetAtomicNum() != 6 or any(b.GetBondType() != Chem.BondType.SINGLE for b in a.GetBonds()):
            return False
        stack += [nb.GetIdx() for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1 and nb.GetIdx() not in seen]
    return True


def _ring_layout(mol, form='structural'):
    """Single ring as a skeletal polygon (no C/H labels). Aromatic 6-rings get an inscribed circle;
    cycloalkanes are plain polygons. Substituents: condensed labels (structural), bare skeletal stubs
    for alkyl (skeletal), or explicit atoms+H (displayed)."""
    molH = Chem.AddHs(mol)
    ring = max(mol.GetRingInfo().AtomRings(), key=len)
    n = len(ring)
    ringset = set(ring)
    is_arom = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
    subs = []          # (ring_atom, root) single-bond substituents
    carbonyls = []     # (ring_atom, O_idx) exocyclic C=O (e.g. lactone / cyclic ketone)
    for ri_atom in ring:
        for nb in mol.GetAtomWithIdx(ri_atom).GetNeighbors():
            j = nb.GetIdx()
            if j in ringset or nb.GetAtomicNum() == 1:
                continue
            bt = mol.GetBondBetweenAtoms(ri_atom, j).GetBondType()
            if nb.GetAtomicNum() == 8 and bt == Chem.BondType.DOUBLE:
                carbonyls.append((ri_atom, j))
            else:
                subs.append((ri_atom, j))
    start = ring.index(subs[0][0]) if subs else 0
    ordered = ring[start:] + ring[:start]
    edge = 1.55
    R = edge / (2 * math.sin(math.pi / n))
    upright = {3: 90, 4: 45, 5: 90, 6: 0}.get(n, 90)
    base = 0.0 if subs else upright
    angles = {ordered[k]: math.radians(base + 360.0 / n * k) for k in range(n)}
    atoms, bonds, circles = {}, [], []
    for i in ring:                                   # ring vertices: C blank, heteroatoms (O/N) labelled
        th = angles[i]
        a = mol.GetAtomWithIdx(i); sym = a.GetSymbol()
        if sym == "C":
            lab = ""
        elif sym == "N" and a.GetTotalNumHs() >= 1:
            lab = "NH"                               # ring N-H (e.g. piperidine, lactam)
        else:
            lab = sym
        atoms[i] = (lab, R * math.cos(th), R * math.sin(th))
    for k in range(n):                               # ring bonds: keep true order unless aromatic (→ circle)
        a, b = ordered[k], ordered[(k + 1) % n]
        order = int(mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble())
        bonds.append((a, b, 1 if (is_arom or order < 2) else order))
    if is_arom:
        circles.append((0.0, 0.0, 0.60 * R))
    for ri_atom, o_idx in carbonyls:                 # exocyclic =O drawn radially outward
        th = angles[ri_atom]
        atoms[o_idx] = ("O", (R + 1.05) * math.cos(th), (R + 1.05) * math.sin(th))
        bonds.append((ri_atom, o_idx, 2))
    pid = 10000
    hid = -1
    reversible = set()
    for ri_atom, root in subs:
        th = angles[ri_atom]
        ra = mol.GetAtomWithIdx(root)
        rx, ry = math.cos(th), math.sin(th); pxr, pyr = -ry, rx
        vx, vy = R * rx, R * ry
        dO = next((nb for nb in ra.GetNeighbors() if nb.GetAtomicNum() == 8 and
                   mol.GetBondBetweenAtoms(root, nb.GetIdx()).GetBondType() == Chem.BondType.DOUBLE), None)
        nN = next((nb for nb in ra.GetNeighbors() if nb.GetAtomicNum() == 7 and nb.GetIdx() not in ringset), None)
        Z, nHr = ra.GetAtomicNum(), ra.GetTotalNumHs()
        heavy_out = [nb.GetIdx() for nb in ra.GetNeighbors() if nb.GetAtomicNum() > 1 and nb.GetIdx() not in ringset]
        # amide substituent → draw out as ring–C(=O)–NH2 (carbonyl explicit)
        if Z == 6 and dO is not None and nN is not None:
            fcx, fcy = vx + rx * 1.4, vy + ry * 1.4
            atoms[root] = ("C", fcx, fcy); bonds.append((ri_atom, root, 1))
            atoms[dO.GetIdx()] = ("O", fcx + pxr * 1.15, fcy + pyr * 1.15); bonds.append((root, dO.GetIdx(), 2))
            nH = nN.GetTotalNumHs()
            nlab = "N" + ("H" if nH >= 1 else "") + (str(nH) if nH >= 2 else "")
            atoms[nN.GetIdx()] = (nlab, fcx + rx * 1.5, fcy + ry * 1.5); bonds.append((root, nN.GetIdx(), 1))
            if rx < -0.3:
                reversible.add(nN.GetIdx())
        elif form == 'displayed' and Z == 6 and dO is not None:
            # carbonyl substituent drawn OUT like the acyclic displayed does (benzaldehyde -CHO,
            # benzoic acid -COOH, acyl chloride -COCl): ring-C(=O)-X, trigonal splay pointing outward
            aR = math.atan2(ry, rx)
            fcx, fcy = vx + rx * 1.4, vy + ry * 1.4
            atoms[root] = ("C", fcx, fcy); bonds.append((ri_atom, root, 1))
            dOx, dOy = math.cos(aR + math.radians(60)), math.sin(aR + math.radians(60))   # =O biased UP (house carbonyl-up)
            dXx, dXy = math.cos(aR - math.radians(60)), math.sin(aR - math.radians(60))   # other bond, splayed down
            atoms[dO.GetIdx()] = ("O", fcx + dOx * 1.3, fcy + dOy * 1.3); bonds.append((root, dO.GetIdx(), 2))
            ex, ey = fcx + dXx * 1.3, fcy + dXy * 1.3
            others = [nb for nb in ra.GetNeighbors() if nb.GetIdx() != dO.GetIdx() and nb.GetIdx() not in ringset]
            if others and others[0].GetAtomicNum() == 8 and others[0].GetTotalNumHs() >= 1:
                o2 = others[0]                                # -OH (acid): O drawn, then the H bent off it
                atoms[o2.GetIdx()] = ("O", ex, ey); bonds.append((root, o2.GetIdx(), 1))
                hx = ex + math.cos(aR - math.radians(115)) * 1.15
                hy = ey + math.sin(aR - math.radians(115)) * 1.15
                atoms[hid] = ("H", hx, hy); bonds.append((o2.GetIdx(), hid, 1)); hid -= 1
            elif others:                                      # -Cl / ester -O-CH3: condensed label on the outer bond
                atoms[pid] = (_sub_label(molH, ringset | {root}, others[0].GetIdx()), ex, ey)
                bonds.append((root, pid, 1)); pid += 1
            elif nHr >= 1:                                    # aldehyde -CHO: the H drawn explicitly
                atoms[hid] = ("H", ex, ey); bonds.append((root, hid, 1)); hid -= 1
        elif form == 'displayed' and not heavy_out and Z in (6, 7, 8):
            # terminal substituent (-CH3 / -OH / -NH2) expanded to explicit atoms + H, radially outward
            cx, cy = vx + rx * 1.4, vy + ry * 1.4
            atoms[root] = (_elem(ra) + _charge_suffix(ra.GetFormalCharge()), cx, cy)
            bonds.append((ri_atom, root, 1))
            if Z == 6:                                    # methyl comb: outward + two perpendicular
                dirs = [(rx, ry), (pxr, pyr), (-pxr, -pyr)]
            else:                                         # O/N: H's BEND off to the sides (never collinear with the ring bond)
                a0 = math.atan2(ry, rx)
                dirs = [(math.cos(a0 + math.radians(60)), math.sin(a0 + math.radians(60))),
                        (math.cos(a0 - math.radians(60)), math.sin(a0 - math.radians(60)))]
            for k in range(min(nHr, len(dirs))):
                dx, dy = dirs[k]
                atoms[hid] = ("H", cx + dx * 1.25, cy + dy * 1.25); bonds.append((root, hid, 1)); hid -= 1
        elif form == 'skeletal' and Z == 6 and _alkyl_branch(mol, ringset, root):
            # alkyl stub: bare skeletal vertices zig-zagging outward (implicit C/H), no labels
            prev, pdir, k = ri_atom, (rx, ry), 0
            stack = [(root, 1.4)]
            placed = {}
            cur = root; ox, oy = vx, vy
            while True:
                a = mol.GetAtomWithIdx(cur)
                s = 0 if k == 0 else (30 if k % 2 == 1 else -30)   # first bond STRAIGHT out (radial); a longer chain then zig-zags
                ddx = rx * math.cos(math.radians(s)) - ry * math.sin(math.radians(s))
                ddy = rx * math.sin(math.radians(s)) + ry * math.cos(math.radians(s))
                ox, oy = ox + ddx * 1.4, oy + ddy * 1.4
                atoms[cur] = ("", ox, oy); bonds.append((prev, cur, 1)); placed[cur] = True
                nxts = [nb.GetIdx() for nb in a.GetNeighbors() if nb.GetIdx() not in placed and nb.GetIdx() not in ringset and nb.GetAtomicNum() == 6]
                if not nxts:
                    break
                prev, cur, k = cur, nxts[0], k + 1
        else:
            lab = _sub_label(molH, ringset, root)
            atoms[pid] = (lab, (R + 1.25) * rx, (R + 1.25) * ry)
            bonds.append((ri_atom, pid, 1)); reversible.add(pid); pid += 1
    return atoms, bonds, circles, reversible


HLEN = 1.25                                              # visual length of a terminal X-H bond (heavy bond == HEAVY_LEN)
_BONDORD = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}


def _fweight(a):
    """Functional weight: heteroatoms and carbonyl carbons pull to the RIGHT (house rule)."""
    w = 2 if a.GetAtomicNum() in (7, 8, 9, 17, 35, 53) else 0
    for nb in a.GetNeighbors():
        if nb.GetAtomicNum() == 8 and a.GetOwningMol().GetBondBetweenAtoms(
                a.GetIdx(), nb.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
            w += 2
    return w


def _carbonyl_Os(mol):
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 8 and a.GetDegree() == 1
            and any(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds())]


def _coords(mol, shorten_h=None):
    """RDKit 2D coords, house-styled: median heavy bond == HEAVY_LEN; the carbon skeleton's principal
    axis laid HORIZONTAL (natural bond angles preserved); functional end on the RIGHT; C=O pointing UP;
    terminal X-H bonds shortened to `shorten_h` so hydrogens tuck in. Returns {idx:(x,y)}."""
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    P = {i: (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(mol.GetNumAtoms())}
    at = mol.GetAtomWithIdx
    hb = [math.hypot(P[b.GetBeginAtomIdx()][0]-P[b.GetEndAtomIdx()][0],
                     P[b.GetBeginAtomIdx()][1]-P[b.GetEndAtomIdx()][1]) for b in mol.GetBonds()
          if at(b.GetBeginAtomIdx()).GetAtomicNum() > 1 and at(b.GetEndAtomIdx()).GetAtomicNum() > 1]
    s = HEAVY_LEN / (statistics.median(hb) if hb else 1.5)
    P = {i: (x*s, y*s) for i, (x, y) in P.items()}
    carbons = [i for i in P if at(i).GetAtomicNum() == 6]
    if len(carbons) >= 2:                                # rotate the carbon principal axis to horizontal
        cx = sum(P[i][0] for i in carbons)/len(carbons); cy = sum(P[i][1] for i in carbons)/len(carbons)
        sxx = sum((P[i][0]-cx)**2 for i in carbons); syy = sum((P[i][1]-cy)**2 for i in carbons)
        sxy = sum((P[i][0]-cx)*(P[i][1]-cy) for i in carbons)
        th = 0.5*math.atan2(2*sxy, sxx-syy); c, sn = math.cos(-th), math.sin(-th)
        P = {i: ((x-cx)*c-(y-cy)*sn, (x-cx)*sn+(y-cy)*c) for i, (x, y) in P.items()}
    heavy = [i for i in P if at(i).GetAtomicNum() > 1]
    xg = sum(P[i][0] for i in heavy)/len(heavy)          # mirror x so the functional end is on the right
    wsum = sum(_fweight(at(i)) for i in heavy)
    if wsum and sum(_fweight(at(i))*P[i][0] for i in heavy)/wsum < xg - 1e-6:
        P = {i: (-x, y) for i, (x, y) in P.items()}
    cO = _carbonyl_Os(mol)                               # flip y so C=O points up
    if cO:
        rise = sum(P[o][1] - P[at(o).GetNeighbors()[0].GetIdx()][1] for o in cO)
        if rise < 0:
            P = {i: (x, -y) for i, (x, y) in P.items()}
    if shorten_h:                                        # tuck terminal H's in toward their heavy neighbour
        for i in [k for k in P if at(k).GetAtomicNum() == 1]:
            nb = at(i).GetNeighbors()[0].GetIdx()
            hx, hy = P[i]; bx, by = P[nb]; d = math.hypot(hx-bx, hy-by) or 1
            P[i] = (bx + (hx-bx)/d*shorten_h, by + (hy-by)/d*shorten_h)
    return P


def _displayed_comb_layout(mol, grid=False, spine_atoms=None, condense=None):
    """Displayed formula in RIGHT-ANGLE 'comb' geometry (authentic A-level, per mined notes):
    horizontal heavy-atom spine with vertical C-H / substituent bonds. Every atom + bond explicit.
    grid=True forces EVERY atom (sp3, sp2, O/N alike) onto 90-deg right angles — the 'natural' form
    that shows all bonds but does NOT follow real bond angles (no trigonal splay, no O-H/N-H bend)."""
    Chem.Kekulize(mol, clearAromaticFlags=True)
    molH = Chem.AddHs(mol)
    at = molH.GetAtomWithIdx
    A = molH.GetNumAtoms()
    carbonylO = {i for i in range(A) if _is_carbonyl_O(at(i))}
    hydroxylO = {i for i in range(A) if at(i).GetAtomicNum() == 8 and i not in carbonylO
                 and sum(1 for nb in at(i).GetNeighbors() if nb.GetAtomicNum() > 1) == 1
                 and any(nb.GetAtomicNum() == 1 for nb in at(i).GetNeighbors())}
    halogenX = {i for i in range(A) if at(i).GetAtomicNum() in (9, 17, 35, 53)
                and sum(1 for nb in at(i).GetNeighbors() if nb.GetAtomicNum() > 1) == 1
                and at(i).GetNeighbors()[0].GetAtomicNum() == 6}   # terminal C-X halogen -> substituent (NJC: C-chain horizontal, X drawn up/down)
    _off = carbonylO | hydroxylO | halogenX              # keep -OH / -X off the backbone -> drawn as substituents
    adj = {i: [] for i in range(A) if at(i).GetAtomicNum() > 1 and i not in _off}
    for i in adj:
        adj[i] = [nb.GetIdx() for nb in at(i).GetNeighbors() if nb.GetIdx() in adj]
    spine = _longest_path(adj, weight=lambda n: 2 if any(nb.GetIdx() in carbonylO for nb in at(n).GetNeighbors())
                                          else (1 if at(n).GetAtomicNum() in (7, 8) else 0),
                          edgewt=lambda a, b: 1 if molH.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() >= 2 else 0)
    def fw(idx):                                         # functional end pulls right
        a = at(idx); w = 2 if a.GetAtomicNum() in (7, 8, 9, 17, 35, 53) else 0
        return w + (2 if any(nb.GetIdx() in carbonylO for nb in a.GetNeighbors()) else 0)
    tot = sum(fw(i) for i in spine)
    if tot and spine and sum(k*fw(i) for k, i in enumerate(spine))/tot < (len(spine)-1)/2.0:
        spine = spine[::-1]
    if spine_atoms is not None:                          # mechanism form: force the backbone to the functional-group bond
        spine = [i for i in spine_atoms if i in adj] or spine
    condense = condense or set()                          # heavy substituents to draw as ONE condensed label (CH3, C2H5...)
    condlab = {}
    # UNIFIED BFS placement. Every atom places its own neighbours at angles set by ITS hybridisation,
    # measured off the bond it arrived on (its parent). sp3 -> 90 deg comb (straight-ahead spine, vertical
    # H); sp2 carbon / O / N -> 120 deg trigonal-or-bent. Nothing is ever pre-forced onto a flat line, so
    # mid-chain carbonyls, the ester/ether -O-, and secondary-amine N all splay instead of going collinear.
    R = math.radians
    spine_pos = {v: k for k, v in enumerate(spine)}
    pos = {}
    def blen(i, j):                                      # uniform bond length (C-H = C-C = C-Br), per NJC displayed
        return HEAVY_LEN
    def collides(x, y):
        return any(abs(x-px) < 0.9 and abs(y-py) < 0.72 for px, py in pos.values())
    def _is_dbl(i, j):
        return molH.GetBondBetweenAtoms(i, j).GetBondType() == Chem.BondType.DOUBLE
    def _angular(i):                                     # sp2 carbon (any =) OR any O / N: bonds splay ~120 deg, never collinear
        if grid:
            return False                                 # grid/natural form: everything on 90-deg right angles
        return at(i).GetAtomicNum() in (7, 8) or any(_is_dbl(i, nb.GetIdx()) for nb in at(i).GetNeighbors())
    def _norm(a):
        return math.atan2(math.sin(a), math.cos(a))
    def _terminal(i):                                    # carbon/atom with a single heavy neighbour (a cap, e.g. -CH3)
        return sum(1 for nb in at(i).GetNeighbors() if nb.GetAtomicNum() > 1) == 1
    def _slots(u, occ):                                  # candidate bond directions (abs radians) given the parent bond
        if _angular(u):
            return [occ[0] + R(120), occ[0] - R(120)] if occ else [R(90), R(330), R(210)]
        if not occ:
            return [R(0), R(90), R(270), R(180)]
        if _terminal(u):                                 # a cap (-CH3): snap its H's to the page-axis comb; only fan if that would collide
            ux, uy = pos[u]
            comb = sorted([R(0), R(90), R(180), R(270)], key=lambda s: -abs(_norm(s - occ[0])))[:3]
            if all(not collides(ux + HEAVY_LEN*math.cos(s), uy + HEAVY_LEN*math.sin(s)) for s in comb):
                return comb
            return [occ[0] + R(180), occ[0] + R(120), occ[0] - R(120)]   # crowded hanging methyl -> fan away from the parent
        return [occ[0] + R(180), occ[0] + R(90), occ[0] - R(90)]
    q = collections.deque()
    root = spine[0] if spine else next(iter(adj), 0)
    pos[root] = (0.0, 0.0); q.append(root)
    while q:
        u = q.popleft(); ux, uy = pos[u]
        occ = [math.atan2(pos[n.GetIdx()][1]-uy, pos[n.GetIdx()][0]-ux)
               for n in at(u).GetNeighbors() if n.GetIdx() in pos]
        unplaced = [nb.GetIdx() for nb in at(u).GetNeighbors() if nb.GetIdx() not in pos]
        if not unplaced:
            continue
        slots = _slots(u, occ)
        def take(key):                                   # pop the free, non-crowding slot that best matches key(); returns (angle, x, y)
            L = HEAVY_LEN
            for s in sorted(slots, key=key):
                if all(abs(_norm(s-t)) > R(48) for t in occ) and not collides(ux+L*math.cos(s), uy+L*math.sin(s)):
                    slots.remove(s); return s, ux+L*math.cos(s), uy+L*math.sin(s)
            for s in sorted(slots, key=key):             # relax crowding, still avoid a hard overlap
                if not collides(ux+L*math.cos(s), uy+L*math.sin(s)):
                    slots.remove(s); return s, ux+L*math.cos(s), uy+L*math.sin(s)
            for r in (L, L*1.25, L*1.5):                 # last resort: sweep any direction/length so two atoms never stack
                for ad in range(0, 360, 15):
                    s = R(ad); x, y = ux+r*math.cos(s), uy+r*math.sin(s)
                    if not collides(x, y):
                        if slots:
                            slots.pop(0)
                        return s, x, y
            s = sorted(slots, key=key)[0]; slots.remove(s); return s, ux+L*math.cos(s), uy+L*math.sin(s)
        def prio(j):
            return 0 if j in carbonylO else (1 if j in spine_pos else (2 if at(j).GetAtomicNum() > 1 else 3))
        for j in sorted(unplaced, key=prio):
            if j in carbonylO:                                       # =O points up
                s, x, y = take(lambda s: -math.sin(s))
            elif j in spine_pos:                                     # backbone heads right (then up)
                s, x, y = take(lambda s: (-round(math.cos(s), 3), -math.sin(s)))
            elif at(j).GetAtomicNum() in (9, 17, 35, 53) or j in hydroxylO:   # halogen / -OH / -OH2+ substituent -> UP (NJC: C-chain horizontal, heteroatom vertical)
                s, x, y = take(lambda s: -math.sin(s))
            elif at(j).GetAtomicNum() > 1:                           # heavy substituent
                heavy_nb_u = sum(1 for nb in at(u).GetNeighbors() if nb.GetAtomicNum() > 1)
                if heavy_nb_u == 2 and occ and not _angular(u):      # chain link: continue the chain STRAIGHT AHEAD
                    straight = occ[0] + R(180)                       # (keeps an ethyl/propyl branch linear; leaves the 90-deg spot for a small H, so a CH3 never turns into a sibling group)
                    s, x, y = take(lambda s: abs(_norm(s - straight)))
                else:                                                # branch point: outward-right, else hangs down
                    s, x, y = take(lambda s: (-round(math.cos(s), 3), math.sin(s)))
            elif at(u).GetFormalCharge() != 0 and at(u).GetAtomicNum() == 6 \
                    and sum(1 for nb in at(u).GetNeighbors() if nb.GetAtomicNum() > 1) == 1:
                s, x, y = take(lambda s: -abs(math.sin(s)))           # charged terminus (carbocation): H's go UP/DOWN, leaving the reactive face open
            else:                                                    # H: fill the widest remaining gap
                s, x, y = take(lambda s: -min((abs(_norm(s-t)) for t in occ), default=math.pi))
            pos[j] = (x, y); occ.append(s)
            if j in condense:                                        # alkyl substituent -> ONE label (CH3, C2H5...), don't draw its atoms
                condlab[j] = _sub_label(molH, set(spine), j)
            else:
                q.append(j)
    atoms = {i: (condlab.get(i, _elem(at(i)) + _charge_suffix(at(i).GetFormalCharge())), pos[i][0], pos[i][1])
             for i in pos}
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), _BONDORD.get(b.GetBondType(), 1))
             for b in molH.GetBonds() if b.GetBeginAtomIdx() in pos and b.GetEndAtomIdx() in pos]
    return atoms, bonds, [], set()


def _displayed_natural_layout(mol):
    """Displayed formula at natural tetrahedral/trigonal bond angles, all bonds one uniform length
    (C-H drawn the same length as C-C / C-Br)."""
    Chem.Kekulize(mol, clearAromaticFlags=True)
    molH = Chem.AddHs(mol)
    P = _coords(molH)                                    # no H-shortening -> uniform bond lengths
    at = molH.GetAtomWithIdx
    atoms = {i: (_elem(at(i)) + _charge_suffix(at(i).GetFormalCharge()), P[i][0], P[i][1])
             for i in range(molH.GetNumAtoms())}
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), _BONDORD.get(b.GetBondType(), 1))
             for b in molH.GetBonds()]
    return atoms, bonds, [], set()


def _natural_grid_layout(mol, spine_atoms=None):
    """Natural form: the CARBON BACKBONE drawn out on a 90-deg grid (each backbone C shown with its H
    at right angles), while every SUBSTITUENT hanging off the backbone is a CONDENSED label (CH3, OH,
    NH2, Br, OCH3...). C=O / C=C drawn. 'Every bond shown but not at the real bond angle', per the NJC
    naming-diagram style (main chain out, substituents condensed).
    spine_atoms: force the backbone to exactly these atoms (e.g. just the C=C for a mechanism, so EVERY
    alkyl group -- both methyls of 2-methylpropene -- condenses to a CH3 label)."""
    Chem.Kekulize(mol, clearAromaticFlags=True)
    molH = Chem.AddHs(mol)
    at = molH.GetAtomWithIdx
    A = molH.GetNumAtoms()
    carbonylO = {i for i in range(A) if _is_carbonyl_O(at(i))}
    def _backbone(i):                                    # carbon, or an in-chain O/N (>=2 heavy neighbours, e.g. ester/ether -O-, secondary amine N)
        a = at(i)
        return a.GetAtomicNum() == 6 or (a.GetAtomicNum() in (7, 8) and i not in carbonylO
                and sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1) >= 2)
    badj = {i: [] for i in range(A) if _backbone(i)}
    for i in badj:
        badj[i] = [nb.GetIdx() for nb in at(i).GetNeighbors() if nb.GetIdx() in badj]
    def _spwt(n):                                        # keep carbonyls, heteroatoms AND C=C / C#C carbons on the backbone
        if any(nb.GetIdx() in carbonylO for nb in at(n).GetNeighbors()):
            return 2
        if at(n).GetAtomicNum() in (7, 8):
            return 1
        if any(molH.GetBondBetweenAtoms(n, nb.GetIdx()).GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)
               for nb in at(n).GetNeighbors()):
            return 1
        return 0
    if spine_atoms is not None:                          # caller forces the backbone (mechanism: just the functional-group bond)
        spine = [i for i in spine_atoms if i in badj]
    else:
        spine = _longest_path(badj, weight=_spwt,
                              edgewt=lambda a, b: 1 if molH.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() >= 2 else 0)
    if not spine:                                        # no backbone -> fall back to the 90-deg grid
        return _displayed_comb_layout(mol, grid=True)
    spineset = set(spine)
    def fw(i):                                           # a heteroatom/functional end pulls that end to the right
        return sum(2 for nb in at(i).GetNeighbors()
                   if nb.GetIdx() not in spineset and nb.GetAtomicNum() in (7, 8, 9, 17, 35, 53))
    tot = sum(fw(i) for i in spine)
    if tot and sum(k*fw(i) for k, i in enumerate(spine))/tot < (len(spine)-1)/2.0:
        spine = spine[::-1]
    pos = {c: (k*HEAVY_LEN, 0.0) for k, c in enumerate(spine)}
    atoms = {c: (_elem(at(c)) + _charge_suffix(at(c).GetFormalCharge()), pos[c][0], pos[c][1]) for c in spine}
    bonds = [(a, b, _BONDORD.get(molH.GetBondBetweenAtoms(a, b).GetBondType(), 1))
             for a, b in zip(spine, spine[1:])]
    hid, pid, reversible = -1, 20000, set()
    L = HEAVY_LEN                                        # ONE uniform bond length for backbone, H and substituents
    def _place(label, cx, cy, d):                        # push a wider label further out so the VISIBLE bond stays uniform
        a, s = _parse_label(label)
        hw = (_txt_w(a, FONT) + _suffix_w(s)) / 2.0 / U
        r = L + max(0.0, hw - 0.30)
        return cx + DIRS[d][0] * r, cy + DIRS[d][1] * r
    for k, c in enumerate(spine):
        cx, cy = pos[c]
        border = {nb.GetIdx(): _BONDORD.get(molH.GetBondBetweenAtoms(c, nb.GetIdx()).GetBondType(), 1)
                  for nb in at(c).GetNeighbors()}
        is_sp = 3 in border.values()                     # triple bond -> sp -> LINEAR (substituent continues horizontally)
        is_sp2 = (not is_sp) and 2 in border.values()    # double bond -> sp2 -> LEFT / TOP / RIGHT (=O up, never down)
        nH = sum(1 for nb in at(c).GetNeighbors() if nb.GetAtomicNum() == 1)
        outer = (['right'] if k == len(spine) - 1 else []) + (['left'] if k == 0 else [])   # 'right' before 'left' -> functional group to the RIGHT
        if is_sp:
            avail = list(outer); fb = 'right'
        elif is_sp2:                                      # sp2: =O straight UP (incl. aldehydes), other groups left/top/right
            avail = ['up'] + outer
            fb = 'right'
        else:                                            # sp3: substituents up & down
            avail = ['up', 'down'] + outer
            fb = 'down'
        subs = [nb.GetIdx() for nb in at(c).GetNeighbors()
                if nb.GetAtomicNum() > 1 and nb.GetIdx() not in spineset]
        subs.sort(key=lambda j: (border[j] < 2, j))      # multiple-bonded (=O, #N) placed first (chain terminus / up / linear)
        for j in subs:
            d = avail.pop(0) if avail else fb
            if border[j] >= 2:                           # double/triple-bonded substituent -> draw the atom (=O, #N)
                lab = _elem(at(j)) + _charge_suffix(at(j).GetFormalCharge())
                x, y = _place(lab, cx, cy, d); atoms[j] = (lab, x, y); bonds.append((c, j, border[j]))
            else:                                        # single-bonded branch -> CONDENSED label (CH3, OH, Br, OCH3...)
                lab = _sub_label(molH, spineset, j)
                x, y = _place(lab, cx, cy, d); atoms[pid] = (lab, x, y); bonds.append((c, pid, 1))
                if DIRS[d][0] < -0.3:
                    reversible.add(pid)
                pid += 1
        for _ in range(nH):                              # explicit H on the backbone atom (uniform visible length)
            d = avail.pop(0) if avail else fb
            x, y = _place("H", cx, cy, d); atoms[hid] = ("H", x, y); bonds.append((c, hid, 1)); hid -= 1
    return atoms, bonds, [], reversible


def _displayed_layout(mol, angles='comb'):
    if mol.GetRingInfo().NumRings() > 0:                  # ring: hexagon+circle; natural condenses substituents, displayed expands them
        return _ring_layout(mol, 'structural' if angles == 'grid' else 'displayed')
    if angles == 'grid':                                 # 'natural' form: backbone out on a 90-deg grid, substituents condensed
        return _natural_grid_layout(mol)
    if angles == 'tetrahedral':                          # RDKit natural tetrahedral fan (kept, not the default 'natural')
        return _displayed_natural_layout(mol)
    return _displayed_comb_layout(mol)


def _skeletal_layout(mol):
    """Skeletal formula: carbons implicit at vertices (blank labels), C-H implicit, heteroatoms
    labelled with their hydrogens (OH, NH2), carbonyl =O drawn."""
    if mol.GetRingInfo().NumRings() > 0:
        return _ring_layout(mol, 'skeletal')              # skeletal ring == hexagon+circle; alkyl substituents = bare stubs
    Chem.Kekulize(mol, clearAromaticFlags=True)
    at = mol.GetAtomWithIdx
    A = mol.GetNumAtoms()
    carbonylO = set(_carbonyl_Os(mol))
    def _bb(i):                                           # backbone = carbons + in-chain (>=2 heavy) O/N
        a = at(i)
        return a.GetAtomicNum() == 6 or (a.GetAtomicNum() in (7, 8) and i not in carbonylO
                and sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1) >= 2)
    badj = {i: [nb.GetIdx() for nb in at(i).GetNeighbors() if _bb(nb.GetIdx())] for i in range(A) if _bb(i)}
    spine = _longest_path(badj,                          # main chain: keep unsaturation + functional group ON it
                          weight=lambda n: 2 if any(nb.GetIdx() in carbonylO for nb in at(n).GetNeighbors())
                                          else (1 if at(n).GetAtomicNum() in (7, 8) else 0),
                          edgewt=lambda a, b: 1 if mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() >= 2 else 0)
    from rdkit.Geometry import Point2D
    rdDepictor.SetPreferCoordGen(False)
    if len(spine) >= 2:                                   # TEMPLATE the main chain onto a horizontal zig-zag; RDKit fills the branches around it
        spineset = set(spine)
        def _ew(i):                                       # order the chain so the functional end is on the RIGHT
            w = _fweight(at(i))
            for nb in at(i).GetNeighbors():
                if nb.GetIdx() not in spineset and nb.GetAtomicNum() in (7, 8, 9, 17, 35, 53):
                    w += 2
            return w
        fwl = [_ew(i) for i in spine]
        tot = sum(fwl)
        if tot and sum(k * w for k, w in enumerate(fwl)) / tot < (len(spine) - 1) / 2.0:
            spine = spine[::-1]
        SL = 1.5
        cmap = {c: Point2D(k * SL * math.cos(math.radians(30)), (k % 2) * SL * math.sin(math.radians(30)))
                for k, c in enumerate(spine)}
        try:
            rdDepictor.Compute2DCoords(mol, canonOrient=False, coordMap=cmap)
        except Exception:
            rdDepictor.Compute2DCoords(mol)
    else:
        rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    P = {i: (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(A)}
    hb = [math.hypot(P[b.GetBeginAtomIdx()][0] - P[b.GetEndAtomIdx()][0], P[b.GetBeginAtomIdx()][1] - P[b.GetEndAtomIdx()][1])
          for b in mol.GetBonds() if at(b.GetBeginAtomIdx()).GetAtomicNum() > 1 and at(b.GetEndAtomIdx()).GetAtomicNum() > 1]
    s = HEAVY_LEN / (statistics.median(hb) if hb else 1.5)
    P = {i: (x * s, y * s) for i, (x, y) in P.items()}
    cO = _carbonyl_Os(mol)                               # C=O points UP
    if cO and sum(P[o][1] - P[at(o).GetNeighbors()[0].GetIdx()][1] for o in cO) < 0:
        P = {i: (x, -y) for i, (x, y) in P.items()}
    atoms = {}
    reversible = set()
    for i in range(mol.GetNumAtoms()):
        a = mol.GetAtomWithIdx(i); Z = a.GetAtomicNum(); nH = a.GetTotalNumHs()
        ch = _charge_suffix(a.GetFormalCharge())
        heavy_nb = sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1)
        if Z == 6:
            lab = ""                                     # carbon vertex: implicit, no C/H drawn
        elif Z == 8:
            lab = ("OH" if nH == 1 else "O") + ch
        elif Z == 7:
            if nH >= 1 and heavy_nb >= 2:
                lab = "N" + ch                           # in-chain N-H: bonds meet N; the H is drawn ABOVE (below)
            else:
                lab = "N" + ("H" if nH >= 1 else "") + (str(nH) if nH >= 2 else "") + ch
        else:
            lab = _elem(a) + ch
        atoms[i] = (lab, P[i][0], P[i][1])
        if heavy_nb == 1 and nH >= 1 and Z in (7, 8):    # terminal -OH / -NH2: flip to HO-/H2N- when the bond enters from the right
            reversible.add(i)
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), _BONDORD.get(b.GetBondType(), 1))
             for b in mol.GetBonds()]
    hid = -1                                              # secondary-amine N-H: put the H ON TOP of the N (not inline)
    for i in range(mol.GetNumAtoms()):
        a = mol.GetAtomWithIdx(i)
        if a.GetAtomicNum() == 7 and a.GetTotalNumHs() >= 1 and \
                sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1) >= 2:
            nHn = a.GetTotalNumHs()
            for h in range(nHn):
                dx = (h - (nHn - 1) / 2.0) * 0.75
                atoms[hid] = ("H", P[i][0] + dx, P[i][1] + 1.25); bonds.append((i, hid, 1)); hid -= 1
    return atoms, bonds, [], reversible


def _stereo_layout(mol, mirror=False, star=False):
    """One enantiomer at a chiral centre in FLYING-WEDGE (natta) style, tuned by overlaying the NJC
    figure: letter 'C' centre; the LARGEST group points straight UP (90 deg, on the mirror axis); the
    next carbon is a plain bond DOWN-RIGHT (-42 deg); the HASHED DASH (usually H) points DOWN-LEFT
    (197 deg, into the plane); the bold WEDGE (the heteroatom / functional group) points DOWN-LEFT
    (-122 deg, out of the plane). Left-side simple group labels reverse so the bonded atom faces the
    bond (OH->HO, CH3->H3C). `mirror` reflects to the other enantiomer. Optional * centre (star=True)."""
    centres = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False, includeUnassigned=True)
    if not centres:
        raise ValueError("no chiral centre: the stereochemistry form needs a molecule with a stereocentre")
    c = centres[0][0]
    molH = Chem.AddHs(mol)
    atH = molH.GetAtomWithIdx

    def branch_size(n):                                  # heavy-atom count of the branch rooted at n
        seen = {c, n}; stack = [n]; cnt = 0
        while stack:
            u = stack.pop(); cnt += 1
            for nb in atH(u).GetNeighbors():
                if nb.GetIdx() not in seen and nb.GetAtomicNum() > 1:
                    seen.add(nb.GetIdx()); stack.append(nb.GetIdx())
        return cnt

    nbrs = [n.GetIdx() for n in atH(c).GetNeighbors()]
    carbons = sorted((n for n in nbrs if atH(n).GetAtomicNum() == 6), key=branch_size, reverse=True)
    hetero = sorted((n for n in nbrs if atH(n).GetAtomicNum() not in (1, 6)),
                    key=lambda n: atH(n).GetAtomicNum(), reverse=True)
    hydro = [n for n in nbrs if atH(n).GetAtomicNum() == 1]
    ordered = carbons + hetero + hydro
    up, dr = ordered[0], ordered[1]                      # LARGEST group straight UP; next carbon plain down-right
    rest = sorted(ordered[2:], key=lambda n: atH(n).GetAtomicNum(), reverse=True)
    wedge, dash = rest[0], rest[1]                        # wedge = heavier group (heteroatom), dash = lighter (usually H)
    L = 1.7
    ang = {up: 90, dr: -42, dash: 197, wedge: -122}      # tuned by overlaying the NJC figure (main group UP)
    P = {c: (0.0, 0.0)}
    for n, a in ang.items():
        LL = L * 0.85 if n == dash else (L * 0.82 if n == up else L)   # top bond a touch shorter; dash closer in
        P[n] = (LL * math.cos(math.radians(a)), LL * math.sin(math.radians(a)))
    if mirror:
        P = {k: (-x, y) for k, (x, y) in P.items()}      # single inversion -> the other enantiomer
    def lbl(n):
        return "H" if atH(n).GetAtomicNum() == 1 else _sub_label(molH, {c}, n).replace("COOH", "CO2H")
    atoms = {c: ("C^*" if star else "C", P[c][0], P[c][1])}
    for n in (up, dr, dash, wedge):
        atoms[n] = (lbl(n), P[n][0], P[n][1])
    bonds = [(c, up, 1), (c, dr, 1), (c, dash, 'dash'), (c, wedge, 'wedge')]
    reversible = {n for n in (up, dr, dash, wedge)        # left-side simple groups face their bond (CH3->H3C, OH->HO)
                  if P[n][0] < -0.30 and branch_size(n) <= 1}
    return atoms, bonds, [], reversible


def _geometric_layout(mol, trans=False):
    """cis/trans (geometric) isomerism across a stereogenic C=C, per the mined notes: horizontal C=C
    (explicit letters, double line), four PLAIN substituent bonds splayed ~120 deg. Reference groups
    same-side = cis, diagonal = trans. No wedge/dash (the alkene unit is planar); never E/Z."""
    Chem.Kekulize(mol, clearAromaticFlags=True)
    db = next(((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
               if b.GetBondType() == Chem.BondType.DOUBLE and not b.IsInRing()
               and b.GetBeginAtom().GetAtomicNum() == 6 and b.GetEndAtom().GetAtomicNum() == 6), None)
    if db is None:
        raise ValueError("no open-chain C=C: cis/trans needs a stereogenic double bond")
    c1, c2 = db
    molH = Chem.AddHs(mol)
    atH = molH.GetAtomWithIdx
    def subs(cx, other):
        return sorted((n.GetIdx() for n in atH(cx).GetNeighbors() if n.GetIdx() != other),
                      key=lambda n: atH(n).GetAtomicNum(), reverse=True)
    s1, s2 = subs(c1, c2), subs(c2, c1)
    def slbl(n):
        return "H" if atH(n).GetAtomicNum() == 1 else _sub_label(molH, {c1, c2}, n)
    if slbl(s1[0]) == slbl(s1[1]) or slbl(s2[0]) == slbl(s2[1]):
        raise ValueError("no cis/trans: an alkene carbon bears two identical groups")
    ref1, oth1 = s1[0], s1[1]
    ref2, oth2 = s2[0], s2[1]
    P = {c1: (-0.75, 0.0), c2: (0.75, 0.0),
         ref1: (-1.75, 1.0), oth1: (-1.75, -1.0)}
    P[ref2], P[oth2] = ((1.75, -1.0), (1.75, 1.0)) if trans else ((1.75, 1.0), (1.75, -1.0))
    atoms = {c1: ("C", P[c1][0], P[c1][1]), c2: ("C", P[c2][0], P[c2][1])}
    for n in (ref1, oth1, ref2, oth2):
        atoms[n] = (slbl(n), P[n][0], P[n][1])
    bonds = [(c1, c2, 2), (c1, ref1, 1), (c1, oth1, 1), (c2, ref2, 1), (c2, oth2, 1)]
    return atoms, bonds, [], set()


def layout(smiles, form='structural', angles='comb'):
    """Layout for one of three representations: 'structural' (condensed, NJC house default),
    'displayed' (every atom + bond at 'natural' bond angles by default, or right-angle 'comb'),
    or 'skeletal' (zig-zag, implicit C/H)."""
    mol = Chem.MolFromSmiles(smiles)
    if form == 'displayed':
        return _displayed_layout(mol, angles)
    if form in ('natural', 'expanded', 'full'):          # every atom + bond drawn on a 90-deg right-angle GRID (no trigonal/bent)
        return _displayed_layout(mol, 'grid')
    if form == 'mechanism':                              # functional group drawn out at TRIGONAL angles (like ethene) + EVERY alkyl condensed
        fg = next(([b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in mol.GetBonds()
                   if b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE) and not b.GetIsAromatic()), None)
        if fg:
            condense = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1 and a.GetIdx() not in fg}
            return _displayed_comb_layout(mol, spine_atoms=fg, condense=condense)
        return _displayed_layout(mol, angles)
    if form == 'skeletal':
        return _skeletal_layout(mol)
    if form == 'stereo':
        return _stereo_layout(mol, angles == 'mirror')
    if form == 'geometric':
        return _geometric_layout(mol, angles == 'trans')
    if mol.GetRingInfo().NumRings() > 0:
        return _ring_layout(mol)
    # CONDENSED structural formula: heavy atoms only, H absorbed into labels (CH3, CH2, OH, NH2…)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    A = mol.GetNumAtoms()
    at = mol.GetAtomWithIdx
    carbonylO = {i for i in range(A) if _is_carbonyl_O(at(i))}
    # fold a terminal halogen on a carbon into that carbon's label (condensed: CH2Br, CHCl2, CCl3)
    absorbed_on = collections.defaultdict(list)          # carbon idx -> [halogen idx]
    absorbedX = set()
    for i in range(A):
        a = at(i)
        if a.GetAtomicNum() in (9, 17, 35, 53) and a.GetDegree() == 1 \
                and a.GetNeighbors()[0].GetAtomicNum() == 6:
            cnb = a.GetNeighbors()[0]
            if any(nb.GetIdx() in carbonylO for nb in cnb.GetNeighbors()):
                continue                                 # don't fold a halogen onto an acyl carbon -> draw -C(=O)-Cl, not "CCl"
            absorbed_on[cnb.GetIdx()].append(i); absorbedX.add(i)
    _skip = carbonylO | absorbedX

    def _halstr(c_idx):                                  # folded halogens on a carbon, ordered F, Cl, Br, I
        z = {'F': 9, 'Cl': 17, 'Br': 35, 'I': 53}
        cnt = collections.Counter(at(x).GetSymbol() for x in absorbed_on.get(c_idx, []))
        return "".join(s + (str(n) if n > 1 else "")
                       for s, n in sorted(cnt.items(), key=lambda kv: z.get(kv[0], 99)))

    adj = {i: [] for i in range(A) if i not in _skip}
    for i in adj:
        for nb in at(i).GetNeighbors():
            if nb.GetIdx() not in _skip:
                adj[i].append(nb.GetIdx())
    spine = _longest_path(adj, weight=lambda n: 2 if any(nb.GetIdx() in carbonylO for nb in at(n).GetNeighbors())
                                          else (1 if at(n).GetAtomicNum() in (7, 8) else 0),
                          edgewt=lambda a, b: 1 if mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() >= 2 else 0)

    # orient: functional end (heteroatoms / carbonyl) on the RIGHT
    def fweight(idx):
        a = at(idx)
        w = 2 if a.GetAtomicNum() in (7, 8, 9, 17, 35, 53) else 0
        if any(nb.GetIdx() in carbonylO for nb in a.GetNeighbors()):
            w += 2
        if idx in absorbed_on:               # carbon bearing a folded halogen also pulls right
            w += 2
        return w
    tot = sum(fweight(i) for i in spine)
    if tot and spine:
        cw = sum(k * fweight(i) for k, i in enumerate(spine)) / tot
        if cw < (len(spine) - 1) / 2.0:
            spine = spine[::-1]

    def blen(i, j):       # carbonyl =O (and other off-backbone double bonds) drawn a little shorter;
        b = mol.GetBondBetweenAtoms(i, j)   # chain C=C stays full length (wide CH2 labels need it)
        return HEAVY_LEN * (0.82 if b and b.GetBondType() == Chem.BondType.DOUBLE else 1.0)

    pos = {}; occ = collections.defaultdict(set)
    xpos = 0.0
    for i in spine:                          # variable spacing: widen where a folded halogen makes a label wide
        pos[i] = (xpos, 0.0)
        xpos += HEAVY_LEN + (_suffix_w(_halstr(i)) / U if i in absorbed_on else 0.0)
    for a, b in zip(spine, spine[1:]):
        occ[a].add('right'); occ[b].add('left')

    def has_carbonyl(idx):
        return any(nb.GetIdx() in carbonylO for nb in at(idx).GetNeighbors())

    def collides(x, y):
        return any(abs(x - px2) < 0.95 and abs(y - py2) < 0.85 for px2, py2 in pos.values())

    def dir_order(u, nb_idx):
        if nb_idx in carbonylO:
            return ['up', 'down', 'right', 'left']
        if has_carbonyl(u):
            return ['right', 'left', 'down', 'up']
        return ['up', 'down', 'right', 'left']

    q = collections.deque(spine)
    while q:
        u = q.popleft(); ux, uy = pos[u]
        for nb in at(u).GetNeighbors():
            j = nb.GetIdx()
            if j in pos or j in absorbedX:
                continue
            placed = False
            L = blen(u, j)
            for d in dir_order(u, j):
                if d in occ[u]:
                    continue
                x, y = ux + DIRS[d][0] * L, uy + DIRS[d][1] * L
                if collides(x, y):
                    continue
                pos[j] = (x, y); occ[u].add(d); occ[j].add(OPP[d]); placed = True
                break
            if not placed:
                for dx, dy in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
                    x, y = ux + dx * HEAVY_LEN * 0.85, uy + dy * HEAVY_LEN * 0.85
                    if not collides(x, y):
                        pos[j] = (x, y); break
                else:
                    pos[j] = (ux + HEAVY_LEN, uy + HEAVY_LEN)
            q.append(j)

    # show H explicitly on carbonyl carbons (aldehydes incl. methanal / methanoic acid)
    carbonylC = {i for i in pos if at(i).GetAtomicNum() == 6
                 and any(nb.GetIdx() in carbonylO for nb in at(i).GetNeighbors())}
    extraH = []          # (hid, parent)
    hid = -1
    for c in list(carbonylC):
        cx, cy = pos[c]
        for _ in range(at(c).GetTotalNumHs()):
            for d in ('right', 'left', 'down', 'up'):
                if d in occ[c]:
                    continue
                x, y = cx + DIRS[d][0] * 1.15, cy + DIRS[d][1] * 1.15
                if collides(x, y):
                    continue
                pos[hid] = (x, y); occ[c].add(d); extraH.append((hid, c)); hid -= 1
                break

    def lab(i):
        if i < 0:
            return "H"
        a = at(i); Z = a.GetAtomicNum(); nH = a.GetTotalNumHs()
        ch = _charge_suffix(a.GetFormalCharge())
        if Z == 6:
            xstr = _halstr(i)
            if i in carbonylC:
                return "C" + xstr + ch           # H drawn explicitly, =O up
            return "C" + ("H" if nH >= 1 else "") + (str(nH) if nH >= 2 else "") + xstr + ch
        if Z == 8:
            return ("OH" if nH == 1 else "O") + ch
        if Z == 7:
            return "N" + ("H" if nH >= 1 else "") + (str(nH) if nH >= 2 else "") + ch
        if Z == 0:                               # dummy = R-group placeholder
            return _RGROUP.get(a.GetAtomMapNum(), "R")
        return a.GetSymbol() + ch

    atoms = {i: (lab(i), pos[i][0], pos[i][1]) for i in pos}
    bonds = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in pos and j in pos:
            o = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2,
                 Chem.BondType.TRIPLE: 3}.get(b.GetBondType(), 1)
            bonds.append((i, j, o))
    for hid_, c in extraH:
        bonds.append((c, hid_, 1))

    # merge runs of fully-folded alkyl carbons (H/halogen only) into one condensed text label:
    # CH3 + CHBr -> "CH3CHBr", so the skeleton reads as condensed text and only functional groups are drawn.
    def _foldable(i):
        aa = at(i)
        return (aa.GetAtomicNum() == 6 and i not in carbonylC
                and not any(b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE) for b in aa.GetBonds())
                and all(nb.GetAtomicNum() in (1, 6, 9, 17, 35, 53) for nb in aa.GetNeighbors()))
    alk = {i for i in pos if i >= 0 and _foldable(i)}
    hadj = collections.defaultdict(set)
    for i, j, o in bonds:
        if o == 1 and i >= 0 and j >= 0:
            hadj[i].add(j); hadj[j].add(i)
    done = set()
    for start in sorted(alk, key=lambda i: pos[i][0]):
        if start in done:
            continue
        run, stack = set(), [start]
        while stack:
            u = stack.pop()
            if u in run:
                continue
            run.add(u)
            stack += [v for v in hadj[u] if v in alk and v not in run]
        done |= run
        if len(run) < 2 or any(len([v for v in hadj[u] if v in run]) > 2 for u in run):
            continue                                          # only simple chains of >=2 fold together
        ends = {u for u in run if len([v for v in hadj[u] if v in run]) <= 1}
        if any(u not in ends and any(v not in run for v in hadj[u]) for u in run):
            continue                                          # an INTERNAL atom has an external bond -> folding would misplace it (e.g. (CH3)2CH-CHO)
        order = sorted(run, key=lambda i: pos[i][0])          # left-to-right reading order
        keep, drop = order[0], set(order[1:])
        atoms[keep] = ("".join(atoms[i][0] for i in order), pos[keep][0], pos[keep][1])
        for i in drop:
            atoms.pop(i, None)
        merged = []
        for i, j, o in bonds:
            if i in drop and j in drop:
                continue
            i = keep if i in drop else i
            j = keep if j in drop else j
            if i != j:
                merged.append((i, j, o))
        bonds = merged
    # terminal -OH / -NH2 flips to HO- / H2N- when the bond enters from the right (renderer gates on direction)
    reversible = {i for i in atoms if i >= 0 and at(i).GetAtomicNum() in (7, 8)
                  and at(i).GetTotalNumHs() >= 1
                  and sum(1 for nb in at(i).GetNeighbors() if nb.GetAtomicNum() > 1) == 1}
    return atoms, bonds, [], reversible


SUBFONT = ImageFont.truetype(FONT_PATH, 23)
_SCRATCH = ImageDraw.Draw(Image.new("RGB", (8, 8)))
_BASE_H = _SCRATCH.textbbox((0, 0), "H", font=FONT)[3] - _SCRATCH.textbbox((0, 0), "H", font=FONT)[1]

def _label_parts(text):
    return re.findall(r'[^0-9]+|[0-9]+', text)

_RGROUP = {0: "R", 1: "R", 2: "R'", 3: "R''", 4: "R'''", 9: "X", 8: "X'"}

def _charge_suffix(fc):
    """Encode a formal charge as a superscript token: +1->'^+', -1->'^-', -2->'^2-'."""
    if not fc:
        return ""
    return "^" + (str(abs(fc)) if abs(fc) > 1 else "") + ("+" if fc > 0 else "-")

def _label_runs(suffix):
    """Split a label suffix into (text, kind) runs; kind in 'n','sub','sup'.
    Bare digits -> subscript; '^…' -> superscript (digits then +/- , or one char)."""
    runs, i, n = [], 0, len(suffix)
    while i < n:
        c = suffix[i]
        if c == '^':
            i += 1
            m = re.match(r'[0-9]*[+\-−]|[+\-−]|.', suffix[i:]) if i < n else None
            s = m.group(0) if m else ''
            if s:
                runs.append((s, 'sup')); i += len(s)
            continue
        if c in '0123456789':
            m = re.match(r'[0-9]+', suffix[i:]); s = m.group(0)
            runs.append((s, 'sub')); i += len(s); continue
        m = re.match(r'[^0-9^]+', suffix[i:]); s = m.group(0)
        runs.append((s, 'n')); i += len(s)
    return runs

def _parse_label(text):
    """Split into (anchor element, suffix). Anchor is the bonded atom (touches the bond)."""
    if text[:2] in ('Cl', 'Br'):
        return text[:2], text[2:]
    return text[:1], text[1:]

def _txt_w(text, font):
    tb = _SCRATCH.textbbox((0, 0), text, font=font); return tb[2] - tb[0]

def _suffix_w(suffix):
    return sum(_txt_w(s, SUBFONT if k in ('sub', 'sup') else FONT) for s, k in _label_runs(suffix))

def _draw_group(dr, anchor, suffix, side, cx, cy):
    aw = _txt_w(anchor, FONT)
    top = cy - _BASE_H / 2 - _SCRATCH.textbbox((0, 0), "H", font=FONT)[1]
    dr.text((cx - aw / 2, top), anchor, fill=(0, 0, 0), font=FONT)   # anchor centred on (cx,cy)
    sw = _suffix_w(suffix)
    x = cx + aw / 2 if side == 'right' else cx - aw / 2 - sw
    for s, kind in _label_runs(suffix):
        if kind == 'sub':
            dr.text((x, top + _BASE_H * 0.46), s, fill=(0, 0, 0), font=SUBFONT); x += _txt_w(s, SUBFONT)
        elif kind == 'sup':
            dr.text((x, top - _BASE_H * 0.30), s, fill=(0, 0, 0), font=SUBFONT); x += _txt_w(s, SUBFONT)
        else:
            dr.text((x, top), s, fill=(0, 0, 0), font=FONT); x += _txt_w(s, FONT)


def draw(atoms, bonds, circles=None, reversible=None, fname=None, return_map=False):
    circles = circles or []
    reversible = reversible or set()
    nbx = collections.defaultdict(list)
    for i, j, o in bonds:
        nbx[i].append(atoms[j][1] - atoms[i][1]); nbx[j].append(atoms[i][1] - atoms[j][1])
    parsed, side, ext = {}, {}, {}   # ext: (left_px, right_px, vert_px) from anchor centre
    for i, (label, x, y) in atoms.items():
        if not label:
            continue
        a, s = _parse_label(label); parsed[i] = (a, s)
        has_r = any(dx > 0.3 for dx in nbx[i]); has_l = any(dx < -0.3 for dx in nbx[i])
        # only reverse (anchor→right, e.g. O2N–) for explicitly reversible atoms (aromatic substituents)
        sd = 'left' if (i in reversible and has_r and not has_l) else 'right'
        side[i] = sd
        aw = _txt_w(a, FONT); sw = _suffix_w(s)
        ext[i] = (aw / 2 + (sw if sd == 'left' else 0),
                  aw / 2 + (sw if sd == 'right' else 0), _BASE_H / 2)

    def exu(i, k):  # extent in coord units
        return (ext[i][k] / U + 0.05) if i in ext else 0
    minx = min(x - exu(i, 0) for i, (l, x, y) in atoms.items())
    maxx = max(x + exu(i, 1) for i, (l, x, y) in atoms.items())
    miny = min(y - exu(i, 2) for i, (l, x, y) in atoms.items())
    maxy = max(y + exu(i, 2) for i, (l, x, y) in atoms.items())
    M = 40
    W = int((maxx - minx) * U + 2 * M); H = int((maxy - miny) * U + 2 * M)
    img = Image.new("RGB", (W, H), "white"); dr = ImageDraw.Draw(img)

    def px(i):
        _, x, y = atoms[i]; return ((x - minx) * U + M, (maxy - y) * U + M)
    def pxc(x, y):
        return ((x - minx) * U + M, (maxy - y) * U + M)
    def gap(i, ux, uy):
        if i not in ext:
            return 0          # empty-label atoms (ring vertices) → no gap, so the hexagon closes
        hx = ext[i][1] if ux > 0 else ext[i][0]
        return abs(ux) * hx + abs(uy) * ext[i][2] + 4

    def edge_gap(i, ux, uy):   # distance to the label's box edge along (ux,uy): a wedge/dash tip reaches
        if i not in ext:       # the label instead of being eaten by a wide suffix on a mostly-vertical bond
            return 0
        hx = ext[i][1] if ux > 0 else ext[i][0]
        tx = hx / abs(ux) if abs(ux) > 1e-3 else 1e9
        ty = ext[i][2] / abs(uy) if abs(uy) > 1e-3 else 1e9
        return min(tx, ty) + 3

    for i, j, order in bonds:
        p, q = px(i), px(j)
        dx, dy = q[0] - p[0], q[1] - p[1]; L = math.hypot(dx, dy) or 1
        ux, uy = dx / L, dy / L
        gj = edge_gap(j, -ux, -uy) if order in ('wedge', 'dash') else gap(j, -ux, -uy)
        p = (p[0] + ux * gap(i, ux, uy), p[1] + uy * gap(i, ux, uy))
        q = (q[0] - ux * gj, q[1] - uy * gj)
        ox, oy = -uy, ux
        if order == 'wedge':                             # bond out of the plane: solid filled triangle
            wl = 9
            dr.polygon([(p[0], p[1]), (q[0] + ox * wl, q[1] + oy * wl),
                        (q[0] - ox * wl, q[1] - oy * wl)], fill=(0, 0, 0))
            continue
        if order == 'dash':                              # bond into the plane: hashed lines widening outward
            dx2, dy2, nd = q[0] - p[0], q[1] - p[1], 6
            for t in range(1, nd + 1):
                f = t / (nd + 1); cx, cy = p[0] + dx2 * f, p[1] + dy2 * f; wl = 1.5 + 8 * f
                dr.line([(cx + ox * wl, cy + oy * wl), (cx - ox * wl, cy - oy * wl)],
                        fill=(0, 0, 0), width=STROKE - 1)
            continue
        offs = {1: [0], 2: [-DBL, DBL], 3: [-TRP, 0, TRP]}[order]
        for s in offs:
            dr.line([(p[0] + ox * s, p[1] + oy * s), (q[0] + ox * s, q[1] + oy * s)],
                    fill=(0, 0, 0), width=STROKE)

    for cx, cy, r in circles:
        c = pxc(cx, cy); rr = r * U
        dr.ellipse([c[0] - rr, c[1] - rr, c[0] + rr, c[1] + rr], outline=(0, 0, 0), width=STROKE)

    for i in parsed:
        cx, cy = px(i)
        _draw_group(dr, parsed[i][0], parsed[i][1], side[i], cx, cy)
    if fname:
        img.save(fname)
    if return_map:
        return img, pxc, W, H          # pxc(mol_x, mol_y) -> (px, py); lets an overlay (mechanism arrows) register to atoms/bonds
    return img


def render(smiles, fname=None, form='structural', angles='comb'):
    out = layout(smiles, form, angles)
    a, b, c = out[0], out[1], out[2]
    rev = out[3] if len(out) > 3 else set()
    return draw(a, b, c, rev, fname)


def render_stereo_pair(smiles, fname=None):
    """The two enantiomers, drawn in ONE shared coordinate frame (identical scale, symmetric about
    the mirror line), split by a vertical dash-dot MIRROR-PLANE line and captioned per the NJC figure."""
    mol = Chem.MolFromSmiles(smiles)
    left = draw(*_stereo_layout(mol, mirror=False))      # equal scale (same U); centre each in an equal cell
    right = draw(*_stereo_layout(mol, mirror=True))
    cellW = max(left.width, right.width) + 24
    cellH = max(left.height, right.height)
    GAP, TOP, BOT = 60, 52, 58
    W, H = cellW * 2 + GAP, cellH + TOP + BOT
    canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)
    canvas.paste(left, ((cellW - left.width) // 2, TOP + (cellH - left.height) // 2))
    canvas.paste(right, (cellW + GAP + (cellW - right.width) // 2, TOP + (cellH - right.height) // 2))
    mx = cellW + GAP // 2
    y = TOP - 4                                          # dash-dot vertical mirror-plane line
    while y < TOP + cellH + 4:
        d.line([(mx, y), (mx, min(y + 8, TOP + cellH + 4))], fill=(0, 0, 0), width=1)
        d.ellipse([mx - 1, y + 11, mx + 1, y + 13], fill=(0, 0, 0)); y += 16
    f = ImageFont.truetype(FONT_PATH, 20)
    def centred(text, yy, cx=mx):
        tb = d.textbbox((0, 0), text, font=f); d.text((cx - (tb[2] - tb[0]) // 2, yy), text, fill=(0, 0, 0), font=f)
    centred("mirror", 3); centred("plane", 24)           # stacked label above the line
    centred("enantiomers /", H - BOT + 8, W // 2)
    centred("non-superimposable mirror images", H - BOT + 30, W // 2)
    if fname:
        canvas.save(fname)
    return canvas


def render_geometric_pair(smiles, fname=None):
    """The cis and trans geometric isomers side by side, labelled 'cis' / 'trans' (never E/Z)."""
    mol = Chem.MolFromSmiles(smiles)
    imgs = [("cis", draw(*_geometric_layout(mol, trans=False))),
            ("trans", draw(*_geometric_layout(mol, trans=True)))]
    if fname:
        tile(imgs, fname, cols=2)
    return imgs


def tile(named_imgs, fname, cols=4):
    f = ImageFont.truetype(FONT_PATH, 24)
    CW = max(im.width for _, im in named_imgs) + 40
    CH = max(im.height for _, im in named_imgs) + 56
    rows = (len(named_imgs) + cols - 1) // cols
    grid = Image.new("RGB", (cols * CW, rows * CH), "white"); d = ImageDraw.Draw(grid)
    for i, (name, im) in enumerate(named_imgs):
        ox, oy = (i % cols) * CW, (i // cols) * CH
        grid.paste(im, (ox + (CW - im.width)//2, oy + (CH - 36 - im.height)//2))
        tb = d.textbbox((0, 0), name, font=f)
        d.text((ox + (CW - (tb[2]-tb[0]))//2, oy + CH - 32), name, fill=(0,0,0), font=f)
    grid.save(fname); print("saved", fname, grid.size)


if __name__ == "__main__":
    tests = [
        ("ethanol", "CCO"), ("propan-2-ol", "CC(O)C"), ("propanone", "CC(C)=O"),
        ("ethanal", "CC=O"), ("ethanoic acid", "CC(=O)O"), ("propanoic acid", "CCC(=O)O"),
        ("chloroethane", "CCCl"), ("ethylamine", "CCN"), ("ethanenitrile", "CC#N"),
        ("methyl ethanoate", "COC(C)=O"), ("ethanamide", "CC(N)=O"),
        ("2-aminopropanoic acid", "CC(N)C(=O)O"),
    ]
    imgs = [(n, render(s)) for n, s in tests]
    tile(imgs, "/tmp/structure_panel.png")
