"""Displayed-formula generator.
Auto-layout from SMILES: horizontal backbone, C=O up, all H shown, Arial, monochrome.
Acyclic molecules (chains + one functional group + simple branches). Aromatic handled separately."""
import math, collections, re
from rdkit import Chem
from PIL import Image, ImageDraw, ImageFont

ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
CALIBRI = "/Applications/Microsoft Word.app/Contents/Resources/DFonts/Calibri.ttf"
FONT_PATH = ARIAL              # label font — Arial for clean, even structure labels
FONT_NAME = "Arial"            # font-family used in SVG output
U = 46; STROKE = 3; GAP = 15; DBL = 6; TRP = 7
HEAVY_LEN = 1.95; H_LEN = 1.0
FONT = ImageFont.truetype(FONT_PATH, 34)
DIRS = {'up': (0, 1), 'down': (0, -1), 'left': (-1, 0), 'right': (1, 0)}
OPP = {'up': 'down', 'down': 'up', 'left': 'right', 'right': 'left'}


def _is_carbonyl_O(a):
    return (a.GetAtomicNum() == 8 and a.GetDegree() == 1
            and any(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds()))


def _longest_path(adj):
    if not adj:
        return []
    def bfs(s):
        seen = {s: [s]}; q = collections.deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in seen:
                    seen[v] = seen[u] + [v]; q.append(v)
        far = max(seen, key=lambda k: len(seen[k]))
        return far, seen[far]
    a, _ = bfs(next(iter(adj)))
    _, path = bfs(a)
    return path


def _sub_label(molH, ringset, root):
    """Condensed label for a benzene substituent rooted at atom `root`."""
    a = molH.GetAtomWithIdx(root)
    Z = a.GetAtomicNum()
    nbrs = [n for n in a.GetNeighbors() if n.GetIdx() not in ringset]
    heavy = [n for n in nbrs if n.GetAtomicNum() > 1]
    nH = sum(1 for n in nbrs if n.GetAtomicNum() == 1)
    if Z == 17: return "Cl"
    if Z == 35: return "Br"
    if Z == 8:  return ("OH" if nH else "O") + _charge_suffix(a.GetFormalCharge())
    if Z == 7:
        nO = sum(1 for n in heavy if n.GetAtomicNum() == 8)
        if nO >= 2: return "NO2"
        return "NH2" if nH == 2 else ("NHR" if heavy else "NH2")
    if Z == 6:
        # what is C bonded to (besides ring)?
        dbO = any(b.GetBondType() == Chem.BondType.DOUBLE and
                  b.GetOtherAtom(a).GetAtomicNum() == 8 for b in a.GetBonds())
        trN = any(b.GetBondType() == Chem.BondType.TRIPLE for b in a.GetBonds())
        if trN: return "CN"
        if dbO:
            oth = [n for n in heavy if n.GetAtomicNum() == 8 and
                   molH.GetBondBetweenAtoms(root, n.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
            if oth: return "COOH"
            if any(n.GetAtomicNum() == 7 for n in heavy): return "CONH2"   # amide
            if any(n.GetAtomicNum() == 6 for n in heavy): return "COCH3"
            return "CHO"
        if nH == 3: return "CH3"
        if heavy and heavy[0].GetAtomicNum() == 6:
            return "CH2CH3"
        return "CH3"
    if Z == 0:
        return _RGROUP.get(a.GetAtomMapNum(), "R")
    return a.GetSymbol()


def _ring_layout(mol):
    """Single ring as a skeletal polygon (no C/H labels). Aromatic 6-rings get an inscribed circle;
    cycloalkanes are plain polygons. Substituents drawn as condensed labels radially outward."""
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
    reversible = set()
    for ri_atom, root in subs:
        th = angles[ri_atom]
        ra = mol.GetAtomWithIdx(root)
        dO = next((nb for nb in ra.GetNeighbors() if nb.GetAtomicNum() == 8 and
                   mol.GetBondBetweenAtoms(root, nb.GetIdx()).GetBondType() == Chem.BondType.DOUBLE), None)
        nN = next((nb for nb in ra.GetNeighbors() if nb.GetAtomicNum() == 7 and nb.GetIdx() not in ringset), None)
        # amide substituent → draw out as ring–C(=O)–NH2 (carbonyl explicit)
        if ra.GetAtomicNum() == 6 and dO is not None and nN is not None:
            rx, ry = math.cos(th), math.sin(th); pxr, pyr = -ry, rx
            vx, vy = R * math.cos(th), R * math.sin(th)
            fcx, fcy = vx + rx * 1.4, vy + ry * 1.4
            atoms[root] = ("C", fcx, fcy); bonds.append((ri_atom, root, 1))
            atoms[dO.GetIdx()] = ("O", fcx + pxr * 1.15, fcy + pyr * 1.15); bonds.append((root, dO.GetIdx(), 2))
            nH = nN.GetTotalNumHs()
            nlab = "N" + ("H" if nH >= 1 else "") + (str(nH) if nH >= 2 else "")
            atoms[nN.GetIdx()] = (nlab, fcx + rx * 1.5, fcy + ry * 1.5); bonds.append((root, nN.GetIdx(), 1))
            if rx < -0.3:
                reversible.add(nN.GetIdx())
        else:
            lab = _sub_label(molH, ringset, root)
            atoms[pid] = (lab, (R + 1.25) * math.cos(th), (R + 1.25) * math.sin(th))
            bonds.append((ri_atom, pid, 1)); reversible.add(pid); pid += 1
    return atoms, bonds, circles, reversible


def layout(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol.GetRingInfo().NumRings() > 0:
        return _ring_layout(mol)
    # CONDENSED structural formula: heavy atoms only, H absorbed into labels (CH3, CH2, OH, NH2…)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    A = mol.GetNumAtoms()
    at = mol.GetAtomWithIdx
    carbonylO = {i for i in range(A) if _is_carbonyl_O(at(i))}

    adj = {i: [] for i in range(A) if i not in carbonylO}
    for i in adj:
        for nb in at(i).GetNeighbors():
            if nb.GetIdx() not in carbonylO:
                adj[i].append(nb.GetIdx())
    spine = _longest_path(adj)

    # orient: functional end (heteroatoms / carbonyl) on the RIGHT
    def fweight(idx):
        a = at(idx)
        w = 2 if a.GetAtomicNum() in (7, 8, 9, 17, 35, 53) else 0
        if any(nb.GetIdx() in carbonylO for nb in a.GetNeighbors()):
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
    for k, i in enumerate(spine):
        pos[i] = (k * HEAVY_LEN, 0.0)
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
            if j in pos:
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
            if i in carbonylC:
                return "C" + ch                  # H drawn explicitly, =O up
            return "C" + ("H" if nH >= 1 else "") + (str(nH) if nH >= 2 else "") + ch
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
    return atoms, bonds, [], set()


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


def draw(atoms, bonds, circles=None, reversible=None, fname=None):
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

    for i, j, order in bonds:
        p, q = px(i), px(j)
        dx, dy = q[0] - p[0], q[1] - p[1]; L = math.hypot(dx, dy) or 1
        ux, uy = dx / L, dy / L
        p = (p[0] + ux * gap(i, ux, uy), p[1] + uy * gap(i, ux, uy))
        q = (q[0] - ux * gap(j, -ux, -uy), q[1] - uy * gap(j, -ux, -uy))
        ox, oy = -uy, ux
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
    return img


def render(smiles, fname=None):
    out = layout(smiles)
    a, b, c = out[0], out[1], out[2]
    rev = out[3] if len(out) > 3 else set()
    return draw(a, b, c, rev, fname)


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
