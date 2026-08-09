"""SVG mechanism renderer: electron-pushing curly arrows, lone pairs, partial charges, composed frames.

Sits ON TOP of structure_svg / structure_draw (which draw the compounds) and overlays the mechanism
annotations in the SAME vector space, addressed by atom/bond index. Output is sharp SVG (rasterise
with `render_frame_png` for preview / PNG fallback). A-level (9729) house style, arrow-pushing.

The curly-arrow construction is reverse-engineered from NJC's own PDF vector paths:
  * shaft  : a cubic Bezier that LEAVES the source bond perpendicular and ARRIVES aiming at the target;
             the tail starts just BELOW the bond it springs from (never inside a double bond).
  * head   : a FILLED, notched harpoon (a closed path) -- not two open barbs.
  * gap    : the head POINTS AT the target atom and stops short of its label (never overlaps it).
  * depth  : auto-computed by the router to clear the source compound's outermost substituent.

Frame spec (all locators are (side, species_index, locator); side in {'L','R'}):
    species  = {'smiles', 'form', 'orient':('bond',i,j)|None, 'label':str|None}
    locator  = ('atom', i) | ('bond', i, j) | ('pt', dx, dy)
    arrow    = {'from':loc_ref, 'to':loc_ref, 'depth':px|None, 'route':bool, 'bow':+/-1,
                'lift':float, 'stop':px, 'kind':'full'|'half', 'startperp':px}
    charge   = {'at':(side,idx,atom), 'text':str, 'dx':px, 'dy':px}
    lone_pair= {'at':(side,idx,atom), 'where':'up'|'down'|..., 'dist':px}
"""
import math, os, tempfile
from PIL import ImageFont
import structure_draw as J
import structure_svg as S
from structure_draw import U, HEAVY_LEN, FONT_PATH, _txt_w

ARROW_W = 3.0
HEAD_LEN = 17.0
HEAD_W = 6.5
NOTCH = 0.34

# ----------------------------------------------------------------------------- arrow primitive
def svg_arrow(p0, p1, depth, lift=0.30, stop=0.0, stroke=ARROW_W, kind='full',
              head_len=HEAD_LEN, head_w=HEAD_W, notch=NOTCH):
    """Curly electron-flow arrow from p0 to the TARGET centre p1.
    depth : signed arc depth (px); + bows to the LEFT of p0->p1, - to the right.
    lift  : slides the control points along the chord so the U isn't a hairpin.
    stop  : pull the tip back this many px from p1 (so the head points AT the atom with a gap).
    kind  : 'full' = 2-electron (double-barb harpoon); 'half' = 1-electron fish-hook (single barb)."""
    dx, dy = p1[0]-p0[0], p1[1]-p0[1]
    L = math.hypot(dx, dy) or 1.0
    ux, uy = dx/L, dy/L
    nx, ny = -uy, ux
    c1 = (p0[0]+nx*depth+ux*L*lift, p0[1]+ny*depth+uy*L*lift)
    c2 = (p1[0]+nx*depth-ux*L*lift, p1[1]+ny*depth-uy*L*lift)
    tx, ty = p1[0]-c2[0], p1[1]-c2[1]; tl = math.hypot(tx, ty) or 1.0; tx, ty = tx/tl, ty/tl
    tip = (p1[0]-tx*stop, p1[1]-ty*stop)
    base = (tip[0]-tx*head_len*0.82, tip[1]-ty*head_len*0.82)
    shaft = (f'<path d="M {p0[0]:.1f} {p0[1]:.1f} C {c1[0]:.1f} {c1[1]:.1f} '
             f'{c2[0]:.1f} {c2[1]:.1f} {base[0]:.1f} {base[1]:.1f}" fill="none" '
             f'stroke="black" stroke-width="{stroke}" stroke-linecap="round"/>')
    px, py = -ty, tx
    b1 = (tip[0]-tx*head_len+px*head_w, tip[1]-ty*head_len+py*head_w)
    b2 = (tip[0]-tx*head_len-px*head_w, tip[1]-ty*head_len-py*head_w)
    nk = (tip[0]-tx*head_len*(1-notch), tip[1]-ty*head_len*(1-notch))
    if kind == 'full':
        head = (f'<path d="M {tip[0]:.1f} {tip[1]:.1f} L {b1[0]:.1f} {b1[1]:.1f} '
                f'L {nk[0]:.1f} {nk[1]:.1f} L {b2[0]:.1f} {b2[1]:.1f} Z" fill="black"/>')
    else:
        head = (f'<path d="M {tip[0]:.1f} {tip[1]:.1f} L {b1[0]:.1f} {b1[1]:.1f} '
                f'L {nk[0]:.1f} {nk[1]:.1f} Z" fill="black"/>')
    return shaft + head

_DIR = {'up': (0, -1), 'down': (0, 1), 'left': (-1, 0), 'right': (1, 0),
        'upleft': (-0.7, -0.7), 'upright': (0.7, -0.7)}

def svg_lone_pair(cx, cy, where='up', dist=24, gap=9, r=4.0):
    dx, dy = _DIR[where] if isinstance(where, str) else where
    dn = math.hypot(dx, dy) or 1.0; dx, dy = dx/dn, dy/dn
    bx, by = cx + dx*dist, cy + dy*dist
    px, py = -dy, dx
    return "".join(f'<circle cx="{bx+px*gap*s:.1f}" cy="{by+py*gap*s:.1f}" r="{r}" fill="black"/>'
                   for s in (-1, 1))

# ----------------------------------------------------------------------------- geometry helpers
def _orient(a, i, j):
    """Rigid-rotate the atoms dict so bond i->j lies horizontal (angles preserved)."""
    xi, yi = a[i][1], a[i][2]
    ang = math.atan2(a[j][2]-yi, a[j][1]-xi)
    ca, sa = math.cos(-ang), math.sin(-ang)
    return {k: (lab, xi + (x-xi)*ca - (y-yi)*sa, yi + (x-xi)*sa + (y-yi)*ca)
            for k, (lab, x, y) in a.items()}

def _lay(sp):
    a, b, circ, rev = J.layout(sp['smiles'], sp.get('form', 'displayed'))
    o = sp.get('orient')
    if o and o[0] == 'bond':
        a = _orient(a, o[1], o[2])
    if sp.get('vflip'):                                    # mirror vertically (puts the other substituent on the arrow's side)
        ys = [y for _, x, y in a.values()]; y0 = (min(ys) + max(ys)) / 2.0
        a = {k: (l, x, 2*y0 - y) for k, (l, x, y) in a.items()}
    inner, W, H, px = S.draw_svg_inner(a, b, circ, rev)
    return dict(a=a, b=b, inner=inner, W=W, H=H, px=px, axis=px(0.0, 0.0)[1], label=sp.get('label'))

def _inset(p0, p1, d):
    dx, dy = p1[0]-p0[0], p1[1]-p0[1]; L = math.hypot(dx, dy) or 1.0
    return (p1[0]-dx/L*d, p1[1]-dy/L*d)

def _auto_face(sp, cc, tgt, sp_gap, startperp=14):
    """Pick (orient, vflip) for an alkene so its pi arrow to the electrophile `tgt` has the clearest lane.
    Searches both C=C orientations x vertical flip; returns the config whose best routable arc clears most."""
    best = None
    def apx(s, ox, oy, i):
        x, y = s['px'](s['a'][i][1], s['a'][i][2]); return (ox+x, oy+y)
    for orient in [tuple(cc), (cc[1], cc[0])]:
        for vflip in (False, True):
            trial = dict(sp); trial['orient'] = ('bond',) + orient; trial['vflip'] = vflip
            e, t = _lay(trial), _lay(tgt)
            axis = max(e['axis'], t['axis']); pad = 70
            oxE, oyE = pad, pad + (axis - e['axis'])
            oxB, oyB = pad + e['W'] + sp_gap, pad + (axis - t['axis'])
            c0, c1 = apx(e, oxE, oyE, cc[0]), apx(e, oxE, oyE, cc[1])
            p0 = ((c0[0]+c1[0])/2, (c0[1]+c1[1])/2 + startperp)
            nb = apx(t, oxB, oyB, 0)
            atoms = [(*apx(e, oxE, oyE, i), 15 if (not e['a'][i][0] or e['a'][i][0] == 'H') else 22) for i in e['a']]
            segs = [(apx(e, oxE, oyE, bi), apx(e, oxE, oyE, bj)) for (bi, bj, o) in e['b'] if {bi, bj} != set(cc)]
            clr = -1e9
            for dx in (0, -15, 15, -30, 30, -46, 46):
                for k in range(22):
                    depth = 0.5*HEAVY_LEN*U + k*0.13*HEAVY_LEN*U
                    md = _min_clear(_shaft_samples((p0[0]+dx, p0[1]), nb, depth, 0.05, 24), atoms, segs)
                    if md > clr: clr = md
                    if clr >= 13: break
                if clr >= 13: break
            if best is None or clr > best[0]:
                best = (clr, ('bond',)+orient, vflip)
    return best[1], best[2]

def _seg_dist(p, a, b):
    """Distance from point p to segment a-b."""
    ax, ay = a; bx, by = b; px, py = p
    dx, dy = bx-ax, by-ay
    L2 = dx*dx + dy*dy
    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((px-ax)*dx + (py-ay)*dy)/L2))
    qx, qy = ax + t*dx, ay + t*dy
    return math.hypot(px-qx, py-qy)

def _min_clear(pts, atom_pts, segs):
    """Min gap from a set of curve points to any obstacle (atom minus its radius, or bond segment)."""
    m = 1e9
    for (x, y) in pts:
        for (ax, ay, rad) in atom_pts:
            d = math.hypot(x-ax, y-ay) - rad
            if d < m: m = d
        for (a, b) in segs:
            d = _seg_dist((x, y), a, b)
            if d < m: m = d
    return m

def _shaft_samples(p0, p1, depth, lift, stop, head_len=HEAD_LEN, n=44):
    """Sample points along the SAME cubic shaft that svg_arrow() draws (for collision testing)."""
    dx, dy = p1[0]-p0[0], p1[1]-p0[1]; L = math.hypot(dx, dy) or 1.0
    ux, uy = dx/L, dy/L; nx, ny = -uy, ux
    c1 = (p0[0]+nx*depth+ux*L*lift, p0[1]+ny*depth+uy*L*lift)
    c2 = (p1[0]+nx*depth-ux*L*lift, p1[1]+ny*depth-uy*L*lift)
    tx, ty = p1[0]-c2[0], p1[1]-c2[1]; tl = math.hypot(tx, ty) or 1.0; tx, ty = tx/tl, ty/tl
    tip = (p1[0]-tx*stop, p1[1]-ty*stop)
    base = (tip[0]-tx*head_len*0.82, tip[1]-ty*head_len*0.82)
    out = []
    for k in range(n+1):
        t = k/n; m = 1-t
        out.append((m*m*m*p0[0] + 3*m*m*t*c1[0] + 3*m*t*t*c2[0] + t*t*t*base[0],
                    m*m*m*p0[1] + 3*m*m*t*c1[1] + 3*m*t*t*c2[1] + t*t*t*base[1]))
    return out

# ----------------------------------------------------------------------------- text
_font_cache = {}
def _font(sz):
    if sz not in _font_cache:
        _font_cache[sz] = ImageFont.truetype(FONT_PATH, sz)
    return _font_cache[sz]

def _text(x, y, s, sz, anchor='start'):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial" font-size="{sz}" '
            f'text-anchor="{anchor}">{s}</text>')

def _rxn_arrow_svg(x0, x1, y, reagent, conditions, rsz=26):
    out = [f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x1:.1f}" y2="{y:.1f}" stroke="black" stroke-width="3"/>',
           f'<path d="M {x1:.1f} {y:.1f} L {x1-14:.1f} {y-6:.1f} L {x1-14:.1f} {y+6:.1f} Z" fill="black"/>']
    cx = (x0+x1)/2.0
    if reagent:
        out.append(_text(cx, y-12, reagent, rsz, 'middle'))
    if conditions:
        out.append(_text(cx, y+rsz+4, conditions, rsz, 'middle'))
    return "".join(out)

# ----------------------------------------------------------------------------- frame
def render_frame_svg(lhs, rhs=(), arrows=(), charges=(), lone_pairs=(),
                     reagent="", conditions="", pad=70, sp_gap=150, arrow_gap=170):
    """Compose one arrow-pushing step into an SVG string. Returns (svg, W, H)."""
    lhs = list(lhs)
    for si, s in enumerate(lhs):                            # resolve orient=('auto',i,j): pick the clearest alkene face
        o = s.get('orient')
        if o and o[0] == 'auto':
            ar = next((a for a in arrows if a.get('route') and a['from'][0] == 'L'
                       and a['from'][1] == si and a['from'][2][0] == 'bond'), None)
            tgt = lhs[ar['to'][1]] if (ar and ar['to'][0] == 'L') else None
            if tgt:
                orient, vflip = _auto_face(s, (o[1], o[2]), tgt, sp_gap)
                s = dict(s); s['orient'] = orient; s['vflip'] = vflip; lhs[si] = s
    L = [_lay(s) for s in lhs]
    R = [_lay(s) for s in rhs]
    axis = max([s['axis'] for s in L + R] + [0.0])
    rowtop = pad
    plus_w = _txt_w("+", _font(40)) + 30
    Hmax = max([s['H'] for s in L + R] + [1])

    def place(side, x):
        blocks = []
        for k, s in enumerate(side):
            if k:
                x += plus_w
            oy = rowtop + (axis - s['axis'])
            blocks.append((s, x, oy)); x += s['W'] + sp_gap
        return blocks, (x - sp_gap)

    Lb, lx = place(L, pad)
    has_arrow = bool(rhs or reagent or conditions)
    ax0 = lx + 30
    ax1 = ax0 + (arrow_gap if has_arrow else 0)
    Rb, rxend = place(R, ax1 + 30) if R else ([], ax1)
    Wt = int(max(rxend, ax1) + pad)
    Ht = int(rowtop + Hmax + pad + 46)

    body = [f'<rect width="{Wt}" height="{Ht}" fill="white"/>']
    blockmap = {'L': Lb, 'R': Rb}
    for blocks in (Lb, Rb):
        for s, ox, oy in blocks:
            body.append(f'<g transform="translate({ox:.1f},{oy:.1f})">{s["inner"]}</g>')
            if s['label']:
                body.append(_text(ox + s['W']/2, oy + s['H'] + 26, s['label'], 24, 'middle'))

    def draw_plus(blocks):
        for k in range(1, len(blocks)):
            s0, ox0, _ = blocks[k-1]; s1, ox1, _ = blocks[k]
            body.append(_text((ox0 + s0['W'] + ox1)/2.0, rowtop + Hmax/2 + 12, "+", 40, 'middle'))
    draw_plus(Lb); draw_plus(Rb)
    if has_arrow:
        body.append(_rxn_arrow_svg(ax0, ax1, rowtop + Hmax/2, reagent, conditions))

    def apx(side, idx, i):
        s, ox, oy = blockmap[side][idx]
        x, y = s['px'](s['a'][i][1], s['a'][i][2]); return (ox + x, oy + y)

    def gpx(ref):
        side, idx, loc = ref
        if loc[0] == 'atom':
            return apx(side, idx, loc[1])
        if loc[0] == 'atomoff':                            # atom centre + (dx,dy) px offset (e.g. aim below a leaving atom)
            ax, ay = apx(side, idx, loc[1]); return (ax + loc[2], ay + loc[3])
        if loc[0] == 'bond':
            p, q = apx(side, idx, loc[1]), apx(side, idx, loc[2]); return ((p[0]+q[0])/2.0, (p[1]+q[1])/2.0)
        if loc[0] == 'pt':
            s, ox, oy = blockmap[side][idx]; return (ox + loc[1], oy + loc[2])
        raise ValueError(loc)

    # obstacle geometry (frame px): atom points (with a clearance radius) + bond segments, per species
    def _obstacles(exclude_bond, target_ref):
        atom_pts, segs = [], []
        tside, tidx, tloc = target_ref
        tatom = tloc[1] if tloc[0] in ('atom', 'atomoff') else None
        for side in ('L', 'R'):
            for idx, (s, ox, oy) in enumerate(blockmap[side]):
                loc = exclude_bond[2] if (exclude_bond and exclude_bond[0] == side and exclude_bond[1] == idx) else None
                exb = loc if (loc and loc[0] == 'bond') else None
                for i in s['a']:
                    if side == tside and idx == tidx and i == tatom:
                        continue                            # don't fence off the atom we're aiming at
                    lab = s['a'][i][0]
                    rad = 15.0 if (not lab or lab == 'H') else 22.0   # labels claim more room than bare vertices/H
                    atom_pts.append((*apx(side, idx, i), rad))
                for (bi, bj, o) in s['b']:
                    if exb and {bi, bj} == {exb[1], exb[2]}:
                        continue                            # the bond the arrow springs from is not an obstacle
                    segs.append((apx(side, idx, bi), apx(side, idx, bj)))
        return atom_pts, segs

    def _min_clear(pts, atom_pts, segs):
        m = 1e9
        for (x, y) in pts:
            for (ax, ay, rad) in atom_pts:
                d = math.hypot(x-ax, y-ay) - rad
                if d < m: m = d
            for (a, b) in segs:
                d = _seg_dist((x, y), a, b)
                if d < m: m = d
        return m

    def route(arr, p0, p1):
        """Search (sideways start-offset, dip depth) for the LEAST-shifted arc that clears everything
        by a comfortable margin. Deepening alone can't dodge a substituent in the descent corridor;
        a small sideways shift can."""
        sgn = 1.0 if arr.get('bow', 1) >= 0 else -1.0
        lift, stop = arr.get('lift', 0.3), arr.get('stop', 22)
        atom_pts, segs = _obstacles(arr['from'], arr['to'])
        margin = arr.get('clear', 12.0)
        step = 0.13 * HEAVY_LEN * U
        best = None                                         # (min_clear, p0x, depth) fallback if nothing hits margin
        for dx in (0, -15, 15, -30, 30, -46, 46):
            px0 = p0[0] + dx
            depth = 0.5 * HEAVY_LEN * U
            for _ in range(22):
                md = _min_clear(_shaft_samples((px0, p0[1]), p1, sgn*depth, lift, stop), atom_pts, segs)
                if md >= margin:
                    return (px0, p0[1]), sgn*depth          # offsets ordered by preference -> first comfortable clear wins
                if best is None or md > best[0]:
                    best = (md, px0, sgn*depth)
                depth += step
        return (best[1], p0[1]), best[2]

    for arr in arrows:
        p0 = gpx(arr['from'])
        p0 = (p0[0] + arr.get('startdx', 0.0), p0[1] + arr.get('startperp', 14))   # offset the tail to the lone pair / off the source bond
        p1 = gpx(arr['to'])
        if arr.get('route'):
            p0, depth = route(arr, p0, p1)
        else:
            depth = arr.get('depth', 40)
        body.append(svg_arrow(p0, p1, depth=depth, lift=arr.get('lift', 0.3),
                              stop=arr.get('stop', 22), kind=arr.get('kind', 'full')))
    for ch in charges:
        x, y = gpx((ch['at'][0], ch['at'][1], ('atom', ch['at'][2])))
        body.append(_text(x + ch.get('dx', -8), y + ch.get('dy', -34), ch['text'], 26))
    for lp in lone_pairs:
        x, y = gpx((lp['at'][0], lp['at'][1], ('atom', lp['at'][2])))
        body.append(svg_lone_pair(x, y, where=lp.get('where', 'up'), dist=lp.get('dist', 26)))

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{Wt}" height="{Ht}" '
           f'viewBox="0 0 {Wt} {Ht}">' + "".join(body) + '</svg>')
    return svg, Wt, Ht

def render_frame_png(scale=3, **kw):
    """Render a frame and rasterise to a PIL image (via PyMuPDF)."""
    import fitz
    svg, W, H = render_frame_svg(**kw)
    f = os.path.join(tempfile.gettempdir(), "_mech.svg"); open(f, 'w').write(svg)
    doc = fitz.open(f); pdf = fitz.open("pdf", doc.convert_to_pdf())
    from PIL import Image
    import io
    pix = pdf[0].get_pixmap(matrix=fitz.Matrix(scale, scale))
    return Image.open(io.BytesIO(pix.tobytes("png")))
