"""Mechanism annotation layer: curly arrows (electron flow), lone pairs, partial charges/bonds.

Sits ON TOP of the per-molecule engine (structure_draw / structure_svg): a species is drawn
by that engine, then electron-flow annotations are overlaid, addressed by atom / bond index in
the SAME pixel space. A-level (9729) house style, arrow-pushing convention.

This module is being built incrementally. Right now: the curly-arrow PRIMITIVE only.
  * full  arrow  (double-barb head)  -> a PAIR of electrons  (heterolysis, most mechanisms)
  * half  arrow  (single-barb 'fish-hook') -> ONE electron    (homolysis, free-radical)
"""
import math
from PIL import Image, ImageDraw, ImageFont
import structure_draw as J
from structure_draw import STROKE, FONT_PATH

ARROW_W = STROKE                 # arrow shaft weight, matches bond stroke
HEAD_LEN = 14                    # arrowhead barb length (px)
HEAD_ANG = math.radians(24)      # half-angle of the arrowhead


def _cubic(p0, c1, c2, p1, n=60):
    out = []
    for k in range(n + 1):
        t = k / n; mt = 1 - t
        out.append((mt*mt*mt*p0[0] + 3*mt*mt*t*c1[0] + 3*mt*t*t*c2[0] + t*t*t*p1[0],
                    mt*mt*mt*p0[1] + 3*mt*mt*t*c1[1] + 3*mt*t*t*c2[1] + t*t*t*p1[1]))
    return out


def _rot(vx, vy, ang):
    ca, sa = math.cos(ang), math.sin(ang)
    return (vx*ca - vy*sa, vx*sa + vy*ca)


def draw_curly(dr, p0, p1, bow=0.4, kind='full', width=ARROW_W, color=(0, 0, 0), lift=0.28, depth=None):
    """Curved electron-flow arrow from p0 to p1, drawn as a U that LEAVES and ARRIVES perpendicular
    to the p0->p1 chord (a real curly arrow) so it clears the bonds it springs from / points at.
    bow: signed depth of the arc as a fraction of |p0 p1| (+ bows LEFT of p0->p1, - bows right).
    depth: absolute arc depth in px (overrides bow*L); use for SHORT arrows (a bond-fission hook)
           where a relative depth would collapse to nothing. Sign follows `bow`.
    lift: how far the control points slide along the chord (rounds the U so it isn't a hairpin).
    kind: 'full' = 2-electron (double barb); 'half' = 1-electron fish-hook (single barb)."""
    dx, dy = p1[0]-p0[0], p1[1]-p0[1]
    L = math.hypot(dx, dy) or 1.0
    ux, uy = dx/L, dy/L                                   # unit chord
    nx, ny = -uy, ux                                      # unit perpendicular (left of p0->p1)
    d = (math.copysign(depth, bow) if depth is not None else bow * L)
    c1 = (p0[0] + nx*d + ux*L*lift, p0[1] + ny*d + uy*L*lift)   # leave p0 perpendicular (toward the bow side)
    c2 = (p1[0] + nx*d - ux*L*lift, p1[1] + ny*d - uy*L*lift)   # arrive p1 perpendicular (from the bow side)
    pts = _cubic(p0, c1, c2, p1)
    dr.line(pts, fill=color, width=width, joint='curve')
    tx, ty = p1[0]-c2[0], p1[1]-c2[1]; tl = math.hypot(tx, ty) or 1.0   # arrival tangent
    tx, ty = tx/tl, ty/tl
    back = (-tx, -ty)
    b1 = _rot(back[0], back[1],  HEAD_ANG)
    b2 = _rot(back[0], back[1], -HEAD_ANG)
    dr.line([p1, (p1[0]+b1[0]*HEAD_LEN, p1[1]+b1[1]*HEAD_LEN)], fill=color, width=width)
    if kind == 'full':                                   # 2-electron: the second barb
        dr.line([p1, (p1[0]+b2[0]*HEAD_LEN, p1[1]+b2[1]*HEAD_LEN)], fill=color, width=width)


def _resolve(spec, AP, BMID, pad):
    """spec: ('atom', i) | ('bond', i, j) | ('pt', px, py)  ->  (x, y) in the padded canvas."""
    t = spec[0]
    if t == 'atom':
        return AP(spec[1])
    if t == 'bond':
        return BMID(spec[1], spec[2])
    if t == 'pt':
        return (spec[1] + pad, spec[2] + pad)
    raise ValueError(f"bad locator {spec!r}")


def _orient(a, i, j):
    """Rigid-rotate the atoms dict so bond i->j lies horizontal (i left, j right). Angles are
    preserved, so a trigonal alkene stays trigonal but its C=C reads horizontally, as in the notes."""
    xi, yi = a[i][1], a[i][2]
    ang = math.atan2(a[j][2] - yi, a[j][1] - xi)
    ca, sa = math.cos(-ang), math.sin(-ang)
    out = {}
    for k, (lab, x, y) in a.items():
        dx, dy = x - xi, y - yi
        out[k] = (lab, xi + dx*ca - dy*sa, yi + dx*sa + dy*ca)
    return out


def render_species(smiles, arrows=(), form='displayed', pad=80, orient=None):
    """Draw one species and overlay curly arrows addressed by atom/bond index.
    arrows: list of {'from': locator, 'to': locator, 'bow': float, 'kind': 'full'|'half'}.
    orient: ('bond', i, j) to rotate that bond horizontal before drawing."""
    a, b, circ, rev = J.layout(smiles, form)
    if orient and orient[0] == 'bond':
        a = _orient(a, orient[1], orient[2])
    img, pxc, W, H = J.draw(a, b, circ, rev, return_map=True)
    big = Image.new("RGB", (W + 2*pad, H + 2*pad), "white"); big.paste(img, (pad, pad))
    dr = ImageDraw.Draw(big)
    def AP(i):
        x, y = pxc(a[i][1], a[i][2]); return (x + pad, y + pad)
    def BMID(i, j):
        p, q = AP(i), AP(j); return ((p[0]+q[0])/2.0, (p[1]+q[1])/2.0)
    for arr in arrows:
        p0 = _resolve(arr['from'], AP, BMID, pad)
        p1 = _resolve(arr['to'], AP, BMID, pad)
        draw_curly(dr, p0, p1, bow=arr.get('bow', 0.4), kind=arr.get('kind', 'full'),
                   lift=arr.get('lift', 0.28), depth=arr.get('depth'))
    return big


_DIR = {'up': (0, -1), 'down': (0, 1), 'left': (-1, 0), 'right': (1, 0),
        'upleft': (-0.7, -0.7), 'upright': (0.7, -0.7)}

def draw_lone_pair(dr, cx, cy, where='up', dist=24, gap=9, r=4.0, color=(0, 0, 0)):
    """Two dots (a lone pair) sitting off an atom centre, in direction `where`."""
    dx, dy = _DIR[where] if isinstance(where, str) else where
    dn = math.hypot(dx, dy) or 1.0; dx, dy = dx/dn, dy/dn
    bx, by = cx + dx*dist, cy + dy*dist                  # pair centre, clear of the label
    px, py = -dy, dx                                      # perpendicular offset between the two dots
    for s in (-1, 1):
        x, y = bx + px*gap*s, by + py*gap*s
        dr.ellipse([x-r, y-r, x+r, y+r], fill=color)


def _lay_species(sp):
    """Lay out one species -> its own PIL image + an atom->local-px mapper + the y of its backbone axis."""
    a, b, circ, rev = J.layout(sp['smiles'], sp.get('form', 'displayed'))
    o = sp.get('orient')
    if o and o[0] == 'bond':
        a = _orient(a, o[1], o[2])
    img, pxc, W, H = J.draw(a, b, circ, rev, return_map=True)
    return dict(img=img, pxc=pxc, a=a, W=W, H=H, axis=pxc(0.0, 0.0)[1], label=sp.get('label'))


def _rxn_arrow(dr, x0, x1, y, reagent="", conditions="", font=None):
    dr.line([(x0, y), (x1, y)], fill=(0, 0, 0), width=STROKE)
    dr.polygon([(x1, y), (x1-13, y-6), (x1-13, y+6)], fill=(0, 0, 0))
    if font:
        if reagent:
            w = dr.textlength(reagent, font=font); dr.text(((x0+x1)/2 - w/2, y-30), reagent, fill=(0,0,0), font=font)
        if conditions:
            w = dr.textlength(conditions, font=font); dr.text(((x0+x1)/2 - w/2, y+8), conditions, fill=(0,0,0), font=font)


# ------- the frame: [ reactants, with curly arrows ]  -->  [ products ] --------------------
def render_frame(lhs, rhs=(), arrows=(), charges=(), lone_pairs=(), reagent="", conditions="",
                 pad=70, sp_gap=54, arrow_gap=150, plus=True):
    """One arrow-pushing step. lhs/rhs: list of species dicts {smiles, form, orient, label}.
    arrows/charges reference atoms by (side, sp_idx, locator); side in {'L','R'}.
      arrow:  {'from': (side, i, loc), 'to': (side, i, loc), 'bow', 'kind'}
      charge: {'at': (side, i, atom), 'text': 'δ+', 'dx':0, 'dy':-26}
    loc: ('atom', i) | ('bond', i, j) | ('pt', dx, dy) relative to that species' top-left."""
    fbig = ImageFont.truetype(FONT_PATH, 30); fsm = ImageFont.truetype(FONT_PATH, 24); frg = ImageFont.truetype(FONT_PATH, 22)
    L = [_lay_species(s) for s in lhs]
    R = [_lay_species(s) for s in rhs]
    axis = max([s['axis'] for s in L + R] + [0])            # common backbone axis (max so nothing clips above)
    rowtop = pad
    # measure the "+" width on a throwaway canvas
    dr0 = ImageDraw.Draw(Image.new("RGB", (4, 4)))
    plus_w = dr0.textlength("+", font=fbig) + 24 if plus else 0
    def place2(side, x):
        blocks = []
        for k, s in enumerate(side):
            if k:
                x += plus_w
            oy = rowtop + (axis - s['axis'])
            blocks.append((s, x, oy)); x += s['W'] + sp_gap
        return blocks, (x - sp_gap)
    Lblocks, lx = place2(L, pad)
    ax0 = lx + 24
    ax1 = ax0 + (arrow_gap if (rhs or reagent or conditions) else 0)
    rx = ax1 + 24 if (rhs or reagent or conditions) else lx + sp_gap
    Rblocks, rxend = place2(R, rx) if R else ([], rx)
    Wt = int(max(rxend, ax1) + pad)
    Ht = int(rowtop + max([s['H'] for s in L + R] + [1]) + pad + 40)
    canvas = Image.new("RGB", (Wt, Ht), "white"); dr = ImageDraw.Draw(canvas)
    blockmap = {'L': Lblocks, 'R': Rblocks}
    for side_blocks in (Lblocks, Rblocks):
        for s, ox, oy in side_blocks:
            canvas.paste(s['img'], (int(ox), int(oy)))
            if s['label']:
                w = dr.textlength(s['label'], font=fsm)
                dr.text((ox + s['W']/2 - w/2, oy + s['H'] + 4), s['label'], fill=(0,0,0), font=fsm)
    # plus signs
    def draw_plus(blocks):
        for k in range(1, len(blocks)):
            s0, ox0, _ = blocks[k-1]; s1, ox1, _ = blocks[k]
            mx = (ox0 + s0['W'] + ox1) / 2.0
            dr.text((mx - dr.textlength("+", font=fbig)/2, rowtop + max(s['H'] for s in L+R)/2 - 18), "+", fill=(0,0,0), font=fbig)
    draw_plus(Lblocks); draw_plus(Rblocks)
    if rhs or reagent or conditions:
        _rxn_arrow(dr, ax0, ax1, rowtop + max(s['H'] for s in L+R)/2, reagent, conditions, frg)

    def gpx(ref):                                          # (side, sp_idx, locator) -> composed-canvas px
        side, idx, loc = ref
        s, ox, oy = blockmap[side][idx]
        def apx(i):
            x, y = s['pxc'](s['a'][i][1], s['a'][i][2]); return (ox + x, oy + y)
        if loc[0] == 'atom': return apx(loc[1])
        if loc[0] == 'bond':
            p, q = apx(loc[1]), apx(loc[2]); return ((p[0]+q[0])/2.0, (p[1]+q[1])/2.0)
        if loc[0] == 'pt':   return (ox + loc[1], oy + loc[2])
        raise ValueError(loc)

    def _clear_depth(arr, p0, p1):
        """Router: dip depth (px) so the arc clears the OUTERMOST substituent of the source compound on
        the bow side, + a bond-length margin. Uses kekule's exact atom coords, so it never guesses."""
        side, idx, loc = arr['from']
        s = blockmap[side][idx][0]
        dx, dy = p1[0]-p0[0], p1[1]-p0[1]; Lc = math.hypot(dx, dy) or 1.0
        ux, uy = dx/Lc, dy/Lc
        sgn = 1.0 if arr.get('bow', 1) >= 0 else -1.0
        nx, ny = -uy*sgn, ux*sgn                            # unit perpendicular on the bow side
        maxproj = 0.0
        for i in s['a']:
            ax, ay = gpx((side, idx, ('atom', i)))
            along = (ax-p0[0])*ux + (ay-p0[1])*uy
            proj = (ax-p0[0])*nx + (ay-p0[1])*ny            # how far this atom sticks out on the bow side
            if -0.25*Lc < along < 1.15*Lc and proj > maxproj:
                maxproj = proj
        return maxproj + arr.get('margin', 0.85) * J.HEAVY_LEN * J.U

    for arr in arrows:
        p0, p1 = gpx(arr['from']), gpx(arr['to'])
        if arr['to'][2][0] == 'atom':                      # stop just OUTSIDE the target atom's label so the head points at it
            dx, dy = p1[0]-p0[0], p1[1]-p0[1]; d = math.hypot(dx, dy) or 1.0
            ins = arr.get('inset', 20)
            p1 = (p1[0] - dx/d*ins, p1[1] - dy/d*ins)
        depth = _clear_depth(arr, p0, p1) if arr.get('route') else arr.get('depth')
        draw_curly(dr, p0, p1, bow=arr.get('bow', 0.4), kind=arr.get('kind', 'full'), depth=depth,
                   lift=arr.get('lift', 0.28))
    for ch in charges:
        x, y = gpx((ch['at'][0], ch['at'][1], ('atom', ch['at'][2])))
        dr.text((x + ch.get('dx', 6), y + ch.get('dy', -30)), ch['text'], fill=(0,0,0), font=fsm)
    for lp in lone_pairs:
        x, y = gpx((lp['at'][0], lp['at'][1], ('atom', lp['at'][2])))
        draw_lone_pair(dr, x, y, where=lp.get('where', 'up'), dist=lp.get('dist', 24))
    return canvas


# --------------------------------------------------------------------------- dev test card
def _test_card(path):
    W, H = 900, 620
    img = Image.new("RGB", (W, H), "white"); dr = ImageDraw.Draw(img)
    f = ImageFont.truetype(FONT_PATH, 20)
    cases = [
        ("full, bow +0.4  (arcs left)",  (80, 120),  (320, 120),  0.4,  'full'),
        ("full, bow -0.4  (arcs right)", (420, 120), (660, 120), -0.4,  'full'),
        ("full, vertical",               (120, 300), (120, 500),  0.35, 'full'),
        ("full, diagonal, big bow",      (330, 500), (540, 300),  0.6,  'full'),
        ("half (fish-hook), bow +0.4",   (610, 480), (840, 320),  0.45, 'half'),
        ("full, shallow bow +0.15",      (620, 560), (860, 560),  0.15, 'full'),
    ]
    for label, p0, p1, bow, kind in cases:
        dr.ellipse([p0[0]-3, p0[1]-3, p0[0]+3, p0[1]+3], fill=(200, 0, 0))   # source dot (red)
        draw_curly(dr, p0, p1, bow=bow, kind=kind)
        dr.text((min(p0[0], p1[0]), min(p0[1], p1[1]) - 26), label, fill=(90, 90, 90), font=f)
    img.save(path)
    return img


if __name__ == "__main__":
    import sys
    _test_card(sys.argv[1] if len(sys.argv) > 1 else "/tmp/curly_test.png")
    print("wrote test card")
