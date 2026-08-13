"""SVG (vector) output for structure_draw: sharp at any zoom, for embedding in Word."""
import math, collections
import structure_draw as J
from structure_draw import (_parse_label, _txt_w, _suffix_w, _label_parts, _label_runs,
                      _suffix_for_side, _anchor_w, _inline_w, _inline_runs,
                      U, STROKE, DBL, TRP, _BASE_H, FONT, SUBFONT, FONT_NAME,
                      DATIVE_HEAD_LEN, DATIVE_HEAD_HALF)

FS = 34        # main glyph px (Arial)
SUBFS = 23     # subscript px
BASE = FS * 0.35
SUB_DROP = 9
SUP_RISE = 11  # superscript (charge) rise above the baseline

def _esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def _inline_svg(text):
    return ''.join(
        f'<tspan font-style="italic">{_esc(run)}</tspan>' if is_italic else _esc(run)
        for run, is_italic in _inline_runs(text)
    )


def _geom_and_elems(atoms, bonds, circles, reversible, font):
    circles = circles or []; reversible = reversible or set()
    nbx = collections.defaultdict(list)
    for i, j, o in bonds:
        nbx[i].append(atoms[j][1] - atoms[i][1]); nbx[j].append(atoms[i][1] - atoms[j][1])
    parsed, side, ext = {}, {}, {}
    for i, (label, x, y) in atoms.items():
        if not label:
            continue
        a, s = _parse_label(label)
        has_r = any(dx > 0.3 for dx in nbx[i]); has_l = any(dx < -0.3 for dx in nbx[i])
        sd = 'left' if (i in reversible and has_r and not has_l) else 'right'
        s = _suffix_for_side(a, s, sd)
        parsed[i] = (a, s)
        side[i] = sd
        aw = _anchor_w(a); sw = _suffix_w(s)
        ext[i] = (aw / 2 + (sw if sd == 'left' else 0),
                  aw / 2 + (sw if sd == 'right' else 0), _BASE_H / 2)

    def exu(i, k): return (ext[i][k] / U + 0.05) if i in ext else 0
    minx = min(x - exu(i, 0) for i, (l, x, y) in atoms.items())
    maxx = max(x + exu(i, 1) for i, (l, x, y) in atoms.items())
    miny = min(y - exu(i, 2) for i, (l, x, y) in atoms.items())
    maxy = max(y + exu(i, 2) for i, (l, x, y) in atoms.items())
    M = 40
    W = int((maxx - minx) * U + 2 * M); H = int((maxy - miny) * U + 2 * M)

    def px(i):
        _, x, y = atoms[i]; return ((x - minx) * U + M, (maxy - y) * U + M)
    def pxc(x, y): return ((x - minx) * U + M, (maxy - y) * U + M)
    def gap(i, ux, uy):
        if i not in ext: return 0
        hx = ext[i][1] if ux > 0 else ext[i][0]
        return abs(ux) * hx + abs(uy) * ext[i][2] + 4

    out = []
    for i, j, order in bonds:
        p = px(i); q = px(j); dx = q[0] - p[0]; dy = q[1] - p[1]; L = math.hypot(dx, dy) or 1
        ux, uy = dx / L, dy / L
        p = (p[0] + ux * gap(i, ux, uy), p[1] + uy * gap(i, ux, uy))
        q = (q[0] - ux * gap(j, -ux, -uy), q[1] - uy * gap(j, -ux, -uy))
        ox, oy = -uy, ux
        if order == 'dative':  # coordinate bond: arrow with the filled head at the acceptor
            bx, by = q[0] - ux * DATIVE_HEAD_LEN, q[1] - uy * DATIVE_HEAD_LEN
            out.append(f'<line x1="{p[0]:.1f}" y1="{p[1]:.1f}" x2="{bx:.1f}" y2="{by:.1f}" '
                       f'stroke="black" stroke-width="{STROKE}" stroke-linecap="round"/>')
            out.append(f'<path d="M {q[0]:.1f} {q[1]:.1f} '
                       f'L {bx + ox * DATIVE_HEAD_HALF:.1f} {by + oy * DATIVE_HEAD_HALF:.1f} '
                       f'L {bx - ox * DATIVE_HEAD_HALF:.1f} {by - oy * DATIVE_HEAD_HALF:.1f} Z" '
                       'fill="black"/>')
            continue
        for s in {1: [0], 2: [-DBL, DBL], 3: [-TRP, 0, TRP]}[order]:
            out.append(f'<line x1="{p[0]+ox*s:.1f}" y1="{p[1]+oy*s:.1f}" '
                       f'x2="{q[0]+ox*s:.1f}" y2="{q[1]+oy*s:.1f}" '
                       f'stroke="black" stroke-width="{STROKE}" stroke-linecap="round"/>')
    for cx, cy, r in circles:
        c = pxc(cx, cy)
        out.append(f'<circle cx="{c[0]:.1f}" cy="{c[1]:.1f}" r="{r*U:.1f}" '
                   f'fill="none" stroke="black" stroke-width="{STROKE}"/>')
    for i in parsed:
        cx, cy = px(i); anchor, suffix = parsed[i]; sd = side[i]
        by = cy + BASE
        aw = _anchor_w(anchor)
        if anchor == "Cl":
            x0 = cx - aw / 2
            out.append(f'<text x="{x0:.1f}" y="{by:.1f}" font-family="{font}" '
                       f'font-size="{FS}" text-anchor="start">{_inline_svg(anchor)}</text>')
        else:
            out.append(f'<text x="{cx:.1f}" y="{by:.1f}" font-family="{font}" font-size="{FS}" '
                       f'text-anchor="middle">{_esc(anchor)}</text>')
        if suffix:
            x = cx + aw / 2 if sd == 'right' else cx - aw / 2 - _suffix_w(suffix)
            for s, kind in _label_runs(suffix):
                if kind == 'sub':
                    out.append(f'<text x="{x:.1f}" y="{by+SUB_DROP:.1f}" font-family="{font}" '
                               f'font-size="{SUBFS}" text-anchor="start">{_esc(s)}</text>')
                    x += _txt_w(s, SUBFONT)
                elif kind == 'sup':
                    out.append(f'<text x="{x:.1f}" y="{by-SUP_RISE:.1f}" font-family="{font}" '
                               f'font-size="{SUBFS}" text-anchor="start">{_esc(s)}</text>')
                    x += _txt_w(s, SUBFONT)
                else:
                    out.append(f'<text x="{x:.1f}" y="{by:.1f}" font-family="{font}" '
                               f'font-size="{FS}" text-anchor="start">{_inline_svg(s)}</text>')
                    x += _inline_w(s, FONT)
    return out, W, H, pxc

def draw_svg(atoms, bonds, circles=None, reversible=None, font=FONT_NAME):
    elems, W, H, _ = _geom_and_elems(atoms, bonds, circles, reversible, font)
    body = f'<rect width="{W}" height="{H}" fill="white"/>' + ''.join(elems)
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
            f'viewBox="0 0 {W} {H}">{body}</svg>'), W, H

def draw_svg_inner(atoms, bonds, circles=None, reversible=None, font=FONT_NAME):
    """Return (inner_svg_elements, W, H, px_func) for composing into a larger SVG."""
    elems, W, H, px = _geom_and_elems(atoms, bonds, circles, reversible, font)
    return ''.join(elems), W, H, px

def _render_svg_mol(mol, font=FONT_NAME, form='structural', angles='comb'):
    """Render an already validated RDKit molecule without parsing it again."""
    a, b, c, rev = J._layout_mol(mol, form, angles)
    return draw_svg(a, b, c, rev, font)


def _render_svg_inner_mol(mol, font=FONT_NAME, form='structural', angles='comb'):
    """Render inner SVG elements for an already validated RDKit molecule."""
    a, b, c, rev = J._layout_mol(mol, form, angles)
    return draw_svg_inner(a, b, c, rev, font)


def render_svg(smiles, font=FONT_NAME, form='structural', angles='comb', dative=False):
    """Render one SMILES input as an SVG string.

    dative=True draws a coordinate bond written in SMILES dative notation as a
    donor-to-acceptor arrow; the default refuses such input rather than drawing
    a plain line silently."""
    mol = J.validate_input(smiles, form, dative)
    if form == 'stereo':
        raise J.UnsupportedStereochemistryError(
            "Stereo wedge/dash rendering is not supported by the SVG backend; "
            "use structure_draw.render_stereo_pair for the supported raster pair."
        )
    return _render_svg_mol(mol, font, form, angles)


def render_svg_inner(smiles, font=FONT_NAME, form='structural', angles='comb', dative=False):
    """Render inner SVG elements for one SMILES input; dative as in render_svg."""
    mol = J.validate_input(smiles, form, dative)
    if form == 'stereo':
        raise J.UnsupportedStereochemistryError(
            "Stereo wedge/dash rendering is not supported by the SVG backend; "
            "use structure_draw.render_stereo_pair for the supported raster pair."
        )
    return _render_svg_inner_mol(mol, font, form, angles)
