#!/usr/bin/env python3
"""Phy CE: add the physics letters to the math faces of src/ui_fontdata.c, in place, and lengthen
every face's index for their codes: omega, alpha, tau, phi italic (ui_mi*, STIX Two Text Italic),
Delta and the degree sign upright (ui_mu*, STIX Two Text). The math layout (mathinput.cc) draws
the names omega, alpha, tau, phi, Delta as these letters and the unit _deg as the degree sign.
Only the new glyphs are rasterized, with mkfonts.py's own functions and settings.
  python3 tools/fonts/add_math_greek.py
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mkfonts as M

P = os.path.join(M.ROOT, 'src', 'ui_fontdata.c')
NEW_IT = ['ω', 'α', 'τ', 'φ']   # italic (ui_mi*)
NEW_UP = ['Δ', '°']             # upright (ui_mu*)
ALL = NEW_IT + NEW_UP
EXTRA = M.EXTRA + [c for c in ALL if c not in M.EXTRA]
N = 0x80 + len(EXTRA)
code = lambda ch: 0x80 + EXTRA.index(ch)
s = open(P).read()
newindex = {}
for size in M.MATH_SIZES:
    for kind, chars, path in (('mi', NEW_IT, 'STIXTwoText-Italic-VF.ttf'), ('mu', NEW_UP, 'STIXTwoText-VF.ttf')):
        name = 'ui_%s%d' % (kind, size)
        bm = re.search(r'static const unsigned char %s_bits\[(\d+)\] = \{\n(.*?)\n\};' % name, s, re.S)
        gm = re.search(r'static const ui_glyph %s_glyphs\[(\d+)\] = \{\n(.*?)\n\};' % name, s, re.S)
        if not bm or not gm:
            print('no face', name); continue
        bits = [int(x) for x in re.findall(r'\d+', bm.group(2))]
        glyphs = [tuple(int(v) for v in g) for g in re.findall(r'\{(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+)\}', gm.group(2))]
        assert len(bits) == int(bm.group(1)) and len(glyphs) == int(gm.group(1)), name
        for ch in chars:
            c = code(ch)
            if any(g[0] == c for g in glyphs): continue
            f = M.font(path, size, 400)
            adv, ox, oy, w, h, pix = M.raster(f, ch)
            glyphs.append((c, len(bits), adv, w, h, ox, oy))
            bits += M.encode(pix)
        if chars: bits.append(0)  # (the blitter reads one byte past a glyph's last pixels)
        glyphs.sort()
        assert len(bits) < 65536
        bl = ['  ' + ','.join(str(b) for b in bits[i:i + 24]) + ',' for i in range(0, len(bits), 24)]
        gl = ['{%d,%d,%d,%d,%d,%d,%d}' % g for g in glyphs]
        gl = ['  ' + ','.join(gl[i:i + 6]) + ',' for i in range(0, len(gl), 6)]
        s = s[:bm.start()] + 'static const unsigned char %s_bits[%d] = {\n%s\n};' % (name, len(bits), '\n'.join(bl)) + s[bm.end():]
        gm = re.search(r'static const ui_glyph %s_glyphs\[(\d+)\] = \{\n(.*?)\n\};' % name, s, re.S)
        s = s[:gm.start()] + 'static const ui_glyph %s_glyphs[%d] = {\n%s\n};' % (name, len(glyphs), '\n'.join(gl)) + s[gm.end():]
        s = re.sub(r'const ui_face %s = \{(\d+),(\d+),(\d+),(\d+),' % name, lambda m: 'const ui_face %s = {%s,%s,%s,%d,' % (name, m.group(1), m.group(2), m.group(3), len(glyphs)), s)
        idx = [255] * N
        for i, g in enumerate(glyphs): idx[g[0]] = i
        assert kind not in newindex or newindex[kind] == idx, 'the %s faces differ' % kind
        newindex[kind] = idx
        print(name, 'glyphs', len(glyphs), 'bits', len(bits))
# every index N long: the math faces' from their glyphs (shared by size), the others 255 for the new codes
for im in list(re.finditer(r'static const unsigned char (ui_\w+_index)\[(\d+)\] = \{\n(.*?)\n\};', s, re.S))[::-1]:
    nm, vals = im.group(1), [int(x) for x in re.findall(r'\d+', im.group(3))]
    vals = newindex['mi'] if nm == 'ui_mi34_index' else newindex['mu'] if nm == 'ui_mu34_index' else (vals + [255] * N)[:N]
    il = ['  ' + ','.join(str(b) for b in vals[i:i + 24]) + ',' for i in range(0, len(vals), 24)]
    s = s[:im.start()] + 'static const unsigned char %s[%d] = {\n%s\n};' % (nm, len(vals), '\n'.join(il)) + s[im.end():]
    print(nm, len(vals))
open(P, 'w').write(s)
print('codes:', ', '.join('%s 0x%x' % (c, code(c)) for c in ALL))
