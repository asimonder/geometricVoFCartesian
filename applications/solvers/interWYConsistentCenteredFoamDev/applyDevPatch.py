#!/usr/bin/env python3
# Regenerates interWYConsistentCenteredFoamDev.C from the pristine copy
# and applies the single-reference correctPhi patch. Every replacement
# is asserted; any mismatch aborts before writing.
import re, sys, pathlib

d = pathlib.Path(__file__).resolve().parent
src = (d/'interWYConsistentCenteredFoamDev.C.pristine').read_text()

def replace_once(text, old, new, tag):
    n = text.count(old)
    if n != 1:
        sys.exit(f"PATCH ABORT: literal '{tag}' matched {n} times, need 1")
    return text.replace(old, new)

def replace_once_re(text, pattern, new, tag):
    ms = list(pattern.finditer(text))
    if len(ms) != 1:
        sys.exit(f"PATCH ABORT: regex '{tag}' matched {len(ms)} times, need 1")
    m = ms[0]
    return text[:m.start()] + new + text[m.end():]

# 1. includes for the inline projection
src = replace_once(
    src,
    '#include "CorrectPhi.H"',
    '#include "CorrectPhi.H"\n#include "adjustPhi.H"\n'
    '#include "fixedValueFvPatchFields.H"',
    'includes')

# 2. startup banner identifying the build
old_banner = 'Info<< "\\nStarting time loop\\n" << endl;'
new_banner = (
    'Info<< "interWYConsistentCenteredFoamDev: correctPhi rebuilt with"\n'
    '        << " the single global reference cell (pRefCell), conditional flux rebuild (v3); library"\n'
    '        << " CorrectPhi rank-local setReference(0,0) removed" << endl;\n\n'
    '    ' + old_banner)
src = replace_once(src, old_banner, new_banner, 'banner')

# 3. the CorrectPhi call block -> local correctPhi.H
block = re.compile(
    r'const surfaceScalarField rhofCorr\s*\(\s*'
    r'1\.0/fvc::interpolate\(1\.0/rho\)\s*\);\s*'
    r'CorrectPhi\s*\(\s*U, phi, p_rgh,\s*'
    r'surfaceScalarField\("rAUf", dtDim/rhofCorr\),\s*'
    r'geometricZeroField\(\),\s*pimple\s*\);',
    re.S)
src = replace_once_re(
    src, block,
    '// single-global-reference projection; see correctPhi.H\n'
    '                        #include "correctPhi.H"',
    'correctPhi block')

(d/'interWYConsistentCenteredFoamDev.C').write_text(src)
print('patch applied cleanly')
