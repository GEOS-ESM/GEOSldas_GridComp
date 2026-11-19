#!/usr/bin/env python3
"""
Compare two GEOSldas experiment directories created by ldas_setup
(e.g., one from the reference code and one from refactored code).

usage: compare_exp.py  REF_EXPDIR  NEW_EXPDIR  [--same REF_PATH NEW_PATH ...]

Compares:
  - directory tree (relative paths of all files, dirs, symlinks)
  - symlink targets (with REF/NEW expdir prefixes normalized)
  - contents of text files in run/ (rc, nml, lenkf.j, ldas_batchrun.j, cap_restart)
    and output/*/rc_out/BCS_info.txt, with expdir prefixes normalized
  - byte equality of other regular files (restarts, BCs copied into output/)

"--same REF_PATH NEW_PATH" declares additional path prefixes that are
equivalent (e.g., the two installation dirs used by the regression harness).

Exits 0 if identical, 1 otherwise.
"""

import os
import sys
import difflib
import filecmp

TEXT_SUFFIXES = ('.rc', '.nml', '.j', '.txt', 'cap_restart')


def walk(root):
    entries = {}
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        for name in dirnames + filenames:
            full = os.path.join(dirpath, name)
            rel  = os.path.relpath(full, root)
            if os.path.islink(full):
                entries[rel] = ('link', os.readlink(full))
            elif os.path.isdir(full):
                entries[rel] = ('dir', None)
            else:
                entries[rel] = ('file', None)
    return entries


def make_norm(expdir, extra):
    pairs = [(p, '<SAME%d>' % i) for i, p in enumerate(extra)]
    pairs += [(expdir, '<EXPDIR>'), (os.path.dirname(expdir), '<EXPHOME>')]
    # longest prefixes first
    pairs.sort(key=lambda x: -len(x[0]))
    def norm(text):
        for old, new in pairs:
            text = text.replace(old, new)
        return text
    return norm


def is_text(rel):
    return rel.startswith('run' + os.sep) and rel.endswith(TEXT_SUFFIXES) \
        or rel.endswith('BCS_info.txt')


def compare(ref, new, same=(), out=sys.stdout):
    ref, new = os.path.abspath(ref), os.path.abspath(new)
    norm_a = make_norm(ref, [os.path.abspath(p) for p, _ in same])
    norm_b = make_norm(new, [os.path.abspath(q) for _, q in same])
    a, b = walk(ref), walk(new)
    ndiff = 0

    def p(msg):
        out.write(msg + '\n')

    for rel in sorted(set(a) - set(b)):
        p('ONLY IN REF: ' + rel); ndiff += 1
    for rel in sorted(set(b) - set(a)):
        p('ONLY IN NEW: ' + rel); ndiff += 1

    for rel in sorted(set(a) & set(b)):
        ka, va = a[rel]
        kb, vb = b[rel]
        if ka != kb:
            p('TYPE DIFFERS: %s (%s vs %s)' % (rel, ka, kb)); ndiff += 1
            continue
        if ka == 'link':
            if norm_a(va) != norm_b(vb):
                p('LINK DIFFERS: %s\n   ref -> %s\n   new -> %s' % (rel, va, vb)); ndiff += 1
        elif ka == 'file':
            fa, fb = os.path.join(ref, rel), os.path.join(new, rel)
            if is_text(rel):
                with open(fa, errors='replace') as f: ta = norm_a(f.read()).splitlines()
                with open(fb, errors='replace') as f: tb = norm_b(f.read()).splitlines()
                if ta != tb:
                    ndiff += 1
                    p('TEXT DIFFERS: ' + rel)
                    for line in difflib.unified_diff(ta, tb, 'ref/' + rel, 'new/' + rel, lineterm=''):
                        p('   ' + line)
            elif not filecmp.cmp(fa, fb, shallow=False):
                p('BINARY DIFFERS: ' + rel); ndiff += 1
    return ndiff


def main(argv):
    args = list(argv)
    same = []
    while '--same' in args:
        i = args.index('--same')
        same.append((args[i+1], args[i+2]))
        del args[i:i+3]
    if len(args) != 2:
        sys.exit(__doc__)
    ndiff = compare(args[0], args[1], same)
    print('\n%d difference(s)' % ndiff)
    return 0 if ndiff == 0 else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
