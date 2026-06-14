#!/usr/bin/env python3
"""
OWL correctness unit tests – master runner.

Runs all six tests and produces a summary figure.

Usage:
    cd test/
    python run_all_tests.py

Each sub-test writes its own plot to test/plots/.  This script additionally
produces a single summary figure (test/plots/summary.png).
"""

import sys, os, time, traceback
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

HERE = os.path.dirname(os.path.abspath(__file__))
PLOT = os.path.join(HERE, 'plots')
os.makedirs(PLOT, exist_ok=True)

# All test scripts live alongside this runner in a single directory.
sys.path.insert(0, HERE)

# ── Import test modules ────────────────────────────────────────────────────────
import importlib

TEST_DEFS = [
    ('test_analytic',          'Test 1',  'Analytic wholespace\nresponse'),
    ('test_linearity_src',     'Test 2',  'Modeling linearity\nw.r.t. source'),
    ('test_adjoint_src',       'Test 3',  'Modeling adjoint\n(reciprocity)'),
    ('test_linearization',     'Test 4',  'Nonlinear operator\nlinearization'),
    ('test_jacobian_linearity','Test 5',  'Jacobian linearity\nw.r.t. model'),
    ('test_jacobian_adjoint',  'Test 6',  'Jacobian adjoint\ntest'),
    ('test_elastic_lamb',      'Test 7',  'Elastic Lamb problem\non a tilted surface'),
]

results = []   # list of (name, short, passed, elapsed, detail)

for mod_name, label, short in TEST_DEFS:
    print()
    print('=' * 70)
    print(f'  Running {label}: {short.replace(chr(10), " ")}')
    print('=' * 70)
    t0  = time.time()
    try:
        spec = importlib.util.spec_from_file_location(
            mod_name, os.path.join(HERE, mod_name + '.py'))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        ret    = mod.run()
        passed = ret[0]
        detail = ''
    except Exception:
        passed = False
        detail = traceback.format_exc()[-500:]
        print(detail)

    elapsed = time.time() - t0
    results.append((label, short, passed, elapsed, detail))
    status = 'PASSED' if passed else 'FAILED'
    print(f'\n  → {label} {status}  ({elapsed:.1f} s)')

# ── Summary table ──────────────────────────────────────────────────────────────
print()
print('=' * 70)
print('  SUMMARY')
print('=' * 70)
all_passed = True
for label, short, passed, elapsed, detail in results:
    status = '✓ PASSED' if passed else '✗ FAILED'
    print(f'  {label}: {status}  ({elapsed:.1f} s)')
    if not passed:
        all_passed = False
print('=' * 70)
print(f'  Overall: {"ALL PASSED" if all_passed else "SOME TESTS FAILED"}')
print()

# ── Summary figure ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 4))
ax.set_xlim(-0.5, len(results) - 0.5)
ax.set_ylim(-0.5, 1.5)
ax.axis('off')

for i, (label, short, passed, elapsed, _) in enumerate(results):
    color = '#2ecc71' if passed else '#e74c3c'
    rect  = mpatches.FancyBboxPatch(
        (i - 0.4, 0.2), 0.8, 0.8,
        boxstyle='round,pad=0.05',
        facecolor=color, edgecolor='white', lw=2
    )
    ax.add_patch(rect)
    ax.text(i, 0.95, label, ha='center', va='center', fontsize=10,
            fontweight='bold', color='white')
    ax.text(i, 0.60, short, ha='center', va='center', fontsize=7.5,
            color='white', multialignment='center')
    status_txt = 'PASSED' if passed else 'FAILED'
    ax.text(i, 0.30, f'{status_txt}\n({elapsed:.0f} s)', ha='center', va='center',
            fontsize=8, color='white', fontweight='bold')

title_color = '#27ae60' if all_passed else '#c0392b'
fig.suptitle('OWL Correctness Unit Tests – '
             + ('All Tests Passed ✓' if all_passed else 'Some Tests Failed ✗'),
             fontsize=13, fontweight='bold', color=title_color, y=1.0)
plt.tight_layout()
fig.savefig(os.path.join(PLOT, 'summary.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f'Summary figure saved to {PLOT}/summary.png')

sys.exit(0 if all_passed else 1)
