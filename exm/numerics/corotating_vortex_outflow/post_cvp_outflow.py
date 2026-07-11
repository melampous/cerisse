#!/usr/bin/env python3
"""
CVP OUTFLOW variant post-processing: corruption metrics with radiating lateral BC.
Compare foextrap vs characteristic lateral BCs.
Output: comparison figure showing whether open-boundary BCs reveal true reflection differences.
"""
import sys, os, glob, json
import numpy as np
try:
    import yt
    yt.set_log_level(50)
except:
    print("yt not available; skipping data processing")
    sys.exit(1)

case_dir = "/shared/cerisse/exm/numerics/corotating_vortex_outflow" if len(sys.argv)<2 else sys.argv[1]

print("CVP OUTFLOW Post-Processing")
print("=" * 60)
print(f"Case dir: {case_dir}")

# Load latest plotfiles for each variant
for variant in ["outflow_foe", "outflow_char"]:
    pdir = f"{case_dir}/plot_{variant}"
    plts = sorted([p for p in glob.glob(f"{pdir}/plt*") if os.path.isfile(f"{p}/Header")])
    if plts:
        print(f"  {variant}: {len(plts)} frames, latest t={float(yt.load(plts[-1]).current_time):.2f}")
    else:
        print(f"  {variant}: NO FRAMES")

print("")
print("Ready for corruption metric extraction (once yt sampling is available).")
print("Expected output: outflow_vs_quiescent_comparison.png")
