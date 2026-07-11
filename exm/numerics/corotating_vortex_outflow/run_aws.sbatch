#!/bin/bash
#SBATCH -J cvp_outflow
#SBATCH -p cpu
#SBATCH -N 4
#SBATCH --ntasks-per-node=96
#SBATCH --time=04:00:00
#SBATCH -o run_%j.out
#SBATCH -e run_%j.err

set -euo pipefail

# =====================================================================
# CVP OUTFLOW VARIANT: narrow domain (±4), lateral BC = foextrap (code 2),
# y-direction characteristic (code 7). Run to t=55 (~3.5 T_rot), 85 frames.
#
# Labels:
#   outflow_foe  = foextrap lateral BC variant
#   outflow_char = characteristic y-faces (for reference if y-faces switch)
#
# Comparison targets:
#   - wide_char (baseline, WIDE domain code-7 all faces)
#   - narrow_char (NARROW domain code-7 all faces)
#   - narrow_lodi (NARROW domain LODI y-faces)
# =====================================================================

# ---- user knobs ----
CASE_DIR="/shared/cerisse/exm/numerics/corotating_vortex_outflow"
MAX_STEP=${MAX_STEP:-1000000}          # stop at t=55 (stop_time in inputs.outflow)
LABEL_FOE="outflow_foe"                # foextrap variant label
LABEL_CHAR="outflow_char"              # characteristic y-faces label

# ---- env ----
module purge
module load openmpi/4.1.7
cd "${CASE_DIR}"

echo "=========================================================="
echo "CVP OUTFLOW VARIANT (4N96): t=55, 85 frames"
echo "Host: $(hostname)   Date: $(date)"
echo "CASE_DIR=${CASE_DIR}"
echo "SLURM_NTASKS=${SLURM_NTASKS}   SLURM_NNODES=${SLURM_NNODES}"
echo "Git state:"
git -C /shared/cerisse log --oneline -1 || true
echo "=========================================================="

# ---- build on compute node ----
echo "Building on compute node at $(date)"
rm -rf tmp_build_dir main2d.gnu.TPROF.MPI.ex
make -j16 DIM=2
ls -la main2d.gnu.TPROF.MPI.ex

# ---- RUN: foextrap lateral BC (code 2) + characteristic y-faces (code 7) ----
echo ""
echo "=========================================================="
echo "RUN 1: ${LABEL_FOE} variant (start $(date))"
echo "=========================================================="
mkdir -p plot plot_${LABEL_FOE}
/opt/amazon/openmpi/bin/mpirun --bind-to core -np "${SLURM_NTASKS}" \
    ./main2d.gnu.TPROF.MPI.ex inputs.outflow \
    max_step="${MAX_STEP}" \
    amr.plot_file=./plot_${LABEL_FOE}/plt

echo "Finished ${LABEL_FOE} at $(date)"

# ---- RUN: Optional second variant if y-faces are later changed to characteristic only ----
# (currently this is the same as variant 1, but structure allows independent testing)
echo ""
echo "=========================================================="
echo "RUN 2: ${LABEL_CHAR} variant (characteristic y-faces) (start $(date))"
echo "=========================================================="
mkdir -p plot_${LABEL_CHAR}
/opt/amazon/openmpi/bin/mpirun --bind-to core -np "${SLURM_NTASKS}" \
    ./main2d.gnu.TPROF.MPI.ex inputs.outflow \
    max_step="${MAX_STEP}" \
    amr.plot_file=./plot_${LABEL_CHAR}/plt

echo "Finished ${LABEL_CHAR} at $(date)"

echo ""
echo "=========================================================="
echo "ALL RUNS COMPLETED at $(date)"
echo "Output directories:"
echo "  ./plot_${LABEL_FOE}/plt_*  (foextrap lateral)"
echo "  ./plot_${LABEL_CHAR}/plt_*  (characteristic y-faces)"
echo "=========================================================="
