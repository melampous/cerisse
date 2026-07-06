# Cerisse on 8× A100 — ready-to-run package (2026-07-05)

**Bottom line:** the code fully supports 8× A100 (AMReX CUDA + MPI, 39 GPU-ported src
files, CUDA-aware MPI). It has been strong-scaling benchmarked to 32 GPUs. The only
mandatory change vs the local build is `CUDA_ARCH = 8.0` (A100 = sm_80), already set in
the files below. Everything marked `__FILL__` is machine-specific (scheduler/account/modules).

---

## Files prepared (all in `exm/scaling/`)

| File | Purpose |
|------|---------|
| `exm/scaling/GNUmakefile.a100` | GPU build config, `USE_CUDA=TRUE`, `CUDA_ARCH=8.0` |
| `exm/scaling/input_bench256`   | 256³ periodic TGV, 500 steps, I/O off — matches the reference scaling table |
| `exm/scaling/run_a100.slurm`   | SLURM launcher, 8 ranks / 8 GPUs / 1 node |
| `exm/scaling/run_a100.pbs`     | PBS launcher (Imperial hx1 style), 8 GPUs |

Absolute path on this machine: `/home/qiaoj/testcerisse/cerisse/exm/scaling/`
AMReX: `/home/qiaoj/testcerisse/cerisse/lib/amrex` · solver src: `.../cerisse/src`

---

## STEP 0 — before you build, fill in 3 machine-specific things
Edit `run_a100.slurm` **or** `run_a100.pbs` (whichever your cluster uses):
1. `__ACCOUNT__` / project code
2. `__GPU_PARTITION__` (SLURM) or the PBS `select=...` line (GPU node syntax)
3. `__CUDA_MODULE__ __CUDA_AWARE_MPI_MODULE__` — a CUDA toolkit **and a CUDA-aware MPI**
   (e.g. `cuda/12.2` + `openmpi/4.1.x-cuda`, or Intel `impi/2021.15`).
In `GNUmakefile.a100`, set `NVCC_CCBIN` to the `g++` that pairs with your CUDA
(CUDA 12.x → `g++-11` or `g++-12`).

## STEP 1 — build (login/compile node with the CUDA module loaded)
```bash
cd /home/qiaoj/testcerisse/cerisse/exm/scaling      # or copy the case to the cluster first
module load __CUDA_MODULE__ __CUDA_AWARE_MPI_MODULE__
make -f GNUmakefile.a100 -j 8                        # ~10-20 min first build (AMReX+chemistry)
# -> produces  main3d.gnu.MPI.CUDA.ex
```

## STEP 2 — smoke test (1 GPU, ~seconds) to confirm the binary runs on a GPU
```bash
# interactive GPU or a 1-task job:
./main3d.gnu.MPI.CUDA.ex inputs max_step=20 amr.n_cell="128 128 128" amr.plot_files_output=0
# expect: it starts, prints steps, and 'Device: NVIDIA A100' in the AMReX banner.
```

## STEP 3 — 8-GPU scaling run (the verification)
```bash
sbatch run_a100.slurm        # SLURM
# or
qsub   run_a100.pbs          # PBS
```
Then compare your per-step time against the reference (below). Run 1, 2, 4, 8 GPUs by
editing `--ntasks`/`--gres` (SLURM) or `ngpus`/`-n` (PBS) to measure your own scaling.

## Reference strong scaling (256³, 500 steps, from `exm/scaling/numbers_marenostrum`)
| GPUs | Nodes | s/step-set | Efficiency |
|------|-------|-----------|-----------|
| 1 | 1 | 49.5 | 1.00 |
| 2 | 1 | 30.6 | 0.81 |
| 4 | 1 | 17.8 | 0.70 |
| **8** | 2 | **9.33** | **0.66** |
| 16 | 4 | 4.87 | 0.64 |
| 32 | 8 | 2.69 | 0.58 |
256³ / 8 GPUs = 2M cells/GPU (comm-bound regime). For higher efficiency use a bigger
grid (512³ → 16M cells/GPU).

---

## Production run (after Step 3 passes)
Copy `GNUmakefile.a100` into your target case directory (adjust the `AMR_SOLVER`/
`AMREX_HOME` relative paths for the new depth), set `USE_GPIBM=TRUE` only if the case
uses immersed-boundary geometry, then `make -f GNUmakefile.a100 -j 8` and launch with the
same script pointed at the case's `inputs`. Always keep `amrex.use_gpu_aware_mpi=1` and
add `amrex.abort_on_out_of_gpu_memory=1` to catch OOM early.

**Sweet spot (where A100 pays off):** large, uniform-ish 3D runs with many cells/GPU and
infrequent regridding — e.g. a high-resolution 3D underexpanded jet. A100 full-rate FP64
is exactly what the local consumer GPU lacked.

**Avoid on multi-GPU:** deep-AMR + moving-immersed-boundary (SRP) — documented as
overhead-bound on A100 (2 GPUs ≈ 1 GPU) with regrid-driven memory growth → OOM. That is a
latency-bound regime; more GPUs do not help it.

---

## One-line summary of the config
`USE_MPI=TRUE  USE_CUDA=TRUE  USE_OMP=FALSE  CUDA_ARCH=8.0`, 1 MPI rank per GPU (8 ranks),
`amrex.use_gpu_aware_mpi=1`, CUDA-aware MPI in the module env.
