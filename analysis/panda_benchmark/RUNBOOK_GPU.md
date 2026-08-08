# AWS GPU runbook for Cerisse

Last audited: 2026-07-26 (us-west-2 p4d bootstrap added, sections 10-11)

This document separates live AWS facts, campaign-specific requirements, and
historical workarounds. Commands that purchase or launch capacity remain
explicit and fail closed.

## 1. Scope

The staged H100 campaign is:

```text
/shared/h100_gpibm_srp_20260726
```

Its method identity is:

```text
pure_shared_gp_full_cartesian
```

The executable uses the frozen pure GP-IBM production path:

```text
shared GP
+ cell-average-aware BI-CWLS
+ curvature-compatible Euler-slip pressure closure
+ one-sided entropy extension
+ full-cell Cartesian LLF-WENO-Z5
+ Cartesian shared-face positivity limiter
+ SSPRK(4,3)
```

It contains no cut-cell, EB, cut-control, open-subface, or aggregate update.
Do not use `/shared/h100_launch/bootstrap_run.sh`; it belongs to an unrelated
older campaign.

## 2. Before purchasing capacity

Run on the HeadNode:

```bash
CONFIRM=DRY bash /shared/h100_launch/buy_launch_h100.sh
```

The AWS API currently accepts `24` as the capacity-duration query window and
may return offerings of different actual lengths. Select using the returned
`Hours`, `Start`, `AZ`, and `Fee`. Offering IDs expire quickly and must never
be copied from an old log.

Purchase remains a separate explicit action:

```bash
OFFERING=cb-... CONFIRM=BUY \
  bash /shared/h100_launch/buy_launch_h100.sh
```

Record the returned reservation ID. Do not launch until the reservation state
is `active`; a bank charge or `scheduled` state is not sufficient.

The launch script checks:

- reservation state and available count;
- AZ-to-subnet mapping;
- an existing pending/running `p5.48xlarge`;
- AMI, key, security groups, and private subnet;
- capacity-block targeting.

Launch:

```bash
CR=cr-... CONFIRM=LAUNCH \
  bash /shared/h100_launch/buy_launch_h100.sh
```

The launcher writes:

```text
/shared/h100_launch/active_h100.env
```

Run scripts use its reservation end time and refuse work that would consume
the final 30-minute safety margin.

## 3. SSH and ownership

The HeadNode does not contain the private key used by manually launched GPU
instances. From the local workstation:

```bash
ssh -i ~/.ssh/aws-hpc-us.pem \
  -o StrictHostKeyChecking=no \
  -o ProxyJump=aws-hpc \
  ec2-user@PRIVATE_IP
```

Run Cerisse as `ec2-user`, not root. User-data runs as root and can leave
root-owned NFS output because `/shared` is exported with `no_root_squash`.
The campaign preflight verifies that `ec2-user` can write its output tree.

Do not use broad commands such as `pkill -f make` or `pkill -f main`; they can
match the invoking shell. Kill a verified PID only.

## 4. Mandatory H100 preflight

On the p5 node:

```bash
bash /shared/h100_gpibm_srp_20260726/bootstrap_h100.sh
```

This mounts `/shared`, starts Fabric Manager, and executes
`preflight_h100.sh`. It starts no Cerisse simulation.

The preflight fails unless all of the following hold:

- exactly eight H100 GPUs are visible;
- Fabric Manager is active;
- `/shared` is mounted, writable, and has at least 250 GiB free;
- the executable SHA-256 matches its manifest;
- all dynamic libraries resolve;
- `libcudart` resolves to `/shared/cuda-12.6/lib64`;
- the executable contains `sm_90` device code;
- the pure GP-IBM method identity is embedded;
- both geometry hashes match;
- inputs remain level-0, `max_grid_size=256`,
  `blocking_factor=64`, and reflux-off;
- MPI local ranks map one-to-one to devices 0 through 7.

The preflight report is saved in:

```text
/shared/h100_gpibm_srp_20260726/manifests/
```

## 5. CUDA and MPI facts

Every campaign script sources:

```text
/shared/h100_gpibm_srp_20260726/source/IBM/cases/
srp_planar_70deg_nozzle_2d/aws_h100_20260726/h100_env.sh
```

This selects CUDA 12.6 and Amazon OpenMPI 4.1.7. Without it, the AMI may select
CUDA 13 and load a different runtime.

The installed Amazon OpenMPI reports:

```text
mpi_built_with_cuda_support = false
```

AMReX therefore uses host-staged MPI. This is valid for correctness but can
limit multi-GPU scaling and prevent NVSwitch from delivering its full
potential. Do not set a CUDA-aware MPI option against this library. Select the
rank count from the measured benchmark.

The staged AMReX commit already contains a CCCL-version guard and the exact
CUDA 12.6 build succeeded. Do not apply the historical
`cuda::minimum -> cub::Min` source patch to this campaign. Any source change
invalidates the executable hash and requires a clean rebuild and regression.

## 6. Grid and memory risks

Prepared grids:

| Grid | Cells | Boxes | Lip-radius resolution |
|---|---:|---:|---:|
| r4h | 1600 x 768 | 21 | 4 cells |
| r8h | 3200 x 1536 | 78 | 8 cells |

The box counts follow `max_grid_size=256` and `blocking_factor=64`.

The r8h problem has 4.92 million cells. A one-rank r8h run concentrates all
temporary WENO, GP, and SSPRK storage on one GPU and may exceed 80 GiB even
when the eight-rank run fits. Consequently:

- the scaling sequence starts at 8 ranks, then 4 and 2;
- the one-rank r8h point is disabled by default;
- enable it only with `BENCH_INCLUDE_NP1=1` after checking memory.

The current MPI is not CUDA-aware, so more ranks are not automatically faster.
Use measured coarse-step time, not nominal GPU count.

## 7. Scientific limitations of this campaign

The H100 inputs use:

```text
x in [-2, 48], y in [-12, 12], level-0 uniform grid
```

This is a resolved-lip, nozzle and near-field qualification domain. It is not
the earlier 16-body-width far-field domain and must not be cited as a
far-field-independent SRP solution.

The smooth quintic expansion has `Ae/A*=4`. The quasi-one-dimensional
area-Mach relation gives a design value near Mach 2.94, but the finite-angle
two-dimensional wall generates characteristic waves. The actual exit Mach,
mass flow, wall pressure, and shock-cell spacing must be measured from the
solution; Mach 2.94 is not an exact two-dimensional oracle.

The imposed `NPR=60` makes the jet underexpanded relative to the ambient state.
The startup is initially quiescent and ramps over four time units. During the
ramp, monitor:

- every SSPRK Forward-Euler bracket for admissibility;
- limiter locations and minimum theta;
- throat and lip GP/BIC fallback counts;
- inlet mass flow and exit-plane Mach profile;
- characteristic far-field failures;
- top/bottom/right boundary pressure histories.

The rounded-lip topology gate is necessary but not sufficient. A failure at
the throat before the causal front reaches the lip must be diagnosed as an
inlet/throat ownership problem, not repaired by changing lip support.

## 8. Required execution order

```bash
# 1. One step, one H100, no plotfile
CONFIRM=RUN \
  bash /shared/h100_gpibm_srp_20260726/run_topology_r4h.sh

# 2. r8h scaling, high rank count first
CONFIRM=RUN \
  bash /shared/h100_gpibm_srp_20260726/benchmark_r8h_scaling.sh

# 3. r4h transient pilot
CONFIRM=RUN \
  bash /shared/h100_gpibm_srp_20260726/run_r4h_pilot.sh

# 4. Reviewed r8h run
RUN_TIMEOUT_SEC=21600 CONFIRM=RUN \
  bash /shared/h100_gpibm_srp_20260726/run_r8h.sh
```

`RUN_TIMEOUT_SEC` must fit inside the remaining paid window plus the enforced
30-minute drain margin. A seven-hour block cannot use the script's twelve-hour
default.

Stop immediately if:

- the topology audit reports a cross-feature support conflict;
- the bow/nozzle baseline loses its expected wave system;
- density or internal energy becomes non-positive after a stage;
- global-theta, all-low, or repeated step rejection occurs unexpectedly;
- all eight GPUs do not map uniquely;
- the shared filesystem approaches the configured free-space floor.

Do not respond by tuning WENO weights, BI-CWLS order, ghost layers, pressure
floors, or dissipation.

## 9. Output and restart

Periodic output can dominate wall time. The topology and scaling gates disable
plotfiles and checkpoints. The pilot and r8h run keep periodic checkpoints.

The current campaign is level-0 with `cns.do_reflux=0`; the historical
`regrid_on_restart`/FluxRegister failure does not apply. If a future AMR run
enables reflux, restart/regrid/FluxRegister reconstruction needs its own
validated input pair and must not inherit this level-0 qualification.

Do not run destructive wildcard cleanup automatically. Before deleting data,
record the newest valid checkpoint and verify its header. Preserve:

- final and latest restartable checkpoint;
- run log;
- input and executable hashes;
- preflight and topology reports;
- positivity and momentum-budget diagnostics.

## 10. A100 differences

For `p4d.24xlarge`:

- eight A100 GPUs, normally 40 GiB each;
- 96 CPU cores rather than 192;
- compile and run a separate `sm_80` executable;
- use a separate build directory and executable manifest;
- never run the H100 `sm_90` binary as the A100 production artifact;
- remeasure memory and scaling because r8h has less memory headroom.

The p4d on-demand capacity observation is time-dependent. Query live capacity
and quota; do not treat a previous `InsufficientInstanceCapacity` result as a
permanent AWS guarantee.

An `sm_90` executable does not execute on A100 and cannot be salvaged. On
2026-07-26 eighteen `sm_90` binaries built successfully in us-east-2 were
discarded when the campaign moved to a p4d Capacity Block. Decide the target
architecture before compiling, not after.

## 11. Standalone instance outside ParallelCluster

Measured record, us-west-2 Capacity Block `cr-098160f9e4c211ffa`,
`p4d.24xlarge`, 8 x A100-40G, us-west-2c, AMI `ami-0649a2303f6276012`
(Deep Learning AMI, user `ubuntu`), instance `i-0936ebc858384e2be`.

Instance running 15:04:24Z; first successful `sm_80` executable 15:23:56Z.
Nineteen minutes, of which the compile itself was 70 s. The remaining time was
the six items below. Every one is avoidable by preflight.

| # | Symptom | Cause | Action | Cost |
|---|---|---|---|---:|
| 1 | `Unable to locate credentials` | `run-instances` outside ParallelCluster attaches no IAM role | `aws ec2 associate-iam-instance-profile --instance-id i-... --iam-instance-profile Name=qiaoj-ec2-ssm-core`; verify with `aws sts get-caller-identity` | 3 min |
| 2 | `s3 ls` / `s3 cp` AccessDenied with the role attached | the SSM-core role carries no S3 permission | presigned URLs, below | 4 min |
| 3 | `/scratch` absent, `mdadm` absent, all eight NVMe devices busy | the DL AMI has already assembled the instance store as LVM `vg.01-lv_ephemeral` on `/opt/dlami/nvme` | `df -h \| grep nvme` first, then `sudo ln -sfn /opt/dlami/nvme /scratch` | 4 min |
| 4 | `mpicxx: command not found` | this AMI ships CUDA and NCCL but no MPI development package; there is no `/opt/amazon/openmpi` and no module system | `sudo apt-get install -y libopenmpi-dev openmpi-bin` | 1 min |
| 5 | `namespace "cuda" has no member "minimum"` | `/usr/local/cuda` points at 13.2, whose CCCL has no `cuda/functional` | select a version whose CCCL provides it, below | 2 min |
| 6 | `No rule to make target '/shared/cerisse_reflux/src/Make.CNS'` | the case `GNUmakefile` carries the absolute solver path of the previous region | rewrite `AMR_SOLVER` in every case, below | 3 min |

### Moving files in without S3 permission

Generate the URLs where credentials exist, fetch them where they do not. The
region must be given explicitly or SigV4 signs for the wrong endpoint:

```bash
# on the HeadNode, which has credentials
aws s3 presign s3://cerisse-usw2-results-332677055650/staging/FILE \
  --expires-in 43200 --region us-west-2
# on the GPU instance
curl -sSfL -o FILE 'URL'
```

166 MiB of case files transferred in seconds within the region. Presigned URLs
carry the signature in the query string: keep them out of shell history that is
shared, and use the shortest lifetime that covers the window.

### Selecting a CUDA version instead of patching AMReX

Test the candidate before building:

```bash
grep -c minimum /usr/local/cuda-12.9/include/cuda/functional   # 1 -> usable
```

On this AMI 12.8, 12.9 and 13.2 were installed; 12.9 provides
`cuda::minimum`, so `CUDA_HOME=/usr/local/cuda-12.9` builds unmodified AMReX.
Prefer this to the historical `cuda::minimum -> cub::Min` source patch: the
patch changes the source, invalidates the executable hash, and forces a clean
rebuild and regression. Diagnosing this from scratch cost 20 minutes in
us-east-2 on 2026-07-18; the header test costs seconds.

### Absolute paths do not survive a region change

Cerisse cases reference the solver by absolute path:

```text
AMR_SOLVER = /shared/cerisse_reflux
AMREX_HOME = $(AMR_SOLVER)/lib/amrex
include $(AMR_SOLVER)/src/Make.CNS
```

On a new machine the tree is elsewhere and `make` aborts in under a second.
Before the first build:

```bash
grep -rl '/shared' */GNUmakefile
sed -i 's#AMR_SOLVER = /shared/cerisse_reflux#AMR_SOLVER = /scratch/cerisse#' */GNUmakefile
```

Audit `inputs`, restart paths, and any staging script for the same reason. A
near-instant `make` failure is almost always a path, not a compiler.

### One case per GPU

Round-robin assignment by launch order, `g=$(( (g+1) % 8 ))`, silently places
the ninth case on a device that is already loaded, and it fails with an
out-of-memory error minutes later. Launch in waves of eight and wait for the
whole wave. Verify with `pgrep -af main2d` and `nvidia-smi`; `ps -eo args |
grep -c` has misreported live processes as dead in both directions.

### Preflight for a standalone GPU instance

```bash
nvidia-smi --query-gpu=index,name,memory.total --format=csv   # eight devices
df -h | grep -E 'nvme|dlami'                                  # scratch space
aws sts get-caller-identity                                   # role attached
which mpicxx nvcc                                             # toolchain
grep -c minimum $CUDA_HOME/include/cuda/functional            # must be 1
grep -h AMR_SOLVER */GNUmakefile | sort -u                    # path resolves
```

All six pass in under a minute and cover every failure recorded above.

## 12. End of paid window

Before termination:

1. stop launching new stages at the safety deadline;
2. verify the latest checkpoint;
3. copy final logs and manifests to the campaign root;
4. confirm no `mpirun` or Cerisse process remains;
5. record the instance ID, reservation ID, end time, and completed step;
6. terminate the instance explicitly.

From the HeadNode:

```bash
aws ec2 terminate-instances --region us-east-2 --instance-ids i-...
```

Capacity Block payment is not refunded by early termination, but termination
prevents accidental work or storage activity on a stale instance.
