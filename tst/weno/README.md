# WENO and TENO reconstruction tests

This test locks the reconstruction methods in `src/rhs/Weno.h`. It checks
WENO-Z5 and TENO constant preservation, smooth periodic derivative order,
TENO6 discontinuity selection, inactive-candidate NaN isolation, and the
ghost-point face-reconstruction plan. It also compiles `Weno.h` together with
`Weno_old.h` to guard the formal and archived public interfaces against
symbol collisions.

Run the CPU tests with:

```bash
make -j4
./main3d.gnu.ex
```

Compile the reconstruction methods with NVCC using:

```bash
make USE_CUDA=TRUE USE_MPI=FALSE \
  tmp_build_dir/o/3d.gnu.CUDA.EXE/cuda_compile.o
```

The TENO5 test targets fifth-order smooth reconstruction.  The TENO6 test
targets sixth-order smooth reconstruction.  These direct spatial tests do not
replace the Cartesian MMS and Shu--Osher solver studies.
