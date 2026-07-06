---
icon: table-cells-large
---

# Boundaries

## Boundary Conditions

There are generally two ways to define boundary conditions in the code:

1. **Through the input**: This involves specifying the boundary conditions for the domain directly in the input file, typically for simpler conditions.
2. **Through the `prob.h` file**: This allows for more customized or complex boundary conditions.

For internal geometries, the immersed-boundary and embedded-boundary treatments are documented separately in [`IBM/ibm_docs/`](../IBM/ibm_docs/) (see `ibm_technical.pdf` and `ibm_thesis.pdf`).

#### domain boundary conditions

They are defined in the `input` file and apply to all faces in that direction (see details in input).

```ini
# 0 = Interior                               3 = Symmetry
# 1 = Inflow / UserBC                        4 = SlipWall =3
# 2 = Outflow (First Order Extrapolation)    5 = NoSlipWall (adiabatic)
# 6 = UserBC / ext_dir                       7 = characteristic far-field
cns.lo_bc = 1 5 0
cns.hi_bc = 2 5 0
```

#### user-derived boundaries

During run-time, complex boundary conditions can be set-up by defining the `bcnormal` function in `prob.h`

```cpp
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[5],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[5], const int idir,
         const int sgn, const Real time, GeometryData const &,
         ProbClosures const &closures, ProbParm const &prob_parm)
         {
            ...
         }
```

This function is called for each face on the boundary. It defines the array `s_ext[5]`, which specifies the values at the ghost points required to build the flux across domain boundaries.

| Option      | Type         | Dimensions | Description                             |
| ----------- | ------------ | :--------: | --------------------------------------- |
| `x`         | Real         |  SPACEDIM  | ghost cell cooordinates                 |
| `dratio`    | Real         |      1     | ghost/first internal distance ratio     |
| `s_int`     | Real         |    NCONS   | flow state inside domain                |
| `s_ext`     | Real         |    NCONS   | flow state to be filled                 |
| `idir`      | int          |      1     | direction (0: x, 1: y, 2: z)            |
| `sign`      | int          |      1     | high or low boundary (1: low, -1: high) |
| `time`      | Real         |      1     | time                                    |
| `geomdata`  | GeometryData |      -     | domain geometry data                    |
| `prob_parm` | ProbParm     |      -     | problem parameter data                  |

For example, a slip wall in the bottom wall (y-direction) would be defined as

```cpp
    if (idir == 1 && sgn == 1) {
        s_ext[URHO] = s_int[URHO];
        s_ext[UMX]  = s_int[UMX];
        s_ext[UMY]  = -s_int[UMY];
        s_ext[UMZ]  = s_int[UMZ];
        s_ext[UET]  = s_int[UET];
    }
```

If the boundary is not defined, the code will use the one specified in the input for that particular face.

NOTE: Boundary conditions are specified for the _conserved_ variables so appropiate conservations may be required.

#### multidimensional LODI / NSCBC hook

For open boundaries that need tangential derivatives, a `prob.h` file may also
define the optional hook

```cpp
bool bcnormal_lodi(const IntVect& iv, Array4<Real> const& state,
                   const Real* x, Real dratio,
                   const Real* s_int, const Real* s_refl, Real* s_ext,
                   int idir, int sgn, Real time,
                   GeometryData const& geomdata,
                   ProbClosures const& closures,
                   ProbParm const& prob_parm);
```

If this function exists and returns `true`, its `s_ext` value is used.  If it is
absent or returns `false`, Cerisse falls back to `bcnormal`.  The reusable
characteristic-amplitude LODI helper is `manual_bc_t::bc_nscbc_lodi_farfield`
in `src/set/bc_types.h`; see `docs/bc_multid_lodi.md` for the model and
current validation status.

## Internal solid boundaries

The immersed-boundary (ghost-cell) and embedded-boundary (cut-cell) treatments
of internal solid geometry are documented in [`IBM/ibm_docs/`](../IBM/ibm_docs/). See
`ibm_technical.pdf` for the implementation reference and `ibm_thesis.pdf` for
the academic-style write-up. Setup steps for `GNUmakefile`, `input`, and
`prob.h`, the wall-model templates, and the `ibfab` marker semantics are
covered there.
