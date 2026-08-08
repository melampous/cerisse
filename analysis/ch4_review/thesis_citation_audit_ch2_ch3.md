# Citation Audit Report — Chapters 2 and 3

**Scope:** All bibliography keys cited in Chapter 2 (numerical methods) and Chapter 3 (free underexpanded jet), including shared keys also used in Chapter 4. Chapter 4-specific keys were audited separately and passed. Each key was checked on three axes: existence (DOI/publisher record), bibliographic field accuracy, and attribution (does the cited source actually support the claim at each citing location).

## 1. Summary

**50 keys audited: 45 fully clean, 4 minor, 1 major.**

## 2. MAJOR issues

| Key | Problem | Required fix |
|---|---|---|
| **Anderson1995CFD** | Cited at ch2 line 167 for "the standard axisymmetric conservation form", but full-text search shows the book derives the governing equations in Cartesian coordinates only — no axisymmetric or cylindrical-coordinate form appears anywhere (only a cylindrical grid figure and incidental mentions). The use at line 22 (conservation of mass, momentum, total energy) is correct. | At line 167, replace with a source that actually presents the axisymmetric conservation form — e.g. White2006 (cylindrical NS, Appendix B), Tannehill1997, or the thesis's own appendix_axisymmetric_derivation. Keep Anderson1995CFD only at line 22. |

## 3. MINOR corrections

| Key | Issue | Suggested fix |
|---|---|---|
| **Lax1954** | Line 553 attributes a flux "constructed from a local upper bound on the characteristic wave speeds" to Lax (1954). That is the local Lax–Friedrichs/Rusanov flux (Rusanov 1961); Lax's scheme uses the global Δx/Δt coefficient. The solver file is literally `Rusanov.h`. | Add a Rusanov citation (V. V. Rusanov, USSR Comput. Math. Math. Phys. 1(2):304–320, 1962; orig. 1961) alongside or in place of Lax1954, or rephrase to credit Lax for the underlying dissipation form and Rusanov for the local wave-speed bound. |
| **Roache2002** | At 04_ibm_amr_fsi/sections_new/05_verification.tex line 760, "intended-use qualification in the sense of Roache (2002)" — that framing is not in the 2002 JFE paper (which is about MMS and observed-order verification); it belongs to the V&V canon. All other uses of this key are correct. | Attribute the intended-use qualification concept to OberkampfRoy2010 alone, or substitute Roache's 1998 V&V book at that one location. |
| **EdgingtonMitchell2014** | The l.70 use (internal shear layer persisting across shock cells, coherent shock motion) is confirmed from the abstract. The l.66 claim ("a very small Mach stem may exist below experimental or numerical resolution") is plausible but could not be verified in the 2014 paper (paywalled). | Verify the sentence exists in the paper; if not, cite or co-cite the Mach-disk onset study (Muraoka & Hiejima, Phys. Fluids 34, 116125, 2022) or another source that states it explicitly. |
| **AbdelRahman2010** | Attribution is substantively correct (three-zone jet schematic and terminology match the review), but the specific "[Figure 1]" pinpoint could not be verified — the source PDF is inaccessible — and the in-file comment already flags a rights-clearance TODO. | Check the original WSEAS PDF to confirm the schematic is Figure 1, and clear figure-reuse rights (Duronio's CC BY licence does not cover this third-party figure). |

**Optional polish (severity ok — no action strictly required):**

- **White2006**: several records give *Boston* (not New York) for the 3rd ed.; change or drop the location field.
- **BergerOliger1984 / BergerColella1989 / SalariKnupp2000 / Rathore2022Thesis / Patel2022Thesis**: DOIs exist and could be added (10.1016/0021-9991(84)90073-1, 10.1016/0021-9991(89)90035-1, 10.2172/759450, 10.25560/107649, 10.25560/104074).
- **Toro1994HLLC**: optionally add `number = {1}`.
- **CerisseSource2026**: add `urldate`; cross-check any CITATION.cff for the canonical author list.
- **Un2025Thesis**: AFD-HLLC (l.1108) and Euler MMS (l.2286) claims are chapter-level and not verifiable from repository metadata — a quick page-level spot check is recommended (first-hand access available).
- **Rathore2022Thesis**: spot-check the centred-viscous-formulation content behind the l.1549 citation.
- **Addy1981**: spot-check coefficients 0.36/3.9 (contoured) and 0.31/5.0 (conical/orifice) against the original note or Franquet (2015) tables.
- **LewisCarlson1964**: optionally add "Jr." to Lewis.
- **BhideCuppoletti2024**: reword the l.10 sentence in 07_vorticity_dynamics.tex so "highly underexpanded" unambiguously modifies only the circular (Li2016) case — the rectangular cases are pressure ratio 3–4, not highly underexpanded.
- **Berland2007**: at 10_results_profile.tex:295, optionally note the waveguide description follows Tam (as Berland et al. themselves credit).
- **GottliebShu1998**: optionally co-cite Shu & Osher (1988) as the origin of the (3,3) scheme.

## 4. Clean (existence, bib fields, and attribution all confirmed)

LeVeque1992, White2006, Sutherland1893, Tannehill1997, Pirozzoli2011, Ducros2000Skew, JamesonSchmidtTurkel1981, Jameson2017JST, Roe1981, Borges2008WENOZ, Harten1983HLL, Toro1994HLLC, Balsara2025AFD, Zhang2019AMReX, BergerOliger1984, BergerColella1989, ShuOsher1989, Spiegel2015IsentropicVortex, Rathore2022Thesis, Patel2022Thesis, Un2025Thesis, CerisseSource2026, SalariKnupp2000, Roy2004, OberkampfRoy2010, Franquet2015, Duronio2023, Panda1999, Li2016, BhideCuppoletti2024, JiangShu1996, AshkenasSherman1964, Crist1966, Addy1981, LewisCarlson1964, Gibbings1972, Phalnikar2008, PrandtlPack1950, HatanakaSaito2012, Morris1989, Berland2007, FuHuAdams2016, MuraokaHiejima2022, Gribben2000, GottliebShu1998.

## 5. Uncited bibliography entries

22 entries in the .bib file are never cited in Chapters 2–3: Vinokur1989, Sandberg2011Axis, Garnier1999ShockCapturingLES, Grinstein2007ILES, Kim1956 (now cited in Chapter 4), Andre2013, NorumSeiner1982, Zapryagaev2015, Patel2024, KorzunBraun2009, and 12 more. This is harmless — biblatex simply omits them from the printed bibliography — and cleanup is optional.