# NASA Langley UPWT Test 1853 engineering-drawing extracts

Official source:

- NASA/TP-2014-218256, *Supersonic Retropropulsion Test 1853 in NASA LaRC Unitary Plan Wind Tunnel Test Section 2*, Scott A. Berry and Matthew N. Rhode, April 2014.
- NTRS record: <https://ntrs.nasa.gov/citations/20140006403>
- Official PDF: <https://ntrs.nasa.gov/api/citations/20140006403/downloads/20140006403.pdf>
- NTRS document ID: `20140006403`.

The files in this directory are lossless PNG rotations of the engineering-drawing raster images extracted from the official NASA PDF. The approximately 222 MB source PDF is intentionally not stored in the repository. `pNNN` means the printed report page number, not a viewer-dependent PDF index.

| File | Printed page | Drawing | Sheet | Content |
|---|---:|---|---:|---|
| `report_p092_DWG1168288_sheet1_aftbody.png` | 92 | 1168288 | 1 | Port aftbody cover and axial/outer dimensions |
| `report_p093_DWG1168288_sheet2_aftbody_sections.png` | 93 | 1168288 | 2 | Port aftbody pressure-port and section details |
| `report_p094_DWG1168288_sheet3_aftbody_starboard.png` | 94 | 1168288 | 3 | Starboard aftbody cover |
| `report_p097_DWG1168291_sheet1_sting.png` | 97 | 1168291 | 1 | Model sting and internal supply passage |
| `report_p098_DWG1168292_sheet1_front_plenum.png` | 98 | 1168292 | 1 | Laser-sintered front plenum |
| `report_p099_DWG1168293_sheet1_interface_plate.png` | 99 | 1168293 | 1 | Nozzle interface plate |
| `report_p100_DWG1168294_sheet1_center_nozzle_plug.png` | 100 | 1168294 | 1 | Center plug and instrumented center nozzle assembly |
| `report_p101_DWG1168294_sheet2_center_nozzle_contour.png` | 101 | 1168294 | 2 | Center nozzle nominal internal contour |
| `report_p102_DWG1168295_sheet1_half_radius_nozzles_plugs.png` | 102 | 1168295 | 1 | Half-radius plugs and nozzle assemblies |
| `report_p103_DWG1168295_sheet2_half_radius_nozzle_contour.png` | 103 | 1168295 | 2 | Uninstrumented half-radius nozzle nominal contour |
| `report_p104_DWG1168295_sheet3_instrumented_half_radius_nozzle.png` | 104 | 1168295 | 3 | Instrumented half-radius nozzle |
| `report_p105_DWG1168296_sheet1_forebody_shell.png` | 105 | 1168296 | 1 | 70-degree sphere-cone forebody shell, coordinate convention, pressure ports |
| `report_p106_DWG1168296_sheet2_forebody_section.png` | 106 | 1168296 | 2 | Forebody rear section and attachment details |
| `report_p108_DWG1168298_sheet1_manifold_machining.png` | 108 | 1168298 | 1 | Brazed manifold assembly, as-machined details |
| `report_p109_DWG1168298_sheet2_manifold_weldment.png` | 109 | 1168298 | 2 | Front-plenum/interface-plate assembly and braze details |

## Interpretation cautions

- The drawings mark several exposed forebody, plug, nozzle-lip, and plenum surfaces as defined by an external CAD file. That CAD file is not included in the public report; the drawings alone do not uniquely recover those surfaces.
- The nozzle dimensions in drawings 1168294 and 1168295 are nominal. Use the as-built throat diameter, virtual-exit diameter/area, area ratio, and divergence angle from Table A-2 for test-specific nozzle flow geometry.
- The `x = 1.0066 in` station on drawing 1168288 is the forward station of the aftbody cover. It is not proof of a forebody cone/cylinder tangency station; the forebody and cover overlap.
- A nominal nose radius near `1.000 in` is strongly supported by the rounded pressure-port coordinates and the 70-degree sphere-cone construction, but no explicit as-built `R_n = 1.000 in` dimension was found in the public drawing set. Treat an analytic one-inch sphere as a documented reconstruction, not a directly stated as-built dimension.
- Table A-2 calls the nozzle exit a **virtual exit**. For off-axis nozzles the externally visible opening is the oblique intersection of the axial passage with the CAD-defined 70-degree forebody surface; it is not necessarily a physical circular plane equal to the virtual-exit disk.
- The part lists reference final assembly drawing `1168299` as `NEXT ASSY`, but that drawing is not reproduced in the public Appendix A drawing set. Consequently the public report does not supply a verified TSP-to-manifold/sting assembly datum chain; global nozzle virtual-exit stations and the global sting translation must remain unspecified for a strict as-built model.
