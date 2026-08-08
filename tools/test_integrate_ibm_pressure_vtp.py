#!/usr/bin/env python3
"""Regression tests for point-quadrature cut-control wall VTP files."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from integrate_ibm_pressure import oriented_measures, read_vtp


def wall_vtp(points_name: str) -> str:
    name = f' Name="{points_name}"' if points_name else ""
    return f"""<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="2" NumberOfVerts="2" NumberOfLines="0"
           NumberOfStrips="0" NumberOfPolys="0">
      <Points>
        <DataArray type="Float64"{name} NumberOfComponents="3" format="ascii">
          0 0 0  1 0 0
        </DataArray>
      </Points>
      <Verts>
        <DataArray type="Int64" Name="connectivity" format="ascii">0 1</DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">1 2</DataArray>
      </Verts>
      <CellData Scalars="Pressure">
        <DataArray type="Float64" Name="Pressure" format="ascii">2 3</DataArray>
        <DataArray type="Float64" Name="QuadratureMeasure" format="ascii">
          0.5 0.25
        </DataArray>
        <DataArray type="Float64" Name="FluidOutwardNormal"
                   NumberOfComponents="3" format="ascii">
          1 0 0  0 1 0
        </DataArray>
      </CellData>
    </Piece>
  </PolyData>
</VTKFile>
"""


class CutControlWallVtpTest(unittest.TestCase):
    def check_contract(self, points_name: str) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "wall.vtp"
            path.write_text(wall_vtp(points_name), encoding="ascii")
            surface = read_vtp(path)
            vectors, centroids = oriented_measures(surface)

        np.testing.assert_array_equal(
            centroids, np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        )
        np.testing.assert_array_equal(
            vectors, np.asarray([[-0.5, 0.0, 0.0], [0.0, -0.25, 0.0]])
        )
        pressure = np.asarray(surface["fields"]["Pressure"])
        force = -np.sum(pressure[:, None] * vectors, axis=0)
        np.testing.assert_array_equal(force, np.asarray([1.0, 0.75, 0.0]))

    def test_named_points_and_verts(self) -> None:
        self.check_contract("Points")

    def test_legacy_unnamed_points_and_verts(self) -> None:
        self.check_contract("")


if __name__ == "__main__":
    unittest.main()
