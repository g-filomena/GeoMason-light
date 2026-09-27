/*
 * Copyright (c) 2023 Gabriele Filomena University of Liverpool, UK
 *
 * This program is free software: it can redistributed and/or modified under the terms of the GNU
 * General Public License 3.0 as published by the Free Software Foundation.
 *
 * See the file "LICENSE" for more information
 */
package sim.field.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Envelope;
import sim.field.grid.IntGrid2D;

class GridLayerTest {

  @Test
  @DisplayName("a cell's polygon contains that cell's centre point")
  void toPolygonAndToPointDescribeTheSameCell() {
    GridLayer layer = new GridLayer(new IntGrid2D(3, 2));
    layer.setMBR(new Envelope(0, 3, 0, 2));

    for (int x = 0; x < 3; x++) {
      for (int y = 0; y < 2; y++) {
        assertTrue(layer.toPolygon(x, y).contains(layer.toPoint(x, y)), "cell " + x + "," + y);
      }
    }
  }

  @Test
  @DisplayName("row 0 is the top row, for polygons as for points")
  void rowZeroIsAtTheTop() {
    GridLayer layer = new GridLayer(new IntGrid2D(2, 2));
    layer.setMBR(new Envelope(0, 2, 0, 2));

    assertEquals(2.0, layer.toPolygon(0, 0).getEnvelopeInternal().getMaxY(), 1e-9);
  }

  @Test
  @DisplayName("setting the grid after the MBR sizes the pixels from the MBR")
  void setGridAfterSetMBR() {
    GridLayer layer = new GridLayer();
    layer.setMBR(new Envelope(0, 10, 0, 20));
    layer.setGrid(new IntGrid2D(5, 4));

    assertEquals(2.0, layer.getPixelWidth(), 1e-9);
    assertEquals(5.0, layer.getPixelHeight(), 1e-9);
    assertEquals(3, layer.toXCoord(7.0));
  }
}
