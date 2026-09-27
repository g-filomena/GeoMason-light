/*
 * Copyright (c) 2023 Gabriele Filomena University of Liverpool, UK
 *
 * This program is free software: it can redistributed and/or modified under the terms of the GNU
 * General Public License 3.0 as published by the Free Software Foundation.
 *
 * See the file "LICENSE" for more information
 */
package sim.io.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Polygon;

class ImporterUtilsTest {

  private static final GeometryFactory FACTORY = new GeometryFactory();

  /** A square ring, clockwise (a Shapefile shell) or counter-clockwise (a hole). */
  private static LinearRing square(double minX, double minY, double side, boolean clockwise) {
    Coordinate a = new Coordinate(minX, minY);
    Coordinate b = new Coordinate(minX, minY + side);
    Coordinate c = new Coordinate(minX + side, minY + side);
    Coordinate d = new Coordinate(minX + side, minY);
    Coordinate[] ring = clockwise ? new Coordinate[] {a, b, c, d, a}
        : new Coordinate[] {a, d, c, b, a};
    return FACTORY.createLinearRing(ring);
  }

  @Test
  @DisplayName("each hole goes to the shell containing it, and no hole becomes a shell")
  void holesGoToTheirOwnShell() {
    LinearRing small = square(0, 0, 1, true);
    LinearRing hole = square(12, 12, 2, false);
    LinearRing large = square(10, 10, 10, true);

    Geometry polygon = ImporterUtils.createPolygon(new Geometry[] {small, hole, large});

    assertTrue(polygon instanceof MultiPolygon);
    assertTrue(polygon.isValid());
    assertEquals(1.0 + 100.0 - 4.0, polygon.getArea(), 1e-9);
    Polygon smallPart = (Polygon) polygon.getGeometryN(0);
    Polygon largePart = (Polygon) polygon.getGeometryN(1);
    assertEquals(0, smallPart.getNumInteriorRing());
    assertEquals(1, largePart.getNumInteriorRing());
  }

  @Test
  @DisplayName("rings all wound the wrong way are read as shells, not dropped")
  void reversedRingsAreShells() {
    Geometry polygon = ImporterUtils.createPolygon(
        new Geometry[] {square(0, 0, 1, false), square(5, 5, 1, false)});

    assertEquals(2.0, polygon.getArea(), 1e-9);
  }

  /** A stream that hands back one byte per read, as network streams may. */
  private static InputStream trickle(byte[] bytes) {
    return new ByteArrayInputStream(bytes) {
      @Override
      public synchronized int read(byte[] b, int off, int len) {
        return super.read(b, off, Math.min(len, 1));
      }
    };
  }

  @Test
  @DisplayName("short reads are completed rather than reported as the end of the stream")
  void shortReadsAreCompleted() throws IOException {
    byte[] bytes = ByteBuffer.allocate(12).order(ByteOrder.LITTLE_ENDIAN).putInt(42)
        .putDouble(2.5).array();
    InputStream in = trickle(bytes);

    assertEquals(42, ImporterUtils.readInt(in, true));
    assertEquals(2.5, ImporterUtils.readDouble(in, true), 0.0);
  }
}
