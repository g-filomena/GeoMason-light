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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import sim.field.geo.VectorLayer;
import sim.testing.CapturedLog;
import sim.util.geo.MasonGeometry;

class ShapeFileImporterTest {

  private static final int ID_WIDTH = 11;

  /** A point record, or a null-shape record when {@code point} is null. */
  private static byte[] record(int number, double[] point) {
    int contentBytes = point == null ? 4 : 20;
    ByteBuffer buffer = ByteBuffer.allocate(8 + contentBytes);
    buffer.order(ByteOrder.BIG_ENDIAN).putInt(number).putInt(contentBytes / 2);
    buffer.order(ByteOrder.LITTLE_ENDIAN).putInt(point == null ? 0 : 1);
    if (point != null) {
      buffer.putDouble(point[0]).putDouble(point[1]);
    }
    return buffer.array();
  }

  /** A .shp of the given records: a 100-byte header the importer skips, then the records. */
  private static byte[] shp(byte[]... records) {
    ByteArrayOutputStream out = new ByteArrayOutputStream();
    out.write(new byte[100], 0, 100);
    for (byte[] record : records) {
      out.write(record, 0, record.length);
    }
    return out.toByteArray();
  }

  /** A .dbf with one numeric field, ID, holding the given values. */
  private static byte[] dbf(String... ids) {
    int headerSize = 32 + 32 + 1;
    int recordSize = 1 + ID_WIDTH;
    ByteBuffer buffer = ByteBuffer.allocate(headerSize + recordSize * ids.length)
        .order(ByteOrder.LITTLE_ENDIAN);
    buffer.put((byte) 3).put(new byte[3]).putInt(ids.length).putShort((short) headerSize)
        .putShort((short) recordSize).put(new byte[20]);
    byte[] name = new byte[11];
    name[0] = 'I';
    name[1] = 'D';
    buffer.put(name).put((byte) 'N').put(new byte[4]).put((byte) ID_WIDTH).put((byte) 0)
        .put(new byte[14]);
    buffer.put((byte) 0x0D);
    for (String id : ids) {
      String padded = String.format("%" + ID_WIDTH + "s", id);
      buffer.put((byte) ' ').put(padded.getBytes(StandardCharsets.US_ASCII));
    }
    return buffer.array();
  }

  @Test
  @DisplayName("a null shape is skipped and logged, and the features after it are still read")
  void nullShapesAreSkipped(@TempDir Path dir) throws Exception {
    File shp = dir.resolve("points.shp").toFile();
    File dbf = dir.resolve("points.dbf").toFile();
    Files.write(shp.toPath(), shp(record(1, new double[] {1, 1}), record(2, null),
        record(3, new double[] {3, 3})));
    Files.write(dbf.toPath(), dbf("1", "2", "12345678901"));
    VectorLayer layer = new VectorLayer();

    try (CapturedLog log = new CapturedLog(ShapeFileImporter.class)) {
      ShapeFileImporter.read(shp.toURI().toURL(), dbf.toURI().toURL(), layer);
      assertEquals(1, log.messages().size());
      assertTrue(log.messages().get(0).startsWith("Skipped 1 records with no geometry"));
    }

    assertEquals(2, layer.size());
    MasonGeometry last = layer.getGeometries().get(1);
    assertEquals(3.0, last.getGeometry().getCoordinate().x, 0.0);
    assertEquals(12345678901L, last.getAttributes().get("ID").getValue());
  }

  @Test
  @DisplayName("numeric fields become the narrowest number that holds them")
  void numbersAreParsedWithoutOverflow() {
    assertEquals(42, ShapeFileImporter.parseNumber("42"));
    assertEquals(12345678901L, ShapeFileImporter.parseNumber("12345678901"));
    assertEquals(1.5, ShapeFileImporter.parseNumber("1.5"));
    assertEquals(100000.0, ShapeFileImporter.parseNumber("1e5"));
    assertEquals("*****", ShapeFileImporter.parseNumber("*****"));
  }
}
