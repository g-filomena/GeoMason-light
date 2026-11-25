package sim.io.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import mil.nga.geopackage.GeoPackage;
import mil.nga.geopackage.GeoPackageManager;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import sim.field.geo.VectorLayer;
import sim.testing.Fixtures;
import sim.util.geo.AttributeValue;
import sim.util.geo.MasonGeometry;

class GeoPackageExporterTest {

  /** A point whose attributes can be read once, for the column types, and not when it is written. */
  private static final class FailsWhenWritten extends MasonGeometry {
    private static final long serialVersionUID = 1L;
    private int reads = 0;

    FailsWhenWritten() {
      super(Fixtures.point(20, 0).getGeometry());
    }

    @Override
    public Map<String, AttributeValue> getAttributes() {
      if (reads++ > 0) {
        throw new IllegalStateException("unwritable");
      }
      return super.getAttributes();
    }
  }

  private static int rowsIn(File file, String table) throws Exception {
    try (GeoPackage geoPackage = GeoPackageManager.open(file)) {
      return geoPackage.getFeatureDao(table).count();
    }
  }

  @Test
  void everyRowIsWritten(@TempDir Path dir) throws Exception {
    List<MasonGeometry> points = new ArrayList<>();
    for (int i = 0; i < 500; i++) {
      MasonGeometry point = Fixtures.point(i, 2 * i);
      point.addIntegerAttribute("rank", i);
      points.add(point);
    }
    File file = dir.resolve("points.gpkg").toFile();

    VectorLayer.writeGPKG(file.getPath(), new VectorLayer(points));

    VectorLayer read = new VectorLayer();
    VectorLayer.readGPKG(file.toURI().toURL(), read);
    assertEquals(500, read.size());
    assertEquals(499, read.getGeometries().get(499).getIntegerAttribute("rank"));
    assertEquals(998.0, read.getGeometries().get(499).getGeometry().getCoordinate().y);
  }

  @Test
  @DisplayName("a write that fails part-way leaves no rows behind")
  void aFailedWriteIsRolledBack(@TempDir Path dir) throws Exception {
    File file = dir.resolve("points.gpkg").toFile();
    VectorLayer layer =
        new VectorLayer(
            new ArrayList<>(
                Arrays.asList(Fixtures.point(0, 0), Fixtures.point(10, 0), new FailsWhenWritten())));

    assertThrows(
        IllegalStateException.class, () -> VectorLayer.writeGPKG(file.getPath(), layer));

    assertEquals(0, rowsIn(file, "points"));
  }
}
