package sim.io.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.net.URL;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import mil.nga.geopackage.GeoPackage;
import mil.nga.geopackage.GeoPackageManager;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import sim.field.geo.VectorLayer;
import sim.testing.Fixtures;

class GeoPackageImporterTest {

  /** A GeoPackage with one feature table, "nodes", holding two points. */
  private static File nodesFile(Path dir) throws Exception {
    File file = dir.resolve("nodes.gpkg").toFile();
    VectorLayer layer =
        new VectorLayer(
            new ArrayList<>(Arrays.asList(Fixtures.point(0, 0), Fixtures.point(10, 0))));
    VectorLayer.writeGPKG(file.getPath(), layer);
    return file;
  }

  /** The same file with a copy of its table under a second name, as a stale layer would sit. */
  private static File twoTableFile(Path dir) throws Exception {
    File file = nodesFile(dir);
    try (GeoPackage geoPackage = GeoPackageManager.open(file)) {
      geoPackage.copyTable("nodes", "old_nodes");
    }
    return file;
  }

  @Test
  void aSingleTableFileIsRead(@TempDir Path dir) throws Exception {
    URL url = nodesFile(dir).toURI().toURL();
    VectorLayer layer = new VectorLayer();

    VectorLayer.readGPKG(url, layer);

    assertEquals(2, layer.size());
  }

  @Test
  @DisplayName("a file with two feature tables is refused, not read as their union")
  void aFileWithTwoTablesIsRefused(@TempDir Path dir) throws Exception {
    URL url = twoTableFile(dir).toURI().toURL();
    VectorLayer layer = new VectorLayer();

    IllegalStateException refused =
        assertThrows(IllegalStateException.class, () -> VectorLayer.readGPKG(url, layer));

    assertTrue(refused.getMessage().contains("nodes"), refused.getMessage());
    assertTrue(refused.getMessage().contains("old_nodes"), refused.getMessage());
    assertEquals(0, layer.size());
  }

  @Test
  void aNamedTableIsReadAlone(@TempDir Path dir) throws Exception {
    URL url = twoTableFile(dir).toURI().toURL();
    VectorLayer layer = new VectorLayer();

    VectorLayer.readGPKG(url, layer, "nodes");

    assertEquals(2, layer.size());
  }

  @Test
  void aMissingTableNameIsRefused(@TempDir Path dir) throws Exception {
    URL url = nodesFile(dir).toURI().toURL();

    assertThrows(
        IllegalArgumentException.class,
        () -> VectorLayer.readGPKG(url, new VectorLayer(), "edges"));
  }
}
