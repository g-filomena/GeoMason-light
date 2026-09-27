package sim.io.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.net.URL;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import mil.nga.geopackage.GeoPackage;
import mil.nga.geopackage.GeoPackageManager;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
import sim.field.geo.VectorLayer;
import sim.testing.Fixtures;
import sim.util.geo.MasonGeometry;

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

  @Test
  @DisplayName("the geometry column is found by its declared name, whatever it is")
  void aGeometryColumnNamedGeomIsNotAnAttribute(@TempDir Path dir) throws Exception {
    File file = nodesFile(dir);
    try (GeoPackage geoPackage = GeoPackageManager.open(file)) {
      geoPackage.execSQL("ALTER TABLE nodes RENAME COLUMN geometry TO geom");
      geoPackage.execSQL(
          "UPDATE gpkg_geometry_columns SET column_name = 'geom' WHERE table_name = 'nodes'");
    }
    VectorLayer layer = new VectorLayer();

    VectorLayer.readGPKG(file.toURI().toURL(), layer);

    assertEquals(2, layer.size());
    assertFalse(layer.getGeometries().get(0).hasAttribute("geom"));
  }

  @Test
  @DisplayName("a feature without a geometry is skipped, not a NullPointerException")
  void aFeatureWithoutGeometryIsSkipped(@TempDir Path dir) throws Exception {
    File file = nodesFile(dir);
    try (GeoPackage geoPackage = GeoPackageManager.open(file)) {
      geoPackage.execSQL("INSERT INTO nodes (geometry) VALUES (NULL)");
    }
    VectorLayer layer = new VectorLayer();

    VectorLayer.readGPKG(file.toURI().toURL(), layer);

    assertEquals(2, layer.size());
  }

  @Test
  @DisplayName("a polygon keeps its holes")
  void polygonHolesSurviveARoundTrip(@TempDir Path dir) throws Exception {
    File file = dir.resolve("blocks.gpkg").toFile();
    LinearRing shell = (LinearRing) ((Polygon) Fixtures.FACTORY.toGeometry(
        new Envelope(0, 10, 0, 10))).getExteriorRing();
    LinearRing hole = (LinearRing) ((Polygon) Fixtures.FACTORY.toGeometry(
        new Envelope(2, 4, 2, 4))).getExteriorRing();
    Polygon block = Fixtures.FACTORY.createPolygon(shell, new LinearRing[] {hole});
    VectorLayer.writeGPKG(file.getPath(),
        new VectorLayer(new ArrayList<>(Arrays.asList(new MasonGeometry(block)))));
    VectorLayer layer = new VectorLayer();

    VectorLayer.readGPKG(file.toURI().toURL(), layer);

    Polygon read = (Polygon) layer.getGeometries().get(0).getGeometry();
    assertEquals(1, read.getNumInteriorRing());
    assertEquals(96.0, read.getArea(), 1e-9);
  }

  @Test
  @DisplayName("text stays text, and a column mixing integers and decimals keeps the decimals")
  void attributeTypesSurviveARoundTrip(@TempDir Path dir) throws Exception {
    File file = dir.resolve("typed.gpkg").toFile();
    MasonGeometry first = Fixtures.point(0, 0);
    first.addStringAttribute("code", "007");
    first.addStringAttribute("ref", "12345678901");
    first.addIntegerAttribute("value", 1);
    MasonGeometry second = Fixtures.point(10, 0);
    second.addStringAttribute("code", "008");
    second.addStringAttribute("ref", "12345678902");
    second.addDoubleAttribute("value", 2.5);
    VectorLayer.writeGPKG(file.getPath(),
        new VectorLayer(new ArrayList<>(Arrays.asList(first, second))));
    VectorLayer layer = new VectorLayer();

    VectorLayer.readGPKG(file.toURI().toURL(), layer);

    List<MasonGeometry> read = layer.getGeometries();
    assertEquals("007", read.get(0).getStringAttribute("code"));
    assertEquals("12345678901", read.get(0).getStringAttribute("ref"));
    assertEquals(1.0, read.get(0).getDoubleAttribute("value"), 0.0);
    assertEquals(2.5, read.get(1).getDoubleAttribute("value"), 0.0);
  }
}
