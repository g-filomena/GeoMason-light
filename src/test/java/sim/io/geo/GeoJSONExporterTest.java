package sim.io.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
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

class GeoJSONExporterTest {

  private static VectorLayer layerOf(MasonGeometry... geometries) {
    return new VectorLayer(new ArrayList<>(Arrays.asList(geometries)));
  }

  private static LinearRing ring(double minX, double minY, double size) {
    return Fixtures.FACTORY.createLinearRing(new Coordinate[] {new Coordinate(minX, minY),
        new Coordinate(minX + size, minY), new Coordinate(minX + size, minY + size),
        new Coordinate(minX, minY + size), new Coordinate(minX, minY)});
  }

  @Test
  void anEmptyLayerIsAnEmptyFeatureCollection() {
    assertEquals("{\"type\":\"FeatureCollection\",\"features\":[]}",
        GeoJSONExporter.toFeatureCollection(new VectorLayer()));
  }

  @Test
  @DisplayName("a point feature carries its coordinates and its properties")
  void pointFeatureIsSerialisedInFull() {
    MasonGeometry point = Fixtures.point(2.5, -3.0);
    point.addStringAttribute("name", "market");

    assertEquals("{\"type\":\"FeatureCollection\",\"features\":["
        + "{\"type\":\"Feature\",\"geometry\":{\"type\":\"Point\",\"coordinates\":[2.5,-3]},"
        + "\"properties\":{\"name\":\"market\"}}]}",
        GeoJSONExporter.toFeatureCollection(layerOf(point)));
  }

  @Test
  @DisplayName("numbers are written plainly, with no trailing zeros and no exponent")
  void numbersAreWrittenPlainly() {
    String json = GeoJSONExporter.toFeatureCollection(layerOf(Fixtures.point(10000000.0, 0.5)));
    assertTrue(json.contains("[10000000,0.5]"), json);
  }

  @Test
  @DisplayName("attribute values keep their JSON type")
  void attributeValuesKeepTheirType() {
    MasonGeometry point = Fixtures.point(0.0, 0.0);
    point.addIntegerAttribute("count", 3);
    String json = GeoJSONExporter.toFeatureCollection(layerOf(point));
    assertTrue(json.contains("\"count\":3"), json);

    MasonGeometry flagged = Fixtures.point(0.0, 0.0);
    flagged.addAttribute("pedestrianised", Boolean.TRUE);
    assertTrue(GeoJSONExporter.toFeatureCollection(layerOf(flagged))
        .contains("\"pedestrianised\":true"));

    MasonGeometry broken = Fixtures.point(0.0, 0.0);
    broken.addDoubleAttribute("slope", Double.NaN);
    assertTrue(GeoJSONExporter.toFeatureCollection(layerOf(broken)).contains("\"slope\":null"));
  }

  @Test
  @DisplayName("attribute names and string values are escaped")
  void stringsAreEscaped() {
    MasonGeometry point = Fixtures.point(0.0, 0.0);
    point.addStringAttribute("label", "a \"quoted\" name\nwith a break\tand a \\ slash");

    String json = GeoJSONExporter.toFeatureCollection(layerOf(point));
    assertTrue(json.contains("a \\\"quoted\\\" name\\nwith a break\\tand a \\\\ slash"), json);
    assertFalse(json.contains("\n"), "raw newlines would break the document");
  }

  @Test
  @DisplayName("control characters become unicode escapes")
  void controlCharactersAreEscaped() {
    MasonGeometry point = Fixtures.point(0.0, 0.0);
    point.addStringAttribute("label", "a" + ((char) 1) + "b");
    assertTrue(GeoJSONExporter.toFeatureCollection(layerOf(point)).contains("a\\u0001b"));
  }

  @Test
  @DisplayName("properties can be left out for geometry-only consumers")
  void propertiesCanBeOmitted() {
    MasonGeometry point = Fixtures.point(1.0, 2.0);
    point.addStringAttribute("name", "market");

    String json = GeoJSONExporter.toFeatureCollection(layerOf(point), false);
    assertTrue(json.contains("\"properties\":{}"), json);
    assertFalse(json.contains("market"));
  }

  @Test
  void lineStringsAndPolygonsKeepTheirType() {
    assertTrue(GeoJSONExporter.toFeatureCollection(layerOf(Fixtures.segment(0.0, 0.0, 1.0, 1.0)))
        .contains("\"type\":\"LineString\",\"coordinates\":[[0,0],[1,1]]"));

    Polygon block = (Polygon) Fixtures.FACTORY.toGeometry(new Envelope(0.0, 1.0, 0.0, 1.0));
    assertTrue(GeoJSONExporter.toFeatureCollection(layerOf(new MasonGeometry(block)))
        .contains("\"type\":\"Polygon\""));
  }

  @Test
  @DisplayName("a polygon with a hole writes the hole as a second ring")
  void polygonHolesAreWrittenAsInnerRings() {
    Polygon withHole = Fixtures.FACTORY.createPolygon(ring(0.0, 0.0, 10.0),
        new LinearRing[] {ring(2.0, 2.0, 2.0)});

    String json = GeoJSONExporter.toFeatureCollection(layerOf(new MasonGeometry(withHole)));
    assertTrue(json.contains("\"type\":\"Polygon\""), json);
    assertTrue(json.contains("[[2,2],[4,2],[4,4],[2,4],[2,2]]"), json);
  }

  @Test
  void multiPartGeometriesKeepTheirType() {
    assertTrue(GeoJSONExporter
        .toFeatureCollection(layerOf(new MasonGeometry(Fixtures.FACTORY
            .createMultiPointFromCoords(
                new Coordinate[] {new Coordinate(0.0, 0.0), new Coordinate(1.0, 1.0)}))))
        .contains("\"type\":\"MultiPoint\""));
  }

  @Test
  @DisplayName("every feature of the layer is written, in order")
  void everyFeatureIsWritten() {
    String json = GeoJSONExporter.toFeatureCollection(
        layerOf(Fixtures.point(0.0, 0.0), Fixtures.point(1.0, 0.0), Fixtures.point(2.0, 0.0)));

    assertEquals(3, json.split("\"type\":\"Feature\"", -1).length - 1);
    assertTrue(json.indexOf("[0,0]") < json.indexOf("[1,0]"));
    assertTrue(json.indexOf("[1,0]") < json.indexOf("[2,0]"));
  }

  @Test
  @DisplayName("a properties provider replaces the feature attributes")
  void propertiesProviderReplacesTheAttributes() {
    MasonGeometry road = Fixtures.segment(0.0, 0.0, 1.0, 0.0);
    road.addIntegerAttribute("edgeID", 7);
    road.addStringAttribute("name", "not exported");

    String json = GeoJSONExporter.toFeatureCollection(layerOf(road), feature -> {
      Map<String, Object> properties = new LinkedHashMap<>();
      properties.put("edgeID", feature.getIntegerAttribute("edgeID"));
      properties.put("volume", 42);
      return properties;
    });

    assertTrue(json.contains("\"properties\":{\"edgeID\":7,\"volume\":42}"), json);
    assertFalse(json.contains("not exported"));
  }

  @Test
  @DisplayName("a provider returning null or an empty map yields empty properties")
  void providerMayDeclineToSupplyProperties() {
    MasonGeometry point = Fixtures.point(1.0, 2.0);
    point.addStringAttribute("name", "market");

    assertTrue(GeoJSONExporter.toFeatureCollection(layerOf(point), feature -> null)
        .contains("\"properties\":{}"));
    assertTrue(
        GeoJSONExporter.toFeatureCollection(layerOf(point), feature -> new LinkedHashMap<>())
            .contains("\"properties\":{}"));
  }

  @Test
  @DisplayName("provided values are typed and escaped like attributes")
  void providedValuesAreTypedAndEscaped() {
    MasonGeometry point = Fixtures.point(0.0, 0.0);

    String json = GeoJSONExporter.toFeatureCollection(layerOf(point), feature -> {
      Map<String, Object> properties = new LinkedHashMap<>();
      properties.put("lux", 12.50);
      properties.put("lit", Boolean.FALSE);
      properties.put("missing", null);
      properties.put("label", "a \"quoted\" street");
      return properties;
    });

    assertTrue(json.contains("\"lux\":12.5"), json);
    assertTrue(json.contains("\"lit\":false"), json);
    assertTrue(json.contains("\"missing\":null"), json);
    assertTrue(json.contains("a \\\"quoted\\\" street"), json);
  }

  @Test
  @DisplayName("the provider form writes the same geometries as the attribute form")
  void providerFormWritesTheSameGeometries() {
    VectorLayer layer = layerOf(Fixtures.segment(0.0, 0.0, 1.0, 1.0), Fixtures.point(2.0, 3.0));

    String withAttributes = GeoJSONExporter.toFeatureCollection(layer, false);
    String withProvider = GeoJSONExporter.toFeatureCollection(layer, feature -> null);

    assertEquals(withAttributes, withProvider);
  }

  @Test
  @DisplayName("write() accepts a properties provider too")
  void writeAcceptsAProvider(@TempDir Path directory) throws IOException {
    MasonGeometry road = Fixtures.segment(0.0, 0.0, 1.0, 0.0);
    Path target = directory.resolve("volumes");

    GeoJSONExporter.write(target.toString(), layerOf(road), feature -> {
      Map<String, Object> properties = new LinkedHashMap<>();
      properties.put("volume", 3);
      return properties;
    });

    Path written = directory.resolve("volumes.geojson");
    assertTrue(Files.exists(written));
    assertTrue(new String(Files.readAllBytes(written), StandardCharsets.UTF_8)
        .contains("\"volume\":3"));
  }

  @Test
  @DisplayName("write() adds the extension when it is missing and keeps it otherwise")
  void writeNamesTheFileSensibly(@TempDir Path directory) throws IOException {
    VectorLayer layer = layerOf(Fixtures.point(1.0, 2.0));

    Path implicit = directory.resolve("junctions");
    GeoJSONExporter.write(implicit.toString(), layer);
    Path written = directory.resolve("junctions.geojson");
    assertTrue(Files.exists(written));
    assertEquals(GeoJSONExporter.toFeatureCollection(layer),
        new String(Files.readAllBytes(written), StandardCharsets.UTF_8));

    Path explicit = directory.resolve("segments.geojson");
    GeoJSONExporter.write(explicit.toString(), layer);
    assertTrue(Files.exists(explicit));
    assertFalse(Files.exists(directory.resolve("segments.geojson.geojson")));

    Path asJson = directory.resolve("nodes.json");
    GeoJSONExporter.write(asJson.toString(), layer);
    assertTrue(Files.exists(asJson));
  }

  @Test
  @DisplayName("VectorLayer.writeGeoJSON() writes the same document")
  void layerConvenienceMethodMatchesTheExporter(@TempDir Path directory) throws Exception {
    VectorLayer layer = layerOf(Fixtures.point(3.0, 4.0));
    Path target = directory.resolve("layer.geojson");

    VectorLayer.writeGeoJSON(target.toString(), layer);

    List<String> lines = Files.readAllLines(target, StandardCharsets.UTF_8);
    assertEquals(1, lines.size());
    assertEquals(GeoJSONExporter.toFeatureCollection(layer), lines.get(0));
  }
}
