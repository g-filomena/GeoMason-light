package sim.field.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import sim.testing.Fixtures;
import sim.util.geo.MasonGeometry;

class VectorLayerTest {

  private VectorLayer layer;
  private MasonGeometry west;
  private MasonGeometry centre;
  private MasonGeometry east;

  private static MasonGeometry junction(double x, double y, int id, String type) {
    MasonGeometry geometry = Fixtures.point(x, y);
    geometry.addIntegerAttribute("nodeID", id);
    geometry.addStringAttribute("type", type);
    return geometry;
  }

  private static Polygon square(double minX, double minY, double size) {
    return (Polygon) Fixtures.FACTORY
        .toGeometry(new Envelope(minX, minX + size, minY, minY + size));
  }

  @BeforeEach
  void buildLayer() {
    west = junction(0.0, 0.0, 1, "primary");
    centre = junction(100.0, 0.0, 2, "secondary");
    east = junction(200.0, 0.0, 3, "secondary");
    layer = new VectorLayer(new ArrayList<>(Arrays.asList(west, centre, east)));
  }

  @Test
  void sizeAndEmptinessTrackTheGeometries() {
    assertEquals(3, layer.size());
    assertFalse(layer.isEmpty());
    assertTrue(layer.isPopulated());

    VectorLayer blank = new VectorLayer();
    assertEquals(0, blank.size());
    assertTrue(blank.isEmpty());
    assertFalse(blank.isPopulated());
  }

  @Test
  void minimumBoundingRectangleSpansTheGeometries() {
    assertEquals(0.0, layer.getMBR().getMinX(), 1e-9);
    assertEquals(200.0, layer.getMBR().getMaxX(), 1e-9);
    assertEquals(200.0, layer.getWidth(), 1e-9);
  }

  @Test
  @DisplayName("getGeometries() copies, geometriesView() does not")
  void geometriesViewIsAnUnmodifiableLiveView() {
    List<MasonGeometry> copy = layer.getGeometries();
    copy.clear();
    assertEquals(3, layer.size());

    List<MasonGeometry> view = layer.geometriesView();
    assertEquals(3, view.size());
    assertThrows(UnsupportedOperationException.class, () -> view.clear());

    layer.addGeometry(junction(300.0, 0.0, 4, "tertiary"));
    assertEquals(4, view.size(), "the view has to follow the layer");
    assertNotSame(view, layer.getGeometries());
  }

  @Test
  @DisplayName("the id lookup resolves the ids set from an attribute")
  void geometriesFromIdsResolveThroughTheIdIndex() {
    layer.setID("nodeID");

    assertEquals(Arrays.asList(1, 2, 3), layer.getIDs());
    List<MasonGeometry> found =
        layer.getGeometriesFromIDs(new HashSet<>(Arrays.asList(1, 3, 99)));

    assertEquals(2, found.size());
    assertTrue(found.contains(west));
    assertTrue(found.contains(east));
  }

  @Test
  @DisplayName("the id lookup notices geometries added or removed after it was built")
  void idIndexIsInvalidatedByStructuralChanges() {
    layer.setID("nodeID");
    assertEquals(1, layer.getGeometriesFromIDs(Collections.singleton(2)).size());

    layer.removeGeometry(centre);
    assertTrue(layer.getGeometriesFromIDs(Collections.singleton(2)).isEmpty());

    MasonGeometry replacement = junction(100.0, 0.0, 2, "secondary");
    layer.addGeometry(replacement);
    layer.setID("nodeID");
    List<MasonGeometry> found = layer.getGeometriesFromIDs(Collections.singleton(2));
    assertEquals(1, found.size());
    assertSame(replacement, found.get(0));
  }

  @Test
  @DisplayName("a removed geometry stops turning up in spatial queries")
  void spatialIndexIsRebuiltAfterARemoval() {
    Envelope around = new Envelope(90.0, 110.0, -10.0, 10.0);
    assertEquals(1, layer.queryField(around).size());

    layer.removeGeometry(centre);

    // No explicit updateSpatialIndex() call: the query has to pick up the change itself.
    assertTrue(layer.queryField(around).isEmpty());
    assertEquals(2, layer.queryField(layer.getMBR()).size());
  }

  @Test
  void distanceQueriesRespectTheirLimits() {
    Geometry origin = west.getGeometry();

    assertEquals(Collections.singletonList(centre),
        layer.featuresBetweenLimits(origin, 50.0, 150.0));
    assertEquals(2, layer.featuresWithinDistance(origin, 150.0).size());
    assertEquals(3, layer.featuresWithinDistance(origin, 500.0).size());
  }

  @Test
  void intersectionQueriesFindTheOverlappingFeatures() {
    Polygon around = square(-10.0, -10.0, 120.0);

    assertTrue(layer.intersects(around));
    assertEquals(2, layer.intersectingFeatures(around).size());
    assertFalse(layer.intersects(square(1000.0, 1000.0, 10.0)));
  }

  @Test
  @DisplayName("containment queries read both ways round")
  void containmentQueriesReadBothWays() {
    Polygon around = square(-10.0, -10.0, 120.0);

    assertEquals(2, layer.containedFeatures(around).size());
    assertTrue(layer.coveringFeatures(west).isEmpty());

    VectorLayer blocks = new VectorLayer(new ArrayList<>(Arrays.asList(
        new MasonGeometry(square(-10.0, -10.0, 120.0)))));
    assertEquals(1, blocks.coveringFeatures(west).size());
    assertTrue(blocks.isCovered(west));
    assertFalse(blocks.isCovered(east));
  }

  @Test
  @DisplayName("filterFeatures() on a string attribute honours equal = false")
  void filterFeaturesOnAStringAttribute() {
    assertEquals(2, layer.filterFeatures("type", "secondary", true).size());
    assertEquals(Collections.singletonList(west),
        layer.filterFeatures("type", "secondary", false));
    assertEquals(3, layer.filterFeatures("type", "motorway", false).size());
  }

  @Test
  @DisplayName("filterFeatures() on an int attribute honours equal = false")
  void filterFeaturesOnAnIntAttribute() {
    assertEquals(Collections.singletonList(centre), layer.filterFeatures("nodeID", 2, true));
    assertEquals(2, layer.filterFeatures("nodeID", 2, false).size());
    assertFalse(layer.filterFeatures("nodeID", 2, false).contains(centre));
  }

  @Test
  void filterFeaturesOnAListOfValues() {
    assertEquals(3,
        layer.filterFeatures("type", Arrays.asList("primary", "secondary"), true).size());
    assertEquals(1, layer.filterFeatures("type", Arrays.asList("primary"), true).size());
    assertEquals(2, layer.filterFeatures("type", Arrays.asList("primary"), false).size());
  }

  @Test
  void selectFeaturesReturnsANewLayer() {
    VectorLayer selected = layer.selectFeatures("nodeID", Arrays.asList(1, 3), true);

    assertEquals(2, selected.size());
    assertEquals(3, layer.size(), "the source layer must not be touched");
    assertEquals(1, layer.selectFeatures("nodeID", Arrays.asList(1, 3), false).size());
  }

  @Test
  @DisplayName("getGeometry() looks a feature up by attribute value")
  void getGeometryFindsAFeatureByAttributeValue() {
    assertSame(centre, layer.getGeometry("nodeID", 2));
    assertSame(west, layer.getGeometry("type", "primary"));
    assertNull(layer.getGeometry("type", "motorway"));
    assertNull(layer.getGeometry("missing", "primary"));
  }

  @Test
  void intColumnReadsAnAttributeAcrossTheLayer() {
    assertEquals(Arrays.asList(1, 2, 3), layer.getIntColumn("nodeID"));
  }

  @Test
  @DisplayName("intersection() reports the features of this layer that meet the other one")
  void intersectionReportsTheOverlappingFeaturesOfThisLayer() {
    VectorLayer blocks = new VectorLayer(new ArrayList<>(Arrays.asList(
        new MasonGeometry(square(-10.0, -10.0, 120.0)))));

    List<MasonGeometry> overlapping = layer.intersection(blocks, true);

    assertEquals(2, overlapping.size());
    assertTrue(overlapping.contains(west));
    assertTrue(overlapping.contains(centre));
    assertFalse(overlapping.contains(east));
  }

  @Test
  @DisplayName("intersection() with inclusive = false is the exact complement")
  void nonInclusiveIntersectionIsTheComplement() {
    VectorLayer blocks = new VectorLayer(new ArrayList<>(Arrays.asList(
        new MasonGeometry(square(-10.0, -10.0, 120.0)))));

    List<MasonGeometry> overlapping = layer.intersection(blocks, true);
    List<MasonGeometry> rest = layer.intersection(blocks, false);

    assertEquals(Collections.singletonList(east), rest);
    assertEquals(layer.size(), overlapping.size() + rest.size());
    for (MasonGeometry feature : rest) {
      assertFalse(overlapping.contains(feature));
    }
  }

  @Test
  @DisplayName("a feature meeting several geometries of the other layer is reported once")
  void intersectionDoesNotRepeatFeatures() {
    VectorLayer blocks = new VectorLayer(new ArrayList<>(Arrays.asList(
        new MasonGeometry(square(-10.0, -10.0, 120.0)),
        new MasonGeometry(square(-20.0, -20.0, 140.0)))));

    List<MasonGeometry> overlapping = layer.intersection(blocks, true);

    assertEquals(2, overlapping.size());
    assertEquals(2, new HashSet<>(overlapping).size());
  }

  @Test
  @DisplayName("intersection() never removes features from the layer it was handed")
  void intersectionDoesNotMutateTheOtherLayer() {
    VectorLayer blocks = new VectorLayer(new ArrayList<>(Arrays.asList(
        new MasonGeometry(square(-10.0, -10.0, 120.0)))));

    layer.intersection(blocks, false);
    blocks.intersection(layer, false);

    assertEquals(1, blocks.size());
    assertEquals(3, layer.size());
  }

  @Test
  void convexHullEnclosesEveryGeometry() {
    VectorLayer corners = new VectorLayer(new ArrayList<>(Arrays.asList(
        Fixtures.point(0.0, 0.0), Fixtures.point(100.0, 0.0), Fixtures.point(100.0, 100.0),
        Fixtures.point(0.0, 100.0))));

    assertEquals(100.0 * 100.0, corners.getConvexHull().getArea(), 1e-6);
    assertTrue(corners.isInsideConvexHull(new Coordinate(50.0, 50.0)));
    assertFalse(corners.isInsideConvexHull(new Coordinate(500.0, 500.0)));
  }

  @Test
  void unionMergesAdjacentPolygons() {
    VectorLayer blocks = new VectorLayer(new ArrayList<>(Arrays.asList(
        new MasonGeometry(square(0.0, 0.0, 10.0)),
        new MasonGeometry(square(10.0, 0.0, 10.0)))));

    assertEquals(200.0, blocks.getUnion().getArea(), 1e-6);
    assertTrue(blocks.isInsideUnion(new Coordinate(15.0, 5.0)));
    assertFalse(blocks.isInsideUnion(new Coordinate(25.0, 5.0)));
  }

  @Test
  void clearEmptiesTheLayerAndItsIndexes() {
    layer.setID("nodeID");
    layer.clear();

    assertTrue(layer.isEmpty());
    assertTrue(layer.queryField(new Envelope(-1000.0, 1000.0, -1000.0, 1000.0)).isEmpty());
    assertTrue(layer.getGeometriesFromIDs(Collections.singleton(1)).isEmpty());
  }

  @Test
  void findGeometryReturnsTheStoredInstance() {
    assertSame(centre, layer.findGeometry(junction(100.0, 0.0, 2, "secondary")));
  }

  @Test
  void geometryLocationIsTheCentroid() {
    assertEquals(0.0, layer.getGeometryLocation(west).getX(), 1e-9);
  }
}
