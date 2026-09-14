package sim.graph;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import sim.field.geo.VectorLayer;
import sim.testing.Fixtures;
import sim.util.geo.MasonGeometry;

/**
 * What {@code fromStreetJunctionsSegments} does with the junction layer it is handed: every node
 * takes its geometry from that layer rather than from a bare generated point.
 */
class GraphJunctionsTest {

  @Test
  @DisplayName("each node takes the supplied junction geometry at its coordinate")
  void nodesAdoptTheSuppliedJunctionGeometries() {
    List<MasonGeometry> segments = Fixtures.gridSegments(3, 3, 100.0);
    VectorLayer junctions = Fixtures.junctionLayer(segments);
    Graph graph = Fixtures.graphOf(segments, junctions);

    assertEquals(9, graph.getNodes().size());
    for (NodeGraph node : graph.getNodes()) {
      MasonGeometry geometry = node.getMasonGeometry();
      assertTrue(junctions.geometriesView().contains(geometry));
      assertNotNull(geometry.getIntegerAttribute("nodeID"), "attributes have to survive");
    }
  }

  @Test
  @DisplayName("the junctions layer the lookups query is the adopted one, one geometry per node")
  void theJunctionsLayerFollowsTheNodes() {
    List<MasonGeometry> segments = Fixtures.gridSegments(3, 3, 100.0);
    Graph graph = Fixtures.graphOf(segments, Fixtures.junctionLayer(segments));

    assertEquals(graph.getNodes().size(), graph.junctions.size());
    for (MasonGeometry junction : graph.junctions.geometriesView()) {
      NodeGraph node = graph.findNode(junction.getGeometry().getCoordinate());
      assertNotNull(node);
      assertSame(junction, node.getMasonGeometry());
      assertNotNull(junction.getIntegerAttribute("district"));
    }
  }

  @Test
  @DisplayName("a junction that matches no node is skipped")
  void unmatchedJunctionsAreSkipped() {
    List<MasonGeometry> segments = Fixtures.gridSegments(3, 3, 100.0);
    VectorLayer junctions = Fixtures.junctionLayer(segments);
    MasonGeometry stray = Fixtures.point(9999.0, 9999.0);
    stray.addIntegerAttribute("nodeID", 999);
    junctions.addGeometry(stray);

    Graph graph = Fixtures.graphOf(segments, junctions);

    assertEquals(9, graph.getNodes().size());
    assertEquals(9, graph.junctions.size());
    assertTrue(!graph.junctions.geometriesView().contains(stray));
    assertNull(graph.findNode(stray.getGeometry().getCoordinate()));
  }

  @Test
  @DisplayName("an empty junction layer leaves the generated points in place")
  void anEmptyJunctionLayerChangesNothing() {
    Graph graph = Fixtures.grid(3, 3, 100.0);

    assertEquals(9, graph.getNodes().size());
    assertEquals(9, graph.junctions.size());
    for (NodeGraph node : graph.getNodes()) {
      assertNotNull(node.getMasonGeometry());
      assertTrue(node.getMasonGeometry().getAttributes().isEmpty());
    }
  }

  @Test
  @DisplayName("the spatial index is built on the adopted geometries")
  void theSpatialIndexUsesTheAdoptedGeometries() {
    List<MasonGeometry> segments = Fixtures.gridSegments(3, 3, 100.0);
    Graph graph = Fixtures.graphOf(segments, Fixtures.junctionLayer(segments));

    List<NodeGraph> inside = graph.getNodesWithinPolygon((org.locationtech.jts.geom.Polygon) Fixtures.FACTORY
        .toGeometry(new org.locationtech.jts.geom.Envelope(-10.0, 110.0, -10.0, 110.0)));

    assertEquals(4, inside.size());
    assertEquals(inside.size(), graph.getContainedNodes(Fixtures.FACTORY
        .toGeometry(new org.locationtech.jts.geom.Envelope(-10.0, 110.0, -10.0, 110.0))).size());
  }

  @Test
  @DisplayName("distance lookups still resolve nodes through the adopted layer")
  void distanceLookupsStillResolveNodes() {
    List<MasonGeometry> segments = Fixtures.gridSegments(5, 5, 100.0);
    Graph graph = Fixtures.graphOf(segments, Fixtures.junctionLayer(segments));
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);

    List<NodeGraph> candidates =
        NodesLookup.getNodesBetweenDistanceInterval(graph, origin, 0.0, 150.0);

    assertEquals(3, candidates.size());
    for (NodeGraph candidate : candidates) {
      assertNotNull(candidate, "a junction with no node would leak a null into the candidates");
    }
  }
}
