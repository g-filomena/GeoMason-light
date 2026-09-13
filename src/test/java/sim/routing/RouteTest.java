package sim.routing;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.planargraph.DirectedEdge;
import sim.graph.EdgeGraph;
import sim.graph.Graph;
import sim.graph.NodeGraph;
import sim.testing.Fixtures;

class RouteTest {

  private Graph graph;
  private Route route;

  @BeforeEach
  void routeAlongAChain() {
    graph = Fixtures.path(5, 100.0);
    route = new Astar().astarRoute(Fixtures.nodeAt(graph, 0.0, 0.0),
        Fixtures.nodeAt(graph, 400.0, 0.0), graph, null);
    assertNotNull(route);
  }

  @Test
  @DisplayName("getLength() is the sum of the edges the route traverses")
  void lengthSumsTheEdgesTraversed() {
    double expected = 0.0;
    for (EdgeGraph edge : route.edgesSequence) {
      expected += edge.getLength();
    }
    assertEquals(expected, route.getLength(), 1e-9);
    assertEquals(400.0, route.getLength(), 1e-9);
  }

  @Test
  @DisplayName("resetRoute() recomputes the length rather than keeping the planned one")
  void resetRouteRecomputesTheLength() {
    List<DirectedEdge> firstTwo =
        new ArrayList<>(route.directedEdgesSequence.subList(0, 2));

    route.resetRoute(firstTwo);

    assertEquals(2, route.edgesSequence.size());
    assertEquals(3, route.nodesSequence.size());
    assertEquals(200.0, route.getLength(), 1e-9);
    assertSame(Fixtures.nodeAt(graph, 200.0, 0.0), route.destinationNode);
  }

  @Test
  void computeRouteSequencesSetsOriginAndDestination() {
    assertSame(Fixtures.nodeAt(graph, 0.0, 0.0), route.originNode);
    assertSame(Fixtures.nodeAt(graph, 400.0, 0.0), route.destinationNode);
    assertSame(route.nodesSequence.get(0), route.originNode);
    assertSame(route.nodesSequence.get(route.nodesSequence.size() - 1), route.destinationNode);
  }

  @Test
  @DisplayName("the route geometry runs from the origin to the destination")
  void lineStringRunsFromOriginToDestination() {
    LineString lineString = route.getLineString();

    assertNotNull(lineString);
    assertEquals(400.0, lineString.getLength(), 1e-9);

    Coordinate[] coordinates = lineString.getCoordinates();
    assertEquals(route.originNode.getCoordinate(), coordinates[0]);
    assertEquals(route.destinationNode.getCoordinate(), coordinates[coordinates.length - 1]);
  }

  @Test
  @DisplayName("the geometry stays continuous when segments are stored back to front")
  void lineStringHandlesReversedSegments() {
    // The chain is walked from x=400 down to x=0, so every segment is traversed against the
    // direction its LineString was digitised in.
    Route backwards = new Astar().astarRoute(Fixtures.nodeAt(graph, 400.0, 0.0),
        Fixtures.nodeAt(graph, 0.0, 0.0), graph, null);

    Coordinate[] coordinates = backwards.getLineString().getCoordinates();
    assertEquals(400.0, coordinates[0].x, 1e-9);
    assertEquals(0.0, coordinates[coordinates.length - 1].x, 1e-9);
    for (int i = 1; i < coordinates.length; i++) {
      assertTrue(coordinates[i].x <= coordinates[i - 1].x, "geometry doubles back at " + i);
    }
  }

  @Test
  void anEmptyRouteHasNoLengthAndNoSequences() {
    Route empty = new Route();
    assertEquals(0.0, empty.getLength(), 0.0);
    assertTrue(empty.edgesSequence.isEmpty());
    assertTrue(empty.nodesSequence.isEmpty());
    assertTrue(empty.getVisitedLocations().isEmpty());
  }

  @Test
  void theDirectedEdgeConstructorStoresTheSequence() {
    Route constructed = new Route(route.directedEdgesSequence);
    assertEquals(route.directedEdgesSequence, constructed.directedEdgesSequence);
    // Nothing is derived until the sequences are computed.
    assertTrue(constructed.edgesSequence.isEmpty());
  }

  @Test
  void visitedLocationsAndAttributesAreCarriedAlong() {
    Set<NodeGraph> visited =
        new HashSet<>(Arrays.asList(Fixtures.nodeAt(graph, 100.0, 0.0)));
    route.setVisitedLocations(visited);
    route.attributes.put("purpose", "work");
    route.social = true;

    assertEquals(visited, route.getVisitedLocations());
    assertEquals("work", route.attributes.get("purpose"));
    assertTrue(route.social);
  }
}
