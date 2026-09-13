package sim.graph;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import sim.testing.Fixtures;
import sim.util.geo.GeometryUtilities;

class GraphUtilsTest {

  private Graph graph;

  @BeforeEach
  void buildGrid() {
    graph = Fixtures.grid(3, 3, 100.0);
  }

  @Test
  @DisplayName("nodesDistance() is symmetric and matches the plain Euclidean distance")
  void nodesDistanceIsSymmetric() {
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph destination = Fixtures.nodeAt(graph, 200.0, 100.0);

    double expected = GeometryUtilities.euclideanDistance(origin.getCoordinate(),
        destination.getCoordinate());
    assertEquals(expected, GraphUtils.nodesDistance(origin, destination), 1e-9);
    // Second call comes off the cache; the answer must not change.
    assertEquals(expected, GraphUtils.nodesDistance(origin, destination), 1e-9);
    assertEquals(expected, GraphUtils.nodesDistance(destination, origin), 1e-9);
  }

  @Test
  void lineStringBetweenNodesJoinsTheTwoCoordinates() {
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph destination = Fixtures.nodeAt(graph, 0.0, 200.0);

    assertEquals(200.0, GraphUtils.LineStringBetweenNodes(origin, destination).getLength(), 1e-9);
  }

  @Test
  void nodesAndEdgesCrossReferenceEachOther() {
    NodeGraph corner = Fixtures.nodeAt(graph, 0.0, 0.0);
    Set<EdgeGraph> edges = new HashSet<>(corner.getEdges());

    Set<NodeGraph> nodes = GraphUtils.nodesFromEdges(edges);
    assertEquals(3, nodes.size());
    assertTrue(nodes.contains(corner));
    assertTrue(GraphUtils.edgesFromNodes(new HashSet<>(Arrays.asList(corner))).containsAll(edges));
  }

  @Test
  void idHelpersReadIdsOffListsAndSets() {
    List<NodeGraph> nodes = graph.getNodes();
    List<EdgeGraph> edges = graph.getEdges();

    assertEquals(nodes.size(), GraphUtils.getNodeIDs(nodes).size());
    assertEquals(nodes.size(), GraphUtils.getNodeIDs(new HashSet<>(nodes)).size());
    assertEquals(edges.size(), GraphUtils.getEdgeIDs(edges).size());
    assertEquals(edges.size(), GraphUtils.getEdgeIDs(new HashSet<>(edges)).size());
  }

  @Test
  @DisplayName("id lookups drop ids the map does not know")
  void idLookupsSkipUnknownIds() {
    NodeGraph node = graph.getNodes().get(0);
    Map<Integer, NodeGraph> nodeMap = new HashMap<>();
    nodeMap.put(node.getID(), node);

    List<NodeGraph> found =
        GraphUtils.getNodesFromNodeIDs(Arrays.asList(node.getID(), 9999), nodeMap);
    assertEquals(1, found.size());
    assertSame(node, found.get(0));

    EdgeGraph edge = graph.getEdges().get(0);
    Map<Integer, EdgeGraph> edgeMap = new HashMap<>();
    edgeMap.put(edge.getID(), edge);
    assertEquals(1,
        GraphUtils.getEdgesFromEdgeIDs(new HashSet<>(Arrays.asList(edge.getID(), 9999)), edgeMap)
            .size());
  }

  @Test
  void findClosestNodePicksTheNearestJunction() {
    NodeGraph closest = GraphUtils.findClosestNode(new Coordinate(95.0, 5.0), graph.getNodes());
    assertSame(Fixtures.nodeAt(graph, 100.0, 0.0), closest);
  }

  @Test
  @DisplayName("smallestEnclosingGeometryBetweenNodes() adapts to the number of nodes")
  void smallestEnclosingGeometryAdaptsToTheInput() {
    NodeGraph one = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph two = Fixtures.nodeAt(graph, 200.0, 0.0);
    NodeGraph three = Fixtures.nodeAt(graph, 200.0, 200.0);
    NodeGraph four = Fixtures.nodeAt(graph, 0.0, 200.0);

    Geometry single = GraphUtils.smallestEnclosingGeometryBetweenNodes(Arrays.asList(one));
    assertEquals(Math.PI * 50.0 * 50.0, single.getArea(), Math.PI * 50.0 * 50.0 * 0.01);

    // Two nodes: the circle centred on their midpoint whose diameter is the distance between
    // them. Its boundary passes through both, so measure the area rather than containment.
    Geometry pair = GraphUtils.smallestEnclosingGeometryBetweenNodes(Arrays.asList(one, two));
    assertEquals(Math.PI * 100.0 * 100.0, pair.getArea(), Math.PI * 100.0 * 100.0 * 0.01);
    assertTrue(pair.covers(Fixtures.FACTORY.createPoint(new Coordinate(100.0, 0.0))));

    Geometry hull = GraphUtils.smallestEnclosingGeometryBetweenNodes(
        new java.util.ArrayList<>(Arrays.asList(one, two, three, four)));
    assertEquals(200.0 * 200.0, hull.getArea(), 1e-6);
  }
}
