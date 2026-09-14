package sim.graph;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTimeoutPreemptively;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ec.util.MersenneTwisterFast;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import sim.testing.Fixtures;

class NodesLookupTest {

  private static Graph districtGrid() {
    Graph graph = Fixtures.grid(15, 15, 100.0);
    for (NodeGraph node : graph.getNodes()) {
      node.dma = "live";
    }
    return graph;
  }

  @Test
  @DisplayName("two generators on the same seed draw the same nodes")
  void drawsAreReproducibleUnderASeededGenerator() {
    Graph graph = Fixtures.grid(5, 5, 100.0);
    MersenneTwisterFast first = new MersenneTwisterFast(42L);
    MersenneTwisterFast second = new MersenneTwisterFast(42L);

    for (int draw = 0; draw < 25; draw++) {
      assertSame(NodesLookup.randomNode(graph, first), NodesLookup.randomNode(graph, second));
    }
  }

  @Test
  @DisplayName("different seeds do not draw the same sequence")
  void differentSeedsDiverge() {
    Graph graph = Fixtures.grid(5, 5, 100.0);
    MersenneTwisterFast first = new MersenneTwisterFast(1L);
    MersenneTwisterFast second = new MersenneTwisterFast(2L);

    boolean diverged = false;
    for (int draw = 0; draw < 25 && !diverged; draw++) {
      diverged = NodesLookup.randomNode(graph, first) != NodesLookup.randomNode(graph, second);
    }
    assertTrue(diverged);
  }

  @Test
  @DisplayName("selectRandomNode() answers null rather than throwing on an empty list")
  void selectRandomNodeIsNullSafe() {
    MersenneTwisterFast random = new MersenneTwisterFast(3L);
    assertNull(NodesLookup.selectRandomNode(null, random));
    assertNull(NodesLookup.selectRandomNode(Collections.<NodeGraph>emptyList(), random));
    assertNull(NodesLookup.selectRandomNode(new ArrayList<NodeGraph>()));
  }

  @Test
  @DisplayName("selectRandomNode() draws only from the list it is given")
  void selectRandomNodeDrawsFromThatList() {
    Graph graph = Fixtures.grid(5, 5, 100.0);
    List<NodeGraph> candidates =
        Arrays.asList(Fixtures.nodeAt(graph, 0.0, 0.0), Fixtures.nodeAt(graph, 100.0, 0.0));
    MersenneTwisterFast random = new MersenneTwisterFast(5L);

    for (int draw = 0; draw < 20; draw++) {
      assertTrue(candidates.contains(NodesLookup.selectRandomNode(candidates, random)));
    }
  }

  @Test
  @DisplayName("getNodesBetweenDistanceInterval() excludes the origin node")
  void nodesBetweenDistanceIntervalExcludesTheOrigin() {
    Graph graph = Fixtures.grid(5, 5, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);

    List<NodeGraph> candidates =
        NodesLookup.getNodesBetweenDistanceInterval(graph, origin, 0.0, 150.0);

    assertFalse(candidates.contains(origin));
    assertTrue(candidates.contains(Fixtures.nodeAt(graph, 100.0, 0.0)));
    assertTrue(candidates.contains(Fixtures.nodeAt(graph, 0.0, 100.0)));
    for (NodeGraph candidate : candidates) {
      assertTrue(GraphUtils.nodesDistance(origin, candidate) <= 150.0);
    }
  }

  @Test
  void randomNodeBetweenDistanceIntervalStaysInsideTheInterval() {
    Graph graph = Fixtures.grid(15, 15, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    MersenneTwisterFast random = new MersenneTwisterFast(11L);

    for (int draw = 0; draw < 20; draw++) {
      NodeGraph node =
          NodesLookup.randomNodeBetweenDistanceInterval(graph, origin, 250.0, 350.0, random);
      assertNotNull(node);
      double distance = GraphUtils.nodesDistance(origin, node);
      assertTrue(distance >= 250.0 && distance <= 350.0, "drawn at " + distance);
    }
  }

  @Test
  @DisplayName("getNodesBetweenDistanceIntervalRegion() keeps only the other regions")
  void nodesBetweenDistanceIntervalRegionSkipsTheOwnRegion() {
    Graph graph = Fixtures.grid(5, 5, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    origin.setRegionID(1);
    Fixtures.nodeAt(graph, 100.0, 0.0).setRegionID(1);
    Fixtures.nodeAt(graph, 0.0, 100.0).setRegionID(2);

    List<NodeGraph> candidates =
        NodesLookup.getNodesBetweenDistanceIntervalRegion(graph, origin, 0.0, 150.0);

    assertFalse(candidates.contains(Fixtures.nodeAt(graph, 100.0, 0.0)));
    assertTrue(candidates.contains(Fixtures.nodeAt(graph, 0.0, 100.0)));
  }

  @Test
  @DisplayName("randomNodeFromDistancesSet() avoids nodes adjacent to the origin")
  void randomNodeFromDistancesSetAvoidsNeighbours() {
    Graph graph = Fixtures.grid(15, 15, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    MersenneTwisterFast random = new MersenneTwisterFast(13L);

    for (int draw = 0; draw < 20; draw++) {
      NodeGraph node = NodesLookup.randomNodeFromDistancesSet(graph, graph.junctions, origin,
          Arrays.asList(300.0f), random);
      assertNotNull(node);
      assertNull(graph.getEdgeBetween(origin, node));
    }
  }

  @Test
  @DisplayName("randomNodeFromDistancesSet() returns null for an empty distance set")
  void randomNodeFromDistancesSetHandlesNoDistances() {
    Graph graph = Fixtures.grid(5, 5, 100.0);
    assertNull(NodesLookup.randomNodeFromDistancesSet(graph, graph.junctions,
        Fixtures.nodeAt(graph, 0.0, 0.0), Collections.<Float>emptyList(),
        new MersenneTwisterFast(17L)));
  }

  @Test
  void candidatesByDmaFilterOnTheRequestedCategory() {
    Graph graph = Fixtures.grid(3, 3, 100.0);
    List<NodeGraph> nodes = graph.getNodes();
    nodes.get(0).dma = "live";
    nodes.get(1).dma = "work";
    nodes.get(2).dma = "visit";
    for (int i = 3; i < nodes.size(); i++) {
      nodes.get(i).dma = "none";
    }

    assertEquals(3, NodesLookup.getCandidatesByDMA(nodes, "random").size());
    assertEquals(2, NodesLookup.getCandidatesByDMA(nodes, "workOrVisit").size());
    assertEquals(1, NodesLookup.getCandidatesByDMA(nodes, "live").size());
    assertEquals(nodes.size() - 3, NodesLookup.getNodesByDMA(graph, "none").size());
  }

  @Test
  @DisplayName("the DMA interval search widens towards the origin as well as away from it")
  void dmaIntervalSearchWidensBothEnds() {
    // The only node of the requested category sits nearer than the lower limit. Widening the
    // upper end alone can never reach it, so the search would fall through to accepting any
    // category and hand back a node that was not the one asked for.
    Graph graph = districtGrid();
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph target = Fixtures.nodeAt(graph, 600.0, 0.0);
    target.dma = "work";

    NodeGraph drawn = NodesLookup.randomNodeBetweenDistanceIntervalDMA(graph, origin, 1000.0,
        1100.0, "work", new MersenneTwisterFast(19L));

    assertSame(target, drawn);
  }

  @Test
  @DisplayName("the DMA interval search gives up instead of widening forever")
  void dmaIntervalSearchTerminatesWhenNoCategoryMatches() {
    Graph graph = districtGrid();
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);

    NodeGraph drawn = assertTimeoutPreemptively(Duration.ofSeconds(20),
        () -> NodesLookup.randomNodeBetweenDistanceIntervalDMA(graph, origin, 200.0, 300.0,
            "nowhere", new MersenneTwisterFast(23L)));

    // No node carries that category, so once the interval has been widened past the cap any
    // candidate is accepted: the search returns something rather than spinning.
    assertNotNull(drawn);
  }

  @Test
  @DisplayName("randomNodeRegion() draws a node of the same region as the origin")
  void randomNodeRegionDrawsFromTheSameRegion() {
    // The region is a node property, read from the node itself: the junction geometries carry no
    // "district" attribute to read it from.
    Graph graph = Fixtures.grid(5, 5, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 200.0, 200.0);
    for (NodeGraph node : graph.getNodes()) {
      node.setRegionID(2);
    }
    origin.setRegionID(1);
    Fixtures.nodeAt(graph, 200.0, 100.0).setRegionID(1);
    MersenneTwisterFast random = new MersenneTwisterFast(31L);

    for (int draw = 0; draw < 10; draw++) {
      NodeGraph drawn = NodesLookup.randomNodeRegion(graph, origin, 150.0, random);
      assertSame(Fixtures.nodeAt(graph, 200.0, 100.0), drawn);
    }
  }

  @Test
  @DisplayName("randomNodeRegion() returns null when the region has nobody else nearby")
  void randomNodeRegionGivesUpWhenAlone() {
    Graph graph = Fixtures.grid(5, 5, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    for (NodeGraph node : graph.getNodes()) {
      node.setRegionID(2);
    }
    origin.setRegionID(1);

    assertNull(NodesLookup.randomNodeRegion(graph, origin, 150.0, new MersenneTwisterFast(37L)));
  }

  @Test
  @DisplayName("the salient lookup lowers the centrality bar until the interval answers")
  void salientLookupRetriesAtALowerPercentile() {
    Graph graph = Fixtures.grid(9, 9, 100.0);
    List<NodeGraph> nodes = graph.getNodes();
    for (int i = 0; i < nodes.size(); i++) {
      nodes.get(i).setCentrality(i);
    }
    graph.generateCentralityMap();
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);

    // A high bar leaves the interval empty; the lookup has to walk the percentile down rather
    // than return null on the first try, and must not discard a node it has already drawn.
    NodeGraph drawn = NodesLookup.randomSalientNodeBetweenDistanceInterval(graph, origin, 100.0,
        250.0, 0.95, new MersenneTwisterFast(41L));

    assertNotNull(drawn);
    double distance = GraphUtils.nodesDistance(origin, drawn);
    assertTrue(distance >= 100.0 && distance <= 250.0, "drawn at " + distance);
  }

  @Test
  @DisplayName("the salient lookup returns null only once the bar cannot be lowered further")
  void salientLookupGivesUpWhenNothingIsInRange() {
    Graph graph = Fixtures.grid(5, 5, 100.0);
    List<NodeGraph> nodes = graph.getNodes();
    for (int i = 0; i < nodes.size(); i++) {
      nodes.get(i).setCentrality(i);
    }
    graph.generateCentralityMap();

    assertNull(NodesLookup.randomSalientNodeBetweenDistanceInterval(graph,
        Fixtures.nodeAt(graph, 0.0, 0.0), 5000.0, 6000.0, 0.5, new MersenneTwisterFast(43L)));
  }

  @Test
  void randomNodeDmaDrawsFromTheRequestedCategory() {
    Graph graph = districtGrid();
    Fixtures.nodeAt(graph, 200.0, 200.0).dma = "work";
    MersenneTwisterFast random = new MersenneTwisterFast(29L);

    for (int draw = 0; draw < 10; draw++) {
      assertSame(Fixtures.nodeAt(graph, 200.0, 200.0),
          NodesLookup.randomNodeDMA(graph, "work", random));
    }
  }
}
