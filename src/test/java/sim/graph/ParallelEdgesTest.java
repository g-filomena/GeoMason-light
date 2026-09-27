/*
 * Copyright (c) 2023 Gabriele Filomena University of Liverpool, UK
 *
 * This program is free software: it can redistributed and/or modified under the terms of the GNU
 * General Public License 3.0 as published by the Free Software Foundation.
 *
 * See the file "LICENSE" for more information
 */
package sim.graph;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.planargraph.DirectedEdge;
import sim.routing.Astar;
import sim.routing.Route;
import sim.testing.Fixtures;

/**
 * Two different streets joining the same pair of junctions: a straight road and a crescent. Every
 * lookup by node pair sees both.
 */
class ParallelEdgesTest {

  private Graph graph;
  private NodeGraph west;
  private NodeGraph east;
  private EdgeGraph straight;
  private EdgeGraph crescent;

  @BeforeEach
  void buildCrescentBesideAStraightRoad() {
    // The crescent is read first, so a last-one-wins map would have kept the straight road and a
    // first-one-wins map the crescent: the lookup must not depend on the order.
    graph = Fixtures.graphOf(Arrays.asList(Fixtures.polyline(0, 0, 50, 60, 100, 0),
        Fixtures.segment(0, 0, 100, 0), Fixtures.segment(-100, 0, 0, 0),
        Fixtures.segment(100, 0, 200, 0)));
    west = Fixtures.nodeAt(graph, 0.0, 0.0);
    east = Fixtures.nodeAt(graph, 100.0, 0.0);
    List<EdgeGraph> between = graph.getEdgesBetween(west, east);
    straight = between.get(0);
    crescent = between.get(1);
  }

  @Test
  @DisplayName("both streets are kept between the pair, shortest first")
  void bothStreetsAreKept() {
    assertEquals(4, graph.getEdges().size());
    assertEquals(2, graph.getEdgesBetween(west, east).size());
    assertEquals(100.0, straight.getLength(), 1e-9);
    assertTrue(crescent.getLength() > straight.getLength());
    assertEquals(graph.getEdgesBetween(west, east), graph.getEdgesBetween(east, west));
  }

  @Test
  @DisplayName("the single-edge lookups answer the shortest street, in either direction")
  void singleLookupsAnswerTheShortest() {
    assertSame(straight, graph.getEdgeBetween(west, east));
    assertSame(straight, graph.getEdgeBetween(east, west));
    assertSame(straight, graph.getDirectedEdgeBetween(west, east).getEdge());
    assertSame(west, graph.getDirectedEdgeBetween(west, east).getFromNode());
    assertSame(east, graph.getDirectedEdgeBetween(east, west).getFromNode());
  }

  @Test
  @DisplayName("every directed edge between the pair runs the requested way")
  void directedLookupsFollowTheDirection() {
    List<DirectedEdge> directed = graph.getDirectedEdgesBetween(west, east);
    assertEquals(2, directed.size());
    for (DirectedEdge directedEdge : directed) {
      assertSame(west, directedEdge.getFromNode());
      assertSame(east, directedEdge.getToNode());
    }
    assertSame(crescent, directed.get(1).getEdge());
  }

  @Test
  @DisplayName("a node joined by two streets is one neighbour, not two")
  void adjacentNodesAreListedOnce() {
    assertEquals(2, west.getAdjacentNodes().size());
    assertEquals(3, west.getEdges().size());
  }

  @Test
  @DisplayName("the shortest route takes the straight road")
  void routeTakesTheShorterStreet() {
    Route route = new Astar().astarRoute(west, east, graph, null);

    assertEquals(Collections.singletonList(straight), route.edgesSequence);
  }

  @Test
  @DisplayName("with the straight road avoided, the route takes the crescent")
  void routeTakesTheCrescentWhenTheRoadIsAvoided() {
    Set<Integer> avoid = new HashSet<>(Collections.singletonList(straight.getID()));

    Route route = new Astar().astarRoute(west, east, graph, avoid);

    assertEquals(Collections.singletonList(crescent), route.edgesSequence);
  }

  @Test
  @DisplayName("an edge set holding only the crescent still joins the pair")
  void islandsSeeAPairJoinedByTheLongerStreetOnly() {
    Set<EdgeGraph> edges = new HashSet<>(graph.getEdges());
    edges.remove(straight);

    assertEquals(1, new Islands(graph).findDisconnectedIslands(edges).size());
  }
}
