/*
 * Copyright (c) 2023 Gabriele Filomena University of Liverpool, UK
 *
 * This program is free software: it can redistributed and/or modified under the terms of the GNU
 * General Public License 3.0 as published by the Free Software Foundation.
 *
 * See the file "LICENSE" for more information
 */
package sim.testing;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import sim.field.geo.VectorLayer;
import sim.graph.EdgeGraph;
import sim.graph.Graph;
import sim.graph.NodeGraph;
import sim.util.geo.MasonGeometry;

/**
 * Builders shared by the test suite. Everything is assembled in memory: no shapefile, no
 * GeoPackage, no display, so the tests run anywhere the library compiles.
 */
public final class Fixtures {

  public static final GeometryFactory FACTORY = new GeometryFactory();

  private Fixtures() {}

  /** A two-point LineString wrapped as a MasonGeometry, the shape a street segment arrives in. */
  public static MasonGeometry segment(double fromX, double fromY, double toX, double toY) {
    LineString line = FACTORY.createLineString(
        new Coordinate[] {new Coordinate(fromX, fromY), new Coordinate(toX, toY)});
    return new MasonGeometry(line);
  }

  /** A point wrapped as a MasonGeometry. */
  public static MasonGeometry point(double x, double y) {
    return new MasonGeometry(FACTORY.createPoint(new Coordinate(x, y)));
  }

  /**
   * Builds a graph from the given segments and numbers its nodes and edges, which importers
   * normally do and which several lookups rely on. No junction layer is supplied, so the nodes
   * keep the bare points the graph generates for them.
   */
  public static Graph graphOf(List<MasonGeometry> segments) {
    return graphOf(segments, new VectorLayer());
  }

  /** As {@link #graphOf(List)}, with a junction layer for the graph to take its node geometries from. */
  public static Graph graphOf(List<MasonGeometry> segments, VectorLayer junctions) {
    Graph graph = new Graph();
    graph.fromStreetJunctionsSegments(junctions, new VectorLayer(segments));

    int nodeID = 0;
    for (NodeGraph node : graph.getNodes()) {
      node.setID(nodeID++);
    }
    int edgeID = 0;
    for (EdgeGraph edge : graph.getEdges()) {
      edge.setID(edgeID++);
    }
    return graph;
  }

  /**
   * The distinct endpoints of the given segments, as the attributed junction layer an importer
   * would hand over: each junction carries a {@code nodeID} and a {@code district}.
   */
  public static VectorLayer junctionLayer(List<MasonGeometry> segments) {
    Map<Coordinate, MasonGeometry> byCoordinate = new LinkedHashMap<>();
    for (MasonGeometry segment : segments) {
      for (Coordinate coordinate : segment.getGeometry().getCoordinates()) {
        byCoordinate.computeIfAbsent(coordinate, c -> point(c.x, c.y));
      }
    }

    VectorLayer layer = new VectorLayer();
    int nodeID = 0;
    for (MasonGeometry junction : byCoordinate.values()) {
      junction.addIntegerAttribute("nodeID", nodeID++);
      junction.addIntegerAttribute("district", 1);
      layer.addGeometry(junction);
    }
    return layer;
  }

  /** The segments of a regular street grid of {@code columns x rows} junctions. */
  public static List<MasonGeometry> gridSegments(int columns, int rows, double spacing) {
    List<MasonGeometry> segments = new ArrayList<>();
    for (int i = 0; i < columns; i++) {
      for (int j = 0; j < rows; j++) {
        if (i < columns - 1) {
          segments.add(segment(i * spacing, j * spacing, (i + 1) * spacing, j * spacing));
        }
        if (j < rows - 1) {
          segments.add(segment(i * spacing, j * spacing, i * spacing, (j + 1) * spacing));
        }
      }
    }
    return segments;
  }

  /** A regular street grid of {@code columns x rows} junctions, {@code spacing} apart. */
  public static Graph grid(int columns, int rows, double spacing) {
    return graphOf(gridSegments(columns, rows, spacing));
  }

  /** As {@link #grid(int, int, double)}, built with an attributed junction layer. */
  public static Graph gridWithJunctions(int columns, int rows, double spacing) {
    List<MasonGeometry> segments = gridSegments(columns, rows, spacing);
    return graphOf(segments, junctionLayer(segments));
  }

  /** A single chain of {@code junctions} nodes, {@code spacing} apart along the x axis. */
  public static Graph path(int junctions, double spacing) {
    List<MasonGeometry> segments = new ArrayList<>();
    for (int i = 0; i < junctions - 1; i++) {
      segments.add(segment(i * spacing, 0.0, (i + 1) * spacing, 0.0));
    }
    return graphOf(segments);
  }

  /** Two chains far enough apart that no edge or path joins them. */
  public static Graph twoDisjointPaths() {
    List<MasonGeometry> segments = new ArrayList<>();
    segments.add(segment(0.0, 0.0, 100.0, 0.0));
    segments.add(segment(100.0, 0.0, 200.0, 0.0));
    segments.add(segment(10000.0, 0.0, 10100.0, 0.0));
    return graphOf(segments);
  }

  public static NodeGraph nodeAt(Graph graph, double x, double y) {
    return graph.findNode(new Coordinate(x, y));
  }
}
