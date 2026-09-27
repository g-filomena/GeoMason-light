# Changelog

## Version 2.3.0

* [Breaking] A pair of nodes can be joined by more than one edge - different streets between the same two junctions, such as a crescent beside a straight road. `Graph`'s adjacency maps held one edge per pair, so whichever street was read last hid the others from every lookup by node pair, and a caller rebuilding a route that way could report a street other than the one walked. The protected `adjacencyMatrix` and `adjacencyMatrixDirected` now map a pair to a list, shortest first.
* [Enhancement] `Graph.getEdgesBetween(from, to)` and `Graph.getDirectedEdgesBetween(from, to)` return every edge between two nodes, shortest first (then lowest id).
* [Note] `getEdgeBetween` and `getDirectedEdgeBetween` keep their signatures and answer the shortest of the parallel edges - the one a shortest path takes - so the answer no longer depends on reading order.
* [Fix] `NodeGraph.getAdjacentNodes()` lists a node joined by parallel edges once; it was listed once per edge.
* [Fix] `Islands` treats two nodes as joined when any of their parallel edges is in the edge set. It checked only the edge `getEdgeBetween` answered, so a set holding only the longer street split the pair into two islands.
* [Enhancement] `RoutingUtils.getPrimalJunction(centroid, otherCentroid, arrivalJunction)` answers the junction a walk crosses from one segment to the next, given the junction it arrived by: the segment's far end. Parallel segments share both ends, and the two-argument form cannot tell which one is meant. `getPreviousJunction` now walks the sequence with it.
* [Enhancement] Added `ParallelEdgesTest` (lookups, direction, adjacency, `Astar`, `Islands` on a crescent beside a straight road) and parallel-street cases in `RoutingUtilsTest`; `Fixtures.polyline(...)` builds a street with bends.
* [Breaking] `GeoPackageImporter` keeps TEXT values as text. It guessed numbers and booleans from them, so "007" came back as 7 and an identifier past the int range threw and aborted the read. Numeric and boolean columns are read as before.
* [Breaking] `CSVUtils.writeLine` quotes a value holding the separator, a double quote or a line break (RFC 4180). Such values were written bare, splitting into extra fields or rows; doubled quotes were written outside any quoted field. With a custom quote character, that character is doubled inside values.
* [Fix] `Route.dualNodesSequence` appended every edge to `edgesSequence` again, doubling it on any route whose sequences had been computed. It now rebuilds both lists.
* [Fix] Shapefile polygons with several outer rings: each polygon was built from the raw ring at its index, so a hole could become an outer boundary, and every hole was given to every polygon. Holes now go to the shell containing them; a file whose rings are all wound the other way is read as shells rather than as an empty geometry.
* [Fix] `ShapeFileImporter` skips null-shape records instead of stopping at the first one, reads to the end of the stream rather than while `available()` is positive (which can be 0 early on network and jar streams), keeps records aligned whatever the shape reader consumes, and closes both files on every path. Numeric fields past the int range are Longs rather than a `NumberFormatException` that aborted the import; unparseable ones are kept as strings. Record and header sizes are read as unsigned.
* [Fix] `ImporterUtils` completes short reads instead of reporting them as the end of the stream.
* [Fix] `VectorLayer.getUnion()` and `getConvexHull()` cast their result to `Polygon`, throwing for disjoint polygons (a MultiPolygon), a single point or collinear points; they return whatever geometry results, and an empty polygon for an empty layer. Both were cached forever: they are now recomputed after any geometry is added, removed, moved or cleared.
* [Fix] `VectorLayer.setGeometryLocation` discards the moved geometry's prepared form, so `coveringFeatures`, `coveredFeatures` and `isCovered` no longer answer for its old position.
* [Fix] `VectorLayer.addGeometry` ignores a geometry that is missing or empty instead of throwing inside the spatial index on its null envelope; such a geometry has no location and is not kept in the layer.
* [Fix] `GridLayer.setGrid` set the pixel size to the number of cells; it now divides the MBR by them when an MBR is set. `GridLayer.toPolygon(x, y)` did not flip the y axis as `toPoint(x, y)` does, so the two described mirrored cells.
* [Fix] `GeoPackageImporter` finds the geometry column by its declared name (QGIS and GDAL use "geom", which was imported as an attribute), skips features with no geometry instead of throwing, and keeps polygon holes.
* [Fix] `GeoPackageExporter` types a column from all its values: integers mixed with decimals make a DOUBLE column (decimals were truncated to the first value's INTEGER), and numbers mixed with other values make a TEXT column.
* [Fix] `GeoJSONExporter` writes an empty geometry as `null`; an empty Point threw.
* [Fix] `GeometryUtilities.screenToWorldPointTransform` called `System.exit(-1)` on a transform it could not invert; it throws `IllegalStateException`. `worldToScreenTransform` pads a zero-width or zero-height extent (one point, points on one line, an empty layer), so drawing such a layer no longer produces that transform.
* [Fix] `GeomPortrayal.draw` set the shared portrayal's `filled` to false on drawing any line, leaving every polygon drawn afterwards unfilled. `hitObject` pads the hit box by `SLOP` rather than half of it.
* [Fix] `NodeGraph.getDualNode` and `getDualNodes` with region-based navigation looked at the far end of each edge from the origin rather than from the node itself, throwing on any node but the origin. `getDualNodes` skips edges with no dual node, as `getDualNode` does.
* [Fix] `NodeGraph.getAdjacentRegion()` rebuilds `adjacentRegionEntries` instead of appending to it on every call. `EdgeGraph.setNodes` replaces the node list instead of appending.
* [Fix] `NodesLookup.getCandidatesByDMA` leaves out nodes with no DMA label instead of throwing. `randomNodeFromDistancesSet` searches the junction layer it is given, which it ignored. Distance lookups leave out junctions matching no graph node instead of returning null entries.
* [Fix] `GraphUtils.smallestEnclosingGeometryBetweenNodes` returns the circle through the two outer nodes for collinear nodes (their hull has no area and contained nothing) and a buffer for nodes at one point; it no longer reorders the caller's list.
* [Fix] `GraphUtils.nodesDistance` no longer caches distances in a static map that grew for the life of the JVM.
* [Fix] A `SubGraph` builds its spatial index, so `getNodesWithinPolygon` no longer always answers empty; `SubGraph(Graph)` registers the graph's nodes and edges, so `findNode` finds them. Populating a `Graph` twice no longer throws from the spatial index.
* [Fix] `Astar` from a node that is itself a target returns a route starting and ending there, not one with null origin and destination.
* [Fix] `Angles.angle` returned NaN for about one destination in twenty-five due north or south, from rounding past the domain of `acos`.
* [Fix] `Utilities.filterMapByPercentile` with a percentile of 0 or an empty map no longer indexes out of bounds.
* [Note] `MasonGeometry.hashCode` documents that equality follows the geometry, so a geometry that moves must not be kept in a hash-based collection.
* [Enhancement] Regression tests for each of the above; new `GridLayerTest`, `NodeGraphTest`, `ImporterUtilsTest`, `ShapeFileImporterTest` and `GeomPortrayalTest`.

## Version 2.2.3

* [Breaking] `VectorLayer.readGPKG(url, layer)` refuses a GeoPackage holding more than one feature table, with an `IllegalStateException` naming them. It used to read every feature table into the one layer, so a file carrying two versions of the same data - an old layer left beside a rewritten one, which writing a GeoPackage layer by name does - loaded as their union, with no error. A file with a single feature table reads exactly as before, and one with none still adds nothing.
* [Enhancement] `VectorLayer.readGPKG(url, layer, tableName)` reads one named feature table, for a file that is meant to hold several. A name the file does not contain is an `IllegalArgumentException` listing the tables it does.
* [Fix] `GeoPackageImporter.read` closes the GeoPackage when reading fails, and closes the stream it copies a jar resource from; neither was closed before.
* [Enhancement] Added `GeoPackageImporterTest`: a single-table file is read, a two-table file is refused and adds nothing, a named table is read alone, and an unknown table name is refused.

## Version 2.2.2

* [Enhancement] `Astar.astarRouteAllowing(origin, destination, graph, Predicate<EdgeGraph>)` admits
  edges by a test applied as the search reaches them, instead of by a set enumerated over the whole
  graph beforehand. A search settles a few hundred edges, so a caller whose exclusion rule is derived
  rather than listed no longer has to materialise it across the network before it can ask.
* [Enhancement] `Astar.astarRouteAllowing(origin, Collection<NodeGraph> targets, graph, predicate)`
  searches for several targets at once and returns the route to whichever is settled first, with
  `reachedTarget()` naming it. The heuristic becomes the distance to the nearest target, which keeps
  it admissible. Scanning candidates previously cost one search each - and when none was reachable,
  one exhaustive failure each.
* [Note] The predicate methods carry their own name rather than overloading `astarRoute`: as
  overloads, `astarRoute(a, b, graph, null)` is ambiguous between `Set<Integer>` and the predicate,
  which stopped this library's own `Islands` compiling and would silently break any caller passing a
  literal null.
* [Fix] `Astar`'s inner loop no longer copies the avoid-set into a fresh `HashSet` per call, resolves
  each neighbour through three `Pair`-allocating map lookups, scans the open set with `contains()`
  and `remove()`, or reconstructs the path by head-insertion into an `ArrayList`. It walks the node's
  own outgoing directed edges, defers deletion, and reverses once. Three of these `Dijkstra` had
  fixed long ago; `Astar` never received them.
* [Enhancement] `AstarTest` covers the new surface: predicate and avoid-set admit the same edges, a
  multi-target search names the nearest target and agrees with the best of the separate searches, an
  unreachable target names nothing, and a predicate admitting nothing reaches nothing.

## Version 2.2.1

* [Fix] `NodeGraph` and `EdgeGraph` now override `hashCode()`, deriving it from the node's coordinate and the edge's centroid. Neither declared one before, so both inherited `Object`'s identity hash, which HotSpot draws from a per-JVM generator - and every `HashMap` or `HashSet` keyed on one iterated in an order that differed between JVM builds. Any decision taken by walking such a collection in order, such as accumulating weights into a cumulative distribution and picking by position, therefore reached a different answer on a different machine from the same seed. Version 2.2.0 addressed this one collection at a time with `LinkedHashSet`; this addresses the cause, and covers callers' collections too. **This changes results**, on every machine, wherever a hash-ordered collection of graph objects was walked - a one-time change, after which the same seed gives the same run anywhere.
* [Note] `equals` is deliberately unchanged and remains identity-based, so nothing about which objects count as equal is affected; only iteration order is. The hash is taken from the coordinate rather than from `getID()` on purpose: `setID(...)` is called after a graph is built, by which time `generateAdjacencyMatrix()` has already stored `Pair<NodeGraph, NodeGraph>` keys whose hash delegates to those nodes. Hashing on a field that changes after insertion would move every such key to a bucket it is not stored in, and adjacency lookups would begin returning `null`. A regression test pins this.
* [Fix] `VectorLayer.intersectingFeatures(Geometry)` returned its results in whatever order the threads happened to finish in: it collected a `parallelStream` into a `ConcurrentLinkedQueue`. That is non-deterministic between two runs on one machine, not merely between machines. It is now sequential and returns the spatial index's own order; the candidate set for a single envelope is small, so there was little to parallelise.
* [Fix] `Route.computeRouteSequences()` threw a `NullPointerException` from inside a stream in `nodesSequence()` when the directed edge sequence contained a `null`, and an `IndexOutOfBoundsException` when it was empty - both several frames from the caller that built the sequence, which is where the mistake is. It now fails with an `IllegalStateException` naming the position of the offending entry.
* [Fix] `Route.computeRouteSequences()` cleared neither `nodesSequence` nor `edgesSequence` before rebuilding them, while `edgeSequence()` appends - so a second call on the same `Route` doubled the edge list, and the length summed from it doubled with it. `resetRoute(...)` cleared them first, so the two entry points disagreed. Both now clear.
* [Enhancement] Added `GraphHashOrderTest`, covering the hashing contract above: the hash is coordinate-derived, survives `setID(...)` and `setRegionID(...)`, keeps adjacency lookups resolving after IDs are assigned, leaves `equals` as identity so that two distinct nodes sharing a coordinate stay distinct in a `HashSet`, and gives identical iteration order for identically built graphs.

## Version 2.2.0

* [Breaking] `NodesLookup` draws through MASON's `MersenneTwisterFast` instead of `java.util.Random`. Every lookup comes in two forms: the one without a generator draws from a per-thread `MersenneTwisterFast` - contention-free under parallel simulation and not reproducible, which is what these lookups have always been - and the overload taking a `MersenneTwisterFast` draws from the caller's own generator, so a model that owns a seeded generator gets the same origins and destinations on every run of the same seed, concurrency included. Callers passing a `java.util.Random` have to switch.
* [Enhancement] `Graph.fromStreetJunctionsSegments(...)` now uses the junction layer it is handed. Each node takes the imported geometry sitting at its coordinate, and the internal `junctions` layer the lookups query is rebuilt from the nodes, so both carry the imported attributes and stay in one-to-one correspondence with the graph. The parameter was previously accepted and discarded, which left every node on a bare generated point and forced callers to re-attach the geometries themselves after the call. Junctions matching no node are skipped, and an empty layer still leaves the generated points in place.
* [Enhancement] Added `Islands.incompleteMerges()`, a count of the island merges that gave up with the islands still separate.
* [Enhancement] Added `GeoJSONExporter.toFeatureCollection(layer, properties)` and `write(fileName, layer, properties)`, which take each feature's properties from a supplied function rather than from its attributes - for writing a layer alongside values the model computed, such as a pedestrian volume per street. Without it, callers rebuilt the document by hand and re-implemented the escaping, the number formatting and the geometry serialisation.
* [Breaking] Removed `NodesLookup.randomNodeFromList(List)` and `randomNodeFromList(List, MersenneTwisterFast)`. Both were character-for-character the same as `selectRandomNode(...)`, which they delegated to; call that instead.
* [Fix] `Route.getLength()` returned 0 for every route: the field behind it was declared and never assigned. It is now summed from the edge sequence whenever that sequence is built, including by `resetRoute(...)`, so a route cut back to the edges actually walked reports the walked length rather than the planned one.
* [Fix] `Islands.mergeConnectedIslands(...)` could fail to terminate. When A* found no route between the closest pair of islands, no edges were added, the island count never fell, and the loop kept asking the same question - a hang rather than a crash, so a run simply stopped making progress. Each pass must now strictly reduce the island count or the merge stops and returns the edges it has, and the give-up is counted rather than silent: a known network left in pieces is one that some origin-destination pairs have no route through.
* [Fix] `NodesLookup.randomNodeBetweenDistanceIntervalDMA(...)` widened only the upper end of its search interval, so a graph too sparse to answer the first interval could only ever be answered by a node further away than asked for - a bias produced by the search rather than by the data. Both ends now widen.
* [Fix] `Angles.isInDirection(...)` rejected directions on the clockwise side of a cone straddling 0 degrees, accepting only its anticlockwise half.
* [Fix] `AttributeValue.getInteger()` and `getDouble()` threw a `NullPointerException` for an attribute carrying no value, because the conditional expression they were written as unboxed the null. They return `null` again, as their boxed return types say they may.
* [Fix] `VectorLayer.filterFeatures(String, String, boolean)` and `filterFeatures(String, int, boolean)` ignored `equal = false` and returned the whole layer instead of its complement; the `List` overload was already correct.
* [Fix] `VectorLayer.getGeometry(String, Object)` compared the `AttributeValue` wrapper against the value asked for, so the lookup always returned `null`.
* [Fix] `VectorLayer.coveringFeatures(MasonGeometry)` threw a `NullPointerException` on any layer whose geometries had not been prepared by an earlier `isCovered(...)` call. It now prepares them on demand, as the sibling relation queries do.
* [Fix] `VectorLayer.intersection(VectorLayer, boolean)` answers about the receiving layer in both directions, so `inclusive = false` is now the exact complement of `inclusive = true`. The exclusive branch used to return the other layer's features minus this one's, which is the complement of nothing, and the inclusive branch repeated a feature once per geometry of the other layer it happened to meet.
* [Fix] `NodesLookup.randomNodeRegion(...)` read the region from a `"district"` string attribute on the junction geometries. Nothing ever put attributes there, so the lookup threw instead of returning a node. It now reads `NodeGraph.getRegionID()`, as `getNodesBetweenDistanceIntervalRegion(...)` already did, and excludes the origin from its own result.
* [Fix] `NodesLookup.randomSalientNodeBetweenDistanceInterval(...)` did the opposite of what its own comment claimed: an empty candidate set returned immediately instead of retrying at a lower percentile, and a draw that did succeed was thrown away whenever the percentile had reached zero. It now lowers the centrality bar until the interval answers.
* [Fix] `Graph.getSalientNodes(...)`, `Graph.getSalientNodesWithinSpace(...)` and `SubGraph.getSubGraphSalientNodes(...)` threw `IndexOutOfBoundsException` at a percentile of 1.0, and on an empty centrality map. The index derived from the percentile is now clamped, and an empty map yields an empty result.
* [Fix] `Islands.findDisconnectedIslands(...)` discarded the ordering of the set it was given, and `GraphUtils.nodesFromEdges(...)`/`edgesFromNodes(...)` never had one. `NodeGraph` and `EdgeGraph` override neither `hashCode` nor `equals`, so a `HashSet` of them iterates in identity-hash order, which HotSpot derives from a per-JVM generator whose values differ between JVM builds. The island list took its order from that iteration, and `mergeConnectedIslands(...)` adds the first bridging edge it finds scanning in that order - so the same model, on the same seed, joined an agent's known network through a different street on a different machine, and every route planned on it followed. These collections are insertion-ordered now, so a caller that supplies a deterministic order keeps it. **This changes results** wherever islands were merged.
* [Performance] `Islands` resolves the closest cross-island pair through one STRtree per island instead of comparing every node of every island against every node of every other, and no longer materialises a tuple and a map entry per cross-island node pair. It also stops enumerating each unordered island pair twice.
* [Performance] `Islands` looks for a bridging edge among adjacent nodes only - the only place an edge can exist - rather than over the full cross product of the islands.
* [Performance] `Islands.findDisconnectedIslands(...)` runs sequentially. The previous `parallelStream` held a lock around its entire body, so every traversal serialised anyway and only the fork-join overhead was left. Its depth-first search also no longer builds and intersects a neighbour list that was never read.
* [Enhancement] Added a JUnit 5 test suite, run by `mvn test`, covering graph construction and adjacency, subgraph mapping, island detection and merging, node lookups, A* routing and route geometry, `VectorLayer` queries, filters and indexes, the GeoJSON exporter, and the angle, attribute, CSV and collection utilities. Fixtures are assembled in memory, so the suite needs no shapefile, GeoPackage or display.
* [Maintenance] Added the `junit-jupiter` test dependency and the Surefire plugin to the build.

## Version 2.1.0

* [Enhancement] Added `GeoPackageExporter`: writes a `VectorLayer` to a single-file OGC GeoPackage (`.gpkg`), the write counterpart of `GeoPackageImporter`. Uses typed columns and imposes no field-name or field-length limits.
* [Enhancement] Added `GeoJSONExporter`: writes a `VectorLayer` to a GeoJSON `FeatureCollection`, with no external JSON dependency. `toFeatureCollection(layer[, includeProperties])` returns the same document as a string.
* [Enhancement] Added `VectorLayer.writeGPKG(...)` and `VectorLayer.writeGeoJSON(...)` convenience methods, mirroring the existing `readGPKG(...)`.
* [Breaking] Removed `ShapeFileExporter`. Vector features are now written in the open GeoPackage and GeoJSON formats; `ShapeFileImporter` is retained, so existing shapefiles can still be read.
* [Enhancement] Added non-copying `VectorLayer` accessors — `isEmpty()`, `isPopulated()`, `size()` and `geometriesView()` (an unmodifiable view) — so callers that only read the geometries avoid the defensive copy made by `getGeometries()`.
* [Performance] `VectorLayer.getGeometriesFromIDs(Set)` now resolves through a lazily-built id index, so its cost scales with the size of the requested set rather than with the number of geometries in the layer.
* [Fix] `NodesLookup` selection methods return `null` on an empty candidate set instead of throwing, so callers no longer have to wrap every lookup in a try/catch. (Recorded late; the change shipped in this release.)
* [Fix] Bounded the search loops in `NodesLookup` with `MAX_EXPANSIONS` and `MAX_DRAW_ATTEMPTS`, so a graph that cannot answer a lookup makes it give up rather than widen or re-draw forever. (Recorded late; the change shipped in this release.)

## Version 2.0.0

* [Breaking] Started a new versioning line after the legacy `1.x` numbering scheme.
* [Breaking] Updated the project Java baseline from Java 8 to Java 11 to align with current GeoPackage dependencies.
* [Breaking] Declared the `sim:mason:21` dependency as `provided`, since MASON is not available from Maven Central and must be supplied by downstream projects.
* [Fix] Updated Maven metadata, including SCM links and project encoding configuration.
* [Fix] Normalised README dependency examples to use the current release version.
* [Fix] Clarified installation instructions for Maven, Eclipse, local builds, and the external MASON dependency.
* [Enhancement] Added package-level Javadocs for the main API packages:

  * `sim`
  * `sim.field.geo`
  * `sim.graph`
  * `sim.io.geo`
  * `sim.portrayal.geo`
  * `sim.routing`
  * `sim.util.geo`
* [Enhancement] Added a Javadoc style guide under `docs/`.
* [Enhancement] Added GitHub Actions validation for Javadoc generation.
* [Maintenance] Updated GitHub Actions configuration to use Java 11.
* [Maintenance] Improved repository hygiene by excluding local Maven settings, signing material, deployment helpers, generated build output, and local IDE files from version control.

## Version 1.1.7, 1.1.8, 1.1.9
Small fixes and new versioning.

## Version 1.16

* Last release under the legacy `1.x` numbering scheme.
* No detailed changelog entry was recorded for this release.

## Version 1.15

* [Enhancement] Allows loading layers from GeoPackage files.
* [Enhancement] Added A* class.
* [Enhancement] Added further utilities in `GraphUtils.java`.

## Version 1.14

* [Enhancement] Cleaned code; removed or simplified redundant functions.
* [Enhancement] Simplified the `SubGraph` class and added utility methods.
* [Enhancement] Replaced concrete `ArrayList` and `HashMap` declarations with generic `List` and `Map` interfaces where appropriate.

## Version 1.13

* No changelog entry was recorded for this release.

## Version 1.12

* [Fix] Replaced `Bag` objects with `ArrayList` objects, usually for `MasonGeometry` objects.
* [Fix] Renamed and consolidated core classes:

  * `GeomVectorField` merged into `VectorLayer`, previously a subclass of `GeomVectorField`.
  * `GeomGridField` renamed to `GridLayer`.
  * `Field` renamed to `Layer`.
* [Enhancement] Cleaned code inherited from GeoMason.
* [Enhancement] Replaced inefficient loops inherited from GeoMason.
* [Enhancement] Simplified partly redundant functions.
* [Enhancement] Removed specific attributes in the `sim.graph` package for better generalisation.

## Version 1.11

* First release.
