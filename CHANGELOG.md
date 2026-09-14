# Changelog

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
