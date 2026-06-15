/**
 * GeoMason-light extension packages for geospatial agent-based modelling with MASON.
 *
 * <p>The {@code sim} namespace contains the public API for geospatial fields, import/export
 * utilities, graph abstractions, routing helpers, portrayals, and geometry utilities. The library is
 * designed to be used from MASON simulations that need vector, raster, or street-network-backed
 * spatial behaviour.</p>
 *
 * <h2>Main package groups</h2>
 * <ul>
 *   <li>{@link sim.field.geo}: raster and vector layer abstractions.</li>
 *   <li>{@link sim.io.geo}: geospatial data importers and exporters.</li>
 *   <li>{@link sim.graph}: network, node, edge, and subgraph structures.</li>
 *   <li>{@link sim.routing}: shortest-path and route construction utilities.</li>
 *   <li>{@link sim.portrayal.geo}: visualisation support for geospatial fields.</li>
 *   <li>{@link sim.util.geo}: geometry, attribute, CSV, and coordinate helper classes.</li>
 * </ul>
 *
 * @see sim.field.geo.VectorLayer
 * @see sim.graph.Graph
 * @see sim.util.geo.MasonGeometry
 */
package sim;
