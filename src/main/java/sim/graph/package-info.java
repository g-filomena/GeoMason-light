/**
 * Graph primitives for street-network and component-based geospatial models.
 *
 * <p>The package provides mutable graph structures centred on {@link sim.graph.Graph},
 * {@link sim.graph.NodeGraph}, and {@link sim.graph.EdgeGraph}. These classes are intended for
 * network-based agent movement, route choice, accessibility, and urban mobility simulations.</p>
 *
 * <h2>Core concepts</h2>
 * <ul>
 *   <li>{@link sim.graph.NodeGraph}: graph vertices backed by spatial coordinates.</li>
 *   <li>{@link sim.graph.EdgeGraph}: graph edges backed by JTS line geometries and attributes.</li>
 *   <li>{@link sim.graph.SubGraph}: connected or analytical subsets of a larger graph.</li>
 *   <li>{@link sim.graph.GraphUtils}: shared graph-processing helpers.</li>
 * </ul>
 *
 * <p>When graph objects are derived from GIS data, preserve the source attributes needed by routing
 * and behaviour models before simplifying or segmenting the network.</p>
 *
 * @see sim.routing.Astar
 * @see sim.routing.Route
 * @see sim.field.geo.VectorLayer
 */
package sim.graph;
