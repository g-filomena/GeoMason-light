/**
 * Portrayal classes for drawing geospatial data in MASON displays.
 *
 * <p>The package provides the bridge between simulation-side geospatial objects and MASON's 2D
 * portrayal system. It contains field and geometry portrayals that can render features held in
 * {@link sim.field.geo.VectorLayer} instances.</p>
 *
 * <h2>Usage notes</h2>
 * <ul>
 *   <li>Use {@link sim.portrayal.geo.GeomVectorFieldPortrayal} to display vector fields.</li>
 *   <li>Use {@link sim.portrayal.geo.GeomPortrayal} to customise individual geometry rendering.</li>
 *   <li>Use {@link sim.portrayal.geo.GeomInfo2D} when display-space geometry metadata is required.</li>
 * </ul>
 *
 * @see sim.field.geo.VectorLayer
 * @see sim.util.geo.MasonGeometry
 */
package sim.portrayal.geo;
