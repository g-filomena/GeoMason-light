/**
 * Import and export utilities for geospatial file formats used by GeoMason-light.
 *
 * <p>This package contains readers and writers for raster and vector formats, including ArcInfo ASCII
 * grids, ESRI Shapefiles, and GeoPackages. Importers normally populate {@link sim.field.geo.GridLayer}
 * or {@link sim.field.geo.VectorLayer} instances with geometry and attribute data that can then be
 * used by simulations, graphs, and portrayals.</p>
 *
 * <h2>Supported workflows</h2>
 * <ul>
 *   <li>Read raster grids with {@link sim.io.geo.ArcInfoASCGridImporter}.</li>
 *   <li>Write raster grids with {@link sim.io.geo.ArcInfoASCGridExporter}.</li>
 *   <li>Read vector features with {@link sim.io.geo.ShapeFileImporter} and
 *       {@link sim.io.geo.GeoPackageImporter}.</li>
 *   <li>Write vector features with {@link sim.io.geo.ShapeFileExporter}.</li>
 * </ul>
 *
 * <p>Callers should validate coordinate reference system assumptions outside these importers when
 * mixing layers from different sources.</p>
 *
 * @see sim.field.geo.VectorLayer
 * @see sim.field.geo.GridLayer
 * @see sim.util.geo.MasonGeometry
 */
package sim.io.geo;
