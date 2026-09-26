package sim.io.geo;

import java.io.File;
import java.io.InputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import mil.nga.geopackage.GeoPackage;
import mil.nga.geopackage.GeoPackageManager;
import mil.nga.geopackage.features.user.FeatureDao;
import mil.nga.geopackage.features.user.FeatureRow;
import mil.nga.geopackage.geom.GeoPackageGeometryData;
import mil.nga.sf.Geometry;
import mil.nga.sf.GeometryType;
import mil.nga.sf.Point;
import sim.field.geo.VectorLayer;
import sim.util.geo.AttributeValue;
import sim.util.geo.MasonGeometry;

public class GeoPackageImporter {

  /**
   * Reads the feature table of a GeoPackage holding exactly one, and populates the provided
   * VectorLayer with its geometries and attributes. A GeoPackage without feature tables adds
   * nothing.
   *
   * @param gpkgURL The URL to the GeoPackage file.
   * @param vectorLayer The VectorLayer to populate with data.
   * @throws IllegalStateException If the GeoPackage holds more than one feature table; use {@link
   *     #read(URL, VectorLayer, String)} to name the one to read.
   * @throws Exception If an error occurs during GeoPackage reading or processing.
   */
  public static void read(URL gpkgURL, VectorLayer vectorLayer) throws Exception {
    read(gpkgURL, vectorLayer, null);
  }

  /**
   * Reads one feature table of a GeoPackage and populates the provided VectorLayer with its
   * geometries and attributes.
   *
   * @param gpkgURL The URL to the GeoPackage file.
   * @param vectorLayer The VectorLayer to populate with data.
   * @param tableName The feature table to read, or null to read the file's only feature table.
   * @throws IllegalArgumentException If the named table is not a feature table of the file.
   * @throws IllegalStateException If no table is named and the file holds more than one.
   * @throws Exception If an error occurs during GeoPackage reading or processing.
   */
  public static void read(URL gpkgURL, VectorLayer vectorLayer, String tableName)
      throws Exception {
    try (GeoPackage geoPackage = GeoPackageManager.open(toFile(gpkgURL))) {
      String table = selectTable(geoPackage.getFeatureTables(), tableName, gpkgURL);
      if (table != null) {
        readTable(geoPackage.getFeatureDao(table), vectorLayer);
      }
    }
  }

  /**
   * The table to read: {@code requested} when named, otherwise the only feature table, or null when
   * there is none. A file with several feature tables and no name is refused rather than read as
   * the union of its tables.
   */
  static String selectTable(List<String> tables, String requested, URL gpkgURL) {
    if (requested != null) {
      if (!tables.contains(requested)) {
        throw new IllegalArgumentException(
            gpkgURL + " has no feature table '" + requested + "'; it holds " + tables);
      }
      return requested;
    }
    if (tables.size() > 1) {
      throw new IllegalStateException(
          gpkgURL
              + " holds "
              + tables.size()
              + " feature tables "
              + tables
              + "; name the one to read with readGPKG(url, layer, tableName)");
    }
    return tables.isEmpty() ? null : tables.get(0);
  }

  /** The file behind {@code gpkgURL}; a resource inside a jar is copied to a temporary file. */
  private static File toFile(URL gpkgURL) throws Exception {
    if ("file".equals(gpkgURL.getProtocol())) {
      return new File(gpkgURL.toURI());
    }
    if ("jar".equals(gpkgURL.getProtocol())) {
      File file = File.createTempFile("temp_", ".gpkg");
      file.deleteOnExit();
      try (InputStream input = gpkgURL.openStream()) {
        Files.copy(input, file.toPath(), StandardCopyOption.REPLACE_EXISTING);
      }
      return file;
    }
    throw new IllegalArgumentException("Unsupported URL protocol: " + gpkgURL.getProtocol());
  }

  /** Adds every feature of one table to {@code vectorLayer}. */
  private static void readTable(FeatureDao featureDao, VectorLayer vectorLayer) {
    for (FeatureRow row : featureDao.queryForAll()) {
      // Parse geometry using GeoPackage-Java's GeometryReader
      GeoPackageGeometryData geometryData = row.getGeometry();

      Geometry sfGeometry = null;
      if (geometryData != null && !geometryData.isEmpty()) {
        sfGeometry = geometryData.getGeometry();
      }

      // Convert to JTS Geometry
      org.locationtech.jts.geom.Geometry jtsGeometry = convertToJTSGeometry(sfGeometry);
      // Extract attributes
      Map<String, AttributeValue> attributes = new HashMap<>();
      for (String columnName : featureDao.getTable().getColumnNames()) {

        if (!columnName.equalsIgnoreCase("geometry")) {
          Object value = row.getValue(columnName);
          attributes.put(columnName, parseAttributeValue(value));
        }
      }

      // Add to VectorLayer
      MasonGeometry masonGeometry = new MasonGeometry();
      masonGeometry.geometry = jtsGeometry;
      masonGeometry.addAttributes(attributes);
      vectorLayer.addGeometry(masonGeometry);
    }
  }

  /**
   * Parses an attribute value and converts it into an AttributeValue object.
   *
   * @param value The raw value to parse.
   * @return An AttributeValue representing the parsed data.
   */
  private static AttributeValue parseAttributeValue(Object value) {
    if (value instanceof String) {
      String rawAttributeValue = ((String) value).trim();
      AttributeValue attributeValue = new AttributeValue();

      if (rawAttributeValue.isEmpty()) {
        attributeValue.setString(rawAttributeValue);
      } else {
        switch (determineType(rawAttributeValue)) {
          case "double":
            attributeValue.setDouble(Double.valueOf(rawAttributeValue));
            break;
          case "integer":
            attributeValue.setInteger(Integer.valueOf(rawAttributeValue));
            break;
          case "boolean":
            attributeValue.setValue(Boolean.valueOf(rawAttributeValue));
            break;
          default:
            attributeValue.setString(rawAttributeValue);
            break;
        }
      }

      return attributeValue;
    } else if (value instanceof Long) {
      // Handle Long values explicitly
      long longValue = (Long) value;
      if (longValue >= Integer.MIN_VALUE && longValue <= Integer.MAX_VALUE) {
        return new AttributeValue((int) longValue); // Fits in Integer range
      }
      return new AttributeValue(longValue); // Store as Long

    } else if (value instanceof Integer) {
      return new AttributeValue(value);
    } else if (value instanceof Double) {
      return new AttributeValue(value);
    } else if (value instanceof Boolean) {
      return new AttributeValue(value);
    } else {
      return new AttributeValue(value);
    }

  }

  /**
   * Determines the type of a string value (double, integer, boolean, or string).
   *
   * @param rawAttributeValue The raw string value to analyze.
   * @return A string representing the determined type.
   */
  private static String determineType(String rawAttributeValue) {
    if (rawAttributeValue.matches("^-?\\d+\\.\\d+$")) {
      return "double";
    } else if (rawAttributeValue.matches("^-?\\d+$")) {
      return "integer";
    } else if (rawAttributeValue.equalsIgnoreCase("true")
        || rawAttributeValue.equalsIgnoreCase("false")) {
      return "boolean";
    } else {
      return "string";
    }
  }

  /**
   * Converts a GeoPackage Geometry into a JTS Geometry.
   *
   * @param sfGeometry The GeoPackage Geometry to convert.
   * @return The corresponding JTS Geometry.
   */
  private static org.locationtech.jts.geom.Geometry convertToJTSGeometry(Geometry sfGeometry) {

    GeometryType geometryType = sfGeometry.getGeometryType();

    if (geometryType.equals(GeometryType.POINT)) {
      Point point = (Point) sfGeometry;
      return new GeometryFactory().createPoint(new Coordinate(point.getX(), point.getY()));
    } else if (geometryType.equals(GeometryType.LINESTRING)) {
      mil.nga.sf.LineString ls = (mil.nga.sf.LineString) sfGeometry;
      return new GeometryFactory().createLineString(ls.getPoints().stream()
          .map(point -> new Coordinate(point.getX(), point.getY())).toArray(Coordinate[]::new));
    } else if (geometryType.equals(GeometryType.MULTILINESTRING)) {
      mil.nga.sf.MultiLineString mls = (mil.nga.sf.MultiLineString) sfGeometry;
      GeometryFactory gf = new GeometryFactory();
      if (mls.getLineStrings().size() == 1) {
        mil.nga.sf.LineString singleLineString = mls.getLineStrings().get(0);
        return gf.createLineString(singleLineString.getPoints().stream()
            .map(point -> new Coordinate(point.getX(), point.getY())).toArray(Coordinate[]::new));
      }
      org.locationtech.jts.geom.LineString[] lineStrings = mls.getLineStrings().stream()
          .map(ls -> gf.createLineString(ls.getPoints().stream()
              .map(point -> new Coordinate(point.getX(), point.getY())).toArray(Coordinate[]::new)))
          .toArray(org.locationtech.jts.geom.LineString[]::new);
      return gf.createMultiLineString(lineStrings);
    } else if (geometryType.equals(GeometryType.POLYGON)) {
      mil.nga.sf.Polygon pg = (mil.nga.sf.Polygon) sfGeometry;
      return new GeometryFactory().createPolygon(pg.getExteriorRing().getPoints().stream()
          .map(point -> new Coordinate(point.getX(), point.getY())).toArray(Coordinate[]::new));
    } else if (geometryType.equals(GeometryType.MULTIPOLYGON)) {
      mil.nga.sf.MultiPolygon mpg = (mil.nga.sf.MultiPolygon) sfGeometry;
      GeometryFactory gf = new GeometryFactory();
      if (mpg.getPolygons().size() == 1) {
        mil.nga.sf.Polygon singlePolygon = mpg.getPolygons().get(0);
        return gf.createPolygon(singlePolygon.getExteriorRing().getPoints()
            .stream().map(point -> new Coordinate(point.getX(), point.getY()))
            .toArray(Coordinate[]::new));
      }
      org.locationtech.jts.geom.Polygon[] polygons = mpg.getPolygons().stream()
          .map(pg -> gf.createPolygon(pg.getExteriorRing().getPoints().stream()
              .map(point -> new Coordinate(point.getX(), point.getY())).toArray(Coordinate[]::new)))
          .toArray(org.locationtech.jts.geom.Polygon[]::new);
      return gf.createMultiPolygon(polygons);
    } else {
      throw new IllegalArgumentException(
          "Unsupported GeoPackage geometry type: " + sfGeometry.getGeometryType());
    }
  }
}
