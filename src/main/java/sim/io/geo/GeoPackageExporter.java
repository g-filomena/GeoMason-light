/*
 * Copyright 2011 by Mark Coletti, Keith Sullivan, Sean Luke, and George Mason University Mason
 * University Licensed under the Academic Free License version 3.0
 *
 * See the file "GEOMASON-LICENSE" for more information
 *
 */
package sim.io.geo;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import mil.nga.geopackage.BoundingBox;
import mil.nga.geopackage.GeoPackage;
import mil.nga.geopackage.GeoPackageManager;
import mil.nga.geopackage.db.GeoPackageDataType;
import mil.nga.geopackage.db.TableColumnKey;
import mil.nga.geopackage.features.columns.GeometryColumns;
import mil.nga.geopackage.features.user.FeatureColumn;
import mil.nga.geopackage.features.user.FeatureDao;
import mil.nga.geopackage.features.user.FeatureRow;
import mil.nga.geopackage.features.user.FeatureTableMetadata;
import mil.nga.geopackage.geom.GeoPackageGeometryData;
import mil.nga.sf.GeometryType;
import sim.field.geo.VectorLayer;
import sim.util.geo.AttributeValue;
import sim.util.geo.MasonGeometry;

/**
 * Writes a {@link VectorLayer} to a single-file OGC GeoPackage (.gpkg), the read counterpart of
 * {@link GeoPackageImporter}. Replaces the ESRI shapefile exporter: no 10-character column-name
 * limit, no 254-character field limit, typed columns, and one file instead of three.
 *
 * <p>The feature table schema is the union of the layer's attribute names; each column's type is
 * inferred from the first non-null value seen for it (text / integer / double / boolean). Geometry
 * is stored in a generic {@code geometry} column, so mixed geometry types in one layer are allowed.
 *
 * <p>Coordinates are written verbatim (no reprojection). The spatial-reference id defaults to
 * {@link #UNDEFINED_CARTESIAN_SRS_ID} (0 = "undefined cartesian"), matching the projected metre
 * coordinates typical of a MASON layer; pass a real EPSG code via {@link #write(String, VectorLayer,
 * long)} when the CRS is known.
 */
public class GeoPackageExporter {

  /** GeoPackage-spec srs_id for "undefined cartesian coordinate reference system". */
  public static final long UNDEFINED_CARTESIAN_SRS_ID = 0L;

  private static final String GEOMETRY_COLUMN = "geometry";
  private static final String ID_COLUMN = "id";

  private GeoPackageExporter() {}

  /**
   * Writes the layer to a GeoPackage with an undefined (projected) coordinate reference system.
   *
   * @param fileName output path; a {@code .gpkg} extension is added when none is present
   * @param vectorLayer the layer to export
   * @throws IOException if the GeoPackage cannot be created or written
   */
  public static void write(String fileName, VectorLayer vectorLayer) throws IOException {
    write(fileName, vectorLayer, UNDEFINED_CARTESIAN_SRS_ID);
  }

  /**
   * Writes the layer to a GeoPackage tagged with the given spatial-reference id.
   *
   * @param fileName output path; a {@code .gpkg} extension is added when none is present
   * @param vectorLayer the layer to export
   * @param srsId EPSG code of the layer's coordinates, or {@link #UNDEFINED_CARTESIAN_SRS_ID}
   * @throws IOException if the GeoPackage cannot be created or written
   */
  public static void write(String fileName, VectorLayer vectorLayer, long srsId)
      throws IOException {
    String path = fileName.toLowerCase().endsWith(".gpkg") ? fileName : fileName + ".gpkg";
    List<MasonGeometry> geometries = vectorLayer.getGeometries();
    if (geometries.isEmpty()) {
      throw new IOException("Cannot export an empty layer to " + path);
    }

    File file = new File(path);
    if (file.exists() && !file.delete()) {
      throw new IOException("Could not overwrite existing GeoPackage: " + path);
    }
    GeoPackageManager.create(file);

    try (GeoPackage geoPackage = GeoPackageManager.open(file)) {
      if (srsId > 0) {
        // Ensure the referenced EPSG row exists in gpkg_spatial_ref_sys.
        geoPackage.getSpatialReferenceSystemDao().getOrCreateFromEpsg(srsId);
      }

      String tableName = sanitiseTableName(file.getName());
      Map<String, GeoPackageDataType> columnTypes = inferColumnTypes(geometries);

      List<FeatureColumn> columns = new ArrayList<>();
      for (Map.Entry<String, GeoPackageDataType> entry : columnTypes.entrySet()) {
        columns.add(FeatureColumn.createColumn(entry.getKey(), entry.getValue()));
      }

      GeometryColumns geometryColumns = new GeometryColumns();
      geometryColumns.setId(new TableColumnKey(tableName, GEOMETRY_COLUMN));
      geometryColumns.setGeometryType(GeometryType.GEOMETRY);
      geometryColumns.setZ((byte) 0);
      geometryColumns.setM((byte) 0);
      geometryColumns.setSrsId(srsId);

      Envelope mbr = vectorLayer.getMBR();
      BoundingBox boundingBox =
          new BoundingBox(mbr.getMinX(), mbr.getMinY(), mbr.getMaxX(), mbr.getMaxY());

      geoPackage.createFeatureTable(
          FeatureTableMetadata.create(geometryColumns, columns, boundingBox));

      FeatureDao featureDao = geoPackage.getFeatureDao(tableName);
      for (MasonGeometry masonGeometry : geometries) {
        Geometry geometry = masonGeometry.getGeometry();
        if (geometry == null) {
          continue;
        }
        FeatureRow row = featureDao.newRow();
        row.setGeometry(GeoPackageGeometryData.create(srsId, toSfGeometry(geometry)));
        writeAttributes(row, masonGeometry, columnTypes);
        featureDao.insert(row);
      }
    } catch (SQLException e) {
      throw new IOException("Failed to write GeoPackage " + path, e);
    }
  }

  // ----------------------------------------------------------------
  // Schema
  // ----------------------------------------------------------------

  /**
   * Builds the column schema: the union of attribute names (insertion order preserved), each typed
   * from its first non-null value. Names colliding with the reserved id / geometry columns are
   * suffixed so no data is lost.
   */
  private static Map<String, GeoPackageDataType> inferColumnTypes(List<MasonGeometry> geometries) {
    Map<String, GeoPackageDataType> types = new LinkedHashMap<>();
    for (MasonGeometry masonGeometry : geometries) {
      Map<String, AttributeValue> attributes = masonGeometry.getAttributes();
      if (attributes == null) {
        continue;
      }
      for (Map.Entry<String, AttributeValue> entry : attributes.entrySet()) {
        String column = safeColumnName(entry.getKey());
        Object value = entry.getValue() == null ? null : entry.getValue().getValue();
        GeoPackageDataType inferred = dataType(value);
        GeoPackageDataType existing = types.get(column);
        // Keep the first concrete type; only upgrade away from the TEXT default when we finally
        // see a typed value.
        if (existing == null || (existing == GeoPackageDataType.TEXT && value != null)) {
          types.put(column, inferred);
        }
      }
    }
    return types;
  }

  private static GeoPackageDataType dataType(Object value) {
    if (value instanceof Boolean) {
      return GeoPackageDataType.BOOLEAN;
    } else if (value instanceof Double || value instanceof Float) {
      return GeoPackageDataType.DOUBLE;
    } else if (value instanceof Integer || value instanceof Long || value instanceof Short) {
      return GeoPackageDataType.INTEGER;
    }
    return GeoPackageDataType.TEXT;
  }

  private static void writeAttributes(
      FeatureRow row, MasonGeometry masonGeometry, Map<String, GeoPackageDataType> columnTypes) {
    Map<String, AttributeValue> attributes = masonGeometry.getAttributes();
    if (attributes == null) {
      return;
    }
    for (Map.Entry<String, AttributeValue> entry : attributes.entrySet()) {
      String column = safeColumnName(entry.getKey());
      Object value = entry.getValue() == null ? null : entry.getValue().getValue();
      if (value == null) {
        continue;
      }
      row.setValue(column, coerce(value, columnTypes.get(column)));
    }
  }

  /** Coerces the raw attribute value to the column's declared type. */
  private static Object coerce(Object value, GeoPackageDataType type) {
    if (type == GeoPackageDataType.INTEGER && value instanceof Number) {
      return ((Number) value).longValue();
    } else if (type == GeoPackageDataType.DOUBLE && value instanceof Number) {
      return ((Number) value).doubleValue();
    } else if (type == GeoPackageDataType.BOOLEAN && value instanceof Boolean) {
      return value;
    } else if (type == GeoPackageDataType.TEXT) {
      return value.toString();
    }
    return value;
  }

  private static String safeColumnName(String name) {
    if (name.equalsIgnoreCase(ID_COLUMN) || name.equalsIgnoreCase(GEOMETRY_COLUMN)) {
      return name + "_attr";
    }
    return name;
  }

  /** Turns a file base name into a valid SQL table identifier. */
  private static String sanitiseTableName(String fileName) {
    String base = fileName;
    int dot = base.lastIndexOf('.');
    if (dot > 0) {
      base = base.substring(0, dot);
    }
    base = base.replaceAll("[^A-Za-z0-9_]", "_");
    if (base.isEmpty()) {
      base = "features";
    }
    if (Character.isDigit(base.charAt(0))) {
      base = "t_" + base;
    }
    return base;
  }

  // ----------------------------------------------------------------
  // JTS -> mil.nga.sf geometry conversion
  // ----------------------------------------------------------------

  private static mil.nga.sf.Geometry toSfGeometry(Geometry geometry) {
    if (geometry instanceof Point) {
      return toSfPoint(((Point) geometry).getCoordinate());
    } else if (geometry instanceof LineString) {
      return toSfLineString((LineString) geometry);
    } else if (geometry instanceof Polygon) {
      return toSfPolygon((Polygon) geometry);
    } else if (geometry instanceof MultiPoint) {
      mil.nga.sf.MultiPoint multiPoint = new mil.nga.sf.MultiPoint();
      for (Coordinate coordinate : geometry.getCoordinates()) {
        multiPoint.addPoint(toSfPoint(coordinate));
      }
      return multiPoint;
    } else if (geometry instanceof MultiLineString) {
      mil.nga.sf.MultiLineString multiLineString = new mil.nga.sf.MultiLineString();
      MultiLineString mls = (MultiLineString) geometry;
      for (int i = 0; i < mls.getNumGeometries(); i++) {
        multiLineString.addLineString(toSfLineString((LineString) mls.getGeometryN(i)));
      }
      return multiLineString;
    } else if (geometry instanceof MultiPolygon) {
      mil.nga.sf.MultiPolygon multiPolygon = new mil.nga.sf.MultiPolygon();
      MultiPolygon mp = (MultiPolygon) geometry;
      for (int i = 0; i < mp.getNumGeometries(); i++) {
        multiPolygon.addPolygon(toSfPolygon((Polygon) mp.getGeometryN(i)));
      }
      return multiPolygon;
    } else if (geometry instanceof GeometryCollection) {
      mil.nga.sf.GeometryCollection<mil.nga.sf.Geometry> collection =
          new mil.nga.sf.GeometryCollection<>();
      GeometryCollection gc = (GeometryCollection) geometry;
      for (int i = 0; i < gc.getNumGeometries(); i++) {
        collection.addGeometry(toSfGeometry(gc.getGeometryN(i)));
      }
      return collection;
    }
    throw new IllegalArgumentException(
        "Unsupported geometry type for GeoPackage export: " + geometry.getGeometryType());
  }

  private static mil.nga.sf.Point toSfPoint(Coordinate coordinate) {
    return new mil.nga.sf.Point(coordinate.x, coordinate.y);
  }

  private static mil.nga.sf.LineString toSfLineString(LineString lineString) {
    mil.nga.sf.LineString sfLineString = new mil.nga.sf.LineString();
    for (Coordinate coordinate : lineString.getCoordinates()) {
      sfLineString.addPoint(toSfPoint(coordinate));
    }
    return sfLineString;
  }

  private static mil.nga.sf.Polygon toSfPolygon(Polygon polygon) {
    mil.nga.sf.Polygon sfPolygon = new mil.nga.sf.Polygon();
    sfPolygon.addRing(ringToSfLineString(polygon.getExteriorRing()));
    for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
      sfPolygon.addRing(ringToSfLineString(polygon.getInteriorRingN(i)));
    }
    return sfPolygon;
  }

  private static mil.nga.sf.LineString ringToSfLineString(LineString ring) {
    mil.nga.sf.LineString sfRing = new mil.nga.sf.LineString();
    for (Coordinate coordinate : ring.getCoordinates()) {
      sfRing.addPoint(toSfPoint(coordinate));
    }
    return sfRing;
  }
}
