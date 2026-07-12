/*
 * Copyright 2011 by Mark Coletti, Keith Sullivan, Sean Luke, and George Mason University Mason
 * University Licensed under the Academic Free License version 3.0
 *
 * See the file "GEOMASON-LICENSE" for more information
 *
 */
package sim.io.geo;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.io.UncheckedIOException;
import java.io.Writer;
import java.math.BigDecimal;
import java.util.Map;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import sim.field.geo.VectorLayer;
import sim.util.geo.AttributeValue;
import sim.util.geo.MasonGeometry;

/**
 * Writes a {@link VectorLayer} to a GeoJSON {@code FeatureCollection} (RFC 7946 structure).
 *
 * <p>Dependency-free: coordinates and attributes are serialised directly, so there is no field-name
 * or field-length limit. Attribute values keep their type — numbers and booleans are written as
 * JSON literals, everything else as a string.
 *
 * <p>Coordinates are written in {@code [x, y]} (easting/northing or lon/lat) order, exactly as they
 * are stored in the layer; no reprojection is performed. Per RFC 7946 GeoJSON coordinates are
 * expected to be WGS84 lon/lat, so reproject the layer beforehand if strict compliance is needed.
 */
public class GeoJSONExporter {

  private GeoJSONExporter() {}

  /**
   * Writes the layer to a GeoJSON file.
   *
   * @param fileName output path; a {@code .geojson} extension is added when none is present
   * @param vectorLayer the layer to export
   * @throws IOException if the file cannot be written
   */
  public static void write(String fileName, VectorLayer vectorLayer) throws IOException {
    String path = hasGeoExtension(fileName) ? fileName : fileName + ".geojson";
    try (Writer writer = new BufferedWriter(new FileWriter(new File(path)))) {
      writeFeatureCollection(writer, vectorLayer, true);
    }
  }

  /**
   * Serialises the layer to a GeoJSON {@code FeatureCollection} string, including feature
   * properties.
   *
   * @param vectorLayer the layer to serialise
   * @return the GeoJSON FeatureCollection as a string
   */
  public static String toFeatureCollection(VectorLayer vectorLayer) {
    return toFeatureCollection(vectorLayer, true);
  }

  /**
   * Serialises the layer to a GeoJSON {@code FeatureCollection} string.
   *
   * @param vectorLayer the layer to serialise
   * @param includeProperties whether to emit each feature's attributes; when {@code false} every
   *        feature carries an empty {@code properties} object (a lighter payload for geometry-only
   *        consumers such as a live map)
   * @return the GeoJSON FeatureCollection as a string
   */
  public static String toFeatureCollection(VectorLayer vectorLayer, boolean includeProperties) {
    StringWriter writer = new StringWriter();
    try {
      writeFeatureCollection(writer, vectorLayer, includeProperties);
    } catch (IOException e) {
      throw new UncheckedIOException(e); // StringWriter does not perform I/O
    }
    return writer.toString();
  }

  private static void writeFeatureCollection(
      Writer writer, VectorLayer vectorLayer, boolean includeProperties) throws IOException {
    writer.write("{\"type\":\"FeatureCollection\",\"features\":[");
    boolean firstFeature = true;
    for (MasonGeometry masonGeometry : vectorLayer.geometriesView()) {
      Geometry geometry = masonGeometry.getGeometry();
      if (geometry == null) {
        continue;
      }
      if (!firstFeature) {
        writer.write(",");
      }
      firstFeature = false;
      writeFeature(writer, masonGeometry, geometry, includeProperties);
    }
    writer.write("]}");
  }

  private static void writeFeature(
      Writer writer, MasonGeometry masonGeometry, Geometry geometry, boolean includeProperties)
      throws IOException {
    writer.write("{\"type\":\"Feature\",\"geometry\":");
    writeGeometry(writer, geometry);
    writer.write(",\"properties\":{");
    Map<String, AttributeValue> attributes =
        includeProperties ? masonGeometry.getAttributes() : null;
    if (attributes != null) {
      boolean firstAttribute = true;
      for (Map.Entry<String, AttributeValue> entry : attributes.entrySet()) {
        if (!firstAttribute) {
          writer.write(",");
        }
        firstAttribute = false;
        writeString(writer, entry.getKey());
        writer.write(":");
        writeValue(writer, entry.getValue() == null ? null : entry.getValue().getValue());
      }
    }
    writer.write("}}");
  }

  // ----------------------------------------------------------------
  // Geometry
  // ----------------------------------------------------------------

  private static void writeGeometry(Writer writer, Geometry geometry) throws IOException {
    if (geometry instanceof Point) {
      writer.write("{\"type\":\"Point\",\"coordinates\":");
      writeCoordinate(writer, ((Point) geometry).getCoordinate());
      writer.write("}");
    } else if (geometry instanceof LineString) {
      writer.write("{\"type\":\"LineString\",\"coordinates\":");
      writeCoordinateArray(writer, geometry.getCoordinates());
      writer.write("}");
    } else if (geometry instanceof Polygon) {
      writer.write("{\"type\":\"Polygon\",\"coordinates\":");
      writePolygonRings(writer, (Polygon) geometry);
      writer.write("}");
    } else if (geometry instanceof MultiPoint) {
      writer.write("{\"type\":\"MultiPoint\",\"coordinates\":");
      writeCoordinateArray(writer, geometry.getCoordinates());
      writer.write("}");
    } else if (geometry instanceof MultiLineString) {
      writer.write("{\"type\":\"MultiLineString\",\"coordinates\":[");
      MultiLineString mls = (MultiLineString) geometry;
      for (int i = 0; i < mls.getNumGeometries(); i++) {
        if (i > 0) {
          writer.write(",");
        }
        writeCoordinateArray(writer, mls.getGeometryN(i).getCoordinates());
      }
      writer.write("]}");
    } else if (geometry instanceof MultiPolygon) {
      writer.write("{\"type\":\"MultiPolygon\",\"coordinates\":[");
      MultiPolygon mp = (MultiPolygon) geometry;
      for (int i = 0; i < mp.getNumGeometries(); i++) {
        if (i > 0) {
          writer.write(",");
        }
        writePolygonRings(writer, (Polygon) mp.getGeometryN(i));
      }
      writer.write("]}");
    } else if (geometry instanceof GeometryCollection) {
      writer.write("{\"type\":\"GeometryCollection\",\"geometries\":[");
      GeometryCollection gc = (GeometryCollection) geometry;
      for (int i = 0; i < gc.getNumGeometries(); i++) {
        if (i > 0) {
          writer.write(",");
        }
        writeGeometry(writer, gc.getGeometryN(i));
      }
      writer.write("]}");
    } else {
      // Unknown/empty geometry: emit JSON null so the feature stays well-formed.
      writer.write("null");
    }
  }

  private static void writePolygonRings(Writer writer, Polygon polygon) throws IOException {
    writer.write("[");
    writeCoordinateArray(writer, polygon.getExteriorRing().getCoordinates());
    for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
      writer.write(",");
      writeCoordinateArray(writer, polygon.getInteriorRingN(i).getCoordinates());
    }
    writer.write("]");
  }

  private static void writeCoordinateArray(Writer writer, Coordinate[] coordinates)
      throws IOException {
    writer.write("[");
    for (int i = 0; i < coordinates.length; i++) {
      if (i > 0) {
        writer.write(",");
      }
      writeCoordinate(writer, coordinates[i]);
    }
    writer.write("]");
  }

  private static void writeCoordinate(Writer writer, Coordinate coordinate) throws IOException {
    writer.write("[");
    writer.write(number(coordinate.x));
    writer.write(",");
    writer.write(number(coordinate.y));
    writer.write("]");
  }

  // ----------------------------------------------------------------
  // Values / JSON primitives
  // ----------------------------------------------------------------

  private static void writeValue(Writer writer, Object value) throws IOException {
    if (value == null) {
      writer.write("null");
    } else if (value instanceof Boolean) {
      writer.write(value.toString());
    } else if (value instanceof Double || value instanceof Float) {
      double d = ((Number) value).doubleValue();
      writer.write((Double.isNaN(d) || Double.isInfinite(d)) ? "null" : number(d));
    } else if (value instanceof Number) {
      writer.write(value.toString());
    } else {
      writeString(writer, value.toString());
    }
  }

  /** Formats a double without scientific notation or trailing zeros, locale-independently. */
  private static String number(double value) {
    return BigDecimal.valueOf(value).stripTrailingZeros().toPlainString();
  }

  private static void writeString(Writer writer, String value) throws IOException {
    writer.write('"');
    for (int i = 0; i < value.length(); i++) {
      char c = value.charAt(i);
      switch (c) {
        case '"':
          writer.write("\\\"");
          break;
        case '\\':
          writer.write("\\\\");
          break;
        case '\n':
          writer.write("\\n");
          break;
        case '\r':
          writer.write("\\r");
          break;
        case '\t':
          writer.write("\\t");
          break;
        case '\b':
          writer.write("\\b");
          break;
        case '\f':
          writer.write("\\f");
          break;
        default:
          if (c < 0x20) {
            writer.write(String.format("\\u%04x", (int) c));
          } else {
            writer.write(c);
          }
      }
    }
    writer.write('"');
  }

  private static boolean hasGeoExtension(String fileName) {
    String lower = fileName.toLowerCase();
    return lower.endsWith(".geojson") || lower.endsWith(".json");
  }
}
