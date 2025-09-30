package sim.io.geo;

import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import org.locationtech.jts.algorithm.CGAlgorithms;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;

public class ImporterUtils {

  final static int NULL_SHAPE = 0;
  final static int POINT = 1;
  final static int POLYLINE = 3;
  final static int POLYGON = 5;
  final static int MULTIPOINT = 8;
  final static int POINTZ = 11;
  final static int POLYLINEZ = 13;
  final static int POLYGONZ = 15;
  final static int MULTIPOINTZ = 18;
  final static int POINTM = 21;
  final static int POLYLINEM = 23;
  final static int POLYGONM = 25;
  final static int MULTIPOINTM = 28;
  final static int MULTIPATCH = 31;
  final static GeometryFactory GEOMETRY_FACTORY = new GeometryFactory();

  public static boolean isSupported(final int shapeType) {
    switch (shapeType) {
      case POINT:
      case POLYLINE:
      case POLYGON:
      case POINTZ:
        return true;
      default:
        return false; // no other types are currently supported
    }
  }

  public static void skip(final InputStream in, final int num)
      throws RuntimeException, IOException {
    final byte[] b = new byte[num];
    final int chk = in.read(b);
    if (chk != num || chk == -1) {
      throw new IOException("Bad seek, chk = " + chk);
    }
  }

  /**
   * Wrapper function which creates a new array of LinearRings and calls the other function.
   */
  static Geometry createPolygon(final Geometry[] parts) {
    LinearRing[] rings = new LinearRing[parts.length];
    for (int i = 0; i < parts.length; i++) {
      rings[i] = (LinearRing) parts[i];
    }

    return createPolygon(rings);
  }

  /**
   * Create a polygon from an array of LinearRings.
   *
   * If there is only one ring the function will create and return a simple polygon. If there are
   * multiple rings, the function checks to see if any of them are holes (which are in
   * counter-clockwise order) and if so, it creates a polygon with holes. If there are no holes, it
   * creates and returns a multi-part polygon.
   *
   */
  private static Geometry createPolygon(final LinearRing[] parts) {

    if (parts.length == 1) {
      return GEOMETRY_FACTORY.createPolygon(parts[0], null);
    }

    final ArrayList<LinearRing> shells = new ArrayList<>();
    final ArrayList<LinearRing> holes = new ArrayList<>();

    for (LinearRing part : parts) {
      if (CGAlgorithms.isCCW(part.getCoordinates())) {
        holes.add(part);
      } else {
        shells.add(part);
      }
    }

    // This will contain any holes within a given polygon
    LinearRing[] holesArray = null;

    if (!holes.isEmpty()) {
      holesArray = new LinearRing[holes.size()];
      holes.toArray(holesArray);
    }

    // single polygon
    if (shells.size() == 1) {
      // It's ok if holesArray is null
      return GEOMETRY_FACTORY.createPolygon(shells.get(0), holesArray);
    } else {
      Polygon[] poly = new Polygon[shells.size()];
      for (int i = 0; i < shells.size(); i++) {
        poly[i] = GEOMETRY_FACTORY.createPolygon(parts[i], holesArray);
      }
      return GEOMETRY_FACTORY.createMultiPolygon(poly);
    }
  }

  static String typeToString(final int shapeType) {
    switch (shapeType) {
      case NULL_SHAPE:
        return "NULL_SHAPE";
      case POINT:
        return "POINT";
      case POLYLINE:
        return "POLYLINE";
      case POLYGON:
        return "POLYGON";
      case MULTIPOINT:
        return "MULTIPOINT";
      case POINTZ:
        return "POINTZ";
      case POLYLINEZ:
        return "POLYLINEZ";
      case POLYGONZ:
        return "POLYGONZ";
      case MULTIPOINTZ:
        return "MULTIPOINTZ";
      case POINTM:
        return "POINTM";
      case POLYLINEM:
        return "POLYLINEM";
      case POLYGONM:
        return "POLYGONM";
      case MULTIPOINTM:
        return "MULTIPOINTM";
      case MULTIPATCH:
        return "MULTIPATCH";
      default:
        return "UNKNOWN";
    }
  }

  public static boolean littleEndian =
      java.nio.ByteOrder.nativeOrder().equals(java.nio.ByteOrder.LITTLE_ENDIAN); // for

  public static byte readByte(final InputStream stream, final boolean littleEndian)
      throws RuntimeException, IOException {
    final byte[] b = new byte[1];
    final int chk = stream.read(b);
    if (chk != b.length || chk == -1) {
      throw new IOException("readByte early termination, chk = " + chk);
    }
    return b[0];
  }

  public static short readShort(final InputStream stream, final boolean littleEndian)
      throws RuntimeException, IOException {
    final byte[] b = new byte[2];
    final int chk = stream.read(b);
    if (chk != b.length || chk == -1) {
      throw new IOException("readShort early termination, chk = " + chk);
    }
    return ByteBuffer.wrap(b).order((littleEndian) ? ByteOrder.LITTLE_ENDIAN : ByteOrder.BIG_ENDIAN)
        .getShort();
  }

  public static int readInt(final InputStream stream, final boolean littleEndian)
      throws RuntimeException, IOException {
    final byte[] b = new byte[4];
    final int chk = stream.read(b);
    if (chk != b.length || chk == -1) {
      throw new IOException("readInt early termination, chk = " + chk);
    }
    return ByteBuffer.wrap(b).order((littleEndian) ? ByteOrder.LITTLE_ENDIAN : ByteOrder.BIG_ENDIAN)
        .getInt();
  }

  public static double readDouble(final InputStream stream, final boolean littleEndian)
      throws RuntimeException, IOException {
    final byte[] b = new byte[8];
    final int chk = stream.read(b);
    if (chk != b.length || chk == -1) {
      throw new IOException("readDouble early termination, chk = " + chk);
    }
    return ByteBuffer.wrap(b).order((littleEndian) ? ByteOrder.LITTLE_ENDIAN : ByteOrder.BIG_ENDIAN)
        .getDouble();
  }
}
