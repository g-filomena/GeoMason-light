package sim.io.geo;

import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;
import org.locationtech.jts.algorithm.CGAlgorithms;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Point;
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
    readFully(in, new byte[num]);
  }

  /**
   * Fills {@code b} from the stream. A single {@code read} may legitimately return fewer bytes than
   * asked for - network and jar streams often do - so this loops until the buffer is full, and
   * only a real end of stream is an error.
   *
   * @param in the stream to read from.
   * @param b the buffer to fill.
   * @throws IOException if the stream ends before the buffer is full.
   */
  static void readFully(final InputStream in, final byte[] b) throws IOException {
    int read = 0;
    while (read < b.length) {
      final int chk = in.read(b, read, b.length - read);
      if (chk == -1) {
        throw new IOException("Unexpected end of stream after " + read + " of " + b.length
            + " bytes");
      }
      read += chk;
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
   * If there is only one ring the function will create and return a simple polygon. Otherwise
   * clockwise rings are shells and counter-clockwise rings are holes, as the Shapefile
   * specification has it; each hole is given to the shell that contains it. A single shell gives a
   * polygon, several give a multi-polygon.
   *
   */
  private static Geometry createPolygon(final LinearRing[] parts) {

    if (parts.length == 1) {
      return GEOMETRY_FACTORY.createPolygon(parts[0], null);
    }

    final List<LinearRing> shells = new ArrayList<>();
    final List<LinearRing> holes = new ArrayList<>();

    for (LinearRing part : parts) {
      if (CGAlgorithms.isCCW(part.getCoordinates())) {
        holes.add(part);
      } else {
        shells.add(part);
      }
    }

    // A file written with the orientation reversed has no clockwise ring at all; rather than
    // returning an empty geometry, read every ring as a shell.
    if (shells.isEmpty()) {
      shells.addAll(holes);
      holes.clear();
    }

    // Each hole goes to the shell that contains it. Building the polygons from the raw parts could
    // turn a hole into an outer boundary, and handing every hole to every shell produced polygons
    // with holes lying outside them.
    final List<Polygon> shellPolygons = new ArrayList<>();
    final List<List<LinearRing>> holesOfShell = new ArrayList<>();
    for (LinearRing shell : shells) {
      shellPolygons.add(GEOMETRY_FACTORY.createPolygon(shell));
      holesOfShell.add(new ArrayList<>());
    }
    for (LinearRing hole : holes) {
      final int owner = shellContaining(shellPolygons, hole);
      if (owner >= 0) {
        holesOfShell.get(owner).add(hole);
      }
    }

    final Polygon[] polygons = new Polygon[shells.size()];
    for (int i = 0; i < shells.size(); i++) {
      polygons[i] = GEOMETRY_FACTORY.createPolygon(shells.get(i),
          holesOfShell.get(i).toArray(new LinearRing[0]));
    }
    return polygons.length == 1 ? polygons[0] : GEOMETRY_FACTORY.createMultiPolygon(polygons);
  }

  /**
   * Finds the shell a hole belongs to: the smallest shell containing the hole's first vertex, so
   * that a hole inside an island inside a lake goes to the island.
   *
   * @return the index of the owning shell, or -1 when no shell contains the hole.
   */
  private static int shellContaining(final List<Polygon> shells, final LinearRing hole) {
    final Point probe = GEOMETRY_FACTORY.createPoint(hole.getCoordinateN(0));
    int owner = -1;
    double ownerArea = Double.MAX_VALUE;
    for (int i = 0; i < shells.size(); i++) {
      final Polygon shell = shells.get(i);
      if (shell.getEnvelopeInternal().contains(probe.getCoordinate()) && shell.covers(probe)
          && shell.getArea() < ownerArea) {
        owner = i;
        ownerArea = shell.getArea();
      }
    }
    return owner;
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
    readFully(stream, b);
    return b[0];
  }

  public static short readShort(final InputStream stream, final boolean littleEndian)
      throws RuntimeException, IOException {
    final byte[] b = new byte[2];
    readFully(stream, b);
    return ByteBuffer.wrap(b).order((littleEndian) ? ByteOrder.LITTLE_ENDIAN : ByteOrder.BIG_ENDIAN)
        .getShort();
  }

  public static int readInt(final InputStream stream, final boolean littleEndian)
      throws RuntimeException, IOException {
    final byte[] b = new byte[4];
    readFully(stream, b);
    return ByteBuffer.wrap(b).order((littleEndian) ? ByteOrder.LITTLE_ENDIAN : ByteOrder.BIG_ENDIAN)
        .getInt();
  }

  public static double readDouble(final InputStream stream, final boolean littleEndian)
      throws RuntimeException, IOException {
    final byte[] b = new byte[8];
    readFully(stream, b);
    return ByteBuffer.wrap(b).order((littleEndian) ? ByteOrder.LITTLE_ENDIAN : ByteOrder.BIG_ENDIAN)
        .getDouble();
  }
}
