/*
 * Copyright 2011 by Mark Coletti, Keith Sullivan, Sean Luke, and George Mason University Mason
 * University Licensed under the Academic Free License version 3.0
 *
 * See the file "GEOMASON-LICENSE" for more information
 *
 */
package sim.io.geo;

import java.io.ByteArrayInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.net.URL;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import sim.field.geo.VectorLayer;
import sim.util.Bag;
import sim.util.geo.AttributeValue;
import sim.util.geo.MasonGeometry;

/**
 * A native Java importer to read ERSI shapefile data into the GeomVectorField. We assume the input
 * file follows the standard ESRI shapefile format.
 */
public class ShapeFileImporter {

  private static final Logger LOGGER = Logger.getLogger(ShapeFileImporter.class.getName());

  /**
   * Not meant to be instantiated
   */
  private ShapeFileImporter() {}

  // Shape types included in ESRI Shapefiles. Not all of these are currently
  // supported.

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

  /**
   * Populate field from the shape file given in fileName
   *
   * @param shpFile to be read from
   * @param dbFile to be read from
   * @param vectorLayer to contain read in data
   * @throws FileNotFoundException
   */
  public static void read(final URL shpFile, final URL dbFile, final VectorLayer vectorLayer)
      throws FileNotFoundException, IOException, Exception {
    read(shpFile, dbFile, vectorLayer, null, MasonGeometry.class);
  }

  public static void read(final String shpPath, final String dbPath, final VectorLayer vectorLayer)
      throws FileNotFoundException, IOException, Exception {
    read((new URI(shpPath)).toURL(), (new URI(dbPath)).toURL(), vectorLayer, null,
        MasonGeometry.class);
  }

  /**
   * Populate field from the shape file given in fileName
   *
   * @param shpFile to be read from
   * @param dbFile to be read from
   * @param vectorLayer to contain read in data
   * @param masked dictates the subset of attributes we want
   * @throws FileNotFoundException
   */
  public static void read(final URL shpFile, final URL dbFile, final VectorLayer vectorLayer,
      final Bag masked) throws FileNotFoundException, IOException, Exception {
    read(shpFile, dbFile, vectorLayer, masked, MasonGeometry.class);
  }

  public static void read(final String shpPath, final String dbPath, final VectorLayer vectorLayer,
      final Bag masked) throws FileNotFoundException, IOException, Exception {
    read((new URI(shpPath)).toURL(), (new URI(dbPath)).toURL(), vectorLayer, masked,
        MasonGeometry.class);
  }

  /**
   * Populate field from the shape file given in fileName
   *
   * @param shpFile to be read from
   * @param dbFile to be read from
   * @param vectorLayer to contain read in data
   * @param masonGeometryClass allows us to over-ride the default MasonGeometry wrapper
   * @throws FileNotFoundException
   */
  public static void read(final URL shpFile, final URL dbFile, final VectorLayer vectorLayer,
      final Class<?> masonGeometryClass) throws FileNotFoundException, IOException, Exception {
    read(shpFile, dbFile, vectorLayer, null, masonGeometryClass);
  }

  public static void read(final String shpPath, final String dbPath, final VectorLayer vectorLayer,
      final Class<?> masonGeometryClass) throws FileNotFoundException, IOException, Exception {
    read((new URI(shpPath)).toURL(), (new URI(dbPath)).toURL(), vectorLayer, null,
        masonGeometryClass);
  }

  public static void read(final Class theClass, final String shpFilePathRelativeToClass,
      final String dbFilePathRelativeToClass, final VectorLayer vectorLayer, final Bag masked,
      final Class<?> masonGeometryClass) throws IOException, Exception {
    read(theClass.getResource(shpFilePathRelativeToClass),
        theClass.getResource(dbFilePathRelativeToClass), vectorLayer, masked, masonGeometryClass);
  }

  /**
   * Reads data from the given shapefile and associated database file, populating the provided
   * VectorLayer.
   *
   * @param shpFile the URL of the shapefile to read from
   * @param dbFile the URL of the associated database file to read from
   * @param vectorLayer the VectorLayer to contain the read data
   * @param masked a Bag that dictates the subset of attributes to include, or null to include all
   *        attributes
   * @param masonGeometryClass the class of MasonGeometry or a subclass, allowing over-ride of the
   *        default MasonGeometry wrapper
   * @throws FileNotFoundException if either the shapefile or database file cannot be found
   * @throws IOException if there is an error reading the files
   * @throws Exception if there is a problem instantiating the MasonGeometry class
   * @throws IllegalArgumentException if the provided masonGeometryClass is not a MasonGeometry
   *         class or subclass
   */
  private static void read(final URL shpFile, final URL dbFile, final VectorLayer vectorLayer,
      final Bag masked, final Class<?> masonGeometryClass)
      throws FileNotFoundException, IOException, Exception {
    if (!MasonGeometry.class.isAssignableFrom(masonGeometryClass)) // Not a subclass? No go
    {
      throw new IllegalArgumentException(
          "masonGeometryClass not a MasonGeometry class or subclass");
    }

    try {
      class FieldDirEntry {
        public String name;
        public int fieldSize;
      }
      InputStream shpFileInputStream;
      InputStream dbFileInputStream;

      try {
        shpFileInputStream = DataObjectImporter.open(shpFile);
        dbFileInputStream = DataObjectImporter.open(dbFile);
      } catch (final IllegalArgumentException e) {
        LOGGER.severe("Either your shpFile or dbFile is missing!");
        throw e;
      }

      try {
        // The header size is 8 bytes in, and is little endian
        ImporterUtils.skip(dbFileInputStream, 8);

        // Both are unsigned 16-bit values: a record may be wider than 32767 bytes.
        final int headerSize = ImporterUtils.readShort(dbFileInputStream, true) & 0xFFFF;
        final int recordSize = ImporterUtils.readShort(dbFileInputStream, true) & 0xFFFF;
        final int fieldCnt = (headerSize - 1) / 32 - 1;

        final FieldDirEntry fields[] = new FieldDirEntry[fieldCnt];
        ImporterUtils.skip(dbFileInputStream, 20); // ImporterUtils.skip 20 ahead.

        final byte c[] = new byte[11];
        final char type[] = new char[fieldCnt];
        int length;

        for (int i = 0; i < fieldCnt; i++) {
          ImporterUtils.readFully(dbFileInputStream, c);
          int j = 0;
          for (j = 0; j < c.length && c[j] != 0; j++) {
            ; // ImporterUtils.skip to first unwritten byte
          }
          final String name = new String(c, 0, j);
          type[i] = (char) ImporterUtils.readByte(dbFileInputStream, true);
          fields[i] = new FieldDirEntry();
          fields[i].name = name;
          ImporterUtils.skip(dbFileInputStream, 4); // data address
          final byte b = ImporterUtils.readByte(dbFileInputStream, true);
          length = (b >= 0) ? (int) b : 256 + b; // Allow 0?
          fields[i].fieldSize = length;

          ImporterUtils.skip(dbFileInputStream, 15);
        }
        dbFileInputStream.close();
        dbFileInputStream = DataObjectImporter.open(dbFile); // Reopen for new seekin'
        ImporterUtils.skip(dbFileInputStream, headerSize); // ImporterUtils.skip the initial stuff.

        final GeometryFactory geomFactory = new GeometryFactory();

        ImporterUtils.skip(shpFileInputStream, 100);

        // Read record by record until the stream ends. available() is only an estimate, and on
        // network and jar streams it can report 0 long before the end, truncating the layer.
        final byte[] recordHeader = new byte[8];
        int skipped = 0;
        while (readRecordHeader(shpFileInputStream, recordHeader)) {
          // The content length, in 16-bit words, is big-endian. Reading the whole record first
          // keeps the stream aligned whatever the shape reader below consumes - a PointZ record,
          // for one, may carry an M value after its Z.
          final int contentLength =
              ByteBuffer.wrap(recordHeader, 4, 4).order(ByteOrder.BIG_ENDIAN).getInt() * 2;
          final byte[] content = new byte[contentLength];
          ImporterUtils.readFully(shpFileInputStream, content);
          final InputStream record = new ByteArrayInputStream(content);

          final int recordType = ImporterUtils.readInt(record, true);

          // Read the attributes; every shape record, null ones included, has a row in the .dbf
          final byte r[] = new byte[recordSize];
          ImporterUtils.readFully(dbFileInputStream, r);

          if (recordType == NULL_SHAPE) {
            // A feature with no geometry, which the format allows anywhere in a file: skipped and
            // counted.
            skipped++;
            continue;
          }

          if (!ImporterUtils.isSupported(recordType)) {
            LOGGER.severe("ShapeType " + ImporterUtils.typeToString(recordType)
                + " not supported.");
            return; // all shapes are the same type so don't bother reading any more
          }

          // Why is this start1 = 1?
          int start1 = 1;

          // Contains all the attribute values keyed by name that will eventually
          // be copied over to a corresponding MasonGeometry wrapper.
          final Map<String, AttributeValue> attributes = new HashMap<>(fieldCnt);

          for (int k = 0; k < fieldCnt; k++) {

            // If the user bothered specifying a mask and the current
            // attribute, as indexed by 'k', is NOT in the mask, then
            // merrily ImporterUtils.skip on to the next attribute
            if (masked != null && !masked.contains(fields[k].name)) {
              // But before we ImporterUtils.skip, ensure that we wind the pointer
              // to the start of the next attribute value.
              start1 += fields[k].fieldSize;

              continue;
            }
            String rawAttributeValue = new String(r, start1, fields[k].fieldSize);
            rawAttributeValue = rawAttributeValue.trim();

            final AttributeValue attributeValue = new AttributeValue();

            if (rawAttributeValue.isEmpty()) {
              // If we've gotten no data for this, then just add the
              // empty string.
              attributeValue.setString(rawAttributeValue);
            } else {
              switch (type[k]) {
                case 'N': // Numeric
                  attributeValue.setValue(parseNumber(rawAttributeValue));
                  break;
                case 'F': { // Float
                  final Object number = parseNumber(rawAttributeValue);
                  attributeValue.setValue(number instanceof Number
                      ? Double.valueOf(((Number) number).doubleValue())
                      : number);
                  break;
                }
                case 'L': // Logical
                  attributeValue.setValue(Boolean.valueOf(rawAttributeValue));
                  break;
                default:
                  attributeValue.setString(rawAttributeValue);
                  break;
              }
            }
            attributes.put(fields[k].name, attributeValue);
            start1 += fields[k].fieldSize;
          }

          // Read the shape
          Geometry geom = null;
          Coordinate pt;
          switch (recordType) {
            case POINT:
              pt = new Coordinate(ImporterUtils.readDouble(record, true),
                  ImporterUtils.readDouble(record, true));
              geom = geomFactory.createPoint(pt);
              break;
            case POINTZ:
              pt = new Coordinate(ImporterUtils.readDouble(record, true),
                  ImporterUtils.readDouble(record, true),
                  ImporterUtils.readDouble(record, true));
              geom = geomFactory.createPoint(pt);
              break;
            case POLYLINE:
            case POLYGON:
              // advance past four doubles: minX, minY, maxX, maxY
              ImporterUtils.skip(record, 32);

              final int numParts = ImporterUtils.readInt(record, true);
              final int numPoints = ImporterUtils.readInt(record, true);

              // get the array of part indices
              final int partIndicies[] = new int[numParts];
              for (int i = 0; i < numParts; i++) {
                partIndicies[i] = ImporterUtils.readInt(record, true);
              }

              // get the array of points
              final Coordinate pointsArray[] = new Coordinate[numPoints];
              for (int i = 0; i < numPoints; i++) {
                pointsArray[i] = new Coordinate(ImporterUtils.readDouble(record, true),
                    ImporterUtils.readDouble(record, true));
              }

              final Geometry[] parts = new Geometry[numParts];

              for (int i = 0; i < numParts; i++) {
                final int start = partIndicies[i];
                int end = numPoints;
                if (i < numParts - 1) {
                  end = partIndicies[i + 1];
                }
                final int size = end - start;
                final Coordinate coords[] = new Coordinate[size];

                for (int j = 0; j < size; j++) {
                  coords[j] = new Coordinate(pointsArray[start + j]);
                }

                if (recordType == ShapeFileImporter.POLYLINE) {
                  parts[i] = geomFactory.createLineString(coords);
                } else {
                  parts[i] = geomFactory.createLinearRing(coords);
                }
              }
              if (recordType == ShapeFileImporter.POLYLINE) {
                final LineString[] ls = new LineString[numParts];
                for (int i = 0; i < numParts; i++) {
                  ls[i] = (LineString) parts[i];
                }
                if (numParts == 1) {
                  geom = parts[0];
                } else {
                  geom = geomFactory.createMultiLineString(ls);
                }
              } else // polygon
              {
                geom = ImporterUtils.createPolygon(parts);
              }
              break;
            default:
              LOGGER.warning("Unknown shape type " + recordType);
          }

          if (geom != null) {
            // The user *may* have created their own MasonGeometry
            // class, so use the given masonGeometry class; by
            // default it's MasonGeometry.
            final MasonGeometry masonGeometry = (MasonGeometry) masonGeometryClass.newInstance();
            masonGeometry.geometry = geom;

            if (!attributes.isEmpty()) {
              masonGeometry.addAttributes(attributes);
            }

            vectorLayer.addGeometry(masonGeometry);
          }
        }
        if (skipped > 0) {
          LOGGER.warning(
              "Skipped " + skipped + " records with no geometry in " + shpFile.getPath());
        }
      } finally {
        dbFileInputStream.close();
        shpFileInputStream.close();
      }
    } catch (final IOException e) {
      LOGGER.severe("Could not read SHP file " + shpFile.getPath() + " with DB file "
          + dbFile.getPath());
      throw e;
    }
  }

  /**
   * Reads the next record header, or reports a clean end of file.
   *
   * @return false if the stream ended before the header started; true once the header is read.
   * @throws IOException if the stream ends part-way through the header.
   */
  private static boolean readRecordHeader(final InputStream in, final byte[] header)
      throws IOException {
    final int first = in.read();
    if (first == -1) {
      return false;
    }
    header[0] = (byte) first;
    final byte[] rest = new byte[header.length - 1];
    ImporterUtils.readFully(in, rest);
    System.arraycopy(rest, 0, header, 1, rest.length);
    return true;
  }

  /**
   * Parses a numeric .dbf value. Whole numbers are Integers where they fit and Longs where they do
   * not, so an 11-digit identifier is read whole; anything with a fraction or an
   * exponent is a Double. A value the writer could not fit in the field (dBase fills it with
   * asterisks) is kept as the raw string rather than failing the whole file.
   */
  static Object parseNumber(final String raw) {
    try {
      final long whole = Long.parseLong(raw);
      if (whole >= Integer.MIN_VALUE && whole <= Integer.MAX_VALUE) {
        return Integer.valueOf((int) whole);
      }
      return Long.valueOf(whole);
    } catch (final NumberFormatException notWhole) {
      try {
        return Double.valueOf(raw);
      } catch (final NumberFormatException notNumber) {
        return raw;
      }
    }
  }
}
