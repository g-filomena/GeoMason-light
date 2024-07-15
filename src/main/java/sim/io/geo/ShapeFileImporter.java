/*
 * Copyright 2011 by Mark Coletti, Keith Sullivan, Sean Luke, and
 * George Mason University Mason University Licensed under the Academic
 * Free License version 3.0
 *
 * See the file "GEOMASON-LICENSE" for more information
 * 
 */
package sim.io.geo;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;

import sim.field.geo.VectorLayer;
import sim.util.Bag;
import sim.util.geo.AttributeValue;
import sim.util.geo.MasonGeometry;

/**
 * A native Java importer to read ERSI shapefile data into the GeomVectorField.
 * We assume the input file follows the standard ESRI shapefile format.
 */
public class ShapeFileImporter {

	/**
	 * Not meant to be instantiated
	 */
	private ShapeFileImporter() {
	}

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
	 * @param shpFile     to be read from
	 * @param dbFile      to be read from
	 * @param vectorLayer to contain read in data
	 * @throws FileNotFoundException
	 */
	public static void read(final URL shpFile, final URL dbFile, final VectorLayer vectorLayer)
			throws FileNotFoundException, IOException, Exception {
		read(shpFile, dbFile, vectorLayer, null, MasonGeometry.class);
	}

	public static void read(final String shpPath, final String dbPath, final VectorLayer vectorLayer)
			throws FileNotFoundException, IOException, Exception {
		read((new URI(shpPath)).toURL(), (new URI(dbPath)).toURL(), vectorLayer, null, MasonGeometry.class);
	}

	/**
	 * Populate field from the shape file given in fileName
	 *
	 * @param shpFile     to be read from
	 * @param dbFile      to be read from
	 * @param vectorLayer to contain read in data
	 * @param masked      dictates the subset of attributes we want
	 * @throws FileNotFoundException
	 */
	public static void read(final URL shpFile, final URL dbFile, final VectorLayer vectorLayer, final Bag masked)
			throws FileNotFoundException, IOException, Exception {
		read(shpFile, dbFile, vectorLayer, masked, MasonGeometry.class);
	}

	public static void read(final String shpPath, final String dbPath, final VectorLayer vectorLayer, final Bag masked)
			throws FileNotFoundException, IOException, Exception {
		read((new URI(shpPath)).toURL(), (new URI(dbPath)).toURL(), vectorLayer, masked, MasonGeometry.class);
	}

	/**
	 * Populate field from the shape file given in fileName
	 *
	 * @param shpFile            to be read from
	 * @param dbFile             to be read from
	 * @param vectorLayer        to contain read in data
	 * @param masonGeometryClass allows us to over-ride the default MasonGeometry
	 *                           wrapper
	 * @throws FileNotFoundException
	 */
	public static void read(final URL shpFile, final URL dbFile, final VectorLayer vectorLayer,
			final Class<?> masonGeometryClass) throws FileNotFoundException, IOException, Exception {
		read(shpFile, dbFile, vectorLayer, null, masonGeometryClass);
	}

	public static void read(final String shpPath, final String dbPath, final VectorLayer vectorLayer,
			final Class<?> masonGeometryClass) throws FileNotFoundException, IOException, Exception {
		read((new URI(shpPath)).toURL(), (new URI(dbPath)).toURL(), vectorLayer, null, masonGeometryClass);
	}

	public static void read(final Class theClass, final String shpFilePathRelativeToClass,
			final String dbFilePathRelativeToClass, final VectorLayer vectorLayer, final Bag masked,
			final Class<?> masonGeometryClass) throws IOException, Exception {
		read(theClass.getResource(shpFilePathRelativeToClass), theClass.getResource(dbFilePathRelativeToClass),
				vectorLayer, masked, masonGeometryClass);
	}

	/**
	 * Reads data from the given shapefile and associated database file, populating
	 * the provided VectorLayer.
	 *
	 * @param shpFile            the URL of the shapefile to read from
	 * @param dbFile             the URL of the associated database file to read
	 *                           from
	 * @param vectorLayer        the VectorLayer to contain the read data
	 * @param masked             a Bag that dictates the subset of attributes to
	 *                           include, or null to include all attributes
	 * @param masonGeometryClass the class of MasonGeometry or a subclass, allowing
	 *                           over-ride of the default MasonGeometry wrapper
	 * @throws FileNotFoundException    if either the shapefile or database file
	 *                                  cannot be found
	 * @throws IOException              if there is an error reading the files
	 * @throws Exception                if there is a problem instantiating the
	 *                                  MasonGeometry class
	 * @throws IllegalArgumentException if the provided masonGeometryClass is not a
	 *                                  MasonGeometry class or subclass
	 */
	private static void read(final URL shpFile, final URL dbFile, final VectorLayer vectorLayer, final Bag masked,
			final Class<?> masonGeometryClass) throws FileNotFoundException, IOException, Exception {
		if (!MasonGeometry.class.isAssignableFrom(masonGeometryClass)) // Not a subclass? No go
		{
			throw new IllegalArgumentException("masonGeometryClass not a MasonGeometry class or subclass");
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
				System.err.println("Either your shpFile or dbFile is missing!");
				throw e;
			}

			// The header size is 8 bytes in, and is little endian
			ImporterUtils.skip(dbFileInputStream, 8);

			final int headerSize = ImporterUtils.readShort(dbFileInputStream, true);
			final int recordSize = ImporterUtils.readShort(dbFileInputStream, true);
			final int fieldCnt = (short) ((headerSize - 1) / 32 - 1);

			final FieldDirEntry fields[] = new FieldDirEntry[fieldCnt];
			ImporterUtils.skip(dbFileInputStream, 20); // ImporterUtils.skip 20 ahead.

			final byte c[] = new byte[20];
			final char type[] = new char[fieldCnt];
			int length;

			for (int i = 0; i < fieldCnt; i++) {
				dbFileInputStream.read(c, 0, 11);
				int j = 0;
				for (j = 0; j < 12 && c[j] != 0; j++)
					; // ImporterUtils.skip to first unwritten byte
				final String name = new String(c, 0, j);
				type[i] = (char) ImporterUtils.readByte(dbFileInputStream, true);
				fields[i] = new FieldDirEntry();
				fields[i].name = name;
				dbFileInputStream.read(c, 0, 4); // data address
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

			while (shpFileInputStream.available() > 0) {
				// advance past two int: recordNumber and recordLength
				// byteBuf.position(byteBuf.position() + 8);

				// byteBuf.order(ByteOrder.LITTLE_ENDIAN);
				ImporterUtils.skip(shpFileInputStream, 8);

				final int recordType = ImporterUtils.readInt(shpFileInputStream, true);

				if (!ImporterUtils.isSupported(recordType)) {
					System.err.println("Error: ShapeFileImporter.ingest(...): ShapeType "
							+ ImporterUtils.typeToString(recordType) + " not supported.");
					return; // all shapes are the same type so don't bother reading any more
				}

				// Read the attributes

				final byte r[] = new byte[recordSize];
				final int chk = dbFileInputStream.read(r);

				// Why is this start1 = 1?
				int start1 = 1;

				// Contains all the attribute values keyed by name that will eventually
				// be copied over to a corresponding MasonGeometry wrapper.
				final Map<String, AttributeValue> attributes = new HashMap<String, AttributeValue>(fieldCnt);

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
						switch (type[k]) { // Numeric case
						case 'N':
							if (rawAttributeValue.length() == 0)
								attributeValue.setString("0");
							if (rawAttributeValue.indexOf('.') != -1)
								attributeValue.setDouble(Double.valueOf(rawAttributeValue));
							else
								attributeValue.setInteger(Integer.valueOf(rawAttributeValue));
							break;
						case 'L': // Logical
							attributeValue.setValue(Boolean.valueOf(rawAttributeValue));
							break;
						case 'F': // Float
							attributeValue.setValue(Double.valueOf(rawAttributeValue));
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
					pt = new Coordinate(ImporterUtils.readDouble(shpFileInputStream, true),
							ImporterUtils.readDouble(shpFileInputStream, true));
					geom = geomFactory.createPoint(pt);
					break;
				case POINTZ:
					pt = new Coordinate(ImporterUtils.readDouble(shpFileInputStream, true),
							ImporterUtils.readDouble(shpFileInputStream, true),
							ImporterUtils.readDouble(shpFileInputStream, true));
					geom = geomFactory.createPoint(pt);
					break;
				case POLYLINE:
				case POLYGON:
					// advance past four doubles: minX, minY, maxX, maxY
					ImporterUtils.skip(shpFileInputStream, 32);

					final int numParts = ImporterUtils.readInt(shpFileInputStream, true);
					final int numPoints = ImporterUtils.readInt(shpFileInputStream, true);

					// get the array of part indices
					final int partIndicies[] = new int[numParts];
					for (int i = 0; i < numParts; i++) {
						partIndicies[i] = ImporterUtils.readInt(shpFileInputStream, true);
					}

					// get the array of points
					final Coordinate pointsArray[] = new Coordinate[numPoints];
					for (int i = 0; i < numPoints; i++) {
						pointsArray[i] = new Coordinate(ImporterUtils.readDouble(shpFileInputStream, true),
								ImporterUtils.readDouble(shpFileInputStream, true));
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

						if (recordType == ShapeFileImporter.POLYLINE)
							parts[i] = geomFactory.createLineString(coords);
						else
							parts[i] = geomFactory.createLinearRing(coords);
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
					System.err.println("Unknown shape type in " + recordType);
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
			dbFileInputStream.close();
			shpFileInputStream.close();
		} catch (final IOException e) {
			System.err.println("Error in ShapeFileImporter!!");
			System.err.println("SHP filename: " + shpFile.getPath() + "; DB filename: " + dbFile.getPath());
			throw e;
		}
	}
}