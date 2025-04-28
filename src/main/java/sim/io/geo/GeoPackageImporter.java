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
	 * Reads GeoPackage data and populates the provided VectorLayer with geometries and attributes.
	 *
	 * @param gpkgURL     The URL to the GeoPackage file.
	 * @param vectorLayer The VectorLayer to populate with data.
	 * @throws Exception If an error occurs during GeoPackage reading or processing.
	 */
	public static void read(URL gpkgURL, VectorLayer vectorLayer) throws Exception {

		File file;

		if ("file".equals(gpkgURL.getProtocol())) {
			// Normal filesystem
			file = new File(gpkgURL.toURI());
		} else if ("jar".equals(gpkgURL.getProtocol())) {
			// Inside a jar: extract to temp file
			InputStream input = gpkgURL.openStream();
			file = File.createTempFile("temp_", ".gpkg");
			file.deleteOnExit();
			Files.copy(input, file.toPath(), StandardCopyOption.REPLACE_EXISTING);
		} else {
			throw new IllegalArgumentException("Unsupported URL protocol: " + gpkgURL.getProtocol());
		}

		GeoPackage geoPackage = GeoPackageManager.open(file);
		// Feature and tile tables
		List<String> features = geoPackage.getFeatureTables();
		// Iterate through feature tables
		for (String featureTableName : features) {
			FeatureDao featureDao = geoPackage.getFeatureDao(featureTableName);

			// Iterate through features
			for (FeatureRow row : featureDao.queryForAll()) {
				// Parse geometry using GeoPackage-Java's GeometryReader
				GeoPackageGeometryData geometryData = row.getGeometry();

				Geometry sfGeometry = null;
				if (geometryData != null && !geometryData.isEmpty())
					sfGeometry = geometryData.getGeometry();

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

		geoPackage.close();
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
			if (longValue >= Integer.MIN_VALUE && longValue <= Integer.MAX_VALUE)
				return new AttributeValue((int) longValue); // Fits in Integer range
			else
				return new AttributeValue(longValue); // Store as Long

		} else if (value instanceof Integer)
			return new AttributeValue((Integer) value);
		else if (value instanceof Double)
			return new AttributeValue((Double) value);
		else if (value instanceof Boolean)
			return new AttributeValue((Boolean) value);
		else
			return new AttributeValue(value);

	}

	/**
	 * Determines the type of a string value (double, integer, boolean, or string).
	 *
	 * @param rawAttributeValue The raw string value to analyze.
	 * @return A string representing the determined type.
	 */
	private static String determineType(String rawAttributeValue) {
		if (rawAttributeValue.matches("^-?\\d+\\.\\d+$"))
			return "double";
		else if (rawAttributeValue.matches("^-?\\d+$"))
			return "integer";
		else if (rawAttributeValue.equalsIgnoreCase("true") || rawAttributeValue.equalsIgnoreCase("false"))
			return "boolean";
		else
			return "string";
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
			if (mls.getLineStrings().size() == 1) {
				mil.nga.sf.LineString singleLineString = mls.getLineStrings().get(0);
				return new GeometryFactory().createLineString(singleLineString.getPoints().stream()
						.map(point -> new Coordinate(point.getX(), point.getY())).toArray(Coordinate[]::new));
			}
		} else if (geometryType.equals(GeometryType.POLYGON)) {
			mil.nga.sf.Polygon pg = (mil.nga.sf.Polygon) sfGeometry;
			return new GeometryFactory().createPolygon(pg.getExteriorRing().getPoints().stream()
					.map(point -> new Coordinate(point.getX(), point.getY())).toArray(Coordinate[]::new));
		} else if (geometryType.equals(GeometryType.MULTIPOLYGON)) {
			mil.nga.sf.MultiPolygon mpg = (mil.nga.sf.MultiPolygon) sfGeometry;
			if (mpg.getPolygons().size() == 1) {
				mil.nga.sf.Polygon singlePolygon = mpg.getPolygons().get(0);
				return new GeometryFactory().createPolygon(singlePolygon.getExteriorRing().getPoints().stream()
						.map(point -> new Coordinate(point.getX(), point.getY())).toArray(Coordinate[]::new));
			}
		} else
			throw new IllegalArgumentException("Unsupported GeoPackage geometry type: " + sfGeometry.getGeometryType());

		return null;
	}
}
