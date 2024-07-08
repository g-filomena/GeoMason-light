/*
 * Copyright 2011 by Mark Coletti, Keith Sullivan, Sean Luke, and
 * George Mason University Mason University Licensed under the Academic
 * Free License version 3.0
 *
 * See the file "GEOMASON-LICENSE" for more information
 * 
 */
package sim.util.geo;

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;

import sim.display.Display2D;
import sim.field.geo.GridLayer;
import sim.field.geo.Layer;
import sim.portrayal.DrawInfo2D;

public class GeometryUtilities {

	GeometryFactory geomFactory = new GeometryFactory();

	/**
	 * Determines the affine transform which converts world coordinates into screen
	 * coordinates. Modified from GeoTools RenderUtilities.java.
	 * 
	 * @param mapExtent MBR in world coordinates mapped to viewport defined in
	 *                  'info'
	 * @param viewport
	 * @return AffineTransform suitable for converting from world to screen
	 *         coordinates
	 */
	public static AffineTransform worldToScreenTransform(final Envelope mapExtent,
			final java.awt.geom.Rectangle2D.Double viewport) {
		double scaleX = viewport.width / mapExtent.getWidth();
		double scaleY = viewport.height / mapExtent.getHeight();

		double tx = -mapExtent.getMinX() * scaleX;
		double ty = (mapExtent.getMinY() * scaleY) + viewport.height;

		AffineTransform at = new AffineTransform(scaleX, 0.0d, 0.0d, -scaleY, tx, ty);
		AffineTransform originTranslation = AffineTransform.getTranslateInstance(viewport.x, viewport.y);
		originTranslation.concatenate(at);

		return originTranslation != null ? originTranslation : at;

	}

	public static org.locationtech.jts.geom.util.AffineTransformation getPortrayalTransform(
			final AffineTransform transform, final Layer field, final Rectangle2D.Double view) {
		AffineTransform worldToScreen = transform;

		double m[] = new double[6];
		worldToScreen.getMatrix(m);
		return new org.locationtech.jts.geom.util.AffineTransformation(m[0], m[2], m[4], m[1], m[3], m[5]);
	}

	/**
	 * Determines the affine transform which converts world coordinates into screen
	 * coordinates. Modified from GeoTools RenderUtilities.java.
	 * 
	 * convenience variant of other worldToSceenTransform()
	 * 
	 * @param mapExtent MBR in world coordinates mapped to viewport defined in
	 *                  'info'
	 * @param info      defines the viewport dimensions
	 * @return AffineTransform suitable for converting from world to screen
	 *         coordinates
	 */
	public static AffineTransform worldToScreenTransform(final Envelope mapExtent, final DrawInfo2D info) {
		return worldToScreenTransform(mapExtent, info.draw);
	}

	/**
	 * Uses the worldToScreen transform to transform the point (x,y) in screen
	 * coordinates to world coordinates.
	 */
	public static Point2D screenToWorldPointTransform(final AffineTransform worldToScreen, double x, double y) {
		// code taken from GeoTools and hacked on
		AffineTransform screenToWorld = null;
		try {
			screenToWorld = worldToScreen.createInverse();
		} catch (Exception e) {
			System.out.println(e);
			System.exit(-1);
		}

		Point2D p = new Point2D.Double();
		screenToWorld.transform(new Point2D.Double(x, y), p);
		return p;
	}

	/**
	 * Computes the MBR for the grid field in display coordinates
	 * 
	 * This is used to determine the display bounds for the grid field for
	 * Display2D.attach().
	 * 
	 * @param outer     denotes MBR that maps to display window in world coordinates
	 * @param display   is the display into which the grid will be rendered
	 * @param gridLayer for which we wish to find bounds in display coordinates
	 * 
	 * @return grid field bounds in display coordinates; will return viewport if
	 *         grid does not intersect given outer MBR
	 */
	static public java.awt.geom.Rectangle2D.Double computeBounds(final Envelope outer, final Display2D display,
			final GridLayer gridLayer) {
		Display2D.InnerDisplay2D innerDisplay = display.insideDisplay;

		// Initialize bounds to that of display viewport
		java.awt.geom.Rectangle2D.Double bounds = new java.awt.geom.Rectangle2D.Double(innerDisplay.xOffset,
				innerDisplay.yOffset, innerDisplay.width, innerDisplay.height);

		AffineTransform transform = GeometryUtilities.worldToScreenTransform(outer, bounds);

		if (isWithinBounds(outer, gridLayer)) {
			// Pretty straightforward; just translate all the corners into
			// display coordinates

			Point2D.Double srcMinPoint = new Point2D.Double(gridLayer.MBR.getMinX(), gridLayer.MBR.getMaxY());
			Point2D destMinPoint = transform.transform(srcMinPoint, null);

			Point2D.Double srcMaxPoint = new Point2D.Double(gridLayer.MBR.getMaxX(), gridLayer.MBR.getMinY());
			Point2D destMaxPoint = transform.transform(srcMaxPoint, null);

			bounds.setRect(destMinPoint.getX(), destMinPoint.getY(), destMaxPoint.getX() - destMinPoint.getX(),
					destMaxPoint.getY() - destMinPoint.getY());

		} else {
			// badness happened
			// not good if the grid isn't even within the outer MBR; this likely
			// means that 'outer' and 'gridLayer' are using different spatial
			// reference systems
			System.err.println("Warning: raster not in display");
		}
		return bounds;
	}

	/**
	 * @param outer     denotes MBR that maps to display window in world coordinates
	 * @param gridLayer for which we wish to find bounds in display coordinates
	 * 
	 * @return true iff 'gridLayer' is within, intersects, or covers 'outer', else
	 *         returns false
	 * 
	 *         Can be used as check for computeBounds()
	 */
	static public boolean isWithinBounds(final Envelope outer, final GridLayer gridLayer) {
		if (outer.contains(gridLayer.MBR) || gridLayer.MBR.contains(outer) || outer.intersects(gridLayer.MBR)) {
			return true;
		}
		return false;
	}

	/**
	 * Computes the Euclidean distance between two locations
	 *
	 * @param originCoord      the origin location;
	 * @param destinationCoord the destination;
	 */
	public static double euclideanDistance(Coordinate originCoord, Coordinate destinationCoord) {
		return Math.sqrt(
				Math.pow(originCoord.x - destinationCoord.x, 2) + Math.pow(originCoord.y - destinationCoord.y, 2));
	}

}
