package sim.util.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;

class GeometryUtilitiesTest {

  @Test
  void euclideanDistanceMatchesPythagoras() {
    assertEquals(5.0, GeometryUtilities.euclideanDistance(new Coordinate(0.0, 0.0),
        new Coordinate(3.0, 4.0)), 1e-9);
  }

  @Test
  void euclideanDistanceIsSymmetricAndZeroOnItself() {
    Coordinate origin = new Coordinate(12.5, -7.25);
    Coordinate destination = new Coordinate(-3.0, 41.0);
    assertEquals(GeometryUtilities.euclideanDistance(origin, destination),
        GeometryUtilities.euclideanDistance(destination, origin), 1e-9);
    assertEquals(0.0, GeometryUtilities.euclideanDistance(origin, origin), 1e-9);
  }

  @Test
  @DisplayName("a single-point extent still gives an invertible transform, centred on the point")
  void singlePointExtentIsDrawable() {
    AffineTransform transform = GeometryUtilities.worldToScreenTransform(
        new Envelope(new Coordinate(5, 5)), new Rectangle2D.Double(0, 0, 100, 100));

    Point2D world = GeometryUtilities.screenToWorldPointTransform(transform, 50, 50);

    assertEquals(5.0, world.getX(), 1e-9);
    assertEquals(5.0, world.getY(), 1e-9);
  }

  @Test
  @DisplayName("points along a horizontal line give an invertible transform")
  void flatExtentIsDrawable() {
    AffineTransform transform = GeometryUtilities.worldToScreenTransform(
        new Envelope(0, 10, 3, 3), new Rectangle2D.Double(0, 0, 100, 100));

    Point2D world = GeometryUtilities.screenToWorldPointTransform(transform, 0, 50);

    assertEquals(0.0, world.getX(), 1e-9);
    assertEquals(3.0, world.getY(), 1e-9);
  }
}
