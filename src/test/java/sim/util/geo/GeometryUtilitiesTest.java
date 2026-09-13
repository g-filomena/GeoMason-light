package sim.util.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;

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
}
