package sim.util.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import sim.graph.NodeGraph;

/** Angles are measured clockwise from the positive y axis, in degrees. */
class AnglesTest {

  private static final double TOLERANCE = 1e-9;

  private static NodeGraph node(double x, double y) {
    NodeGraph node = new NodeGraph(new Coordinate(x, y));
    node.setMasonGeometry(new MasonGeometry(sim.testing.Fixtures.FACTORY
        .createPoint(new Coordinate(x, y))));
    return node;
  }

  @Test
  @DisplayName("angle() reads clockwise from north")
  void angleFollowsCompassBearings() {
    Coordinate origin = new Coordinate(0.0, 0.0);
    assertEquals(0.0, Angles.angle(origin, new Coordinate(0.0, 100.0)), TOLERANCE);
    assertEquals(90.0, Angles.angle(origin, new Coordinate(100.0, 0.0)), TOLERANCE);
    assertEquals(180.0, Angles.angle(origin, new Coordinate(0.0, -100.0)), TOLERANCE);
    assertEquals(270.0, Angles.angle(origin, new Coordinate(-100.0, 0.0)), TOLERANCE);
  }

  @Test
  @DisplayName("the node and coordinate overloads agree")
  void nodeOverloadMatchesCoordinateOverload() {
    NodeGraph origin = node(0.0, 0.0);
    NodeGraph destination = node(300.0, 400.0);
    assertEquals(Angles.angle(origin.getCoordinate(), destination.getCoordinate()),
        Angles.angle(origin, destination), TOLERANCE);
  }

  @Test
  @DisplayName("differenceAngles() takes the short way round 0")
  void differenceAnglesWrapsAroundZero() {
    assertEquals(10.0, Angles.differenceAngles(10.0, 20.0), TOLERANCE);
    assertEquals(100.0, Angles.differenceAngles(200.0, 300.0), TOLERANCE);
    assertEquals(20.0, Angles.differenceAngles(350.0, 10.0), TOLERANCE);
    assertEquals(20.0, Angles.differenceAngles(10.0, 350.0), TOLERANCE);
    assertEquals(180.0, Angles.differenceAngles(0.0, 180.0), TOLERANCE);
  }

  @Test
  @DisplayName("isInDirection() accepts both sides of the cone")
  void isInDirectionAcceptsBothSides() {
    assertTrue(Angles.isInDirection(90.0, 80.0, 30.0));
    assertTrue(Angles.isInDirection(90.0, 100.0, 30.0));
    assertFalse(Angles.isInDirection(90.0, 140.0, 30.0));
    assertFalse(Angles.isInDirection(90.0, 40.0, 30.0));
  }

  @Test
  @DisplayName("isInDirection() accepts both sides of a cone straddling 0")
  void isInDirectionHandlesConeAcrossZero() {
    // A cone pointing due north spans, say, 344..16 degrees. Both halves belong to it: the
    // clockwise half used to be rejected, so an agent heading north only ever saw destinations
    // to its left.
    assertTrue(Angles.isInDirection(0.0, 350.0, 30.0));
    assertTrue(Angles.isInDirection(0.0, 5.0, 30.0));
    assertTrue(Angles.isInDirection(350.0, 5.0, 30.0));
    assertFalse(Angles.isInDirection(0.0, 90.0, 30.0));
    assertFalse(Angles.isInDirection(0.0, 180.0, 30.0));
  }

  @Test
  @DisplayName("getCoordAngle() is the inverse of angle()")
  void getCoordAngleInvertsAngle() {
    NodeGraph origin = node(0.0, 0.0);
    for (double bearing : new double[] {0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0}) {
      Coordinate coordinate = Angles.getCoordAngle(origin, 500.0, bearing);
      assertEquals(500.0, GeometryUtilities.euclideanDistance(origin.getCoordinate(), coordinate),
          1e-6, "distance for bearing " + bearing);
      assertEquals(bearing, Angles.angle(origin.getCoordinate(), coordinate), 1e-6,
          "bearing round trip for " + bearing);
    }
  }

  @Test
  @DisplayName("viewField() spans both nodes and has an area")
  void viewFieldCoversBothNodes() {
    NodeGraph origin = node(0.0, 0.0);
    NodeGraph destination = node(0.0, 1000.0);
    Geometry viewField = Angles.viewField(origin, destination, 70.0);

    assertTrue(viewField.getArea() > 0.0);
    assertTrue(viewField.covers(origin.getMasonGeometry().getGeometry()));
    assertTrue(viewField.covers(destination.getMasonGeometry().getGeometry()));
  }

  @Test
  @DisplayName("viewField() clamps a field of view of 180 degrees or more")
  void viewFieldClampsWideFields() {
    NodeGraph origin = node(0.0, 0.0);
    NodeGraph destination = node(0.0, 1000.0);
    assertEquals(Angles.viewField(origin, destination, 140.0).getArea(),
        Angles.viewField(origin, destination, 220.0).getArea(), 1e-6);
  }
}
