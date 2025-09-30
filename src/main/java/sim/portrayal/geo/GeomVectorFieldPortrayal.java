/*
 * Copyright 2011 by Mark Coletti, Keith Sullivan, Sean Luke, and George Mason University Mason
 * University Licensed under the Academic Free License version 3.0
 *
 * See the file "GEOMASON-LICENSE" for more information
 *
 */
package sim.portrayal.geo;

import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.util.HashMap;
import java.util.List;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;
import sim.field.geo.VectorLayer;
import sim.portrayal.DrawInfo2D;
import sim.portrayal.FieldPortrayal2D;
import sim.portrayal.LocationWrapper;
import sim.portrayal.Portrayal;
import sim.portrayal.SimplePortrayal2D;
import sim.util.Bag;
import sim.util.geo.MasonGeometry;

/**
 * Portrayal for MasonGeometry objects. The portrayal handles drawing and hit-testing (for
 * inspectors). GeomVectorFieldPortrayal overrides getPortrayalForObject to do a different thing
 * than normal FieldPortrayals. Specifically: - The object passed in is expected to be a
 * MasonGeometry. From this we extract USER, the MASON user data of the geometry, and GEOMETRY, the
 * JTS Geometry object. - If there is a portrayalForAll, return it. - If user exists and is a
 * Portrayal, return user as the Portrayal. - If a portrayal is registered for user, return it. - If
 * a portrayal is registered for geometry, return it. - If a portrayal is registered for user's
 * class, return it. - If a portrayal is registered for the geometry's class, return it. - If there
 * is a portrayalForRemainder, return it. - Else return the getDefaultPortrayal
 *
 * Note that nowhere do we return portrayals for null objects: there is no PortrayalForNull and no
 * DefaultNullPortrayal. Indeed, the method setPortrayalForNull will throw an error -- you are not
 * permitted to call it.
 */
public class GeomVectorFieldPortrayal extends FieldPortrayal2D {

  /** Throws an exception. Do not call this method. */
  @Override
  public void setPortrayalForNull(Portrayal portrayal) {
    // this bad boy throws an exception
    throw new RuntimeException(
        "setPortrayalForNull(Portrayal) may NOT be called on a GeomVectorFieldPortrayal");
  }

  /**
   * Returns the appropriate Portrayal. See the class header for more information on the
   * implementation of this method.
   */
  @Override
  public Portrayal getPortrayalForObject(Object obj) {
    // return the portrayal-for-all if any
    if (portrayalForAll != null) {
      return portrayalForAll;
    }

    MasonGeometry masonGeometry = (MasonGeometry) obj;
    Geometry geometry = masonGeometry.getGeometry();
    Object user = masonGeometry.getUserData();

    Portrayal tmp;

    // we don't check for null values of obj, so this is simpler than the
    // one in FieldPortrayal

    if (user != null && user instanceof Portrayal) {
      return (Portrayal) user;
    }
    if (portrayalForNonNull != null) {
      return portrayalForNonNull;
    }
    if (((portrayals != null /* && !portrayals.isEmpty() */) && // a little
    // efficiency
    // -- avoid
    // making
    // weak keys
    // etc.
        ((tmp = ((Portrayal) (portrayals.get(user)))) != null))
        || ((portrayals != null /* && !portrayals.isEmpty() */) && // a little
        // efficiency
        // -- avoid
        // making
        // weak keys
        // etc.
            ((tmp = ((Portrayal) (portrayals.get(geometry)))) != null))) {
      return tmp;
    }
    if (user != null && (classPortrayals != null /*
                                                  * && !classPortrayals.isEmpty ()
                                                  */) && // a little
    // efficiency --
    // avoid making weak
    // keys etc.
        ((tmp = ((Portrayal) (classPortrayals.get(user.getClass())))) != null)) {
      return tmp;
    }
    if (geometry != null && (classPortrayals != null /*
                                                      * && !classPortrayals.isEmpty ()
                                                      */) && // a little
    // efficiency --
    // avoid making
    // weak keys
    // etc.
        ((tmp = ((Portrayal) (classPortrayals.get(geometry.getClass())))) != null)) {
      return tmp;
    }
    if (portrayalForRemainder != null) {
      return portrayalForRemainder;
    }

    return getDefaultPortrayal();
  }

  /** The underlying portrayal */
  GeomPortrayal defaultPortrayal = new GeomPortrayal();

  /** Default constructor */
  public GeomVectorFieldPortrayal() {
    super();
    setImmutableField(false);
  }

  /** Constructor which sets the field's immutable flag */
  public GeomVectorFieldPortrayal(boolean immutableField) {
    super();
    setImmutableField(immutableField);
  }

  /** Return the underlying portrayal */
  @Override
  public Portrayal getDefaultPortrayal() {
    return defaultPortrayal;
  }

  /** Caches immutable fields. */
  BufferedImage buffer = null;

  RenderingHints hints = null;

  /** Handles hit-testing and drawing of the underlying geometry objects. */
  @Override
  protected void hitOrDraw(Graphics2D graphics, DrawInfo2D info, Bag putInHere) {
    if (field == null) {
      return;
    }

    // If we're drawing (and not inspecting), re-fresh the buffer if the
    // associated field is immutable.
    if (graphics != null && immutableField && !info.precise) {

      VectorLayer geomField = (VectorLayer) field;
      double x = info.clip.x;
      double y = info.clip.y;
      boolean dirty = false;

      // make a new buffer? or did the user change the zoom? Or change the
      // rendering hints?
      if (buffer == null || buffer.getWidth() != info.clip.width
          || buffer.getHeight() != info.clip.height || hints == null
          || !hints.equals(graphics.getRenderingHints())) {
        hints = graphics.getRenderingHints();
        buffer = new BufferedImage((int) info.clip.width, (int) info.clip.height,
            BufferedImage.TYPE_INT_ARGB);
        dirty = true;
      }

      // handles the case for scrolling
      if (geomField.drawX != x || geomField.drawY != y) {
        dirty = true;
      }

      // save the origin of the drawn region for later
      geomField.drawX = x;
      geomField.drawY = y;

      // re-draw into the buffer
      if (dirty) {
        clearBufferedImage(buffer);
        Graphics2D newGraphics = (Graphics2D) buffer.getGraphics();
        newGraphics.setRenderingHints(hints);
        hitOrDraw2(newGraphics, new DrawInfo2D(info, -x, -y), putInHere);
        newGraphics.dispose();
      }

      // draw buffer on screen
      graphics.drawImage(buffer, (int) x, (int) y, null);
    } else if (graphics == null) // we're just hitting
    {
      hitOrDraw2(graphics, info, putInHere);
    } else
    // might as well clear the buffer -- likely we're doing precise drawing
    {
      // do regular MASON-style drawing
      buffer = null;
      hitOrDraw2(graphics, info, putInHere);
    }
  }

  /** Clears the BufferedImage by setting all the pixels to RGB(0,0,0,0) */
  void clearBufferedImage(BufferedImage image) {
    int len = image.getHeight() * image.getWidth();
    WritableRaster raster = image.getRaster();
    int[] data = new int[len];
    for (int i = 0; i < len; i++) {
      data[i] = 0;
    }
    raster.setDataElements(0, 0, image.getWidth(), image.getHeight(), data);
  }

  /**
   * Helper function which performs the actual hit-testing and drawing for both immutable fields and
   * non-immutable fields.
   *
   * The objects in the field can either use GeomPortrayal or any SimplePortrayal2D for drawing.
   *
   */
  void hitOrDraw2(Graphics2D graphics, DrawInfo2D info, Bag putInHere) {
    VectorLayer geomField = (VectorLayer) field;
    if (geomField == null) {
      return;
    }

    boolean objectSelected = !selectedWrappers.isEmpty();

    geomField.updateTransform(info);

    List<MasonGeometry> geometries;
    geometries = geomField.queryField(geomField.clipEnvelope);

    if (geometries == null || geometries.isEmpty()) {
      // FIXME This is a hack to correct for situation where when
      // doing hit, not drawing, the incorrect geometries are returned.
      geometries = geomField.getGeometries();

      // Sometimes there really *isn't* anything to render.
      if (geometries.isEmpty()) {
        return;
      }
    }

    GeomInfo2D gInfo = new GeomInfo2D(info, geomField.worldToScreen);

    final double xScale = info.draw.width / geomField.getFieldWidth();
    final double yScale = info.draw.height / geomField.getFieldHeight();
    GeomInfo2D newinfo = new GeomInfo2D(new DrawInfo2D(info.gui, info.fieldPortrayal,
        new Rectangle2D.Double(0, 0, xScale, yScale), info.clip), geomField.worldToScreen);
    newinfo.fieldPortrayal = this;

    for (MasonGeometry masonGeometry : geometries) {

      // FIXME: *Why* the hell would this happen?
      if (masonGeometry == null) {
        continue;
      }

      Geometry geom = masonGeometry.getGeometry();

      if (geomField.clipEnvelope.intersects(geom.getEnvelopeInternal())) {
        Portrayal p = getPortrayalForObject(masonGeometry);

        if (!(p instanceof SimplePortrayal2D)) {
          throw new RuntimeException("Unexpected Portrayal " + p + " for object " + masonGeometry
              + " -- expected a SimplePortrayal2D or a GeomPortrayal");
        }

        SimplePortrayal2D portrayal = (SimplePortrayal2D) p;

        if (graphics == null) {
          if (portrayal.hitObject(masonGeometry, info)) {
            // XXX getGeometryLocation merely returns the centroid of
            // the MasonGeometry object once it finds it in the GeomVectorField;
            // however, we *just got it* from that same field, so why
            // do we need to find it again? Just directly get the
            // centroid of the object and be done with it.
            // putInHere.add(new LocationWrapper(gm, geomField.getGeometryLocation(gm), this));
            putInHere.add(new LocationWrapper(masonGeometry,
                masonGeometry.getGeometry().getCentroid(), this));
          }
        } else {
          if (portrayal instanceof GeomPortrayal) {
            portrayal.draw(masonGeometry, graphics, gInfo);
          } else { // have a SimplePortrayal2D,
            Point pt = masonGeometry.getGeometry().getCentroid();
            pt.apply(geomField.jtsTransform);
            pt.geometryChanged();

            newinfo.selected = (objectSelected && selectedWrappers.get(masonGeometry) != null);
            newinfo.draw.x = pt.getX();
            newinfo.draw.y = pt.getY();
            portrayal.draw(masonGeometry, graphics, newinfo);
          }
        }
      }
    }
  }

  /** Sets the underlying field, after ensuring its a GeomVectorField. */
  @Override
  public void setField(Object field) {
    if (field instanceof VectorLayer) {
      super.setField(field);
    } // sets dirty field already
    else {
      throw new RuntimeException("Invalid field for GeomFieldPortrayal: " + field);
    }
  }

  HashMap<Object, LocationWrapper> selectedWrappers = new HashMap<>();

  @Override
  public boolean setSelected(LocationWrapper wrapper, boolean selected) {
    if ((wrapper == null) || (wrapper.getFieldPortrayal() != this)) {
      return true;
    }

    Object obj = wrapper.getObject();
    boolean b = getPortrayalForObject(obj).setSelected(wrapper, selected);
    if (selected) {
      if (!b) {
        return false;
      }
      selectedWrappers.put(obj, wrapper);
    } else {
      selectedWrappers.remove(obj);
    }
    return true;
  }

}
