/*
 * Copyright (c) 2023 Gabriele Filomena University of Liverpool, UK
 *
 * This program is free software: it can redistributed and/or modified under the terms of the GNU
 * General Public License 3.0 as published by the Free Software Foundation.
 *
 * See the file "LICENSE" for more information
 */
package sim.portrayal.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Envelope;
import sim.portrayal.DrawInfo2D;
import sim.testing.Fixtures;
import sim.util.geo.MasonGeometry;

class GeomPortrayalTest {

  @Test
  @DisplayName("drawing a line does not stop the shared portrayal filling polygons")
  void drawingALineLeavesPolygonsFilled() {
    GeomPortrayal portrayal = new GeomPortrayal(Color.RED, true);
    BufferedImage image = new BufferedImage(10, 10, BufferedImage.TYPE_INT_ARGB);
    Graphics2D graphics = image.createGraphics();
    DrawInfo2D info = new DrawInfo2D(null, null, new Rectangle2D.Double(0, 0, 10, 10),
        new Rectangle2D.Double(0, 0, 10, 10));

    portrayal.draw(Fixtures.segment(0, 0, 1, 1), graphics, info);
    portrayal.draw(new MasonGeometry(Fixtures.FACTORY.toGeometry(new Envelope(2, 8, 2, 8))),
        graphics, info);
    graphics.dispose();

    assertTrue(portrayal.filled);
    assertEquals(Color.RED.getRGB(), image.getRGB(5, 5), "the polygon's interior is painted");
  }
}
