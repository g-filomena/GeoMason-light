package sim.util.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class AttributeValueTest {

  @Test
  @DisplayName("numeric getters convert across number types")
  void numericGettersConvertAcrossTypes() {
    // Importers hand back whichever numeric type the source format used; asking for an int should
    // not depend on that choice.
    assertEquals(Integer.valueOf(3), new AttributeValue(3.7).getInteger());
    assertEquals(Integer.valueOf(3), new AttributeValue(3L).getInteger());
    assertEquals(Integer.valueOf(3), new AttributeValue(3).getInteger());
    assertEquals(Double.valueOf(3.0), new AttributeValue(3).getDouble());
    assertEquals(Double.valueOf(3.5), new AttributeValue(3.5f).getDouble());
  }

  @Test
  @DisplayName("numeric getters return null for an absent value")
  void numericGettersTolerateNull() {
    assertNull(new AttributeValue().getInteger());
    assertNull(new AttributeValue().getDouble());
    assertNull(new AttributeValue().getString());
  }

  @Test
  void settersAndGettersRoundTrip() {
    AttributeValue attribute = new AttributeValue();

    attribute.setInteger(7);
    assertEquals(Integer.valueOf(7), attribute.getInteger());

    attribute.setDouble(1.25);
    assertEquals(Double.valueOf(1.25), attribute.getDouble());

    attribute.setBoolean(true);
    assertTrue(attribute.getBoolean());

    attribute.setString("residential");
    assertEquals("residential", attribute.getString());

    attribute.setArray(new ArrayList<>(Arrays.asList("a", "b")));
    assertEquals(Arrays.asList("a", "b"), attribute.getArray());
  }

  @Test
  @DisplayName("equals() and hashCode() consider value and visibility")
  void equalityConsidersValueAndHiddenFlag() {
    AttributeValue visible = new AttributeValue("centre", false);
    AttributeValue sameValue = new AttributeValue("centre", false);
    AttributeValue hidden = new AttributeValue("centre", true);
    AttributeValue otherValue = new AttributeValue("periphery", false);

    assertEquals(visible, sameValue);
    assertEquals(visible.hashCode(), sameValue.hashCode());
    assertFalse(visible.equals(hidden));
    assertFalse(visible.equals(otherValue));
    assertFalse(visible.equals(null));
    assertFalse(visible.equals("centre"));
  }

  @Test
  void cloneCopiesValueAndVisibility() {
    AttributeValue original = new AttributeValue(42, true);
    AttributeValue copy = (AttributeValue) original.clone();

    assertEquals(original, copy);
    copy.setValue(43);
    assertEquals(Integer.valueOf(42), original.getInteger());
  }

  @Test
  void toStringReportsValueAndVisibility() {
    assertEquals("Value: 5 Hidden: false", new AttributeValue(5).toString());
  }

  @Test
  void hiddenFlagIsMutable() {
    AttributeValue attribute = new AttributeValue("x");
    assertFalse(attribute.isHidden());
    attribute.setHidden(true);
    assertTrue(attribute.isHidden());
  }
}
