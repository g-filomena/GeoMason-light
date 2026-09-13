package sim.util.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.io.StringWriter;
import java.util.Arrays;
import java.util.Collections;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class CSVUtilsTest {

  @Test
  void writesCommaSeparatedValuesByDefault() throws IOException {
    StringWriter writer = new StringWriter();
    CSVUtils.writeLine(writer, Arrays.asList("nodeID", "x", "y"));
    assertEquals("nodeID,x,y\n", writer.toString());
  }

  @Test
  void honoursACustomSeparator() throws IOException {
    StringWriter writer = new StringWriter();
    CSVUtils.writeLine(writer, Arrays.asList("a", "b"), ';');
    assertEquals("a;b\n", writer.toString());
  }

  @Test
  @DisplayName("a blank separator falls back to a comma")
  void blankSeparatorFallsBackToComma() throws IOException {
    StringWriter writer = new StringWriter();
    CSVUtils.writeLine(writer, Arrays.asList("a", "b"), ' ');
    assertEquals("a,b\n", writer.toString());
  }

  @Test
  void wrapsEachValueInTheCustomQuote() throws IOException {
    StringWriter writer = new StringWriter();
    CSVUtils.writeLine(writer, Arrays.asList("a", "b"), ',', '"');
    assertEquals("\"a\",\"b\"\n", writer.toString());
  }

  @Test
  @DisplayName("embedded double quotes are doubled")
  void doublesEmbeddedQuotes() throws IOException {
    StringWriter writer = new StringWriter();
    CSVUtils.writeLine(writer, Collections.singletonList("say \"hi\""));
    assertEquals("say \"\"hi\"\"\n", writer.toString());
  }

  @Test
  void writesAnEmptyLineForNoValues() throws IOException {
    StringWriter writer = new StringWriter();
    CSVUtils.writeLine(writer, Collections.<String>emptyList());
    assertEquals("\n", writer.toString());
  }
}
