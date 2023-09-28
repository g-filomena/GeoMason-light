package sim.util.geo;

import java.io.IOException;
import java.io.Writer;
import java.util.List;

/**
 * The class provides utility methods for working with CSV (Comma-Separated Values) data.
 * It allows you to write CSV lines to a Writer, handling various separators and custom quotes.
 */
public class CSVUtils {

	private static final char DEFAULT_SEPARATOR = ',';

    /**
     * Writes a CSV line to the provided Writer using the specified separator and no custom quotes.
     *
     * @param w         The Writer to write the CSV line to.
     * @param values    The list of values to be written as a CSV line.
     * @param separators The separator character to use between values.
     * @throws IOException If an I/O error occurs while writing to the Writer.
     */
    public static void writeLine(Writer w, List<String> values, char separators) throws IOException {
        writeLine(w, values, separators, ' ');
    }

    /**
     * Writes a CSV line to the provided Writer using the specified separator and no custom quotes.
     *
     * @param w         The Writer to write the CSV line to.
     * @param values    The list of values to be written as a CSV line.
     * @param separators The separator character to use between values.
     * @throws IOException If an I/O error occurs while writing to the Writer.
     */
    public static void writeLine(Writer w, List<String> values, char separators) throws IOException {
        writeLine(w, values, separators, ' ');
    }
	
    /**
     * Writes a CSV line to the provided Writer using the specified separator and custom quotes.
     *
     * @param w            The Writer to write the CSV line to.
     * @param values       The list of values to be written as a CSV line.
     * @param separators    The separator character to use between values.
     * @param customQuote  The character to use for custom quoting, or ' ' for no custom quotes.
     * @throws IOException If an I/O error occurs while writing to the Writer.
     */
    public static void writeLine(Writer w, List<String> values, char separators, char customQuote) throws IOException {

        boolean first = true;

        // Default customQuote is empty

        if (separators == ' ')
            separators = DEFAULT_SEPARATOR;

        final StringBuilder sb = new StringBuilder();
        for (final String value : values) {
            if (!first)
                sb.append(separators);
            if (customQuote == ' ')
                sb.append(followCVSformat(value));
            else
                sb.append(customQuote).append(followCVSformat(value)).append(customQuote);

            first = false;
        }
        sb.append("\n");
        w.append(sb.toString());
    }	
	
    /**
     * Helper method to format a value for CSV by handling special characters.
     *
     * @param value The value to be formatted.
     * @return The formatted value suitable for CSV.
     */
    private static String followCVSformat(String value) {
        String result = value;
        if (result.contains("\""))
            result = result.replace("\"", "\"\"");
        return result;
    }
}