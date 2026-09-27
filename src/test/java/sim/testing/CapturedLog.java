/*
 * Copyright (c) 2023 Gabriele Filomena University of Liverpool, UK
 *
 * This program is free software: it can redistributed and/or modified under the terms of the GNU
 * General Public License 3.0 as published by the Free Software Foundation.
 *
 * See the file "LICENSE" for more information
 */
package sim.testing;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Handler;
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * Collects the messages a class logs through {@code java.util.logging} while it is open, for a
 * test to assert on. Use in a try-with-resources block.
 */
public final class CapturedLog extends Handler implements AutoCloseable {

  // Held so the logger, which java.util.logging references weakly, outlives the capture.
  private final Logger logger;
  private final List<LogRecord> records = new ArrayList<>();

  /**
   * Starts capturing the logger named after {@code source}.
   *
   * @param source the class whose logger to capture.
   */
  public CapturedLog(Class<?> source) {
    logger = Logger.getLogger(source.getName());
    logger.addHandler(this);
  }

  /**
   * The messages logged so far, in order.
   *
   * @return the messages.
   */
  public synchronized List<String> messages() {
    return records.stream().map(LogRecord::getMessage).collect(Collectors.toList());
  }

  @Override
  public synchronized void publish(LogRecord record) {
    records.add(record);
  }

  @Override
  public void flush() {}

  @Override
  public void close() {
    logger.removeHandler(this);
  }
}
