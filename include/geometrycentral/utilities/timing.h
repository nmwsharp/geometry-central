#pragma once

#include <chrono>
#include <stdio.h>
#include <string>

#define NOW (std::chrono::steady_clock::now())
#define START_TIMING(name) auto generated_timer_777_##name = NOW;

// pretty-prints to stdout
#define FINISH_TIMING_PRINT(name)                                                                                      \
  auto generated_timer_777_elapsed_##name =                                                                            \
      std::chrono::duration_cast<std::chrono::microseconds>(NOW - generated_timer_777_##name);                         \
  std::cout << "--- TIMER RESULT: section " << #name << " took "                                                       \
            << pretty_time(generated_timer_777_elapsed_##name.count()) << std::endl;

// reports as integer wth microseconds
#define FINISH_TIMING(name)                                                                                            \
  (std::chrono::duration_cast<std::chrono::microseconds>(NOW - generated_timer_777_##name).count())
// reports as double with seconds
#define FINISH_TIMING_SEC(name)                                                                                            \
  (((double)std::chrono::duration_cast<std::chrono::microseconds>(NOW - generated_timer_777_##name).count())/(1e6))

inline std::string pretty_time(long long microsec) {
  // Useful constants
  long long MILLIS = 1000;
  long long SECOND = 1000 * MILLIS;
  long long MINUTE = 60 * SECOND;
  long long HOUR = 60 * MINUTE;

  char buffer[256];

  // Greater than 1 hour
  if (microsec > HOUR) {
    long long hours = microsec / HOUR;
    microsec -= hours * HOUR;
    long long minutes = microsec / MINUTE;

    sprintf(buffer, "%lld hr, %lld min", hours, minutes);

    return std::string(buffer);
  }

  // Greater than 1 minute
  else if (microsec > MINUTE) {
    long long minutes = microsec / MINUTE;
    microsec -= minutes * MINUTE;
    long long seconds = minutes / SECOND;

    sprintf(buffer, "%lld min, %lld sec", minutes, seconds);

    return std::string(buffer);
  }

  // Greater than 1 second
  else if (microsec > SECOND) {
    double seconds = microsec / (double)SECOND;

    sprintf(buffer, "%.2f sec", seconds);

    return std::string(buffer);
  }

  // Greater than 1 millisecond
  else if (microsec > MILLIS) {
    double millis = microsec / (double)MILLIS;

    sprintf(buffer, "%.2f ms", millis);

    return std::string(buffer);
  }

  // Times less than 1 millisecond
  else {
    double micros = microsec;

    sprintf(buffer, "%.2f µs", micros); // note the nifty unicode \mu. I
                                        // apologize when this breaks something
                                        // later.

    return std::string(buffer);
  }
}
