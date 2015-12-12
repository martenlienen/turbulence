#include "MultiTimer.h"

#include <sstream>

void MultiTimer::start(std::string timer) {
  gettimeofday(&this->starts[timer], NULL);
}

void MultiTimer::stop(std::string timer) {
  struct timeval current;
  gettimeofday(&current, NULL);

  struct timeval diff = this->difference(current, this->starts[timer]);

  struct timeval& acc = this->times[timer];
  acc.tv_sec += diff.tv_sec;
  acc.tv_usec += diff.tv_usec;
}

std::string MultiTimer::toString() {
  std::stringstream out;

  for (auto& entry : this->times) {
    double seconds = entry.second.tv_sec + entry.second.tv_usec / 1e6;

    out << entry.first << " " << seconds << std::endl;
  }

  return out.str();
}

MultiTimer* MultiTimer::instance = NULL;

MultiTimer* MultiTimer::get() {
  if (instance == NULL) {
    instance = new MultiTimer();
  }

  return instance;
}

struct timeval MultiTimer::difference(struct timeval a, struct timeval b) {
  struct timeval diff;
  diff.tv_sec = a.tv_sec - b.tv_sec;
  diff.tv_usec = a.tv_usec - b.tv_usec;

  return diff;
}
