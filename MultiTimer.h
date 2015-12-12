#ifndef _MULTI_TIMER_H_
#define _MULTI_TIMER_H_

#include <memory>
#include <string>
#include <unordered_map>

#include <sys/time.h>

class MultiTimer {
 public:
  void start(std::string timer);
  void stop(std::string timer);

  std::string toString();

  static MultiTimer* get();

 private:
  // The accumulated time of each timer
  std::unordered_map<std::string, struct timeval> times;

  // The starting times of currently running timers
  std::unordered_map<std::string, struct timeval> starts;

  static MultiTimer* instance;

  struct timeval difference(struct timeval a, struct timeval b);
};

#endif // _MULTI_TIMER_H_
