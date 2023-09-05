#ifndef CARPETX_CARPETX_TIMER_HXX
#define CARPETX_CARPETX_TIMER_HXX

#include <cctk.h>

#ifdef __CUDACC__
#include <nvToolsExt.h>
#endif

#include <string>

namespace CarpetX {

class Interval;

class Timer {
public:
  std::string name;

private:
  int handle;
#ifdef __CUDACC__
  nvtxRangeId_t range;
#endif

public:
  Timer() = delete;
  Timer(const Timer &) = delete;
  Timer &operator=(const Timer &) = delete;
  Timer(Timer &&) = default;
  Timer &operator=(Timer &&) = default;

  Timer(const std::string &name);

  void start();
  void stop();

  void print() const;
};

class Interval {
  Timer &timer;

public:
  Interval() = delete;
  Interval(const Interval &) = delete;
  Interval(Interval &&) = delete;
  Interval &operator=(const Interval &) = delete;
  Interval &operator=(Interval &&) = delete;

  Interval(Timer &timer) : timer(timer) { timer.start(); }
  ~Interval() { timer.stop(); }
};

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_TIMER_HXX
