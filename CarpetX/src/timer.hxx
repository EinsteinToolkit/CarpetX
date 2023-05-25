#ifndef CARPETX_CARPETX_TIMER_HXX
#define CARPETX_CARPETX_TIMER_HXX

#include <cctk.h>

#ifdef __CUDACC__
#include <nvToolsExt.h>
#endif

#include <string>

namespace CarpetX {
using namespace std;

class Interval;

class Timer {
  friend class Interval;

  string name;
  int handle;

public:
  Timer() = delete;
  Timer(const Timer &) = delete;
  Timer &operator=(const Timer &) = delete;
  Timer(Timer &&) = default;
  Timer &operator=(Timer &&) = default;

  Timer(const string &name);

  string get_name() const { return name; }
  void print() const;
};

class Interval {
  const Timer &timer;
#ifdef __CUDACC__
  nvtxRangeId_t range;
#endif

public:
  Interval() = delete;
  Interval(const Timer &timer);
  ~Interval();
};

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_TIMER_HXX
