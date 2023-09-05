#include "timer.hxx"

#include <cassert>

namespace CarpetX {

Timer::Timer(const std::string &name)
    : name(name), handle([&]() {
        int handle1;
#pragma omp critical
        handle1 = CCTK_TimerCreate(name.c_str());
        return handle1;
      }()) {
}

void Timer::start() {
#pragma omp critical
  {
    CCTK_TimerStartI(handle);
#ifdef __CUDACC__
    range = nvtxRangeStartA(name.c_str());
#endif
  }
}

void Timer::stop() {
#pragma omp critical
  {
#ifdef __CUDACC__
    nvtxRangeEnd(range);
#endif
    CCTK_TimerStopI(handle);
  }
}

void Timer::print() const {
#pragma omp critical
  {
    const int is_running = CCTK_TimerIsRunningI(handle);
    assert(is_running >= 0);
    // We need to stop timers before printing them; running timers
    // cannot be printed
    if (is_running)
      CCTK_TimerStopI(handle);
    const int num_clocks = CCTK_NumClocks();
    for (int clock = 0; clock < num_clocks; ++clock)
      CCTK_TimerPrintDataI(handle, clock);
    if (is_running)
      CCTK_TimerStartI(handle);
  }
}

} // namespace CarpetX
