// TOOD: Look at
// <https://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes>,
// and see whether we can improve the code below

#pragma GCC optimize("O0")

// needed for dladdr, best at the top to avoid inconsistent includes
#define _GNU_SOURCE 1

#include "backtrace.hh"
#include "dist.hh"

#include <cctk.h>
#include <cctk_Parameters.h>

#include <cstring>
#include <fstream>
#include <ostream>

// These are used for backtraces. HAVE_XXX requires cctk.h
#include <sys/types.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif
#ifdef HAVE_DLADDR
#include <dlfcn.h>
#endif
#ifdef HAVE___CXA_DEMANGLE
#include <cxxabi.h>
#endif

using namespace std;

namespace CarpetLib {

// Generate a stack backtrace and send it to the given output
// stream
void generate_backtrace(ostream &stacktrace) {
  const int MAXSTACK = 100;
  static void *addresses[MAXSTACK];

  stacktrace << "Backtrace from rank " << dist::rank() << " pid " << getpid()
             << ":" << endl;

  int n = 0;
#ifdef HAVE_BACKTRACE
  n = backtrace(addresses, MAXSTACK);
#endif
  if (n < 2) {
    stacktrace << "Backtrace not available!\n";
  } else {
    auto oldflags = stacktrace.flags();
    stacktrace.setf(ios::hex);
#ifdef HAVE_BACKTRACE_SYMBOLS
    char **names = backtrace_symbols(addresses, n);
#endif
    for (int i = 2; i < n; i++) {
      char *demangled = NULL;
// Attempt to demangle this if possible
// Get the nearest symbol to feed to demangler
#ifdef HAVE_DLADDR
      Dl_info info;

      if (dladdr(addresses[i], &info) != 0) {
        int stat;
// __cxa_demangle is a naughty obscure backend and no
// self-respecting person would ever call it directly. ;-)
// However it is a convenient glibc way to demangle syms.
#ifdef HAVE___CXA_DEMANGLE
        demangled = abi::__cxa_demangle(info.dli_sname, 0, 0, &stat);
#endif
      }
#endif

      if (demangled != NULL) {

#if 0
          // Chop off the garbage from the raw symbol
          char *loc = strchr(names[i], '(');
          if (loc != NULL) *loc = '\0';
#endif

        stacktrace << i - 1 << ". " << demangled << "   ["
#ifdef HAVE_BACKTRACE_SYMBOLS
                   << names[i]
#else
		   << addresses[i]
#endif
                   << "]" << '\n';
        free(demangled);
      } else { // Just output the raw symbol
        stacktrace << i - 1 << ". "
#ifdef HAVE_BACKTRACE_SYMBOLS
                   << names[i]
#else
		   << addresses[i]
#endif
                   << '\n';
      }
    }
#ifdef HAVE_BACKTRACE_SYMBOLS
    free(names);
#endif
    stacktrace.flags(oldflags);
  }
}

// Output a stack backtrace file backtrace.<rank>.txt
void write_backtrace_file(void) {
  // declaring parameters is "safe" since it really only accesses
  // some hidden global structures, and no fragile linked list or
  // similar
  DECLARE_CCTK_PARAMETERS;

  ofstream myfile;
  stringstream ss;

  ss << out_dir << "/"
     << "backtrace." << dist::rank() << ".txt";
  string filename = ss.str();

  cerr << "Writing backtrace to " << filename << endl;
  myfile.open(filename.c_str());
  generate_backtrace(myfile);
  myfile
      << "\n"
      << "The hexadecimal addresses in this backtrace can also be interpreted\n"
      << "with a debugger (e.g. gdb), or with the 'addr2line' (or "
         "'gaddr2line')\n"
      << "command line tool: 'addr2line -e cactus_sim <address>'.\n";
  myfile.close();
}

} // namespace CarpetLib

//////////////////////////////////////////////////////////////////////////////

#include <signal.h>

namespace CarpetLib {

void signal_handler(int const signum) {
  pid_t const pid = getpid();

  cerr << "Rank " << dist::rank() << " with PID " << pid << " "
       << "received signal " << signum << endl;
  // Restore the default signal handler
  signal(signum, SIG_DFL);

  write_backtrace_file();

  // Re-raise the signal to be caught by the default handler
  kill(pid, signum);
}

void request_backtraces() {
  signal(SIGQUIT, signal_handler);
  signal(SIGILL, signal_handler);
  signal(SIGABRT, signal_handler);
  signal(SIGFPE, signal_handler);
  signal(SIGBUS, signal_handler);
  signal(SIGSEGV, signal_handler);
}

} // namespace CarpetLib

////////////////////////////////////////////////////////////////////////////////

#include <cctk.h>
#include <cctk_Arguments.h>

namespace CarpetLib {

extern "C" void CarpetLib_BacktraceTest(CCTK_ARGUMENTS) {
  CCTK_INFO("Generating backtrace...");
  kill(0, SIGABRT);
  CCTK_WARN(CCTK_WARN_ABORT, "Backtrace test failed");
}
}
