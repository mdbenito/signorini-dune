/******************************************************************************
 * benchmark.hpp                                                              *
 ******************************************************************************/

#ifndef SIGNORINI_BENCHMARK_HPP
#define SIGNORINI_BENCHMARK_HPP

#include <string>
#include <map>
#include <ctime>

/*! Poor man's benchmarking
 
 Having report() here is not especially nice, but not worth more thinking 
 either. Maybe I could attach objects to tasks and use this to correctly report
 and indent, etc.
 */
class Benchmark {
  
  typedef struct {
    clock_t start = 0;
    clock_t  stop = 0;
    int    indent = 0;  // for display
  } BenchmarkData;
  
  std::map<std::string, BenchmarkData> data;

  int running;
  
public:
  Benchmark () : running (0) { }
  
  void   start (std::string name, bool display=true);
  clock_t stop (std::string name, bool display=true);
  
  inline bool isRunning (std::string name) const;
  inline double elapsed (std::string name) const;
  
  void report (std::string name, std::string what, bool time=true) const;
};

Benchmark& bench ();  // returns singleton

#endif // defined (SIGNORINI_BENCHMARK)
