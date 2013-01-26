//
//  benchmark.h
//  Signorini
//
//  Created by Miguel de Benito Delgado on 26/01/13.
//  Copyright (c) 2013 8027. All rights reserved.
//

#ifndef __Signorini__benchmark__
#define __Signorini__benchmark__

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <ctime>

/*! Poor man's benchmarking
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

  /*
  inline bool isRunning (std::string name)
  {
    auto it = data.find (name);
    return it != data.end() && it->second.stop != 0;
  }
  */

  inline double elapsed (std::string name)
  {
    auto it = data.find (name);
    if (it != data.end() && it->second.stop != 0)
      return (double)(it->second.stop - it->second.start) / (double)CLOCKS_PER_SEC;
    return 0;
  }

  void start (std::string name)
  {
    BenchmarkData b;
    b.start = std::clock();
    b.indent = running++;
    data[name] = b;
  }
  
  clock_t stop (std::string name, bool display=true)
  {
    auto it = data.find (name);
    if (it != data.end()) {
      it->second.stop = std::clock();
      --running;
      double lapse = elapsed (name);
      if (display) {
        auto oldPrecision = std::cout.precision();
        std::cout << std::string (it->second.indent, '\t')
                  << name << " done in "
                  << std::setprecision (4) << lapse << " seconds.\n";
        std::cout.precision (oldPrecision);
      }
      return lapse;
    }
    return 0;
  }
};

#endif /* defined(__Signorini__benchmark__) */
