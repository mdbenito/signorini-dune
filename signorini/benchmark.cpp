/******************************************************************************
 * benchmark.cpp                                                              *
 ******************************************************************************/

#include "benchmark.hpp"
#include <iostream>
#include <iomanip>

using std::string;
using std::cout;


Benchmark& bench ()
{
  static Benchmark* b = NULL;
  if (b == NULL) b = new Benchmark;
  return *b;
}


void Benchmark::start (string name, bool display)
{
  BenchmarkData b;
  b.start = std::clock();
  b.indent = running++;
  data[name] = b;
  if (display)
    cout << string (b.indent, '\t') << name << ":\n";
}


clock_t Benchmark::stop (string name, bool display)
{
  auto it = data.find (name);
  if (it != data.end()) {
    it->second.stop = std::clock();
    --running;
    double lapse = elapsed (name);
    if (display) {
      auto oldPrecision = cout.precision();
      cout << string (it->second.indent, '\t')
           << name << " done in "
           << std::setprecision (4) << lapse << " seconds.\n";
      cout.precision (oldPrecision);
    }
    return lapse;
  }
  return 0;
}


double Benchmark::elapsed (string name) const
{
  auto it = data.find (name);
  if (it != data.end())
    return (double)(std::clock() - it->second.start) / (double)CLOCKS_PER_SEC;
  return 0;
}


bool Benchmark::isRunning (string name) const
{
  auto it = data.find (name);
  return it != data.end() && it->second.stop != 0;
}

void Benchmark::report (string name, string what, bool time) const
{
  auto it = data.find (name);
  if (it != data.end()) {
    cout << string (it->second.indent, '\t') << what;
    if (time) {
      double lapse = elapsed (name);
      auto oldPrecision = cout.precision();
      cout << " (+" << std::setprecision (4) << lapse << " s)\n";
      cout.precision (oldPrecision);
    }
  } else {
    cout << string (running, '\t') << "Attempt to report on inexistent task "
         << name << ".\n";
  }
}
