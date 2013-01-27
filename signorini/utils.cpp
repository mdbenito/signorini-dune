/******************************************************************************
 * utils.cpp                                                                  *
 ******************************************************************************/

#include "utils.hpp"

using std::string;

string operator+ (const string& a, const string& b)
{
  return string(a).append(b);
}
