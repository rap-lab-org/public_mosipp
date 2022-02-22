
/*******************************************
 * Author: Zhongqiang Richard Ren. 
 * All Rights Reserved. 
 *******************************************/


#ifndef ZHONGQIANGREN_SEARCH_MOSPPUTIL_H_
#define ZHONGQIANGREN_SEARCH_MOSPPUTIL_H_

// #include <unordered_map>
// #include <vector>
// #include <iostream>
#include "graph.hpp"

namespace rzq{
namespace search{

/**
 * @brief Epsilon-dominance.
 */
template<typename IterData>
bool EpsDominance(IterData v1, IterData v2, double eps=0.0, bool less=true) {
  auto i2 = v2.begin();
  for (auto i1 = v1.begin(); i1 != v1.end(); i1++){
    if (less) {
      if (*i1 > (1.0+eps)*(*i2)) { return false; }
    }else{
      if (*i1 < (1.0+eps)*(*i2)) { return false; }
    }
    i2++;
  }
  return true;
};

/**
 * @brief result data structure.
 */
struct MOSPPResult {
  bool success = false;
  std::unordered_map< long, std::vector<long> > paths;
  std::unordered_map< long, std::vector<long> > times;
  std::unordered_map< long, basic::CostVector > costs;
  long n_generated = 0;
  long n_expanded = 0;
  double rt_initHeu = 0.0;
  double rt_search = 0.0;
};

/**
 *
 */
std::ostream& operator<<(std::ostream& os, const MOSPPResult& res) ;


/**
 * @brief A search label, which identifies a partial solution path.
 */
struct Label {
  Label() {};
  Label(long id0, long v0, basic::CostVector g0, basic::CostVector f0, long t0, long tb0) {
    id = id0; v = v0; g = g0; f = f0; t = t0; tb = tb0;
  };
  long id; // label's id, make it easy to look up.
  long v;
  basic::CostVector g;
  basic::CostVector f;
  // specific for (MO)-SIPP.
  long t; // the arrival time at a state.
  long tb; // the ending time of the safe interval in a state.
};

/**
 * @brief for cout.
 */
std::ostream& operator<<(std::ostream& os, const Label& l) ;

} // end namespace search
} // end namespace rzq


#endif  // ZHONGQIANGREN_SEARCH_MOSPPUTIL_H_
