
/*******************************************
 * Author: Zhongqiang Richard Ren. 
 * All Rights Reserved. 
 *******************************************/

#include "mospp_util.hpp"

namespace rzq{
namespace search{


std::ostream& operator<<(std::ostream& os, const MOSPPResult& res) {
  os << "MOSIPP-Result:{success:" << int(res.success) << ", n_solutions:" << res.costs.size() 
     << ", n_generated:" << res.n_generated 
     << ", n_expanded:" << res.n_expanded << ", rt_initHeu:" << res.rt_initHeu << ", rt_search:"
     << res.rt_search << "; solutions:{";
  for (auto ii : res.paths) {
    os << ",id(" << ii.first << "),cost=" << res.costs.at(ii.first) << ",nodes{";
    for (auto k : ii.second) {
      os << k << ",";
    }
    os << "},time{";
    for (auto k : res.times.at(ii.first)) {
      os << k << ",";
    }
    os << "}";
  }
  os << "}}";
  return os;
};


std::ostream& operator<<(std::ostream& os, const Label& l)
{
  std::string s;
  s = "{id:" + std::to_string(l.id) + ",v:" + std::to_string(l.v) + ",t:" 
    + std::to_string(l.t) + ",tb:" + std::to_string(l.tb) + ",g:" 
    + l.g.ToStr() + ",f:" + l.f.ToStr() + "}";
  os << s;
  return os;
};


} // end namespace search
} // end namespace rzq
