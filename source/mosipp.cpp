
/*******************************************
 * Author: Zhongqiang Richard Ren. 
 * All Rights Reserved. 
 *******************************************/

#include "avltree.hpp"
#include "mosipp.hpp"
#include <set>
#include <memory>
#include <chrono>

namespace rzq{
namespace search{

SIPPStateSpace::SIPPStateSpace() {
  return;
};

SIPPStateSpace::~SIPPStateSpace() {
  return;
};

void SIPPStateSpace::SetGraphPtr(basic::Graph* g) {
  _graph = g;
  return;
};

void SIPPStateSpace::AddNodeCstr(long nid, long t) {
  if ( _avl_node.find(nid) == _avl_node.end() ) {
    _avl_node[nid] = basic::AVLTree<long>();
  }
  _avl_node[nid].Add(t); // add this unsafe interval.
  // std::cout << " add node cstr = " << nid << " t = " << t << std::endl;
  return;
};

void SIPPStateSpace::AddEdgeCstr(long u, long v, long t) {
  if ( _avl_edge.find(u) == _avl_edge.end() ) {
    _avl_edge[u] = std::unordered_map<long, basic::AVLTree<long> >();
  }
  if ( _avl_edge[u].find(v) == _avl_edge[u].end() ) {
    _avl_edge[u][v] = basic::AVLTree<long>();
  }
  _avl_edge[u][v].Add(t);
  return;
};

bool SIPPStateSpace::FindSafeItv(long nid, long t, long* ta, long* tb) {
  *ta = 0;
  *tb = MOSIPP_MAX_TIME;
  if (_avl_node.find(nid) == _avl_node.end()) {
    return true; // [0, Tmax]
  }
  if (_avl_node[nid].Find(t).h != 0) {
    // the input t overlaps with exact a node constraint.
    *tb = 0;
    return false;
  }
  long ub = MOSIPP_MAX_TIME;
  if ( _avl_node[nid].FindMinMore(t, &ub) ) {
    *tb = ub - 1;
  }
  long lb = 0;
  if ( _avl_node[nid].FindMaxLess(t, &lb) ) {
    *ta = lb + 1;
  }
  return true;
};

long SIPPStateSpace::GetDuration(long u, long v) {
  // the first component of the graph must be the traversal time.
  return _graph->GetCost(u,v)[0];
};


std::vector<Label> SIPPStateSpace::GetSuccLabels(long u, long v, long ta, long tb) {
  if (DEBUG_MOSIPP >= 3) {
    std::cout << "[DEBUG] >>>>>>> SIPPStateSpace::_GetSuccItv with input " << u 
      << ", " << v << ", " << ta << ", " << tb << std::endl;
  }
  std::vector<Label> out;
  std::vector<long> tas, tbs;
  _GetSuccItv(u,v,ta,tb, &tas, &tbs);
  for (size_t i = 0; i < tas.size(); i++) {
    if (DEBUG_MOSIPP >= 3) {
      std::cout << "[DEBUG] >>>>>>> SIPPStateSpace::_GetSuccItv " << tas[i] << ", " << tbs[i] << std::endl;
    }
    if (ENABLE_SAFE_CHECK_MOSIPP) {
      if (tas[i] <= ta) {
        std::cout << "[ERROR] SIPPStateSpace::GetSuccLabels tas[i]=" << tas[i] << " < ta=" << ta << std::endl;
        throw std::runtime_error( "[ERROR] SIPPStateSpace::GetSuccLabels wrong interval !" );
      }
      long d = GetDuration(u,v);
      if (tas[i] > tb + d) {
        std::cout << "[ERROR] SIPPStateSpace::GetSuccLabels tas[i]=" << tas[i] << " > tb+d=" << tb+d << std::endl;
        throw std::runtime_error( "[ERROR] SIPPStateSpace::GetSuccLabels wrong interval !" );
      }
    }
    out.emplace_back(-1, v, basic::CostVector(), basic::CostVector(), tas[i], tbs[i]);
  }
  return out;
};


bool SIPPStateSpace::ViolateNC(long u, long t) {
  if (_avl_node[u].Find(t).h != 0) {
    // the input t overlaps with exact a node constraint.
    return true;
  }
  return false;
};

bool SIPPStateSpace::ViolateEC(long u, long v, long t) {
  if (_avl_edge.find(u) == _avl_edge.end() ) {
    return false;
  }
  if (_avl_edge[u].find(v) == _avl_edge[u].end() ) {
    return false;
  }
  if (_avl_edge[u][v].Find(t).h == 0 ) {
    return false;
  }
  return true;
};



void SIPPStateSpace::_GetSuccItv(long u, long v, long ta, long tb, 
  std::vector<long>* out_tas, std::vector<long>* out_tbs) 
{
  if ( (ENABLE_SAFE_CHECK_MOSIPP) && (tb < ta) ) {
    std::cout << "[ERROR] SIPPStateSpace::_GetSuccItv tb !< ta, " << tb << " !< " << ta << std::endl;
    throw std::runtime_error( "[ERROR] SIPPStateSpace::_GetSuccItv tb !< ta !?" );
  }
  long d = GetDuration(u,v); 
  std::vector<long> unsafe_points;
  _FindAllUnsafeTimeStep(u, v, &unsafe_points);

  if (DEBUG_MOSIPP >= 3) {
    for (auto k: unsafe_points) {std::cout << "[DEBUG] >>>>>>>>>>> unsafe k = " << k << std::endl;}
  }

  long t0 = ta + d;
  long tub = (tb >= MOSIPP_MAX_TIME-d) ? MOSIPP_MAX_TIME : tb + d; // upper bound, this is impl in a way to avoid overflow.
  
  if (unsafe_points.size() == 0) {
    out_tas->push_back(t0);
    out_tbs->push_back(MOSIPP_MAX_TIME); 
    return;
  }
  
  for (size_t ii = 0; ii < unsafe_points.size(); ii++) {
    const long& tk = unsafe_points[ii];
    if (tk <= t0) {
      if (tk == t0) {
        t0++;
        if (t0 > tub) {break;} // everytime when t0 changes, need to check
      }
      continue;
    }

    out_tas->push_back(t0);
    out_tbs->push_back(tk-1); // right before the next unsafe point

    t0 = tk + 1; // note, if MOSIPP_MAX_TIME is not set properly, this may OVERFLOW !!!
    if (t0 > tub) {break;} // everytime when t0 changes, need to check
  }

  if (t0 <= tub) {
    out_tas->push_back(t0);
    out_tbs->push_back(MOSIPP_MAX_TIME); //
  }else{
    // nothing.
  }

  if (ENABLE_SAFE_CHECK_MOSIPP) {
    for (size_t ii = 0; ii < unsafe_points.size(); ii++) {
      for (size_t jj = 0; jj < out_tas->size(); jj++) {
        if (unsafe_points[ii] >= out_tas->at(jj) && unsafe_points[ii] <= out_tbs->at(jj)) {
          std::cout << "[ERROR] SIPPStateSpace::_GetSuccItv unsafe point, " << unsafe_points[ii] 
            << " is within [" << out_tas->at(jj) << "," << out_tbs->at(jj) << "]" << std::endl;
          throw std::runtime_error( "[ERROR] SIPPStateSpace::_GetSuccItv computed interval is not safe !?" );
        }
      }
    }
  }
  return;
};


void SIPPStateSpace::_FindAllUnsafeTimeStep(long u, long v, std::vector<long> *out) {

  long d = GetDuration(u,v); 
  std::vector<long> unsafe_points1;
  if ( _avl_node.find(v) != _avl_node.end() ) {
    _avl_node[v].ToSortedVector(&unsafe_points1);
  }
  std::vector<long> unsafe_points2;
  if (_avl_edge.find(u) != _avl_edge.end() ) {
    if (_avl_edge[u].find(v) != _avl_edge[u].end() ){
      _avl_edge[u][v].ToSortedVector(&unsafe_points2);
    }
  }
  size_t i1=0, i2=0;
  std::vector<long>& unsafe_points = *out;
  if (unsafe_points1.size() == 0) {
    for (auto kk : unsafe_points2) {
      unsafe_points.push_back(kk+d); // don't forget this d !!
    }
  }else if (unsafe_points2.size() == 0){
    unsafe_points = unsafe_points1;
  }else{
    while ( i1 != unsafe_points1.size() || i2 != unsafe_points2.size() ) {
      if (i1 == unsafe_points1.size() ) {
        unsafe_points.push_back(unsafe_points2[i2]+d);
        i2++;
        continue;
      }else if (i2 == unsafe_points2.size()) {
        unsafe_points.push_back(unsafe_points1[i1]);
        i1++;
        continue;
      }else{
        long usfp2 = unsafe_points2[i2]+d;
        if ( unsafe_points1[i1] < usfp2 ) {
          unsafe_points.push_back(unsafe_points1[i1]);
          i1++;
          continue;
        }else if ( unsafe_points1[i1] > usfp2 ) {
          unsafe_points.push_back(usfp2);
          i2++;
          continue;
        }else { // unsafe_points1[i1] === usfp2 
          unsafe_points.push_back(usfp2);
          i1++;
          i2++;
          continue;
        }
      }
    }
  }
  return;
};


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


FrontierLinear::FrontierLinear() {
  return;
};

FrontierLinear::~FrontierLinear() {
  return;
};

bool FrontierLinear::Check(basic::CostVector g) {
  if (c_wait.size() == 0) {
    throw std::runtime_error("[ERROR], FrontierLinear::Check, wait_cost not set!");
  }
  for (auto iter : labels) {
    if (_LabelDom(iter.second.g, g)) {return true;}
  }
  return false;
};

void FrontierLinear::Update(Label& l, std::unordered_set<long> *deleted) {
  Filter(l, deleted);
  Add(l);
  return;
};

bool FrontierLinear::Remove(Label& l) {

  for (auto iter: labels) {
    if (l.g==iter.second.g) {
      labels.erase(iter.first);
      return true;
    }
  }
  return false;
};

void FrontierLinear::Add(Label& l) {
  this->labels[l.id] = l;
  return;
};


void FrontierLinear::Filter(Label& l, std::unordered_set<long> *deleted) {
  auto newLabels = labels;
  for (auto iter : labels) {
    if (_LabelDom(l.g, iter.second.g)) {
      newLabels.erase(iter.first);
      if (deleted) {
        deleted->insert(iter.first);
      }
    }
  }
  labels = newLabels;
  return;
};


void FrontierLinear::SetWaitCost(basic::CostVector& c) {
  c_wait = c;
}

bool FrontierLinear::_LabelDom(basic::CostVector& g1, basic::CostVector& g2) {
  if (g2[0] < g1[0]) {return false;}
  long dt = (g2[0] - g1[0]);
  if (EpsDominance(g1 + c_wait*dt, g2)) {return true;}
  return false;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

MOSIPP::MOSIPP() {};

MOSIPP::~MOSIPP() {};

void MOSIPP::SetGraphPtr(basic::Graph* g) {
  _graph = g;
  _sss.SetGraphPtr(g);
};

void MOSIPP::InitHeu(long vd) {
  auto tstart = std::chrono::steady_clock::now();
  _dijks.resize(_graph->GetCostDim());
  for (size_t i = 0; i<_graph->GetCostDim(); i++) {
    _dijks[i].SetGraphPtr(_graph);
    _dijks[i].Search(vd, i);
  }
  auto tend = std::chrono::steady_clock::now();
  _res.rt_search = std::chrono::duration<double>(tend-tstart).count();
  return ;
};

int MOSIPP::Search(long vo, long vo_ta, long vo_tb, long vd, long vd_ta, long vd_tb, double time_limit) {
  // ### init heu ###
  InitHeu(vd);

  // ### init ###
  auto tstart = std::chrono::steady_clock::now(); // heuristic computation time is not counted.
  int retFlag = _InitSearch(vo, vo_ta, vo_tb, vd, vd_ta, vd_tb);
  if (retFlag == 0) {
    // this can happen for MO-CBS-TSD usage.
    _res = search::MOSPPResult();
    return 0;
  }

  if (DEBUG_MOSIPP > 0) {
    std::cout << "[DEBUG] Init, _lo = " << _lo << std::endl;
  }

  bool timeout = false;

  // ### main search loop ###
  while ( !_open.empty() ) {
    // ## select label l, lexicographic order ##
    Label l = _label[ _open.begin()->second ];
    _open.erase(_open.begin());

    if (DEBUG_MOSIPP > 0) {
      std::cout << "[DEBUG] ### Pop l = " << l << std::endl;
    }

    // ## lazy dominance check ##
    if ( _FrontierCheck(l) ) {
      if (DEBUG_MOSIPP > 1) {
        std::cout << "[DEBUG] Frontier-check, dom, cont..." << std::endl;
      }
      continue;
    }
    if ( _SolDomCheck(l.f) ) {
      if (DEBUG_MOSIPP > 1) {
        std::cout << "[DEBUG] Solution-check, dom, cont..." << std::endl;
      }
      continue;
    }
    if (l.v == vd) {
      // reach goal.
      _AddSols(l);
      continue;
    }

    _UpdateFrontier(l);
    if (DEBUG_MOSIPP > 1) {
      std::cout << "[DEBUG] ### Exp. " << l << std::endl;
    }

    auto tnow = std::chrono::steady_clock::now();
    auto time_passed = std::chrono::duration<double>(tnow-tstart).count();
    if (time_passed > time_limit) {
      if (DEBUG_MOSIPP > 1) {
        std::cout << "[DEBUG] ### MOSIPP - TIMEOUT !! " << l << std::endl;
      }
      timeout = true;
      break;
    }

    // ## expand label l ##
    _Expand(l);

  } // end while

  if (DEBUG_MOSIPP > 1) {
    std::cout << "[DEBUG] ### MOSIPP - Exit while loop" << std::endl;
  }

  // timeout, no success.
  _res.success = ! timeout;
  auto tend = std::chrono::steady_clock::now();
  _res.rt_search = std::chrono::duration<double>(tend-tstart).count();

  // ### post-process the results ###
  _PostProcRes();

  if (DEBUG_MOSIPP > 1) {
    std::cout << "[DEBUG] ### MOSIPP - Exit Search()" << std::endl;
  }

  return int(_res.success); // TODO, extend to more return flags.
};


basic::CostVector MOSIPP::_Heuristic(long v) {
  auto out = basic::CostVector(0, _graph->GetCostDim());
  for (size_t cdim = 0; cdim < out.size(); cdim++) {
    out[cdim] = _dijks[cdim].GetCost(v);
    // out[cdim] = 0;
    if (out[cdim] < 0) {
      throw std::runtime_error( "[ERROR], unavailable heuristic !?" );
    }
  }
  return out;
};

MOSPPResult MOSIPP::GetResult() const {
  return _res;
};

void MOSIPP::SetGrid(basic::GridkConn& g) {
  temp_g = g;
  SetGraphPtr(&temp_g);
};


void MOSIPP::AddNodeCstr(long nid, long t) {
  _sss.AddNodeCstr(nid, t);
};
  
void MOSIPP::AddEdgeCstr(long u, long v, long t) {
  _sss.AddEdgeCstr(u, v, t);
};

void MOSIPP::SetWaitCost(std::vector<long>& wait_cost) {
  _wait_cvec.resize(wait_cost.size());
  for (size_t ii = 0; ii < wait_cost.size(); ii++){
    _wait_cvec[ii] = wait_cost[ii];
    if ((0 == ii) && (wait_cost[0] != 1)) {
      std::cout << "[WARNING] MOSIPP::SetWaitCost, the 1st dim in the cost "
        << "vec is transition time, wait cost is set to " 
        << wait_cost[0] << " (non-unit), be careful!" << std::endl;
    }
  }
  
  return;
};

////////////// Protected //////////////

long MOSIPP::_GenLabelId() {
  return _label_id_gen++;
};

bool MOSIPP::_FrontierCheck(Label l) {

  // if l.v does not reach the goal v_d, need to augment the vector to consider time step. 
  auto ss = _L2S(l); // _State2String(_Label2State(l));
  if (_alpha.find( ss ) == _alpha.end()) {return false;}
  // return _alpha[ ss ].Check(AugVec(l));
  return _alpha[ ss ].Check(l.g);
};

bool MOSIPP::_BuildPath(long lid, std::vector<long>* path, std::vector<long>* times) {
  std::vector<long> out, out2;
  if (DEBUG_MOSIPP) {std::cout << "## _BuildPath first label = " << _label[lid] << std::endl;}
  out.push_back(_label[lid].v);
  out2.push_back(_label[lid].t);
  while( _parent.find(lid) != _parent.end() ) {
    if (DEBUG_MOSIPP) {std::cout << "## _BuildPath label = " << _label[_parent[lid]] << std::endl;}
    out.push_back(_label[_parent[lid]].v);
    out2.push_back(_label[_parent[lid]].t);
    lid = _parent[lid];
  }
  path->clear();
  times->clear();
  for (size_t i = 0; i < out.size(); i++) {
    path->push_back(out[out.size()-1-i]);
    times->push_back(out2[out2.size()-1-i]);
  }
  return true;
};


void MOSIPP::_Expand(Label& l) {
  
  _res.n_expanded++;
  std::unordered_set<long> succs = _graph->GetSuccs(l.v);

  for (const auto& u : succs) { // loop over vertices
    if (DEBUG_MOSIPP > 0) {
      std::cout << "[DEBUG] >>>> Loop auto u : succs, u = " << u << std::endl;
    }
    std::vector<long> tas, tbs;
    auto l2all = _sss.GetSuccLabels(l.v, u, l.t, l.tb);
    for (auto l2 : l2all) {
      
      if (DEBUG_MOSIPP > 0) {
        std::cout << "[DEBUG] >>>> Loop v= " << u << " for loop l' = " << l2 << std::endl;
      }

      // get wait time and cost, generate new labels
      long waitTime = l2.t - l.t - _sss.GetDuration(l.v, l2.v); // wait time
      l2.id = _GenLabelId();
      l2.g = l.g + _graph->GetCost(l.v, l2.v) + (_GetWaitCost(l.v) * waitTime);
      l2.f = l2.g + _Heuristic(l2.v);
      _label[l2.id] = l2;
      _parent[l2.id] = l.id;

      if (DEBUG_MOSIPP > 0) {
        std::cout << "[DEBUG] >>>> Loop v= " << u << " gen l' = " << l2 << std::endl;
      }
      if (_FrontierCheck(l2)) {
        if (DEBUG_MOSIPP > 1) {
          std::cout << "[DEBUG] ----- Frontier-Check, dom, cont..." << std::endl;
        }
        continue;
      }
      if (_SolDomCheck(l2.f)) {
        if (DEBUG_MOSIPP > 1) {
          std::cout << "[DEBUG] ----- Solution-Check, dom, cont..." << std::endl;
        }
        continue;
      }
      if (DEBUG_MOSIPP > 0) {
        std::cout << "[DEBUG] >>>>  Add to open..." << std::endl;
      }
      _res.n_generated++;
      _open.insert( std::make_pair(l2.f, l2.id) );

    } // end for l2
  } // end for u
  return ;
};

int MOSIPP::_InitSearch(long vo, long vo_ta, long vo_tb, long vd, long vd_ta, long vd_tb) {

  _vo = vo;
  _vd = vd;

  basic::CostVector zero_vec;
  zero_vec.resize(_graph->GetCostDim(), 0);

  // about starting label.
  if ((vo_ta == -1) && (vo_tb == -1)) {
    long tta, ttb;
    bool ret_flag = _sss.FindSafeItv(vo, 0, &tta, &ttb);
    if (!ret_flag) {
      std::cout << "[ERROR] MOSIPP::_InitSearch, fail to find start interval for input, vo=" << vo 
        << "itv[" << vo_ta << "," << vo_tb << "]" << std::endl;
      throw std::runtime_error("[ERROR] MOSIPP::_InitSearch, fail to find start interval for input !");
    }
    _lo = Label(_GenLabelId(), vo, zero_vec, _Heuristic(_vo), tta, ttb); // initial label.
  } else {
    std::cout << "[ERROR] MOSIPP::_InitSearch, unacceptable input, vo=" << vo 
      << "itv[" << vo_ta << "," << vo_tb << "]" << std::endl;
    throw std::runtime_error("[ERROR] MOSIPP::_InitSearch, unacceptable input!");
  }

  _label[_lo.id] = _lo;
  _res.n_generated++;
  _open.insert( std::make_pair(_lo.f, _lo.id) );
  return 1;
};

basic::CostVector MOSIPP::_Label2State(const Label& l) {
  auto out = basic::CostVector(0,3);
  long ta=-1, tb=-1;
  auto ret_flag = _sss.FindSafeItv(l.v, l.t, &ta, &tb);
  
  if (!ret_flag) {
    std::cout << "[ERROR] MOSIPP::_Label2State, input label " << l << " overlaps with obstacles." << std::endl;
    throw std::runtime_error( "[ERROR] MOSIPP::_Label2State, input label overlaps with obstacles !?" );
  }
  
  if (ENABLE_SAFE_CHECK_MOSIPP && tb < ta) {
    std::cout << "[ERROR] MOSIPP::_Label2State, FindSafeItv " << ta << tb << ", tb < ta !?" << std::endl;
    throw std::runtime_error( "[ERROR] MOSIPP::_Label2State, FindSafeItv, tb < ta !? !?" ); 
  }

  out[0] = l.v;
  out[1] = ta;
  out[2] = tb;
  return out;
};

std::string MOSIPP::_State2String(const basic::CostVector& g) {
  std::string out = std::to_string(g[0]) + "," + std::to_string(g[1]) + "," + std::to_string(g[2]);
  return out;
};

basic::CostVector MOSIPP::_GetWaitCost(long u) {
  return _wait_cvec;
};

void MOSIPP::_UpdateFrontier(Label l) {
  std::string s = _L2S(l);
  if (_alpha.find(s) == _alpha.end()) {
    if (DEBUG_MOSIPP > 2) {
      std::cout << "[DEBUG] new Frontier for label " << l << std::endl;
    }
    // _alpha[s] = new Frontier;
    _alpha[s] = FrontierLinear();
    _alpha[s].SetWaitCost(_wait_cvec);
  }
  _alpha[s].Update(l); //AugVec(l));

  if (DEBUG_MOSIPP > 1) {
    std::cout << "[DEBUG] ----->> UpdateF. label ids = {";
    for (auto iter : _alpha[s].labels) {
      std::cout << iter.first << ",";
    }
    std::cout << "} " << std::endl;
  }
};

bool MOSIPP::_SolDomCheck(basic::CostVector& f) {
  for (auto iter : _sol_label_ids) {
    if (EpsDominance(_label[iter].g, f)) {
      return true;
    }
  }
  return false;
};

void MOSIPP::_AddSols(Label& l) {

  // filter sol, not necessary in this impl.
  std::vector<long> new_sol_id;
  for (auto lid : _sol_label_ids){
    if (EpsDominance(l.g, _label[lid].g)) {
      continue;
    }else{
      new_sol_id.push_back(lid);
    }
  }
  if (new_sol_id.size() < _sol_label_ids.size()) {
    _sol_label_ids = new_sol_id;
  }

  // add sol
  _sol_label_ids.push_back(l.id);
  return ;
};


void MOSIPP::_PostProcRes() {
  // ### post-process the results ###
  if (_sol_label_ids.size() > 0) {
    for (auto lid : _sol_label_ids) {
      _res.paths[lid] = std::vector<long>();
      _res.times[lid] = std::vector<long>();
      bool ret_flag = _BuildPath(lid, &(_res.paths[lid]), &(_res.times[lid]) );
        // when this flag is used, remember to check here for safety.
      _res.costs[lid] = _label[lid].g;
    }
  }
  return ; 
};

/////////////////////////////////////////////////////////////////
////////////////// RunMOSIPP /////////////////////
/////////////////////////////////////////////////////////////////

int RunMOSIPPGrid(
  basic::GridkConn& g, long vo, long vd, double time_limit, 
  basic::CostVector& wait_cost, std::vector< std::vector<long> >& ncs, 
  std::vector< std::vector<long> >& ecs, MOSPPResult* res) 
{
  MOSIPP mosipp;
  mosipp.SetGrid(g); // polymorphism
  mosipp.SetWaitCost(wait_cost);
  for (auto nc: ncs) {
    mosipp.AddNodeCstr(nc[0], nc[1]);
  }
  for (auto ec: ecs) {
    mosipp.AddEdgeCstr(ec[0], ec[1], ec[2]);
  }
  int outFlag = mosipp.Search(vo, -1, -1, vd, -1, -1, time_limit);
  *res = mosipp.GetResult();
  return outFlag;
};

} // end namespace search
} // end namespace rzq
