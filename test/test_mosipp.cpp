
/*******************************************
 * Author: Zhongqiang Richard Ren. 
 * All Rights Reserved. 
 *******************************************/

#include "mosipp.hpp"
#include <iostream>


int ToyEg();

int main() {
  ToyEg();
  return 0;
};


int ToyEg() {

  // static environment
  rzq::basic::Grid static_world; // static obstacles appears in this 2d grid.
  int r = 3; // rows (y)
  int c = 3; // columns (x)
  static_world.Resize(r,c);
  static_world.Set(1,1,1); // set grid[y=1,x=1] = 1, a static obstacle.

  // declare cost structure.
  rzq::basic::Grid cost1;
  cost1.Resize(r,c,1); // each action takes unit time.
  rzq::basic::Grid cost2; // communication quality
  cost2.Resize(r,c,2); // the communication quality at each location is by default 2.
  cost2.Set(1,0,1); // the lower left passage has better comm quality (left comm cost).
  cost2.Set(2,0,1);
  cost2.Set(2,1,1);
  std::vector<rzq::basic::Grid> cost_grids; // cost vectors (implemented at nodes rather than edges), arrival cost at a node.
  cost_grids.push_back(cost1); // cost_grids[0] = e.g. traversal time cost, 
  cost_grids.push_back(cost2); // cost_grids[1] = e.g. communication quality cost, ...

  // workspace graph (static)
  rzq::basic::GridkConn g; // this is the graph (impl as a grid) that represents the workspace
  g.Init(static_world, cost_grids);

  // dynamic obstacles along known trajectories are represented as node/edge constraints
  std::vector< std::vector<long> > node_constraints;
  // each element must be of length two in format: [node_id, time_step]. 
  // node_id = y*num_cols + x; time_step = the time step when the obstacle appears.
  node_constraints.emplace_back( std::vector<long>({6,1}) ); // node id = 6 (i.e. y=2,x=0), t=1
  node_constraints.emplace_back( std::vector<long>({6,2}) ); // t=2
  node_constraints.emplace_back( std::vector<long>({6,3}) ); // 
  // node_constraints.emplace_back( std::vector<long>({6,4}) ); // node id =6 is blocked for a few consecutive time steps.

  std::vector< std::vector<long> > edge_constraints; // This is not massively tested yet... place holder...

  long vo = 0; // start node id.
  long vd = 8; // goal node id.
  double time_limit = 999999; // seconds
  rzq::basic::CostVector wait_cost(1,2); // all values=1, length=2.
  rzq::search::MOSPPResult res;
  rzq::search::RunMOSIPPGrid(g, vo, vd, time_limit, wait_cost, node_constraints, edge_constraints, &res);
  std::cout << res << std::endl;

  // print paths, times and costs
  std::cout << " reprint solutions for more clarity:" << std::endl;
  for (auto iter : res.paths) {
    long k = iter.first; // id of a Pareto-optipmal solution
    // path nodes
    std::cout << " path nodes = ";
    for (auto xx : res.paths[k]) {
      std::cout << xx << ", ";
    }
    std::cout << std::endl;
    // times
    std::cout << " times = ";
    for (auto xx : res.times[k]) {
      std::cout << xx << ", ";
    }
    std::cout << std::endl;
    // cost
    std::cout << " cost = " << res.costs[k] << std::endl;
  }
  return 1;
};