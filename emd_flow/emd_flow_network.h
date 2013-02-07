#ifndef __EMD_FLOW_NETWORK_H__
#define __EMD_FLOW_NETWORK_H__

#include <vector>
#include <lemon/list_graph.h>
#include <lemon/maps.h>
#include <lemon/capacity_scaling.h>

class EMDFlowNetwork {
 public:
  EMDFlowNetwork(const std::vector<std::vector<double> >& amplitudes,
      int k = 0);
  void set_sparsity(int k);
  void run_flow(double lambda);
  int get_EMD_used();
  double get_supported_amplitude_sum();
  void get_support(std::vector<std::vector<bool> >* support);
  int get_num_nodes();
  int get_num_edges();
  int get_num_columns();
  int get_num_rows();
  ~EMDFlowNetwork();

 private:
  typedef lemon::CapacityScaling<lemon::ListDigraph, int, double> AlgType;

  std::vector<std::vector<double> > a_;
  int k_;
  // number of rows
  int r_;
  // number of columns
  int c_;

  // nodes corresponding to the matrix entries
  std::vector<std::vector<lemon::ListDigraph::Node> > innode_;
  std::vector<std::vector<lemon::ListDigraph::Node> > outnode_;
  // arcs from innodes to outnodes
  std::vector<std::vector<lemon::ListDigraph::Arc> > nodearcs_;
  // arcs corresponding to the EMD
  std::vector<std::vector<std::vector<lemon::ListDigraph::Arc> > > colarcs_;
  // graph
  lemon::ListDigraph g_;
  lemon::ListDigraph::ArcMap<int> capacity_;
  lemon::ListDigraph::ArcMap<double> cost_;
  // source
  lemon::ListDigraph::Node s_;
  // sink
  lemon::ListDigraph::Node t_;

  // algorithm
  AlgType* alg_;

  void apply_lambda(double lambda);
  int extract_emd_cost();
  double extract_amp_sum();
};

#endif
