#ifndef __EMD_FLOW_NETWORK_LEMON_H__
#define __EMD_FLOW_NETWORK_LEMON_H__

#include "emd_flow_network.h"

#include <cmath>
#include <cstdio>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

template <typename MCMFAlgorithm>
class EMDFlowNetworkLemon : public EMDFlowNetwork {
 public:
  EMDFlowNetworkLemon(const std::vector<std::vector<double> >& amplitudes) 
      : EMDFlowNetwork(), capacity_(g_), cost_(g_) {
    construct_graph(amplitudes);
  }

  void set_sparsity(int k) {
    k_ = k;
    alg_->stSupply(s_, t_, k_);
  }

  void run_flow(double lambda) {
    apply_lambda(lambda);
    alg_->costMap(cost_);
    alg_->run();
  }

  int get_EMD_used() {
    return extract_emd_cost();
  }

  double get_supported_amplitude_sum() {
    return extract_amp_sum();
  }

  void get_support(std::vector<std::vector<bool> >* support) {
    support->resize(r_);
    for (int row = 0; row < r_; ++row) {
      (*support)[row].resize(c_);
      for (int col = 0; col < c_; ++col) {
        (*support)[row][col] = (alg_->flow(nodearcs_[row][col]) > 0);
      }
    }
  }

  int get_num_nodes() {
    return countNodes(g_);
  }

  int get_num_edges() {
    return countArcs(g_);
  }

  int get_num_columns() {
    return c_;
  }

  int get_num_rows() {
    return r_;
  }

  ~EMDFlowNetworkLemon() {
    delete alg_;
  }

 private:
  // amplitudes
  std::vector<std::vector<double> > a_;
  // sparsity
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
  MCMFAlgorithm* alg_;

  void apply_lambda(double lambda) {
    for (int row = 0; row < r_; ++row) {
      for (int col = 0; col < c_ - 1; ++col) {
        for (int dest = 0; dest < r_; ++dest) {
          cost_[colarcs_[row][col][dest]] = lambda * abs(row - dest);
        }
      }
    }
  }

  int extract_emd_cost() {
    int emd_cost = 0;
    for (int row = 0; row < r_; ++row) {
      for (int col = 0; col < c_ - 1; ++col) {
        for (int dest = 0; dest < r_; ++dest) {
          if (alg_->flow(colarcs_[row][col][dest]) > 0) {
            emd_cost += abs(row - dest);
            if (alg_->flow(colarcs_[row][col][dest]) != 1) {
              fprintf(stderr, "ERROR: nonzero flow on a column edge is not "
                  "1.\n");
            }
          }
        }
      }
    }
    return emd_cost;
  }

  double extract_amp_sum() {
    double amp_sum = 0;
    for (int row = 0; row < r_; ++row) {
      for (int col = 0; col < c_; ++col) {
        if (alg_->flow(nodearcs_[row][col]) > 0) {
          amp_sum += abs(a_[row][col]);
          if (alg_->flow(nodearcs_[row][col]) != 1) {
            fprintf(stderr, "ERROR: nonzero flow on a node edge is not 1.\n");
          }
        }
      }
    }
    return amp_sum;
  }

  void construct_graph(const std::vector<std::vector<double> >& amplitudes) {
    // number of rows
    r_ = amplitudes.size();
    // number of columns
    c_ = amplitudes[0].size();

    a_.resize(r_);
    for (int row = 0; row < r_; ++row) {
      a_[row].resize(c_);
      for (int col = 0; col < c_; ++col) {
        a_[row][col] = amplitudes[row][col];
      }
    }

    // source and sink
    s_ = g_.addNode();
    t_ = g_.addNode();

    // create two nodes for each entry in the matrix
    innode_.resize(r_);
    outnode_.resize(r_);
    for (int ii = 0; ii < r_; ++ii) {
      innode_[ii].resize(c_);
      outnode_[ii].resize(c_);
      for (int jj = 0; jj < c_; ++jj) {
        innode_[ii][jj] = g_.addNode();
        outnode_[ii][jj] = g_.addNode();
      }
    }

    // add arcs from innodes to outnodes
    nodearcs_.resize(r_);
    for (int ii = 0; ii < r_; ++ii) {
      nodearcs_[ii].resize(c_);
      for (int jj = 0; jj < c_; ++jj) {
        nodearcs_[ii][jj] = g_.addArc(innode_[ii][jj], outnode_[ii][jj]);
        cost_[nodearcs_[ii][jj]] = - abs(a_[ii][jj]);
        capacity_[nodearcs_[ii][jj]] = 1;
      }
    }

    // add arcs from source to column 1
    for (int ii = 0; ii < r_; ++ii) {
      lemon::ListDigraph::Arc a = g_.addArc(s_, innode_[ii][0]);
      cost_[a] = 0;
      capacity_[a] = 1;
    }

    // add arcs from column c to sink
    for (int ii = 0; ii < r_; ++ii) {
      lemon::ListDigraph::Arc a = g_.addArc(outnode_[ii][c_ - 1], t_);
      cost_[a] = 0;
      capacity_[a] = 1;
    }

    // add arcs between columns
    colarcs_.resize(r_);
    for (int row = 0; row < r_; ++row) {
      colarcs_[row].resize(c_ - 1);
      for (int col = 0; col < c_ - 1; ++col) {
        colarcs_[row][col].resize(r_);
        for (int dest = 0; dest < r_; ++dest) {
          colarcs_[row][col][dest] =
            g_.addArc(outnode_[row][col], innode_[dest][col + 1]);
          cost_[colarcs_[row][col][dest]] = abs(row - dest);
          capacity_[colarcs_[row][col][dest]] = 1;
        }
      }
    }

    alg_ = new MCMFAlgorithm(g_);
    alg_->upperMap(capacity_);

    set_sparsity(0);
  }

};

#endif
