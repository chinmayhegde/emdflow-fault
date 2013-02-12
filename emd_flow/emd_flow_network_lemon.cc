#include "emd_flow_network_lemon.h"

#include <cstdio>
#include <cmath>

using namespace std;
using namespace lemon;

template <typename MCMFAlgorithm>
EMDFlowNetworkLemon::EMDFlowNetworkLemon(
    const vector<vector<double> >& amplitudes, int k)
      : k_(k), capacity_(g_), cost_(g_) {

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
    ListDigraph::Arc a = g_.addArc(s_, innode_[ii][0]);
    cost_[a] = 0;
    capacity_[a] = 1;
  }

  // add arcs from column c to sink
  for (int ii = 0; ii < r_; ++ii) {
    ListDigraph::Arc a = g_.addArc(outnode_[ii][c_ - 1], t_);
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

  set_sparsity(k);
}

EMDFlowNetworkLemon::~EMDFlowNetworkLemon() {
  delete alg_;
}

void EMDFlowNetworkLemon::set_sparsity(int k) {
  k_ = k;
  alg_->stSupply(s_, t_, k_);
}

void EMDFlowNetworkLemon::run_flow(double lambda) {
  apply_lambda(lambda);
  alg_->costMap(cost_);
  alg_->run();
}

int EMDFlowNetworkLemon::get_EMD_used() {
  return extract_emd_cost();
}

double EMDFlowNetworkLemon::get_supported_amplitude_sum() {
  return extract_amp_sum();
}

void EMDFlowNetworkLemon::get_support(
    std::vector<std::vector<bool> >* support) {
  support->resize(r_);
  for (int row = 0; row < r_; ++row) {
    (*support)[row].resize(c_);
    for (int col = 0; col < c_; ++col) {
      (*support)[row][col] = (alg_->flow(nodearcs_[row][col]) > 0);
    }
  }
}

void EMDFlowNetworkLemon::apply_lambda(double lambda) {
  for (int row = 0; row < r_; ++row) {
    for (int col = 0; col < c_ - 1; ++col) {
      for (int dest = 0; dest < r_; ++dest) {
        cost_[colarcs_[row][col][dest]] = lambda * abs(row - dest);
      }
    }
  }
}

int EMDFlowNetworkLemon::extract_emd_cost() {
  int emd_cost = 0;
  for (int row = 0; row < r_; ++row) {
    for (int col = 0; col < c_ - 1; ++col) {
      for (int dest = 0; dest < r_; ++dest) {
        if (alg_->flow(colarcs_[row][col][dest]) > 0) {
          emd_cost += abs(row - dest);
          if (alg_->flow(colarcs_[row][col][dest]) != 1) {
            fprintf(stderr, "ERROR: nonzero flow on a column edge is not 1.\n");
          }
        }
      }
    }
  }
  return emd_cost;

}

double EMDFlowNetworkLemon::extract_amp_sum() {
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

int EMDFlowNetworkLemon::get_num_nodes() {
  return countNodes(g_);
}

int EMDFlowNetworkLemon::get_num_edges() {
  return countArcs(g_);
}

int EMDFlowNetworkLemon::get_num_columns() {
  return c_;
}

int EMDFlowNetworkLemon::get_num_rows() {
  return r_;
}
