#include <vector>
#include <cmath>
#include <cstdio>
#include <ctime>

#include "emd_flow_network.h"

using namespace lemon;
using namespace std;

void emd_flow(
    const vector<vector<double> >& a,
    int k,
    int emd_bound,
    vector<vector<bool> >* result,
    int* emd_cost,
    double* amp_sum,
    double* final_lambda,
    void (*output_function)(const char*),
    bool verbose) {

  clock_t total_time_begin = clock();

  const int kOutputBufferSize = 1000;
  char output_buffer[kOutputBufferSize];

  int r = a.size();
  int c = a[0].size();

  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize, "r = %d,  c = %d,  k = %d,  "
        "emd_bound = %d\n", r, c, k, emd_bound);
    output_function(output_buffer);
  }

  // build graph
  clock_t graph_construction_time_begin = clock();

  EMDFlowNetwork network(a, k);

  clock_t graph_construction_time = clock() - graph_construction_time_begin;

  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize, "The graph has %d nodes and %d "
        "edges.\n", network.get_num_nodes(), network.get_num_edges());
    output_function(output_buffer);
    snprintf(output_buffer, kOutputBufferSize, "Total construction time: %lf "
        "s\n ", static_cast<double>(graph_construction_time) / CLOCKS_PER_SEC);
    output_function(output_buffer);

  }

  // make lambda larger until we find a solution that fits into the EMD budget
  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize,
        "Finding large enough value of lambda ...\n");
    output_function(output_buffer);
  }

  double lambda_high = 0.01;
  while (true) {
    network.run_flow(lambda_high);
    int cur_emd_cost = network.get_EMD_used();
    double cur_amp_sum = network.get_supported_amplitude_sum();

    if (verbose) {
      snprintf(output_buffer, kOutputBufferSize, "l: %lf  EMD: %d  amp sum: %lf"
          "\n", lambda_high, cur_emd_cost, cur_amp_sum);
      output_function(output_buffer);
    }

    if (cur_emd_cost <= emd_bound) {
      break;
    } else {
      lambda_high = lambda_high * 2;
    }
  }

  // binary search on lambda
  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize, "Binary search on lambda ...\n");
    output_function(output_buffer);
  }

  double lambda_low = 0;
  double lambda_eps = 0.00001;
  while(lambda_high - lambda_low > lambda_eps) {
    double cur_lambda = (lambda_high + lambda_low) / 2;
    network.run_flow(cur_lambda);
    int cur_emd_cost = network.get_EMD_used();
    double cur_amp_sum = network.get_supported_amplitude_sum();

    if (verbose) {
      snprintf(output_buffer, kOutputBufferSize, "l_cur: %lf  (l_low: %lf, "
          "l_high: %lf)  EMD: %d  amp sum: %lf\n", cur_lambda, lambda_low,
          lambda_high, cur_emd_cost, cur_amp_sum);
      output_function(output_buffer);
    }

    if (cur_emd_cost <= emd_bound) {
      lambda_high = cur_lambda;
    } else {
      lambda_low = cur_lambda;
    }
  }

  // run with final lambda
  network.run_flow(lambda_high);
  *emd_cost = network.get_EMD_used();
  *amp_sum = network.get_supported_amplitude_sum();
  *final_lambda = lambda_high;
  network.get_support(result);

  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize, "Final l: %lf, amp sum: %lf, "
        "EMD cost: %d\n", lambda_high, *amp_sum, *emd_cost);
    output_function(output_buffer);
  }


  clock_t total_time = clock() - total_time_begin;
  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize, "Total time %lf s\n",
        static_cast<double>(total_time) / CLOCKS_PER_SEC);
    output_function(output_buffer);
  }

  return;
}
