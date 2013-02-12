#ifndef __EMD_FLOW_H__
#define __EMD_FLOW_H__

#include <vector>

#include "emd_flow_network_factory.h"

void emd_flow(
    const std::vector<std::vector<double> >& a,
    int k,
    int emd_bound_low,
    int emd_bound_high,
    double lambda_high,
    double lambda_eps,
    std::vector<std::vector<bool> >* result,
    int* emd_cost,
    double* amp_sum,
    double* final_lambda,
    EMDFlowNetworkFactory::EMDFlowNetworkType alg_type,
    void (*output_function)(const char*),
    bool verbose);

#endif
