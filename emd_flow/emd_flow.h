#ifndef __EMD_FLOW_H__
#define __EMD_FLOW_H__

#include <vector>

void emd_flow(
    const std::vector<std::vector<double> >& a,
    int k,
    int emd_bound,
    std::vector<std::vector<bool> >* result,
    int* emd_cost,
    double* amp_sum,
    double* final_lambda,
    void (*output_function)(const char*),
    bool verbose);

#endif
