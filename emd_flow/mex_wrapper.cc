#include <vector>
#include <string>
#include <set>
#include <sstream>

#include <math.h>
#include <matrix.h>
#include <mex.h>

#include "mex_helper.h"
#include "emd_flow.h"
#include "emd_flow_network_factory.h"

using namespace std;

void output_function(const char* s) {
  mexPrintf(s);
  mexEvalString("drawnow;");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs < 3) {
    mexErrMsgTxt("At least three input argument required (amplitudes, sparsity,"
        " EMD budget.");
  }
  if (nlhs > 4) {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  vector<vector<double> > a;
  if (!get_double_matrix(prhs[0], &a)) {
    mexErrMsgTxt("Amplitudes need to be a two-dimensional double array.");
  }

  int k = 0;
  if (!get_double_as_int(prhs[1], &k)) {
    mexErrMsgTxt("Sparsity has to be a double scalar.");
  }

  int emd_bound_low = 0;
  int emd_bound_high = 0;
  if (!get_double_as_int(prhs[2], &emd_bound_low)) {
    if (!get_double_interval_as_ints(prhs[2], &emd_bound_low,
        &emd_bound_high)) {
      mexErrMsgTxt("EMD budget has to be a double scalar or a double "
          "interval.");
    }
  } else {
    emd_bound_high = emd_bound_low;
  }

  // optional parameters
  bool verbose = false;
  double lambda_high = 1.0;
  double lambda_eps = 0.0001;
  if (nrhs == 4) {
    set<string> known_options;
    known_options.insert("verbose");
    known_options.insert("lambda_high");
    known_options.insert("lambda_eps");
    vector<string> options;
    if (!get_fields(prhs[3], &options)) {
      mexErrMsgTxt("Cannot get fields from options argument.");
    }
    for (size_t ii = 0; ii < options.size(); ++ii) {
      if (known_options.find(options[ii]) == known_options.end()) {
        ostringstream ss;
        ss << "Unknown option \"" << options[ii].c_str() << "\"";
        string msg = ss.str();
        mexErrMsgTxt(msg.c_str());
      }
    }

    if (has_field(prhs[3], "verbose")
        && !get_bool_field(prhs[3], "verbose", &verbose)) {
      mexErrMsgTxt("verbose flag has to be a boolean scalar.");
    }

    if (has_field(prhs[3], "lambda_high")
        && !get_double_field(prhs[3], "lambda_high", &lambda_high)) {
      mexErrMsgTxt("lambda_high flag has to be a boolean scalar.");
    }

    if (has_field(prhs[3], "lambda_eps")
        && !get_double_field(prhs[3], "lambda_eps", &lambda_eps)) {
      mexErrMsgTxt("lambda_eps flag has to be a boolean scalar.");
    }
  }

  vector<vector<bool> > result;
  int emd_cost;
  double amp_sum;
  double final_lambda;

  emd_flow(a, k, emd_bound_low, emd_bound_high, lambda_high, lambda_eps,
      &result, &emd_cost, &amp_sum, &final_lambda,
      EMDFlowNetworkFactory::kShortestAugmentingPath, output_function, verbose);

  if (nlhs >= 1) {
    set_double_matrix(&(plhs[0]), result);
  }

  if (nlhs >= 2) {
    set_double(&(plhs[1]), emd_cost);
  }

  if (nlhs >= 3) {
    set_double(&(plhs[2]), amp_sum);
  }

  if (nlhs >= 4) {
    set_double(&(plhs[3]), final_lambda);
  }

  /*
  mexPrintf("r = %d, c = %d, k = %d, EMD budget = %d\n", r, c, k, emd_budget);
  for (int ii = 0; ii < r; ++ii) {
    for (int jj = 0; jj < c; ++jj) {
      mexPrintf("%lf ", a[ii][jj]);
    }
    mexPrintf("\n");
  }
  */
}
