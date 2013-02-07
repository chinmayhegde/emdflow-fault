#include <vector>
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include "emd_flow.h"

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

  int numdims = 0;
  const mwSize* dims;
  
  bool verbose = false;
  if (nrhs == 4) {
    numdims = mxGetNumberOfDimensions(prhs[3]);
    dims = mxGetDimensions(prhs[3]);
    if (numdims != 2 || dims[0] != 1 || dims[1] != 1) {
      mexErrMsgTxt("Verbose flag has to be a scalar.");
    }
    verbose = static_cast<bool*>(mxGetData(prhs[3]))[0];
  }
  
  if (nlhs > 4) {
    mexErrMsgTxt("Too many output arguments.");
  }

  int r = 0, c = 0;

  numdims = mxGetNumberOfDimensions(prhs[0]);
  dims = mxGetDimensions(prhs[0]);
  if (numdims != 2) {
    mexErrMsgTxt("Amplitudes need to be a two-dimensional array.");
  }
  if (!mxIsClass(prhs[0], "double")) {
    mexErrMsgTxt("Amplitudes need to be a double array.");
  }
  r = dims[0];
  c = dims[1];
  double* a_linear = mxGetPr(prhs[0]);
  vector<vector<double> > a;
  a.resize(r);
  for (int ir = 0; ir < r; ++ir) {
    a[ir].resize(c);
    for (int ic = 0; ic < c; ++ic) {
      a[ir][ic] = a_linear[ir + ic * r];
    }
  }

  numdims = mxGetNumberOfDimensions(prhs[1]);
  dims = mxGetDimensions(prhs[1]);
  if (numdims != 2 || dims[0] != 1 || dims[1] != 1) {
    mexErrMsgTxt("Sparsity has to be a scalar.");
  }
  int k = mxGetPr(prhs[1])[0];

  numdims = mxGetNumberOfDimensions(prhs[2]);
  dims = mxGetDimensions(prhs[2]);
  if (numdims != 2 || dims[0] != 1 || dims[1] != 1) {
    mexErrMsgTxt("EMD budget has to be a scalar.");
  }
  int emd_budget = mxGetPr(prhs[2])[0];

  vector<vector<bool> > result;
  int emd_cost;
  double amp_sum;
  double final_lambda;

  emd_flow(a, k, emd_budget, &result, &emd_cost, &amp_sum, &final_lambda,
      output_function, verbose);

  if (nlhs >= 1) {
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dims = mxGetDimensions(prhs[0]);
    plhs[0] = mxCreateNumericArray(numdims, dims, mxUINT8_CLASS, mxREAL);
    unsigned char* result_linear = static_cast<unsigned char*>(
        mxGetData(plhs[0]));

    for (int ir = 0; ir < r; ++ir) {
      for (int ic = 0; ic < c; ++ic) {
        result_linear[ir + ic * r] = result[ir][ic];
      }
    }
  }

  if (nlhs >= 2) {
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *(mxGetPr(plhs[1])) = emd_cost;
  }

  if (nlhs >= 3) {
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *(mxGetPr(plhs[2])) = amp_sum;
  }

  if (nlhs >= 4) {
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *(mxGetPr(plhs[3])) = final_lambda;
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
