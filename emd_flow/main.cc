#include <vector>
#include <cstdio>
#include <cmath>
#include <boost/program_options.hpp>

#include "emd_flow.h"

using namespace std;
namespace po = boost::program_options;

// rows, columns
int r, c;
// sparsity
int k;
// EMD bound
int emd_bound;
// amplitudes
std::vector<std::vector<double> > a;
// result
std::vector<std::vector<bool> > result;

void output_function(const char* s) {
  fprintf(stderr, s);
  fflush(stderr);
}

int main(int argc, char** argv)
{
  po::options_description desc("Allowed options");
  desc.add_options()
      ("matrix_output", po::value<string>(), "File for binary output matrix")
      ("square_amplitudes", "Square all input amplitudes");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm); 

  scanf("%d %d %d %d", &r, &c, &k, &emd_bound);

  a.resize(r);
  for (int ii = 0; ii < r; ++ii) {
    a[ii].resize(c);
    for (int jj = 0; jj < c; ++jj) {
      scanf("%lg", &(a[ii][jj]));
      a[ii][jj] = abs(a[ii][jj]);
    }
  }

  if (vm.count("square_amplitudes")) {
    fprintf(stderr, "Squaring all amplitudes ...\n");
    for (int ii = 0; ii < r; ++ii) {
      for (int jj = 0; jj < c; ++jj) {
        a[ii][jj] = a[ii][jj] * a[ii][jj];
      }
    }
  }

  
  int emd_cost = 0;
  double amp_sum = 0.0;
  double final_lambda = 0.0;
  emd_flow(a, k, emd_bound, &result, &emd_cost, &amp_sum, &final_lambda,
      output_function, true);

  for (int jj = 0; jj < c; ++jj) {
    fprintf(stderr, "col %d:\n", jj + 1);
    for (int ii = 0; ii < r; ++ii) {
      if (result[ii][jj]) {
        fprintf(stderr, " row %d, amplitude %lf\n", ii + 1, a[ii][jj]);
      }
    }
  }

  printf("%lf\n", amp_sum);

  if (vm.count("matrix_output")) {
    string output_file_name = vm["matrix_output"].as<string>();
    FILE* output_file = fopen(output_file_name.c_str(), "w");
    for (int ii = 0; ii < r; ++ii) {
      for (int jj = 0; jj < c; ++jj) {
        if (result[ii][jj]) {
          fprintf(output_file, "1 ");
        } else {
          fprintf(output_file, "0 ");
        }
      }
      fprintf(output_file, "\n");
    }
    fclose(output_file);
  }

  return 0;
}
