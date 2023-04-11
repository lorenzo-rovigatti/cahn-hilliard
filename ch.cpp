#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

#define D 1.0
#define H 1
#define DT 0.01

template<int dims>
struct CahnHilliard {
  int N = 64;
  int N_minus_one;
  int bits;
  int size;
  double epsilon;

  std::vector<double> psi;

  CahnHilliard(int mN, double eps, double psi_average) : N(mN), epsilon(eps) {
    double log2N = std::log2(N);
    if(ceil(log2N) != floor(log2N)) {
      fprintf(stderr, "N should be a power of 2\n");
      exit(1);
    }

    N_minus_one = N - 1;
    bits = (int) log2N;
    
    size = N;
    for(int i = 1; i < dims; i++) {
      size *= N;
    }

    psi.resize(size);
    std::generate(psi.begin(), psi.end(), [psi_average]() {return (drand48() - 0.5) + psi_average;});
  }

  void fill_coords(int coords[dims], int idx);
  int cell_idx(int coords[dims]);

  double cell_laplacian(std::vector<double> &field, int idx);

  void evolve();

  void print_state(std::string filename);
};

template<int dims>
void CahnHilliard<dims>::fill_coords(int coords[dims], int idx) {
  for(int d = 0; d < dims; d++) {
    coords[d] = idx & N_minus_one;
    idx >>= bits; // divide by N
  }
}

template<int dims>
int CahnHilliard<dims>::cell_idx(int coords[dims]) {
  int idx = 0;
  int multiply_by = 1;
  for(int d = 0; d < dims; d++) {
    idx += coords[d] * multiply_by;
    multiply_by <<= bits; // multiply by N
  }
  return idx;
}

template<>
double CahnHilliard<2>::cell_laplacian(std::vector<double> &field, int idx) {
  int coords_xy[2];
  fill_coords(coords_xy, idx);

  int coords_xmy[2] = {
			  (coords_xy[0] - 1 + N) & N_minus_one,
			  coords_xy[1]
  };
  
  int coords_xym[2] = {
			  coords_xy[0],
			  (coords_xy[1] - 1 + N) & N_minus_one
  };
  
  int coords_xpy[2] = {
			  (coords_xy[0] + 1) & N_minus_one,
			  coords_xy[1]
  };
  
  int coords_xyp[2] = {
			  coords_xy[0],
			  (coords_xy[1] + 1) & N_minus_one
  };

  return (field[cell_idx(coords_xmy)] + field[cell_idx(coords_xpy)] + field[cell_idx(coords_xym)] + field[cell_idx(coords_xyp)] - 4 * field[idx]) / (H * H);
}

template<int dims>
void CahnHilliard<dims>::evolve() {
  static std::vector<double> free_energy_der(psi.size());
  for(unsigned int idx = 0; idx < psi.size(); idx++) {
    free_energy_der[idx] = cell_laplacian(psi, idx) + epsilon * psi[idx] - psi[idx] * psi[idx] * psi[idx];
  }

  for(unsigned int idx = 0; idx < psi.size(); idx++) {
    psi[idx] += -D * cell_laplacian(free_energy_der, idx) * DT;
  }
}

template<int dims>
void CahnHilliard<dims>::print_state(std::string filename) {
  std::ofstream out(filename.c_str());
  
  for(int idx = 0; idx < size; idx++) {
    if(idx > 0) {
      int modulo = N;
      for(int d = 1; d < dims; d++) {
	if(idx % modulo == 0) {
	  out << std::endl;
	}
	modulo <<= N;
      }
    }
    out << psi[idx] << " ";
  }

  out.close();
}

int main(int argc, char *argv[]) {
  srand48(51328);

  if(argc < 5) {
    fprintf(stderr, "Usage is %s N epsilon psi_average steps\n", argv[0]);
    exit(1);
  }

  int N = std::atoi(argv[1]);
  double epsilon = std::atof(argv[2]);
  double psi_average = std::atof(argv[3]);
  long long int steps = std::atol(argv[4]);

  CahnHilliard<2> system(N, epsilon, psi_average);
  
  for(int t = 0; t < steps; t++) {
    system.evolve();
  }

  system.print_state("last.dat");

  return 0;
}
