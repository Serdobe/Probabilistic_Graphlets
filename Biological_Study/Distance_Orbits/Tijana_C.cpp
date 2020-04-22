  #include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>
#include <array>
#include <fstream>

/*
==================================
  Functions and global variables
==================================
*/

// Functions
int init(int argc, char const *argv[]);
double* get_weights();
double** signature_similarity4(double GDV[][15], int m);
void write_results(double** D, int m);

// Global variables
using std::fstream;
fstream fin, fout;


/*
=================================
              main
=================================
*/
int main(int argc, char const *argv[]) {

// import file
  if (!init(argc, argv)) {
		std::cerr << "Stopping!" << '\n';
		return 0;
	}

// read matrix in
  char c;
  int m, n;
  fin >> c >> m >> n;
// assert c == '#'
// assert n == 15
  double GDV[m][15];
  for (auto i = 0; i < m; i++) {
    for (auto j = 0; j < n; j++) {
      fin >> GDV[i][j];
    }
  }

  double** D = signature_similarity4(GDV, m);
  write_results(D, m);

  return 0;
}


/*
=================
  init function
=================
*/
int init(int argc, char const *argv[]) {
// open input and output files
	if (argc!=3) {
		std::cerr << "Incorrect number of arguments." << '\n';
		std::cerr << "Usage: tijana.exe [GDV matrix - input file] "
                 "[similarity matrix - output file]" << '\n';
		return 0;
	}

	fin.open(argv[1], fstream::in);
	fout.open(argv[2], fstream::out | fstream::binary);
	if (fin.fail()) {
		std::cerr << "Failed to open file " << argv[1] << '\n';
		return 0;
	}
	if (fout.fail()) {
		std::cerr << "Failed to open file " << argv[2] << '\n';
		return 0;
	}
  return 1;
}

/*
======================
  Calculate weights
======================
*/
double* get_weights() {
  static double WEIGHTS[15] = {1,2,2,2,3,4,3,3,4,3,4,4,4,4,3};
  for (auto i = 0; i < 15; i++)
    WEIGHTS[i] = 1 - std::log(WEIGHTS[i]) / std::log(15);
  return WEIGHTS;
}


/*
=========================
  signature_similarity4
=========================
*/

// Calculate log_GDV
double** signature_similarity4(double GDV[][15], int m) {

  double logGDV[m][15];
  for (auto i = 0; i < m; i++) {
    for (auto j = 0; j < 15; j++) {
      logGDV[i][j] = std::log(GDV[i][j]+1);
    }
  }

// Create empty distance matrix
  double **D = new double*[m];
  for(auto i = 0; i < m; ++i)
    D[i] = new double[m];

  double* WEIGHTS = get_weights();
  double WEIGHT_SUM = 0;

  for (auto i = 0; i < 15; i++) {
    WEIGHT_SUM += WEIGHTS[i];
  }

  double sum = 0;
  for (auto u = 0; u < m; u++) {
    for (auto v = u+1; v < m; v++) {
      sum = 0;
      for (auto i = 0; i < 15; i++) {
        sum += WEIGHTS[i] * std::abs(logGDV[u][i] - logGDV[v][i]) /
              std::log(std::max(GDV[u][i], GDV[v][i])+2);
      }
      D[u][v] = sum/WEIGHT_SUM;
      D[v][u] = D[u][v];
    }
  }
  return D;
}


/*
=================
  write_results
=================
*/
void write_results(double** D, int m){
  for (auto i=0;i<m;i++) {
		for (auto j=0;j<m;j++)
			fout << D[i][j] << ' ';
    fout << '\n';
	}
	fout.close();
}

