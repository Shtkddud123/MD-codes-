#include "pteros/pteros.h"
#include "pteros/core/mol_file.h"
#include <Eigen/Dense>

using namespace std;
using namespace pteros;

using Eigen::MatrixXd;
System sys1;

int main() {
  
  sys1.load("dppc-dupc-chol-eq.gro"); // Load data into system
  Selection sel0(sys1,"all");
  cout << sel0 << endl;

  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;
}
