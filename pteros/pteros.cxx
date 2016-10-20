#include "pteros/pteros.h"
#include "pteros/core/mol_file.h"

using namespace std;
using namespace pteros;

System sys1;

int main() {
  
  sys1.load("dppc-dupc-chol-eq.gro"); // Load data into system
  Selection sel0(sys1,"all");
  cout << sel0 << endl;

}
