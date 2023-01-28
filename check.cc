#include <iostream>
#include <vector>
#include <climits>
#include <cassert>
#include <fstream>

using namespace std;

class Option {
public:
  int c;
  int n; // c out of n
  Option():c(0),n(0){}
  Option(int _c,int _n):c(_c),n(_n){}
};

int nClasses, nOptions, nCars; 
vector<Option> options;
vector<vector<bool>> classOption;
vector<int> carsForClass;
int inputUB;


int costSolution;
vector<int> solution;
void readInput(ifstream& in) {
  in >> nCars >> nOptions >> nClasses;
  options = vector<Option>(nOptions);
  classOption = vector<vector<bool>>(nClasses,vector<bool>(nOptions,false));
  carsForClass = vector<int>(nClasses);
  for (int o = 0; o < nOptions; ++o) in >> options[o].c;
  for (int o = 0; o < nOptions; ++o) in >> options[o].n;
  for (int c = 0; c < nClasses; ++c) {
    int aux;
    in >> aux >> carsForClass[c];
    for (int o = 0; o < nOptions; ++o) {
      in >> aux;
      classOption[c][o] = (aux == 1);
    }
  }

  string s; in >> s;
  in >> inputUB;
}

void readSolution(ifstream& in){
  double t;
  in >> costSolution >> t;
  int x;
  while (in >> x) solution.push_back(x);
}

int countViolations (const vector<int>& solution, const vector<vector<bool>>& classOption) {
  int v = 0;
  vector<int> countOption(nOptions,0);
  for (uint i = 0; i < solution.size(); ++i) {
    for (int o = 0; o < nOptions; ++o) {
      if (classOption[solution[i]][o]) 	++countOption[o];
      int ejectingPos = i-options[o].n;
      if (ejectingPos >= 0 and classOption[solution[ejectingPos]][o]) --countOption[o];
      if (countOption[o] > options[o].c) v+=(countOption[o] - options[o].c);
    }
  }


  // Now count ending windows
  for (int o = 0; o < nOptions; ++o) {
    for (int i = int(solution.size()) - options[o].n; i < int(solution.size()) - 1; ++i) {
      if (classOption[solution[i]][o]) --countOption[o];
      if (countOption[o] > options[o].c) v += (countOption[o] - options[o].c);
    }
  }
  return v;
}

void checkSolution( ){

  bool error = false;
  // Check number of cars
  if (int(solution.size()) != nCars) {
    cout << "ERROR: Solution has " << solution.size() << " cars instead of " << nCars << endl;
    error = true;
  }

  // Check number of cars from a concrete class
  for (int c = 0; c < nClasses; ++c) {
    int n = 0;
    for (auto& x:solution) if (x == c) ++n;
    if (n != carsForClass[c]) {
      cout << "ERROR: Solution has " << n << " cars of class " << c << " instead of " << carsForClass[c] << endl;
      error = true;
    }
  }

  // Check bound
  int v = countViolations(solution,classOption);
  if (v != costSolution) {
    cout << "ERROR: Solution has a cost of " << v << ", but it is claimed to have " << costSolution << endl;
    error = true;
  }

  if (not error) cout << "Solution is OK" << endl;
}

int main (int argc, char** argv ){
  // Input Solution
  if (argc != 3) {
    cout << "Syntax is: " << argv[0] << " input_file sol_file" << endl;
    exit(1);
  }

  
  ifstream inputFile(argv[1],ifstream::in);
  ifstream solFile(argv[2],ifstream::in);

  cout << "Checking " << argv[1] << endl;
  readInput(inputFile);
  readSolution(solFile);
  checkSolution();
}
