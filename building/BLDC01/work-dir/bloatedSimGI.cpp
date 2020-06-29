#include <ppl.hh>
#include <iostream>
#include <utility>
#include <vector>
#include <fstream>
#include <typeinfo>

using namespace std;

void print_box(NNC_Polyhedron poly){
  Generator_System gs=poly.minimized_generators();
  Generator_System::const_iterator i;
  double divisor, dividend;
  int dim;
  cout << "POLY: " << endl;
  for(i=gs.begin();i!=gs.end();++i){
    if(i->is_point()){
      divisor=mpz_get_d(i->divisor().get_mpz_t());
      dim=int(i->space_dimension());
      cout << "POINT: ";
      for(int j=0;j<dim;j++){
        dividend=mpz_get_d(i->coefficient(Variable(j)).get_mpz_t());
        cout<<dividend/divisor<<" ";
      }
      cout<<endl;
    }
  }
  cout << endl;
}

double getMultFactor(double *pt){
  int multiplier=0, tmp_mul, str_len;
  char buffer[100];
  char *dot_loc;

  for(int i=0; i<51; i++){
    sprintf(buffer, "%lf", pt[i]);
    str_len = strlen(buffer);
    dot_loc = strchr(buffer,'.');
    if(dot_loc){
      tmp_mul = (str_len-1)-(dot_loc-buffer);
      if(tmp_mul>multiplier){
        multiplier=tmp_mul;
      }
    }
  }

  return pow(10, multiplier);
}

double getMultFactor(double *ptLower, double *ptUpper){
  int lowerMult = getMultFactor(ptLower);
  int upperMult = getMultFactor(ptUpper);
  int multiplier = lowerMult > upperMult ? lowerMult : upperMult;
  return multiplier;
}

extern "C" bool invariantSatisfied(int curMode, double *ptLower, double *ptUpper){
  NNC_Polyhedron box_poly;
  double mult_factor = getMultFactor(ptLower, ptUpper);
  return true;
}

extern "C" vector<pair<NNC_Polyhedron, int> > hitsGuard(int curMode, double *ptLower, double *ptUpper){
  vector<pair<NNC_Polyhedron, int> > toRet;
  NNC_Polyhedron box_poly;
  double mult_factor = getMultFactor(ptLower, ptUpper);
  return toRet;
}

