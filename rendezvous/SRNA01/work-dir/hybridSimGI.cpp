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

  for(int i=0; i<7; i++){
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
  double mult_factor = getMultFactor(ptUpper);
  if(curMode==1){
    Constraint_System cs_box;
    Variable xp(0);
    Variable yp(1);

    cs_box.insert(mult_factor*xp==mult_factor*ptUpper[1]);
    cs_box.insert(mult_factor*yp==mult_factor*ptUpper[2]);
    box_poly = NNC_Polyhedron(cs_box);
    Pointset_Powerset<NNC_Polyhedron> box(box_poly);
    Pointset_Powerset<NNC_Polyhedron> invariant(2,UNIVERSE);
    Pointset_Powerset<NNC_Polyhedron> curInv;
    NNC_Polyhedron curPoly;
    Constraint_System cs;
    curInv = Pointset_Powerset<NNC_Polyhedron>(2,EMPTY);

    cs.set_space_dimension(2);
    cs.insert(10*yp < -1000);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    cs.set_space_dimension(2);
    cs.insert(10*xp + 10*yp < -1411.0);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    cs.set_space_dimension(2);
    cs.insert(10*xp < -1000);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    cs.set_space_dimension(2);
    cs.insert(-10*xp + 10*yp > 1411.0);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    cs.set_space_dimension(2);
    cs.insert(10*yp > 1000);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    cs.set_space_dimension(2);
    cs.insert(10*xp + 10*yp > 1411.0);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    cs.set_space_dimension(2);
    cs.insert(10*xp > 1000);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    cs.set_space_dimension(2);
    cs.insert(-10*xp + 10*yp < -1411.0);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    return invariant.contains(box);
  }
  if(curMode==2){
    Constraint_System cs_box;
    Variable xp(0);
    Variable yp(1);

    cs_box.insert(mult_factor*xp==mult_factor*ptUpper[1]);
    cs_box.insert(mult_factor*yp==mult_factor*ptUpper[2]);
    box_poly = NNC_Polyhedron(cs_box);
    Pointset_Powerset<NNC_Polyhedron> box(box_poly);
    Pointset_Powerset<NNC_Polyhedron> invariant(2,UNIVERSE);
    Pointset_Powerset<NNC_Polyhedron> curInv;
    NNC_Polyhedron curPoly;
    Constraint_System cs;
    curInv = Pointset_Powerset<NNC_Polyhedron>(2,EMPTY);

    cs.set_space_dimension(2);
    cs.insert(10*yp >= -1000);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    curInv = Pointset_Powerset<NNC_Polyhedron>(2,EMPTY);

    cs.set_space_dimension(2);
    cs.insert(10*xp + 10*yp >= -1411.0);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    curInv = Pointset_Powerset<NNC_Polyhedron>(2,EMPTY);

    cs.set_space_dimension(2);
    cs.insert(10*xp >= -1000);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    curInv = Pointset_Powerset<NNC_Polyhedron>(2,EMPTY);

    cs.set_space_dimension(2);
    cs.insert(-10*xp + 10*yp <= 1411.0);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    curInv = Pointset_Powerset<NNC_Polyhedron>(2,EMPTY);

    cs.set_space_dimension(2);
    cs.insert(10*yp <= 1000);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    curInv = Pointset_Powerset<NNC_Polyhedron>(2,EMPTY);

    cs.set_space_dimension(2);
    cs.insert(10*xp + 10*yp <= 1411.0);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    curInv = Pointset_Powerset<NNC_Polyhedron>(2,EMPTY);

    cs.set_space_dimension(2);
    cs.insert(10*xp <= 1000);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    curInv = Pointset_Powerset<NNC_Polyhedron>(2,EMPTY);

    cs.set_space_dimension(2);
    cs.insert(-10*xp + 10*yp >= -1411.0);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    return invariant.contains(box);
  }
  return true;
}

extern "C" vector<pair<NNC_Polyhedron, int> > hitsGuard(int curMode, double *ptLower, double *ptUpper){
  vector<pair<NNC_Polyhedron, int> > toRet;
  NNC_Polyhedron box_poly;
  double mult_factor = getMultFactor(ptUpper);
  if(curMode==1){
    Constraint_System cs_box;
    Variable yp(0);
    Variable xp(1);

    cs_box.insert(mult_factor*xp==mult_factor*ptUpper[1]);
    cs_box.insert(mult_factor*yp==mult_factor*ptUpper[2]);
    box_poly = NNC_Polyhedron(cs_box);
    Constraint_System cs;
    cs.set_space_dimension(2);
    cs.insert(10*yp >= -1000);
    cs.insert(10*xp + 10*yp >= -1411.0);
    cs.insert(10*xp >= -1000);
    cs.insert(-10*xp + 10*yp <= 1411.0);
    cs.insert(10*yp <= 1000);
    cs.insert(10*xp + 10*yp <= 1411.0);
    cs.insert(10*xp <= 1000);
    cs.insert(-10*xp + 10*yp >= -1411.0);
    NNC_Polyhedron guard(cs);
    if(guard.contains(box_poly)){
      Constraint_System cs_box;
      Variable Simu_time_temp(0);
      Variable xp_temp(1);
      Variable yp_temp(2);
      Variable xd_temp(3);
      Variable yd_temp(4);
      Variable t_temp(5);
      Variable loc_temp(6);

      cs_box.insert(mult_factor*Simu_time_temp==mult_factor*ptUpper[0]);
      cs_box.insert(mult_factor*xp_temp==mult_factor*ptUpper[1]);
      cs_box.insert(mult_factor*yp_temp==mult_factor*ptUpper[2]);
      cs_box.insert(mult_factor*xd_temp==mult_factor*ptUpper[3]);
      cs_box.insert(mult_factor*yd_temp==mult_factor*ptUpper[4]);
      cs_box.insert(mult_factor*t_temp==mult_factor*ptUpper[5]);
      cs_box.insert(mult_factor*loc_temp==mult_factor*ptUpper[6]);
      box_poly = NNC_Polyhedron(cs_box);
      Constraint_System cs_gd;
      cs_gd.set_space_dimension(14);
      Variable Simu_time_new(7);
      Variable xp_new(8);
      Variable yp_new(9);
      Variable xd_new(10);
      Variable yd_new(11);
      Variable t_new(12);
      Variable loc_new(13);
      box_poly.add_space_dimensions_and_embed(7);
      cs_gd.insert(10*loc_new==20);
      cs_gd.insert(xp_new==xp_temp);
      cs_gd.insert(yp_new==yp_temp);
      cs_gd.insert(t_new==t_temp);
      cs_gd.insert(Simu_time_new==Simu_time_temp);
      cs_gd.insert(yd_new==yd_temp);
      cs_gd.insert(xd_new==xd_temp);
      NNC_Polyhedron guard_reset(cs_gd);
      guard_reset.intersection_assign(box_poly);
      Variables_Set vars;
      vars.insert(Simu_time_temp);
      vars.insert(xp_temp);
      vars.insert(yp_temp);
      vars.insert(xd_temp);
      vars.insert(yd_temp);
      vars.insert(t_temp);
      vars.insert(loc_temp);
      guard_reset.remove_space_dimensions(vars);
      toRet.push_back(make_pair(guard_reset,2));
    }
  }
  return toRet;
}

