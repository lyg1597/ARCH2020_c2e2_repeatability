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

  for(int i=0; i<8; i++){
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
  if(curMode==1){
    Constraint_System cs_box;
    Variable px(0);
    Variable py(1);
    Variable trans(2);

    if(ptLower[3]<ptUpper[3]){
      cs_box.insert(mult_factor*px>=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px<=mult_factor*ptUpper[3]);
    }
    else{
      cs_box.insert(mult_factor*px<=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px>=mult_factor*ptUpper[3]);
    }

    if(ptLower[4]<ptUpper[4]){
      cs_box.insert(mult_factor*py>=mult_factor*ptLower[4]);
      cs_box.insert(mult_factor*py<=mult_factor*ptUpper[4]);
    }
    else{
      cs_box.insert(mult_factor*py<=mult_factor*ptLower[4]);
      cs_box.insert(mult_factor*py>=mult_factor*ptUpper[4]);
    }

    if(ptLower[6]<ptUpper[6]){
      cs_box.insert(mult_factor*trans>=mult_factor*ptLower[6]);
      cs_box.insert(mult_factor*trans<=mult_factor*ptUpper[6]);
    }
    else{
      cs_box.insert(mult_factor*trans<=mult_factor*ptLower[6]);
      cs_box.insert(mult_factor*trans>=mult_factor*ptUpper[6]);
    }

    box_poly = NNC_Polyhedron(cs_box);
    Pointset_Powerset<NNC_Polyhedron> box(box_poly);
    Pointset_Powerset<NNC_Polyhedron> invariant(3,UNIVERSE);
    Pointset_Powerset<NNC_Polyhedron> curInv;
    NNC_Polyhedron curPoly;
    Constraint_System cs;
    curInv = Pointset_Powerset<NNC_Polyhedron>(3,EMPTY);

    cs.set_space_dimension(3);
    cs.insert(726542528005361.0*px + 1000000000000000*py <= 0);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    curInv = Pointset_Powerset<NNC_Polyhedron>(3,EMPTY);

    cs.set_space_dimension(3);
    cs.insert(-726542528005361.0*px + 1000000000000000*py >= 0);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    curInv = Pointset_Powerset<NNC_Polyhedron>(3,EMPTY);

    cs.set_space_dimension(3);
    cs.insert(10*trans <= 35.0);
    curPoly = NNC_Polyhedron(cs);
    curInv.add_disjunct(curPoly);
    cs.clear();

    invariant.intersection_assign(curInv);

    return !(invariant.is_disjoint_from(box));
  }
  return true;
}

extern "C" vector<pair<NNC_Polyhedron, int> > hitsGuard(int curMode, double *ptLower, double *ptUpper){
  vector<pair<NNC_Polyhedron, int> > toRet;
  NNC_Polyhedron box_poly;
  double mult_factor = getMultFactor(ptLower, ptUpper);
  if(curMode==1){
    Constraint_System cs_box;
    Variable vy(0);
    Variable px(1);
    Variable vx(2);

    if(ptLower[1]<ptUpper[1]){
      cs_box.insert(mult_factor*vx>=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx<=mult_factor*ptUpper[1]);
    }
    else{
      cs_box.insert(mult_factor*vx<=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx>=mult_factor*ptUpper[1]);
    }

    if(ptLower[2]<ptUpper[2]){
      cs_box.insert(mult_factor*vy>=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy<=mult_factor*ptUpper[2]);
    }
    else{
      cs_box.insert(mult_factor*vy<=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy>=mult_factor*ptUpper[2]);
    }

    if(ptLower[3]<ptUpper[3]){
      cs_box.insert(mult_factor*px>=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px<=mult_factor*ptUpper[3]);
    }
    else{
      cs_box.insert(mult_factor*px<=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px>=mult_factor*ptUpper[3]);
    }

    box_poly = NNC_Polyhedron(cs_box);
    Constraint_System cs;
    cs.set_space_dimension(3);
    cs.insert(1000*px >= -2.0);
    cs.insert(10*vx >= 0);
    cs.insert(10*vy <= 0);
    NNC_Polyhedron guard(cs);
    if(!guard.is_disjoint_from(box_poly)){
      Constraint_System cs_box;
      Variable Simu_time_temp(0);
      Variable vx_temp(1);
      Variable vy_temp(2);
      Variable px_temp(3);
      Variable py_temp(4);
      Variable ii_temp(5);
      Variable trans_temp(6);
      Variable one_temp(7);

      if(ptLower[0]<ptUpper[0]){
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptUpper[0]);
      }
      else{
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptUpper[0]);
      }

      if(ptLower[1]<ptUpper[1]){
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptUpper[1]);
      }
      else{
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptUpper[1]);
      }

      if(ptLower[2]<ptUpper[2]){
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptUpper[2]);
      }
      else{
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptUpper[2]);
      }

      if(ptLower[3]<ptUpper[3]){
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptUpper[3]);
      }
      else{
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptUpper[3]);
      }

      if(ptLower[4]<ptUpper[4]){
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptUpper[4]);
      }
      else{
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptUpper[4]);
      }

      if(ptLower[5]<ptUpper[5]){
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptUpper[5]);
      }
      else{
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptUpper[5]);
      }

      if(ptLower[6]<ptUpper[6]){
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptUpper[6]);
      }
      else{
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptUpper[6]);
      }

      if(ptLower[7]<ptUpper[7]){
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptUpper[7]);
      }
      else{
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptUpper[7]);
      }

      box_poly = NNC_Polyhedron(cs_box);
      Constraint_System cs_gd;
      cs_gd.set_space_dimension(16);
      Variable Simu_time_new(8);
      Variable vx_new(9);
      Variable vy_new(10);
      Variable px_new(11);
      Variable py_new(12);
      Variable ii_new(13);
      Variable trans_new(14);
      Variable one_new(15);
      box_poly.add_space_dimensions_and_embed(8);
      cs_gd.insert(10*ii_new==10*ii_temp + 32.0*vx_temp - 32.0*vy_temp);
      cs_gd.insert(10*vx_new==0);
      cs_gd.insert(10*vy_new==0);
      cs_gd.insert(trans_new==one_temp + trans_temp);
      cs_gd.insert(Simu_time_new==Simu_time_temp);
      cs_gd.insert(px_new==px_temp);
      cs_gd.insert(one_new==one_temp);
      cs_gd.insert(py_new==py_temp);
      cs_gd.insert(1000*px_temp >= -2.0);
      cs_gd.insert(10*vx_temp >= 0);
      cs_gd.insert(10*vy_temp <= 0);
      NNC_Polyhedron guard_reset(cs_gd);
      guard_reset.intersection_assign(box_poly);
      Variables_Set vars;
      vars.insert(Simu_time_temp);
      vars.insert(vx_temp);
      vars.insert(vy_temp);
      vars.insert(px_temp);
      vars.insert(py_temp);
      vars.insert(ii_temp);
      vars.insert(trans_temp);
      vars.insert(one_temp);
      guard_reset.remove_space_dimensions(vars);
      toRet.push_back(make_pair(guard_reset,2));
    }
  }
  if(curMode==1){
    Constraint_System cs_box;
    Variable vy(0);
    Variable px(1);
    Variable vx(2);

    if(ptLower[1]<ptUpper[1]){
      cs_box.insert(mult_factor*vx>=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx<=mult_factor*ptUpper[1]);
    }
    else{
      cs_box.insert(mult_factor*vx<=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx>=mult_factor*ptUpper[1]);
    }

    if(ptLower[2]<ptUpper[2]){
      cs_box.insert(mult_factor*vy>=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy<=mult_factor*ptUpper[2]);
    }
    else{
      cs_box.insert(mult_factor*vy<=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy>=mult_factor*ptUpper[2]);
    }

    if(ptLower[3]<ptUpper[3]){
      cs_box.insert(mult_factor*px>=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px<=mult_factor*ptUpper[3]);
    }
    else{
      cs_box.insert(mult_factor*px<=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px>=mult_factor*ptUpper[3]);
    }

    box_poly = NNC_Polyhedron(cs_box);
    Constraint_System cs;
    cs.set_space_dimension(3);
    cs.insert(1000*px >= -2.0);
    cs.insert(10*vx <= 0);
    cs.insert(10*vy >= 0);
    NNC_Polyhedron guard(cs);
    if(!guard.is_disjoint_from(box_poly)){
      Constraint_System cs_box;
      Variable Simu_time_temp(0);
      Variable vx_temp(1);
      Variable vy_temp(2);
      Variable px_temp(3);
      Variable py_temp(4);
      Variable ii_temp(5);
      Variable trans_temp(6);
      Variable one_temp(7);

      if(ptLower[0]<ptUpper[0]){
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptUpper[0]);
      }
      else{
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptUpper[0]);
      }

      if(ptLower[1]<ptUpper[1]){
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptUpper[1]);
      }
      else{
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptUpper[1]);
      }

      if(ptLower[2]<ptUpper[2]){
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptUpper[2]);
      }
      else{
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptUpper[2]);
      }

      if(ptLower[3]<ptUpper[3]){
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptUpper[3]);
      }
      else{
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptUpper[3]);
      }

      if(ptLower[4]<ptUpper[4]){
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptUpper[4]);
      }
      else{
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptUpper[4]);
      }

      if(ptLower[5]<ptUpper[5]){
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptUpper[5]);
      }
      else{
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptUpper[5]);
      }

      if(ptLower[6]<ptUpper[6]){
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptUpper[6]);
      }
      else{
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptUpper[6]);
      }

      if(ptLower[7]<ptUpper[7]){
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptUpper[7]);
      }
      else{
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptUpper[7]);
      }

      box_poly = NNC_Polyhedron(cs_box);
      Constraint_System cs_gd;
      cs_gd.set_space_dimension(16);
      Variable Simu_time_new(8);
      Variable vx_new(9);
      Variable vy_new(10);
      Variable px_new(11);
      Variable py_new(12);
      Variable ii_new(13);
      Variable trans_new(14);
      Variable one_new(15);
      box_poly.add_space_dimensions_and_embed(8);
      cs_gd.insert(10*ii_new==10*ii_temp + 32.0*vx_temp - 32.0*vy_temp);
      cs_gd.insert(10*vx_new==0);
      cs_gd.insert(10*vy_new==0);
      cs_gd.insert(trans_new==one_temp + trans_temp);
      cs_gd.insert(Simu_time_new==Simu_time_temp);
      cs_gd.insert(px_new==px_temp);
      cs_gd.insert(one_new==one_temp);
      cs_gd.insert(py_new==py_temp);
      cs_gd.insert(1000*px_temp >= -2.0);
      cs_gd.insert(10*vx_temp <= 0);
      cs_gd.insert(10*vy_temp >= 0);
      NNC_Polyhedron guard_reset(cs_gd);
      guard_reset.intersection_assign(box_poly);
      Variables_Set vars;
      vars.insert(Simu_time_temp);
      vars.insert(vx_temp);
      vars.insert(vy_temp);
      vars.insert(px_temp);
      vars.insert(py_temp);
      vars.insert(ii_temp);
      vars.insert(trans_temp);
      vars.insert(one_temp);
      guard_reset.remove_space_dimensions(vars);
      toRet.push_back(make_pair(guard_reset,2));
    }
  }
  if(curMode==1){
    Constraint_System cs_box;
    Variable vy(0);
    Variable px(1);
    Variable vx(2);

    if(ptLower[1]<ptUpper[1]){
      cs_box.insert(mult_factor*vx>=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx<=mult_factor*ptUpper[1]);
    }
    else{
      cs_box.insert(mult_factor*vx<=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx>=mult_factor*ptUpper[1]);
    }

    if(ptLower[2]<ptUpper[2]){
      cs_box.insert(mult_factor*vy>=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy<=mult_factor*ptUpper[2]);
    }
    else{
      cs_box.insert(mult_factor*vy<=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy>=mult_factor*ptUpper[2]);
    }

    if(ptLower[3]<ptUpper[3]){
      cs_box.insert(mult_factor*px>=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px<=mult_factor*ptUpper[3]);
    }
    else{
      cs_box.insert(mult_factor*px<=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px>=mult_factor*ptUpper[3]);
    }

    box_poly = NNC_Polyhedron(cs_box);
    Constraint_System cs;
    cs.set_space_dimension(3);
    cs.insert(1000*px >= -2.0);
    cs.insert(10*vx <= 0);
    cs.insert(10*vy <= 0);
    NNC_Polyhedron guard(cs);
    if(!guard.is_disjoint_from(box_poly)){
      Constraint_System cs_box;
      Variable Simu_time_temp(0);
      Variable vx_temp(1);
      Variable vy_temp(2);
      Variable px_temp(3);
      Variable py_temp(4);
      Variable ii_temp(5);
      Variable trans_temp(6);
      Variable one_temp(7);

      if(ptLower[0]<ptUpper[0]){
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptUpper[0]);
      }
      else{
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptUpper[0]);
      }

      if(ptLower[1]<ptUpper[1]){
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptUpper[1]);
      }
      else{
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptUpper[1]);
      }

      if(ptLower[2]<ptUpper[2]){
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptUpper[2]);
      }
      else{
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptUpper[2]);
      }

      if(ptLower[3]<ptUpper[3]){
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptUpper[3]);
      }
      else{
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptUpper[3]);
      }

      if(ptLower[4]<ptUpper[4]){
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptUpper[4]);
      }
      else{
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptUpper[4]);
      }

      if(ptLower[5]<ptUpper[5]){
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptUpper[5]);
      }
      else{
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptUpper[5]);
      }

      if(ptLower[6]<ptUpper[6]){
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptUpper[6]);
      }
      else{
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptUpper[6]);
      }

      if(ptLower[7]<ptUpper[7]){
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptUpper[7]);
      }
      else{
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptUpper[7]);
      }

      box_poly = NNC_Polyhedron(cs_box);
      Constraint_System cs_gd;
      cs_gd.set_space_dimension(16);
      Variable Simu_time_new(8);
      Variable vx_new(9);
      Variable vy_new(10);
      Variable px_new(11);
      Variable py_new(12);
      Variable ii_new(13);
      Variable trans_new(14);
      Variable one_new(15);
      box_poly.add_space_dimensions_and_embed(8);
      cs_gd.insert(10*ii_new==10*ii_temp + 32.0*vx_temp - 32.0*vy_temp);
      cs_gd.insert(10*vx_new==0);
      cs_gd.insert(10*vy_new==0);
      cs_gd.insert(trans_new==one_temp + trans_temp);
      cs_gd.insert(Simu_time_new==Simu_time_temp);
      cs_gd.insert(px_new==px_temp);
      cs_gd.insert(one_new==one_temp);
      cs_gd.insert(py_new==py_temp);
      cs_gd.insert(1000*px_temp >= -2.0);
      cs_gd.insert(10*vx_temp <= 0);
      cs_gd.insert(10*vy_temp <= 0);
      NNC_Polyhedron guard_reset(cs_gd);
      guard_reset.intersection_assign(box_poly);
      Variables_Set vars;
      vars.insert(Simu_time_temp);
      vars.insert(vx_temp);
      vars.insert(vy_temp);
      vars.insert(px_temp);
      vars.insert(py_temp);
      vars.insert(ii_temp);
      vars.insert(trans_temp);
      vars.insert(one_temp);
      guard_reset.remove_space_dimensions(vars);
      toRet.push_back(make_pair(guard_reset,2));
    }
  }
  if(curMode==1){
    Constraint_System cs_box;
    Variable vy(0);
    Variable px(1);
    Variable vx(2);

    if(ptLower[1]<ptUpper[1]){
      cs_box.insert(mult_factor*vx>=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx<=mult_factor*ptUpper[1]);
    }
    else{
      cs_box.insert(mult_factor*vx<=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx>=mult_factor*ptUpper[1]);
    }

    if(ptLower[2]<ptUpper[2]){
      cs_box.insert(mult_factor*vy>=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy<=mult_factor*ptUpper[2]);
    }
    else{
      cs_box.insert(mult_factor*vy<=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy>=mult_factor*ptUpper[2]);
    }

    if(ptLower[3]<ptUpper[3]){
      cs_box.insert(mult_factor*px>=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px<=mult_factor*ptUpper[3]);
    }
    else{
      cs_box.insert(mult_factor*px<=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px>=mult_factor*ptUpper[3]);
    }

    box_poly = NNC_Polyhedron(cs_box);
    Constraint_System cs;
    cs.set_space_dimension(3);
    cs.insert(1000*px >= -2.0);
    cs.insert(10*vx >= 0);
    cs.insert(10*vy >= 0);
    NNC_Polyhedron guard(cs);
    if(!guard.is_disjoint_from(box_poly)){
      Constraint_System cs_box;
      Variable Simu_time_temp(0);
      Variable vx_temp(1);
      Variable vy_temp(2);
      Variable px_temp(3);
      Variable py_temp(4);
      Variable ii_temp(5);
      Variable trans_temp(6);
      Variable one_temp(7);

      if(ptLower[0]<ptUpper[0]){
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptUpper[0]);
      }
      else{
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptUpper[0]);
      }

      if(ptLower[1]<ptUpper[1]){
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptUpper[1]);
      }
      else{
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptUpper[1]);
      }

      if(ptLower[2]<ptUpper[2]){
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptUpper[2]);
      }
      else{
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptUpper[2]);
      }

      if(ptLower[3]<ptUpper[3]){
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptUpper[3]);
      }
      else{
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptUpper[3]);
      }

      if(ptLower[4]<ptUpper[4]){
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptUpper[4]);
      }
      else{
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptUpper[4]);
      }

      if(ptLower[5]<ptUpper[5]){
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptUpper[5]);
      }
      else{
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptUpper[5]);
      }

      if(ptLower[6]<ptUpper[6]){
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptUpper[6]);
      }
      else{
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptUpper[6]);
      }

      if(ptLower[7]<ptUpper[7]){
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptUpper[7]);
      }
      else{
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptUpper[7]);
      }

      box_poly = NNC_Polyhedron(cs_box);
      Constraint_System cs_gd;
      cs_gd.set_space_dimension(16);
      Variable Simu_time_new(8);
      Variable vx_new(9);
      Variable vy_new(10);
      Variable px_new(11);
      Variable py_new(12);
      Variable ii_new(13);
      Variable trans_new(14);
      Variable one_new(15);
      box_poly.add_space_dimensions_and_embed(8);
      cs_gd.insert(10*ii_new==10*ii_temp + 32.0*vx_temp - 32.0*vy_temp);
      cs_gd.insert(10*vx_new==0);
      cs_gd.insert(10*vy_new==0);
      cs_gd.insert(trans_new==one_temp + trans_temp);
      cs_gd.insert(Simu_time_new==Simu_time_temp);
      cs_gd.insert(px_new==px_temp);
      cs_gd.insert(one_new==one_temp);
      cs_gd.insert(py_new==py_temp);
      cs_gd.insert(1000*px_temp >= -2.0);
      cs_gd.insert(10*vx_temp >= 0);
      cs_gd.insert(10*vy_temp >= 0);
      NNC_Polyhedron guard_reset(cs_gd);
      guard_reset.intersection_assign(box_poly);
      Variables_Set vars;
      vars.insert(Simu_time_temp);
      vars.insert(vx_temp);
      vars.insert(vy_temp);
      vars.insert(px_temp);
      vars.insert(py_temp);
      vars.insert(ii_temp);
      vars.insert(trans_temp);
      vars.insert(one_temp);
      guard_reset.remove_space_dimensions(vars);
      toRet.push_back(make_pair(guard_reset,2));
    }
  }
  if(curMode==1){
    Constraint_System cs_box;
    Variable vx(0);
    Variable vy(1);
    Variable px(2);
    Variable py(3);
    Variable trans(4);

    if(ptLower[1]<ptUpper[1]){
      cs_box.insert(mult_factor*vx>=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx<=mult_factor*ptUpper[1]);
    }
    else{
      cs_box.insert(mult_factor*vx<=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx>=mult_factor*ptUpper[1]);
    }

    if(ptLower[2]<ptUpper[2]){
      cs_box.insert(mult_factor*vy>=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy<=mult_factor*ptUpper[2]);
    }
    else{
      cs_box.insert(mult_factor*vy<=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy>=mult_factor*ptUpper[2]);
    }

    if(ptLower[3]<ptUpper[3]){
      cs_box.insert(mult_factor*px>=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px<=mult_factor*ptUpper[3]);
    }
    else{
      cs_box.insert(mult_factor*px<=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px>=mult_factor*ptUpper[3]);
    }

    if(ptLower[4]<ptUpper[4]){
      cs_box.insert(mult_factor*py>=mult_factor*ptLower[4]);
      cs_box.insert(mult_factor*py<=mult_factor*ptUpper[4]);
    }
    else{
      cs_box.insert(mult_factor*py<=mult_factor*ptLower[4]);
      cs_box.insert(mult_factor*py>=mult_factor*ptUpper[4]);
    }

    if(ptLower[6]<ptUpper[6]){
      cs_box.insert(mult_factor*trans>=mult_factor*ptLower[6]);
      cs_box.insert(mult_factor*trans<=mult_factor*ptUpper[6]);
    }
    else{
      cs_box.insert(mult_factor*trans<=mult_factor*ptLower[6]);
      cs_box.insert(mult_factor*trans>=mult_factor*ptUpper[6]);
    }

    box_poly = NNC_Polyhedron(cs_box);
    Constraint_System cs;
    cs.set_space_dimension(5);
    cs.insert(726542528005361.0*px + 1000000000000000*py >= 0);
    cs.insert(587785252292473.0*vx + 809016994374947.0*vy > 0);
    cs.insert(10*trans < 80);
    NNC_Polyhedron guard(cs);
    if(!guard.is_disjoint_from(box_poly)){
      Constraint_System cs_box;
      Variable Simu_time_temp(0);
      Variable vx_temp(1);
      Variable vy_temp(2);
      Variable px_temp(3);
      Variable py_temp(4);
      Variable ii_temp(5);
      Variable trans_temp(6);
      Variable one_temp(7);

      if(ptLower[0]<ptUpper[0]){
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptUpper[0]);
      }
      else{
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptUpper[0]);
      }

      if(ptLower[1]<ptUpper[1]){
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptUpper[1]);
      }
      else{
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptUpper[1]);
      }

      if(ptLower[2]<ptUpper[2]){
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptUpper[2]);
      }
      else{
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptUpper[2]);
      }

      if(ptLower[3]<ptUpper[3]){
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptUpper[3]);
      }
      else{
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptUpper[3]);
      }

      if(ptLower[4]<ptUpper[4]){
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptUpper[4]);
      }
      else{
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptUpper[4]);
      }

      if(ptLower[5]<ptUpper[5]){
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptUpper[5]);
      }
      else{
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptUpper[5]);
      }

      if(ptLower[6]<ptUpper[6]){
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptUpper[6]);
      }
      else{
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptUpper[6]);
      }

      if(ptLower[7]<ptUpper[7]){
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptUpper[7]);
      }
      else{
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptUpper[7]);
      }

      box_poly = NNC_Polyhedron(cs_box);
      Constraint_System cs_gd;
      cs_gd.set_space_dimension(16);
      Variable Simu_time_new(8);
      Variable vx_new(9);
      Variable vy_new(10);
      Variable px_new(11);
      Variable py_new(12);
      Variable ii_new(13);
      Variable trans_new(14);
      Variable one_new(15);
      box_poly.add_space_dimensions_and_embed(8);
      cs_gd.insert(100000000000000*ii_new==100000000000000*ii_temp + 774867751838096.0*vx_temp + 1066513964386099.0*vy_temp);
      cs_gd.insert(100000000000000*vx_new==-42329949064832.0*vx_temp - 195900368634417.0*vy_temp);
      cs_gd.insert(1000000000000*vy_new==-346343193165.0*vx_temp + 523299490640.0*vy_temp);
      cs_gd.insert(trans_new==one_temp + trans_temp);
      cs_gd.insert(Simu_time_new==Simu_time_temp);
      cs_gd.insert(px_new==px_temp);
      cs_gd.insert(one_new==one_temp);
      cs_gd.insert(py_new==py_temp);
      cs_gd.insert(726542528005361.0*px_temp + 1000000000000000*py_temp >= 0);
      cs_gd.insert(587785252292473.0*vx_temp + 809016994374947.0*vy_temp > 0);
      cs_gd.insert(10*trans_temp < 80);
      NNC_Polyhedron guard_reset(cs_gd);
      guard_reset.intersection_assign(box_poly);
      Variables_Set vars;
      vars.insert(Simu_time_temp);
      vars.insert(vx_temp);
      vars.insert(vy_temp);
      vars.insert(px_temp);
      vars.insert(py_temp);
      vars.insert(ii_temp);
      vars.insert(trans_temp);
      vars.insert(one_temp);
      guard_reset.remove_space_dimensions(vars);
      toRet.push_back(make_pair(guard_reset,1));
    }
  }
  if(curMode==1){
    Constraint_System cs_box;
    Variable vx(0);
    Variable vy(1);
    Variable px(2);
    Variable py(3);
    Variable trans(4);

    if(ptLower[1]<ptUpper[1]){
      cs_box.insert(mult_factor*vx>=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx<=mult_factor*ptUpper[1]);
    }
    else{
      cs_box.insert(mult_factor*vx<=mult_factor*ptLower[1]);
      cs_box.insert(mult_factor*vx>=mult_factor*ptUpper[1]);
    }

    if(ptLower[2]<ptUpper[2]){
      cs_box.insert(mult_factor*vy>=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy<=mult_factor*ptUpper[2]);
    }
    else{
      cs_box.insert(mult_factor*vy<=mult_factor*ptLower[2]);
      cs_box.insert(mult_factor*vy>=mult_factor*ptUpper[2]);
    }

    if(ptLower[3]<ptUpper[3]){
      cs_box.insert(mult_factor*px>=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px<=mult_factor*ptUpper[3]);
    }
    else{
      cs_box.insert(mult_factor*px<=mult_factor*ptLower[3]);
      cs_box.insert(mult_factor*px>=mult_factor*ptUpper[3]);
    }

    if(ptLower[4]<ptUpper[4]){
      cs_box.insert(mult_factor*py>=mult_factor*ptLower[4]);
      cs_box.insert(mult_factor*py<=mult_factor*ptUpper[4]);
    }
    else{
      cs_box.insert(mult_factor*py<=mult_factor*ptLower[4]);
      cs_box.insert(mult_factor*py>=mult_factor*ptUpper[4]);
    }

    if(ptLower[6]<ptUpper[6]){
      cs_box.insert(mult_factor*trans>=mult_factor*ptLower[6]);
      cs_box.insert(mult_factor*trans<=mult_factor*ptUpper[6]);
    }
    else{
      cs_box.insert(mult_factor*trans<=mult_factor*ptLower[6]);
      cs_box.insert(mult_factor*trans>=mult_factor*ptUpper[6]);
    }

    box_poly = NNC_Polyhedron(cs_box);
    Constraint_System cs;
    cs.set_space_dimension(5);
    cs.insert(-726542528005361.0*px + 1000000000000000*py <= 0);
    cs.insert(587785252292473.0*vx - 809016994374947.0*vy > 0);
    cs.insert(10*trans < 80);
    NNC_Polyhedron guard(cs);
    if(!guard.is_disjoint_from(box_poly)){
      Constraint_System cs_box;
      Variable Simu_time_temp(0);
      Variable vx_temp(1);
      Variable vy_temp(2);
      Variable px_temp(3);
      Variable py_temp(4);
      Variable ii_temp(5);
      Variable trans_temp(6);
      Variable one_temp(7);

      if(ptLower[0]<ptUpper[0]){
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptUpper[0]);
      }
      else{
        cs_box.insert(mult_factor*Simu_time_temp<=mult_factor*ptLower[0]);
        cs_box.insert(mult_factor*Simu_time_temp>=mult_factor*ptUpper[0]);
      }

      if(ptLower[1]<ptUpper[1]){
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptUpper[1]);
      }
      else{
        cs_box.insert(mult_factor*vx_temp<=mult_factor*ptLower[1]);
        cs_box.insert(mult_factor*vx_temp>=mult_factor*ptUpper[1]);
      }

      if(ptLower[2]<ptUpper[2]){
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptUpper[2]);
      }
      else{
        cs_box.insert(mult_factor*vy_temp<=mult_factor*ptLower[2]);
        cs_box.insert(mult_factor*vy_temp>=mult_factor*ptUpper[2]);
      }

      if(ptLower[3]<ptUpper[3]){
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptUpper[3]);
      }
      else{
        cs_box.insert(mult_factor*px_temp<=mult_factor*ptLower[3]);
        cs_box.insert(mult_factor*px_temp>=mult_factor*ptUpper[3]);
      }

      if(ptLower[4]<ptUpper[4]){
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptUpper[4]);
      }
      else{
        cs_box.insert(mult_factor*py_temp<=mult_factor*ptLower[4]);
        cs_box.insert(mult_factor*py_temp>=mult_factor*ptUpper[4]);
      }

      if(ptLower[5]<ptUpper[5]){
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptUpper[5]);
      }
      else{
        cs_box.insert(mult_factor*ii_temp<=mult_factor*ptLower[5]);
        cs_box.insert(mult_factor*ii_temp>=mult_factor*ptUpper[5]);
      }

      if(ptLower[6]<ptUpper[6]){
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptUpper[6]);
      }
      else{
        cs_box.insert(mult_factor*trans_temp<=mult_factor*ptLower[6]);
        cs_box.insert(mult_factor*trans_temp>=mult_factor*ptUpper[6]);
      }

      if(ptLower[7]<ptUpper[7]){
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptUpper[7]);
      }
      else{
        cs_box.insert(mult_factor*one_temp<=mult_factor*ptLower[7]);
        cs_box.insert(mult_factor*one_temp>=mult_factor*ptUpper[7]);
      }

      box_poly = NNC_Polyhedron(cs_box);
      Constraint_System cs_gd;
      cs_gd.set_space_dimension(16);
      Variable Simu_time_new(8);
      Variable vx_new(9);
      Variable vy_new(10);
      Variable px_new(11);
      Variable py_new(12);
      Variable ii_new(13);
      Variable trans_new(14);
      Variable one_new(15);
      box_poly.add_space_dimensions_and_embed(8);
      cs_gd.insert(100000000000000*ii_new==100000000000000*ii_temp + 774867751838096.0*vx_temp - 1066513964386099.0*vy_temp);
      cs_gd.insert(100000000000000*vx_new==-42329949064832.0*vx_temp + 195900368634417.0*vy_temp);
      cs_gd.insert(1000000000000*vy_new==346343193165.0*vx_temp + 523299490640.0*vy_temp);
      cs_gd.insert(trans_new==one_temp + trans_temp);
      cs_gd.insert(Simu_time_new==Simu_time_temp);
      cs_gd.insert(px_new==px_temp);
      cs_gd.insert(one_new==one_temp);
      cs_gd.insert(py_new==py_temp);
      cs_gd.insert(-726542528005361.0*px_temp + 1000000000000000*py_temp <= 0);
      cs_gd.insert(587785252292473.0*vx_temp - 809016994374947.0*vy_temp > 0);
      cs_gd.insert(10*trans_temp < 80);
      NNC_Polyhedron guard_reset(cs_gd);
      guard_reset.intersection_assign(box_poly);
      Variables_Set vars;
      vars.insert(Simu_time_temp);
      vars.insert(vx_temp);
      vars.insert(vy_temp);
      vars.insert(px_temp);
      vars.insert(py_temp);
      vars.insert(ii_temp);
      vars.insert(trans_temp);
      vars.insert(one_temp);
      guard_reset.remove_space_dimensions(vars);
      toRet.push_back(make_pair(guard_reset,1));
    }
  }
  return toRet;
}

