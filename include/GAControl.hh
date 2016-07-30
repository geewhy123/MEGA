/****************************************
 *file:GAControl.hh
 * header containing function prototypes 
 * and member variables of the GAControl class
 *
 * 
 * Gary Yan, Peter Zhang
 * June 13,2012
 ****************************************/


#ifndef GACONTROL_H
#define GACONTROL_H

#include <iostream>
#include <string>


#include "GAProblem.hh"
#include "parentSelectionParams.hh"
#include "recombinationParams.hh"
#include "mutationParams.hh"

using namespace std;
class GAProblem;
class GAControl{
  int popSize;
  int initialSeed;

  parentSelectionParams params;
  recombinationParams rparams;
  mutationParams mparams;

  string fileName;
  string fitnessAssignmentType;
  string recombinationType;
  string mutationType;
  string parentSelectionType;
  double distIndex;
  int nGens;
  
public:
   
  GAControl();
  ~GAControl();
  void setPopSize(int popSize);
  int getPopSize();
  void setInitialSeed(int seed);
  int getInitialSeed();

  parentSelectionParams getParentSelectionParams();
  void setParentSelectionParams(parentSelectionParams param);

  recombinationParams getRecombinationParams();
  void setRecombinationParams(recombinationParams rparams);

  mutationParams getMutationParams();
  void setMutationParams(mutationParams mparams);

  void showAll();
  void setOutputFile( string s);
   string getOutputFile();
  void setFitnessAssignmentType(string s);
string getFitnessAssignmentType();
  void setRecombinationType(string s);
  string getRecombinationType();
  void setMutationType(string s);
 string getMutationType();
  void setParentSelectionType(string s);
  string getParentSelectionType();
  void setDistIndex(double distIndex);
  double getDistIndex();

  void setnGens(int nGens);// max number of generations
  double getnGens();
};

#endif
