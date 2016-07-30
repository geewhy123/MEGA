/****************************************
 *file:GAProblem.hh
 * header containing function prototypes and 
 * member variables of the GAProblem class
 * 
 *
 * 
 * Gary Yan, Peter Zhang
 * June 13,2012
 ****************************************/


#ifndef GAPROBLEM_H 
#define GAPROBLEM_H 
#include <iostream>
#include <armadillo>
#include "GAControl.hh"
#include "mutationParams.hh"
#include "recombinationParams.hh"
#include "parentSelectionParams.hh"

using namespace std;
using namespace arma;
class GAControl;
class GAProblem
{
  mat pop;
  mat popObj;
  
  mat currentPop;
  int currentGen;
  int nInputs;
  int nObj;
  vec inputLowerLims;
  vec inputUpperLims;
  vec outputLowerLims;
  vec outputUpperLims;

  double d_lim;
  int L; // see initializePopulation for description of converged and d_lim
  vec maxCrowdingDist; // vector of L numbers from past L generations, each number is the max crowding dist in each geneneration.

  mat penaltyCoefficients;// matrix of size (nIneqConstr + nEqConstr ) x nObj, containing the penalty coefficient for each combination.  
  //Ineq constraints are placed in rows first, followed by eq Constraints

  double (**objFuncPointers)(vec x, void* P);
  double (**ineqConstraintPointers)(vec x, void* Pg);
  double (**eqConstraintPointers)(vec x, void* Ph);

  void **objFuncParameterPointers;
  void **ineqConstraintParameterPointers;
  void **eqConstraintParameterPointers;


  string constraintType;
  int nIneqConstr;
  int nEqConstr;
public:
	// d_lim and L used to calculate the crowding distance stabilization
	// d_lim is allowed upper bound for SD of the max. crowding distances in the past L genertions 
  mat initializePopulation(GAControl C);//, int runCount, double d_lim, int L);
  void setnInputs(int nInputs);
  int getnInputs();
  void setInputLowerLims(vec inputLowerLims);
  vec getInputLowerLims();
  void setInputUpperLims(vec inputUpperLims);
  vec getInputUpperLims();
  
  void setOutputLowerLims(vec outputLowerLims);
  vec getOutputLowerLims();
  void setOutputUpperLims(vec outputUpperLims);
  vec getOutputUpperLims();
  void setnObj(int nObj);
  int getnObj();
  
  bool hasConverged();
  
  //	void setDistIndex(double distIndex);
  //	double getDistIndex();
	
  vec evaluateObjectives(vec x, double (**f)(vec x,void* p/*,void** Pg = NULL, void** Ph= NULL*/),void** objParams);
  vec evaluateIneqConstraints(vec x, double (**g)(vec x,void* p),void** constrParams);
  vec evaluateEqConstraints(vec x, double (**h)(vec x,void* p),void** constrParams);


  mat sortByVec(vec sortBy,mat input, int asc);  

  mat constraintOperation(GAControl* C,mat popAndObj,mat ineqConstraints,mat eqConstraints);
  mat penalizeObjectives(GAControl* C,mat popAndObj,mat ineqConstraints,mat eqConstraints);

  mat fitnessAssignment(GAControl* C,mat popInput,mat objVals, bool isPureNewPop = false, int genCount = 0);
  mat fitnessAssignmentNSGAII(mat popInput, mat objVals, bool isPureNewPop = false, int genCount = 0);//, const string type);

  mat parentSelection(GAControl* C, mat popAndFitness);	
  mat parentSelectionBINTOUR(GAControl* C,mat popAndFitness);
  mat recombination(GAControl* C, mat parents);
  mat recombinationSBX(GAControl* C,mat parents);
  mat mutation(GAControl* C, mat offspring);	
  mat mutationPOLY(GAControl* C, mat offspring);
	
  double NDHV(mat points,vec refMin,vec refMax);
  
  void showAll();
  
  void evolveOneGen(GAControl* C,bool isNewRun);
  void solve(GAControl* C);
  
  double getDLim();
  void setDLim(double DLim);
  int getL();
  void setL(int L);
  vec getMaxCrowdingDist();
  void setMaxCrowdingDist(vec mcd);
  
  void setObjFuncPointers(double(**objFuncPointers)(vec x,void* P));//objectives
  double (**getObjFuncPointers())(vec x, void* P);
  void setIneqConstraintPointers(double(**ineqConstraintPointers)(vec x,void* Pg));//equality/inequalityconstraints
  double (**getIneqConstraintPointers())(vec x, void* Pg);
  void setEqConstraintPointers(double(**eqConstraintPointers)(vec x,void* Ph));
  double (**getEqConstraintPointers())(vec x, void* Ph);
  
  void setObjFuncParameterPointers(void** objFuncParameterPointers);
  void **getObjFuncParameterPointers();
  
  void setIneqConstraintParameterPointers(void** ineqConstraintParameterPointers);
  void **getIneqConstraintParameterPointers();
  
  void setEqConstraintParameterPointers(void** eqConstraintParamaterPointers);
  void **getEqConstraintParameterPointers();
  void readFile(const char *s, GAControl* C);
  
  
  void setCurrentPop(mat pop);
  mat getCurrentPop();

  void setConstraintType(string s);
  string getConstraintType();


  void setnIneqConstr(int n);
  int getnIneqConstr();
  void setnEqConstr(int n);
  int getnEqConstr();

  void setCurGen(int iGen);
  int getCurGen();

  mat getPenaltyCoefficients();
  void setPenaltyCoefficients(mat P);
}; 



#endif
