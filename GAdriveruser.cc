/****************************************
 *file:GAdriveruser.cc
 * A sample user input driver code that runs
 * a constrained GA optimization
 *
 * 
 * Gary Yan, Peter Zhang
 * June 13,2012
 ****************************************/


#include <iostream>
#include <armadillo>
#include "WindFarmEvaluation.hh"
 #include "include/GAControl.hh"
#include "include/GAProblem.hh"
#include <fstream>
#include <string>
#include <vector> 
#include <time.h>
#include <math.h>
#include <float.h>
#include "GA.hh"
using namespace std;
using namespace arma;
#define PI 3.1415926535897932 // less expensive than acos(-1)


//----------------define constraint functions----------------

//g1 function
// inequality constraint g(x) <= 0
// input through x,y coordinates of the turbine placements 
// this sample constraint function restricts the placement
// in the square A = {(x,y) |  292.89 < x < 1707.10 and 292.89 < y < 1707.10 } 
// Output:
// - returns values of g(x) in inequality constraint g(x) <= 0
// Input:
// - vec x : vector of input x values for a particular wind farm layout
//           the first half of the vector are all x values, followed by all y values
// - void * inequalityConstraintParams : array of parameters pertaining to this constraint ( not used)
double g1(vec x, void* inequalityConstraintParams){
  int count = 0;// count stores the number of turbines placed in the forbidden region 

  //restrict square A = {(x,y) |  292.89 < x < 1707.10 and 292.89 < y < 1707.10 } 
  for(int i=0;i< x.n_elem/2 ; i++){
    count += ((x(i) <((2000.0+sqrt(2)*1000)/2.0)) && (x(i) > ((2000.0-sqrt(2)*1000)/2.0)) && (x(i+x.n_elem/2)< ((2000.0+sqrt(2)*1000)/2.0)) && (x(i+x.n_elem/2) > ((2000.0-sqrt(2)*1000)/2.0)));
  }
  return count;

}


int main()
{
    
  //set seeds (optional) for 10 runs
  // for random, set seeds = time(NULL)
  vec seeds = ones(10); // assuming a max number of runs to be 10
  seeds(0) = 28;
  seeds(1) = 8574;
  seeds(2) = 1023;
  seeds(3) = 1995;
  seeds(4) = 4077;
  seeds(5) = 1002;
  seeds(6) = 9938;
  seeds(7) = 6760;
  seeds(8) = 203;
  seeds(9) = 1331;
  
  GAProblem *G = new GAProblem;
  GAControl *C = new GAControl();
  


  //Comment below for reading from file
  
  int nTurbines = 15;
  int nInputs = 2*nTurbines;
  int nObj = 2;
  int rootn = (int)round(sqrt(nInputs)+1.5);
  int popSize = 2*rootn*rootn; // nGen equations from Pen 
  G->setnInputs(nInputs);
  G->setnObj(nObj);
  //endComment


  //Uncomment below for reading from file 
  // and Comment out set functions
  // readFile read parameters are incomplete..
  // missing:
  // - setConstraintType
  // - setnIneqConstr
  // - setnEqConstr
  // - setPenaltyCoefficients
  //
  //   G->readFile("files/_input.dat",C);// reads parameters from file (old version) - not used because parameters overwritten for wind farm problem
  //endUncomment
    


  // Setting up turbine characteristics
   int q=1;
   int v = 1;
  double z = 60.0; // turbine height, homogeneous for now
  double drot = 40.0; // turbine diameter
  double z0 = 0.3; // ground roughness factor
  
  cout << "Generating wind speed data..." << endl;
  mat windSpeeds (q, v);
  for (int w = 0; w < q; w ++)
    {
      for (int l = 0; l < v; l ++)
	{
	  if (q == 36)
	    {
	      if (l % 3 == 0)
		windSpeeds(w,l) = 8.0;
	      else if (l % 3 == 1)
		windSpeeds(w,l) = 12.0;
	      else if (l % 3 == 2)
		windSpeeds(w,l) = 17.0;
	    }
	  else
	    windSpeeds(w,l) = 12.0;
	}
		
    }
  cout << "Generating wind prob data..." << endl;
  // new prob
  mat windProb (q,v);
  for (int w = 0; w < q; w ++)
    {
      for (int l = 0; l < v; l ++)
	{
	  if (q == 36)
	    {
	      if (l % 3 == 0)
		windProb(w,l) = 0.00485;
	      else if (l % 3 == 1)
		windProb(w,l) = 0.008;
	      else if (l % 3 == 2)
		windProb(w,l) = 0.011;
	    }
	  else
	    windProb(w,l) = 1.0 / q / v;
	}
    }
       
  if (q == 36)
    {
      //special directions
      windProb(26,1) = 0.010135;
      windProb(27,1) = 0.012135;
      windProb(28,1) = 0.014595;
      windProb(29,1) = 0.014185;
      windProb(30,1) = 0.019255;
      windProb(31,1) = 0.014185;
      windProb(32,1) = 0.014595;
      windProb(33,1) = 0.012135;
      windProb(34,1) = 0.010135;
               
      windProb(26,2) = 0.013;
      windProb(27,2) = 0.01695;
      windProb(28,2) = 0.0185;
      windProb(29,2) = 0.0298;
      windProb(30,2) = 0.035;
      windProb(31,2) = 0.0298;
      windProb(32,2) = 0.0185;
      windProb(33,2) = 0.01695;
      windProb(34,2) = 0.013;
    }
  vec windDirections (q);
  for (int w = 0; w < q; w ++)
    windDirections(w) = 2.0 * PI / q * w; // for example, -PI is wind blowing from north to south
  


  ///-Comment out below if using readFile
  double turbineFullPower = 0.0;
  for (int w = 0; w < q; w ++)
    for (int l = 0; l < v; l ++)
      turbineFullPower += powerCurve(windSpeeds(w,l)) * windProb(w,l);
  
    vec outputLowerLims(nObj); 
    outputLowerLims(0) = -1 * turbineFullPower * nTurbines;
    if(nObj >= 2)
      outputLowerLims(1) = 45;//noise 
    vec outputUpperLims(nObj);
    outputUpperLims(0) = 0.0;// lower limits of objectives
    if(nObj >= 2)
      outputUpperLims(1) = 70;
    
    
        
    vec inputLowerLims = zeros(nInputs);// upper limits for objectives
    vec inputUpperLims = 2000*ones(nInputs);
    
    int nParents = popSize/2;
    int nOffspring = popSize;
    
    int distIndex = 2;
    double mutationProb = 0.05;
    double recombProb = 0.8;
    
    parentSelectionParams pparam;
    pparam.nParents = nParents;
    
    recombinationParams rparam;
    rparam.recombProb = recombProb;
    rparam.nOffspring = nOffspring;
    
    mutationParams mparam;
    mparam.mutationProb = mutationProb;
    double d_lim = 0.005;
    double L = 100;
 
    int nGens = 200*rootn;
    C->setOutputFile("maxEnergy");   
    C->setParentSelectionParams(pparam);
    C->setRecombinationParams(rparam);
    C->setMutationParams(mparam);
    C->setDistIndex(distIndex);
    G->setOutputLowerLims(outputLowerLims);
    G->setOutputUpperLims(outputUpperLims);
    G->setInputLowerLims(inputLowerLims);
    G->setInputUpperLims(inputUpperLims);
    C->setPopSize(popSize);

    C->setInitialSeed(seeds(0));

    G->setDLim(d_lim);
    G->setL(L);

    C->setFitnessAssignmentType("NSGA-II");
    C->setRecombinationType("SBX");
    C->setMutationType("POLYNOMIAL");
    C->setParentSelectionType("BIN-TOUR");
    C->setnGens(nGens);    
    G->setConstraintType("PENALTY");
    

    int nIneqConstr = 1;
    int nEqConstr = 0;
    G->setnIneqConstr(nIneqConstr);
    G->setnEqConstr(nEqConstr);

    mat penalty(nIneqConstr+nEqConstr,nObj);
    penalty(0,0) = 10000;
    penalty(0,1) = 60;
    G->setPenaltyCoefficients(penalty);
    ///---readfileendcomment
    
  



  
  //objFunctionParameters 1
      objP * p = new objP;

      p->windDirections = windDirections;
      p->windSpeeds = windSpeeds;
      p->windProb = windProb;
      p->drot = drot;
      p->z = z;
      p->z0 = z0;
      

      //objFunctionParameters 2
      objP2 * p2 = new objP2;
      
           
      //arrays of function pointers
      double (*objFuncPointers[2])(vec x,void* P);// objectives
      double (*QPtrs[2])(double f, vec g, vec h, int iGen);//penalized objectives
      double (*ineqConstraintPointers[1])(vec x, void* Pg);//ineq constraints
      
      
      //objects for objectives 
      void *objFuncParameterPointers[2];
      objFuncParameterPointers[0] = p;
      objFuncParameterPointers[1] = p2;
     
      objFuncPointers[0] = &f1;
      objFuncPointers[1] =  &f2;
     
      ineqConstraintPointers[0] = &g1;
      void* ineqConstraintParameterPointers[1];
     

      G->setObjFuncParameterPointers(objFuncParameterPointers);// set objective objects
      G->setObjFuncPointers(objFuncPointers); // set objective functions

      G->setIneqConstraintParameterPointers(ineqConstraintParameterPointers);//set ineqConstraint objects
      G->setIneqConstraintPointers(ineqConstraintPointers);//set ineqConstraint functions



      //if there are equalityConstraints
      //double (*eqConstraintPointers[nEqConstr])(vec x, void* Ph);//eq constraints
      //void * eqConstraintParameterPointers[nEqConstr];
      //G->setEqConstraintParameterPointers(EqConstraintParameterPointers);
      //G->setEqConstraintPointers(EqConstraintPointers);



      //--------displaying GAProblem and GAControl class variables-----
      //G->showAll();
      //C->showAll();

      
      //-------------calling the optimization one generation at a time 
      G->evolveOneGen(C,true);
      ofstream ss;
      ss.open("secondGen.txt");
      ss<<G->getCurrentPop();
      ss.close();
      
      G->evolveOneGen(C,false); 
      ss.open("thirdGen.txt");
      ss<<G->getCurrentPop();
      ss.close();
      
      G->evolveOneGen(C,false);
      ss.open("fourthGen.txt");
      ss<<G->getCurrentPop();
      ss.close();
      //...etc


      //--------------calling the entire GA optimization -------------
      G->solve(C);
      
      return 0;
}
