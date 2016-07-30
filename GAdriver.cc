/****************************************
 *file:GAdriver.cc
 * 
 * 
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
using namespace std;
using namespace arma;
#define PI 3.1415926535897932 // less expensive than acos(-1)
#define EPS 1e-4


double powerCurve(double speed)
{
  return 0.3 * pow(speed,3);
}


/* user defined objective functions f1, f2,...
 * inputs : x - vector of inputs
 *          P - pointer to objective parameters struct
 *
 */

// f1 function
// first objective function of the wind farm optimization problem - energy (power) maximization
// receptors are placed uniformly 40m apart along the at a heighe of 1.75m 
// along the edge of the 2000m x 2000m wind farm
// Outputs:
// - double: total power generated by turbine layout
// Inputs: 
// - vec x : vector of coordinate inputs(first half x coordinates, followed by y coordinates)
// - void *objFuncParam: array of parameters pertaining to this objective function
double f1(vec x,void* objFuncParam){
  objP* nP = static_cast<objP*>(objFuncParam);
    
  WindFarmEvaluation e;
  
  double power = e.calculatePower(x, nP->windDirections, nP->windSpeeds, nP->windProb, nP->z, nP->drot, nP->z0, powerCurve);
  
  return power*-1 ;// turn problem into minimization problem

}



// f2 function
// second objective function of the wind farm optimization problem - noise minimization
// receptors are placed uniformly 40m apart along the at a heighe of 1.75m 
// along the edge of the 2000m x 2000m wind farm
// Outputs:
// - double: maximum noise detected by the receptors
// Inputs: 
// - vec x : vector of coordinate inputs(first half x coordinates, followed by y coordinates)
// - void *objFuncParam: array of parameters pertaining to this objective function
double f2(vec x,void* objFuncParam){
  objP2* nP = static_cast<objP2*>(objFuncParam);
      
  WindFarmEvaluation e;
  
  vec rcoordinates(480);
  
  for(int i=0; i < 40; i++){
    rcoordinates(i) = 50*i;
    rcoordinates(i+160) = 0;
    rcoordinates(i+320) = 1.75;
  }
  for (int i=40 ; i < 80;i++){
    rcoordinates(i) = 2000;
    rcoordinates(i+160) = 50*(i-40);
    rcoordinates(i+320) =1.75;
  }
  for(int i=80; i <120;i++){
    rcoordinates(i) = 2000-50*(i-80);
    rcoordinates(i+160) = 2000;
    rcoordinates(i+320) = 1.75;
  }
  for(int i=120; i < 160; i++){
    rcoordinates(i) = 0;
    rcoordinates(i+160) =  2000-50*(i-120);   
    rcoordinates(i+320) =1.75;
  }  
  
  vec noises = e.calculateSoundField(x, 60.0, rcoordinates, 1, 1, 1);  
  
  return noises.max();
  
}