/****************************************
 *file:GA.cc
 * contains functions of the GAProblem class 
 * involved in the general GA optimization algorithm
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
#include "GA.hh"
using namespace std;
using namespace arma;

#define PI 3.1415926535897932 // less expensive than acos(-1)
#define EPS 1e-4

// read parameters from file and writes into GAProblem and GAControl instances.
// called in driver as alternative to defining parameters in application
// Output: 
// - none (calls GAProblem set methods)
// Inputs: 
// - const char* s : file name string
// - GAControl *C : pointer to GAControl instance
//
// Text file named s contains GA problem parameters in the order specified 
// with values following the name separated by a colon ':'
// '*' denotes repeated values for vectors
void GAProblem::readFile(const char* s,GAControl *C){
  ifstream in;
  string str;
  in.open(s);
  if(!in){
    cerr<< "ERROR in readFile(): CANNOT READ INPUT FILE \n" << "run aborted." << endl;
    exit(1);
  }
  
  //-----------nInputs-----------                                                                                                                             
  in >> str;
  size_t found = str.find(":");
  str = str.substr(found+1);
  int nInputs = atoi(str.c_str());
  this->setnInputs(nInputs);
  
  //-----------nObj---------------
  in >> str;
  found = str.find(":");
  str = str.substr(found+1);
  int nObj = atoi(str.c_str());
  this->setnObj(nObj);
  
  size_t end;
  vec lim= zeros(nInputs);
  string stri;
  
  //-----------inputLowerLims----
  
  string stru,strui;
  in >> stru;
  str = stru;  
  found = stru.find(":");
  stru = stru.substr(found+1);
  
  if(stru[stru.size()-1]=='*'){
    end = stru.find("*");
    strui = stru.substr(0,end);
    double value = atof(strui.c_str());
    this->setInputLowerLims(value*ones(nInputs));
  }
  else{

    for(int i=0; i < nInputs; i++){
    found = str.find(":");
    str = str.substr(found+1);
    end = str.find(":");
    stri = str.substr(0,end);
    lim(i) = atof(stri.c_str());
    }
    this->setInputLowerLims(lim);
  }
    
  //-----------inputUpperLims
  in >> stru;
  str = stru;
  found = stru.find(":");
  stru = stru.substr(found+1);

  if(stru[stru.size()-1]=='*'){
    end = stru.find("*");
    strui = stru.substr(0,end);
    double value = atof(strui.c_str());
    this->setInputUpperLims(value*ones(nInputs));
  }
  else{
    
    for(int i=0; i < nInputs; i++){
      found = str.find(":");
      str = str.substr(found+1);
      end = str.find(":");
      stri = str.substr(0,end);
      lim(i) = atof(stri.c_str());
    }
    this->setInputUpperLims(lim);
  }

  
  //-----------outputLowerLims---
  vec outLim(nObj);
  in >> stru;
  str = stru;
  found = stru.find(":");
  stru = stru.substr(found+1);

  if(stru[stru.size()-1]=='*'){
    end = stru.find("*");
    strui = stru.substr(0,end);
    double value = atof(strui.c_str());
    this->setOutputLowerLims(value*ones(nObj));
  }
  else{

    for(int i=0; i < nObj; i++){
      found = str.find(":");
      str = str.substr(found+1);
      end = str.find(":");
      stri = str.substr(0,end);
      outLim(i) = atof(stri.c_str());
    }
    this->setOutputLowerLims(outLim);
  }
  
  
  //-----------outputUpperLims---

  in >> stru;
  str = stru;
  found = stru.find(":");
  stru = stru.substr(found+1);

  if(stru[stru.size()-1]=='*'){
    end = stru.find("*");
    strui = stru.substr(0,end);
    double value = atof(strui.c_str());
    this->setOutputUpperLims(value*ones(nObj));
  }
  else{

    for(int i=0; i < nObj; i++){
      found = str.find(":");
      str = str.substr(found+1);
      end = str.find(":");
      stri = str.substr(0,end);
      outLim(i) = atof(stri.c_str());
    }
    this->setOutputUpperLims(outLim);
  }

  //------fitnessAssignmentType-----
  in >> str;
  found = str.find(":");
  str = str.substr(found+1);
  C->setFitnessAssignmentType(str); 

  //------parentSelectionType-----                                                                                                                  
  in >> str;
  found = str.find(":");
  str = str.substr(found+1);
  C->setParentSelectionType(str);

  //------recombinationType-----                                                                                                                  
  in >> str;
  found = str.find(":");
  str = str.substr(found+1);
  C->setRecombinationType(str);

  //------mutationType-----                                                                                                                  
  in >> str;
  found = str.find(":");
  str = str.substr(found+1);
  C->setMutationType(str);
  //------------nParents--------                                                                  

  parentSelectionParams pparam; 

  in >> str;                                                                                                                                           found = str.find(":");                                                                                                                               str = str.substr(found+1);                                                                                                                           pparam.nParents = atoi(str.c_str());                                                             
  C->setParentSelectionParams(pparam);                                                              

  // ---------recombinationParams
  recombinationParams rparam;
  in >> str;                                                                                                                                           found = str.find(":");                                                                                                                               str = str.substr(found+1);                                                                                                                           rparam.recombProb = atof(str.c_str());                                                                                                             

  //-----------nOffspring-------                                      
   in >> str;                                                                                                                                           found = str.find(":");                                                                                                                               str = str.substr(found+1);                                                                                                                           rparam.nOffspring = atoi(str.c_str());                                                                                                               C->setRecombinationParams(rparam);   


   //-----------mutationParams----   
   mutationParams mparam;    
   in >> str;          
   found = str.find(":");                                                                                                                       
   str = str.substr(found+1);                                                                                                                      
   mparam.mutationProb = atof(str.c_str());                                                                                                             C->setMutationParams(mparam);
   
   //-----------distIndex------                                                                                                                      
   in >> str;                                                                                                                                           found = str.find(":");                                                                                                                               str = str.substr(found+1);           
   C->setDistIndex(atof(str.c_str()));       
   

   //-----------popSize-----------                                                                                                                    
   in >> str;                                                                                                                                       
   found = str.find(":");                                                                                                                            
   str = str.substr(found+1);                                                                                                                           C->setPopSize(atoi(str.c_str()));  
  
   //-----------d_lim-------------
   in >> str;
   found = str.find(":");
   str = str.substr(found+1);
   this->setDLim(atof(str.c_str()));
   //-----------L-----------------
   in >> str;
   found = str.find(":");
   str = str.substr(found+1);
   this->setL(atoi(str.c_str()));
   
   //-----------randomSeed------                                                                                                                       
   in >> str;
   found = str.find(":");
   str = str.substr(found+1);
   
   if(strcmp(str.c_str(),"RANDOM")==0)
     C->setInitialSeed(time(0));
   else
     C->setInitialSeed(atoi(str.c_str()));
   
   //-----------outputFile-------
   in >> str;
   found = str.find(":");
   str = str.substr(found+1);
   C->setOutputFile(str);

   //------constraintType-----                                                                                                                  
   in >> str;
   
   found = str.find(":");
   str = str.substr(found+1);
   this->setConstraintType(str);
   
   //-----------nIneqConstr-----------                                                                                  
   in >> str;
   found = str.find(":");
   str = str.substr(found+1);
   int nIneqConstr = atoi(str.c_str());
   this->setnIneqConstr(nIneqConstr);
   
   //-----------nEqConstr-----------                                                                                  
   in >> str;

   found = str.find(":");
   str = str.substr(found+1);
   int nEqConstr = atoi(str.c_str());

   this->setnEqConstr(nEqConstr);
   


   //-----------nGens-----------                                                                                  
   in >> str;

   found = str.find(":");
   str = str.substr(found+1);
   int nGens = atoi(str.c_str());

   C->setnGens(nGens);
   
   //-----------penaltyCoefficients-----
   mat penaltyCoefficients = zeros(nIneqConstr+nEqConstr,nObj);
   in >> stru;
   str = stru;
   found = stru.find(":");
   stru = stru.substr(found+1);

   /*  if(stru[stru.size()-1]=='*'){
    end = stru.find("*");
    strui = stru.substr(0,end);
    double value = atof(strui.c_str());
    this->setInputUpperLims(value*ones(nInputs));
  }
  else{*/
   
   for(int i = 0; i < nIneqConstr+nEqConstr; i++){
     if(i > 0) found = str.find(";");
     for(int j=0; j < nObj; j++){
       found = str.find(":");
   
       str = str.substr(found+1);
       end = str.find(":");
   
       stri = str.substr(0,end);
       penaltyCoefficients(i,j) = atof(stri.c_str());
       //}
      }
    
   }
   this->setPenaltyCoefficients(penaltyCoefficients);

}


// helper function for appending data outputs to specified file 
// Outputs:
// - none
// Inputs:
// - const char * s : file name string
// - mat data : matrix of data values to be written to file
void writeToFile(const char* s, mat data){
  ofstream myfile;
  myfile.open(s,ios::app);
  myfile<< data;
  myfile.close();
  return; 
}


// Propagates one generation of the genetic algorithm.
// result is saved in member variable : currentPopulation of
// GAProblem class
// Outputs:
// - none
// Inputs:
// - GAControl *C : pointer to GAControl instance containing problem parameters
// - bool isNewRun: input true if it is first generation ( need to randomly initialize),
//   false otherwise
void GAProblem::evolveOneGen(GAControl *C,bool isNewRun){
  
    time_t startTime, endTime;
    startTime = time(NULL);

    // retrieve GA problem parameters         
    int nInputs = this->getnInputs();
    int popSize = C->getPopSize();
    int nGens = C->getnGens();
    
    int nObj = this->getnObj();		
    vec outputLowerLims = this->getOutputLowerLims(); 
    vec outputUpperLims = this->getOutputUpperLims();
    vec inputLowerLims = this->getInputLowerLims();
    vec inputUpperLims = this->getInputUpperLims();
    
    int nParents = C->getParentSelectionParams().nParents;
    int nOffspring =C->getRecombinationParams().nOffspring;
    double recombProb = C->getRecombinationParams().recombProb;    
    int distIndex = C->getDistIndex();
    double mutationProb = C->getMutationParams().mutationProb;
    
    double (**f)(vec x, void* P) = this->getObjFuncPointers();
    void** p = this->getObjFuncParameterPointers();

    int iGen = 1;
    this->setCurGen(iGen);		
        				
    mat archive = 1e20*ones(10*popSize,nInputs+nObj+2);//preallocated matrix size of best solutions
    
    mat initPop;
    mat popObj = zeros<mat>(popSize,nObj);
    mat joined;//initial input and objective values    

    //initialize if it is first generation
    if(isNewRun){
        initPop =  this->initializePopulation(*C);//  we don't need to store all output file, we only need the seed
              
	//evaluate objectives
	for(int i=0; i < popSize;i++){
	    popObj.row(i) = trans(this->evaluateObjectives(trans(initPop.row(i)),f,p));
	}
            
	joined = join_rows(initPop,popObj); // form one matrix of size (popSize x (nInputs + nObj)) 
    }
    else{
        initPop = this->getCurrentPop().rows(0,popSize-1).cols(0,nInputs-1);//get previous population as starting point
	popObj = this->getCurrentPop().rows(0,popSize-1).cols(nInputs,nInputs+nObj-1);//get previous population objectives
	joined = join_rows(initPop,popObj); // form one matrix 
    }

    mat ineqConstraints, eqConstraints; // matrix of inequality and equality constraints of size (popSize x nIneqConstr) and (popSize x nEqConstr)
    
    int nIneqConstraints = this->getnIneqConstr();
    
    //if there are inequality constraints, evaluate them
    if(nIneqConstraints >0){
      ineqConstraints = zeros(popSize,nIneqConstraints); 
      double (** ineqConstraintPointers)(vec x, void* ineqConstraintParameterPointers) = this->getIneqConstraintPointers();
      void** ineqConstraintParameterPointers = this->getIneqConstraintParameterPointers(); 
      for(int i=0; i < popSize; i++){
       	ineqConstraints.row(i) = trans(this->evaluateIneqConstraints(trans(initPop.row(i)),ineqConstraintPointers,ineqConstraintParameterPointers));
      }
    }

    int nEqConstraints = this->getnEqConstr();
    //if there are equality constraints, evaluate them
    if(nEqConstraints >0){
      eqConstraints = zeros(popSize,nEqConstraints);
      double (** eqConstraintPointers)(vec x, void* eqConstraintParameterPointers) = this->getEqConstraintPointers();
      void** eqConstraintParameterPointers = this->getEqConstraintParameterPointers();
      for(int i=0; i < popSize; i++){
	eqConstraints.row(i) = trans(this->evaluateEqConstraints(trans(initPop.row(i)),eqConstraintPointers,eqConstraintParameterPointers));
      }
    }


    //call constraint handling method before assigning fitness (currently only has penalty method)
    joined = this->constraintOperation(C,joined,ineqConstraints,eqConstraints);

    //evaluate fitess of population
    mat fitness = this->fitnessAssignment(C,joined.cols(0,nInputs-1),joined.cols(nInputs,nInputs+nObj-1));
    
    // setup to keep best solutions in archive
    archive.rows(0,popSize-1).cols(0,nInputs+nObj) = fitness.cols(0,nInputs+nObj);
      

    // Convergence criteria: Either total number of generation is reached, or 
    // NDHV+CrowdingDist have both stabilized
      
    time_t genStartTime, genEndTime;
    genStartTime = time(NULL);

    // --------main loop in GA -----------------       
    
    mat parents = this->parentSelection(C,fitness); // select parents
    
    mat recomb = this->recombination(C,parents); // recombination for offspring
    
    mat mutate = this->mutation(C,recomb); // mutation operator
    
    time_t objStartTime, objEndTime;
    objStartTime = time(NULL);

    mat newObj(mutate.n_rows,nObj); // evaluate objectives of new population
    for(int i=0; i < mutate.n_rows;i++){
        newObj.row(i) = trans(this->evaluateObjectives(trans(mutate.row(i)),f,p));
    }
    objEndTime = time(NULL);
    
    mat newPopObj = join_rows(mutate,newObj); // merge old and new populations  
    mat newPopObjJoined(mutate.n_rows+parents.n_rows,newPopObj.n_cols);// join population with objectives
    newPopObjJoined = join_cols(newPopObj,fitness.cols(0,nInputs+nObj-1));

    mat fit = zeros(newPopObjJoined.n_rows,nInputs+nObj+2); // fit will contain sorted merged population

    //re-evaluate objectives (if objectives are function of generation number, cannot reuse old objective values)
    //for improved runtime, only do this for dynamic penalty, and modify for other cases
    //
    if(1){//this->getConstraintType().compare("PENALTY-DYNAMIC")==0{
      mat newPopNewObj = zeros(newPopObjJoined.n_rows,nInputs+nObj);
      newPopNewObj.cols(0,nInputs-1) = newPopObjJoined.cols(0,nInputs-1);
      for(int i=0; i < newPopNewObj.n_rows; i++){
	newPopNewObj.row(i).cols(nInputs,nInputs+nObj-1) = trans(this->evaluateObjectives(trans(newPopObjJoined.cols(0,nInputs-1).row(i)),objFuncPointers,objFuncParameterPointers));	  
	}
      
    //if there are inequality constraints, evaluate them      
      if(nIneqConstraints >0){
	ineqConstraints = zeros(newPopObjJoined.n_rows,nIneqConstraints); 
	
	for(int i=0; i < ineqConstraints.n_rows; i++){
	  ineqConstraints.row(i) = trans(this->evaluateIneqConstraints(trans(newPopObjJoined.cols(0,nInputs-1).row(i)),ineqConstraintPointers,ineqConstraintParameterPointers));
	}
      }
      
    //if there are equality constraints, evaluate them
      if(nEqConstraints >0){
	eqConstraints = zeros(newPopObjJoined.n_rows,nEqConstraints);
	for(int i=0; i < eqConstraints.n_rows; i++){
	  eqConstraints.row(i) = trans(this->evaluateEqConstraints(trans(newPopObjJoined.cols(0,nInputs-1).row(i)),eqConstraintPointers,eqConstraintParameterPointers));
	  }
      }
    
      //constraint handling operator before assigning fitness  
      newPopNewObj = this->constraintOperation(C,newPopNewObj,ineqConstraints,eqConstraints);

      //assign fitness
      fit = this->fitnessAssignment(C,newPopNewObj.cols(0,nInputs-1),newPopNewObj.cols(nInputs,nInputs+nObj-1),true,0);
    }
    
    this->setCurrentPop(fit);// set new population as current
    mat zrs = zeros(newPopObj.n_rows,2);
    mat largeArchive = join_cols(join_rows(newPopObj,zrs), archive); // merged archive of current best solutions
    
    /////re-evaluation slows down optimization but is needed for dynamic penalty
    for(int i=0; i < largeArchive.n_rows && (largeArchive(i,0) < 1e10); i++){
      largeArchive.row(i).cols(nInputs,nInputs+nObj-1) = trans(this->evaluateObjectives(trans(largeArchive.cols(0,nInputs-1).row(i)),objFuncPointers,objFuncParameterPointers));
    }
          
    //if there are inequality constraints, evaluate them
    if(nIneqConstraints >0){
      ineqConstraints = zeros(largeArchive.n_rows,nIneqConstraints); 
      for(int i=0; i < ineqConstraints.n_rows; i++){
	ineqConstraints.row(i) = trans(this->evaluateIneqConstraints(trans(largeArchive.cols(0,nInputs-1).row(i)),ineqConstraintPointers,ineqConstraintParameterPointers));
      }
    }
    //if there are equality constraints, evaluate them
    if(nEqConstraints >0){
      eqConstraints = zeros(largeArchive.n_rows,nEqConstraints);
      for(int i=0; i < eqConstraints.n_rows; i++){
	eqConstraints.row(i) = trans(this->evaluateEqConstraints(trans(largeArchive.cols(0,nInputs-1).row(i)),eqConstraintPointers,eqConstraintParameterPointers));
      }
    }
    
    //constraint operator before assigning fitness
    largeArchive = this->constraintOperation(C,largeArchive,ineqConstraints,eqConstraints);
    //assign fitness
    largeArchive = this->fitnessAssignment(C,largeArchive.cols(0,nInputs-1),largeArchive.cols(nInputs,nInputs+nObj-1));
    archive = largeArchive.rows(0,10*popSize-1);
    for(int i=0; i < archive.n_rows; i++){
      if(archive(i,nInputs+nObj) != 1)//only take rank 1 as best 
	archive.row(i) = 1e20*ones(1,archive.n_cols);
    }
    
    int archiverows = 0;
    for(archiverows = 0; archiverows < archive.n_rows; archiverows++){
      if(archive(archiverows,0) == 1e20){
	break;
      }
    }
    

    //write archive of best solution to file
    string arcString;
    ofstream arc;
    arcString.append("_pareto");
    
    arc.open(arcString.c_str());
    arc << archive.rows(0,archiverows-1) ;
    arc.close();
    
        
    genEndTime = time(NULL);
    
    endTime = time(NULL);
    cout << "Total run time: " << difftime(endTime, startTime) << endl;
    
}


// Runs the entire GA algorithm for at most the number of generations set by nGen 
// Optimization for 1 objective stops when either:
// - nGens is reached, or
// - the normalized change in the objective value between each generation for the past ndhv_tail = 250 generations 
//   have all been less than EPS = 1e-4 (checking starts after 0.1*nGen generations)
//
// Optimization for 2 objectives stops when either:
// - nGens is reached, or
// - the standard deviations of the crowding distances of all the individuals in the last L generations are less than d_lim, proposed by
// "A steady performance criterion for Pareto-based evolutionary algorithms" Rudenko and Schoenauer
//
// Solutions are written to directory named "data" and 2 files per generation with
// root filename specified in GAControl variable "fileName":
// - _pareto : current best solutions found in the optimization 
// - _pop : current population of the optimization
// Outputs
// - none
// Inputs
// - GAControl *C : pointer to GAControl instance with the parameter values for GA optimization
void GAProblem::solve(GAControl* C)
{
  
    time_t startTime, endTime;
    startTime = time(NULL);

    
    int ndhv_tail = 250; // this is desired length of convergence tail at the end of run. should change to dynamic in the future (i.e. 50% of total number of gens)
    
    int iGen = 1;
    this->setCurGen(iGen);
    
    ifstream in;
    string str;
    string rootStr = "data/"; // output path directory to population and archive solutions
        
    rootStr.append(C->getOutputFile());// output root file name for solutions
        
    // -----gets parameters  
    int nInputs = this->getnInputs();
    int popSize = C->getPopSize();
    int nGens = C->getnGens();
   
    int nObj = this->getnObj();		
    vec outputLowerLims = this->getOutputLowerLims(); 
    vec outputUpperLims = this->getOutputUpperLims();
    vec inputLowerLims = this->getInputLowerLims();
    vec inputUpperLims = this->getInputUpperLims();
    
    int nParents = C->getParentSelectionParams().nParents;
    int nOffspring =C->getRecombinationParams().nOffspring;
    double recombProb = C->getRecombinationParams().recombProb;    
    int distIndex = C->getDistIndex();
    double mutationProb = C->getMutationParams().mutationProb;
    
        
    string fileString;
    fileString.append(rootStr);
    
    fileString.append("_gen");
    
    mat archive = 1e20*ones(10*popSize,nInputs+nObj+2);//preallocated matrix size of archived (best) solutions
    mat eqConstraints;// matrix of equality constraints for the population , size (popSize x nEqConstr)
    mat ineqConstraints;// matrix of inequality constraints for the population, size (popSize x nIneqConstr)
    mat initPop =  this->initializePopulation(*C);
    // runCount is used to "control" the random seed, so that we don't need to store all output file, we only need the seed
    
    mat popObj = zeros<mat>(popSize,nObj);
    
    double (**objFuncPointers)(vec x, void* P) = this->getObjFuncPointers();
    void** objFuncParameterPointers = this->getObjFuncParameterPointers();

    for(int i=0; i < popSize;i++){
      popObj.row(i) = trans(this->evaluateObjectives(trans(initPop.row(i)),objFuncPointers,objFuncParameterPointers));//initial population objectives
    }
    

    int nIneqConstraints = this->getnIneqConstr();
    if(nIneqConstraints >0){
      
      ineqConstraints = zeros(popSize,nIneqConstraints); 
      double (** ineqConstraintPointers)(vec x, void* ineqConstraintParameterPointers) = this->getIneqConstraintPointers();
      void** ineqConstraintParameterPointers = this->getIneqConstraintParameterPointers(); 
      for(int i=0; i < popSize; i++){
       	ineqConstraints.row(i) = trans(this->evaluateIneqConstraints(trans(initPop.row(i)),ineqConstraintPointers,ineqConstraintParameterPointers));
      }
    }
    int nEqConstraints = this->getnEqConstr();
    if(nEqConstraints >0){
      eqConstraints = zeros(popSize,nEqConstraints);
      double (** eqConstraintPointers)(vec x, void* eqConstraintParameterPointers) = this->getEqConstraintPointers();
      void** eqConstraintParameterPointers = this->getEqConstraintParameterPointers();
      for(int i=0; i < popSize; i++){
	eqConstraints.row(i) = trans(this->evaluateEqConstraints(trans(initPop.row(i)),eqConstraintPointers,eqConstraintParameterPointers));
      }
    }

    
    mat joined(initPop.n_rows,initPop.n_cols+popObj.n_cols); 
    
    joined = join_rows(initPop,popObj);// initial population input and objectives
    writeToFile("_functionEvals.dat",joined);
    
    
    //string operations for writing files 

    ///6.11
    joined = this->constraintOperation(C,joined,ineqConstraints,eqConstraints);
    ///6.11
    mat fitness = this->fitnessAssignment(C,joined.cols(0,nInputs-1),joined.cols(nInputs,nInputs+nObj-1));
    stringstream genStream;
    genStream << iGen;
    string genString = fileString;
    
    genString.append(genStream.str());
    string popString = genString;
    popString.append("_pop");
    writeToFile(popString .c_str(),fitness);   
    string arcString = genString;
    arcString.append("_pareto");
    writeToFile(arcString.c_str(),fitness);
    
    archive.rows(0,popSize-1).cols(0,nInputs+nObj) = fitness.cols(0,nInputs+nObj);// update archive
    
    int converged = 0;//convergence flag
    
    int startCheck = floor(0.1*nGens);
    vec checkObj = ones(ndhv_tail);
    
    double prev;
    

      // --------main loop in GA -----------------       
    for(int ii=0; ii < nGens && (!converged); ii++){ //loop over number of generations
      cout << endl << "Start generation " << ii << endl;
      
      time_t genStartTime, genEndTime;
      genStartTime = time(NULL);
            
      mat parents = this->parentSelection(C,fitness); // select parents
      
      mat recomb = this->recombination(C,parents); // recombination for offspring
      
      mat mutate = this->mutation(C,recomb); // mutation operator
      
      
      time_t objStartTime, objEndTime;
      objStartTime = time(NULL);
      mat newObj(mutate.n_rows,nObj); // evaluate objectives of new population
      for(int i=0; i < mutate.n_rows;i++){
	newObj.row(i) = trans(this->evaluateObjectives(trans(mutate.row(i)),objFuncPointers,objFuncParameterPointers));
      }
      objEndTime = time(NULL);
      cout << "Objective evaluation time: " << difftime(objEndTime, objStartTime) << " seconds" << endl;
      mat newPopObj = join_rows(mutate,newObj); // objective values for new population
      ///
      mat newPopObjJoined(mutate.n_rows+parents.n_rows,newPopObj.n_cols);
      newPopObjJoined = join_cols(newPopObj,fitness.cols(0,nInputs+nObj-1));//merge old and new populations
      
      mat fit = zeros(newPopObjJoined.n_rows,nInputs+nObj+2); // fit will contain sorted merged population
      // to save computational time: remove reevaluation 
      
      if(1){//this->getConstraintType().compare("PENALTY-DYNAMIC"){
	mat d = zeros(newPopObjJoined.n_rows,nInputs+nObj);
	d.cols(0,nInputs-1) = newPopObjJoined.cols(0,nInputs-1);
	for(int i=0; i < d.n_rows; i++){
	  d.row(i).cols(nInputs,nInputs+nObj-1) = trans(this->evaluateObjectives(trans(newPopObjJoined.cols(0,nInputs-1).row(i)),objFuncPointers,objFuncParameterPointers));	  
	}

	//if there are inequality constraints, evaluate them      	
	if(nIneqConstraints >0){
	  ineqConstraints = zeros(newPopObjJoined.n_rows,nIneqConstraints); 
	  for(int i=0; i < ineqConstraints.n_rows; i++){
	    ineqConstraints.row(i) = trans(this->evaluateIneqConstraints(trans(newPopObjJoined.cols(0,nInputs-1).row(i)),ineqConstraintPointers,ineqConstraintParameterPointers));
	  }
	}
	//if there are equality constraints, evaluate them      
	if(nEqConstraints >0){
	  eqConstraints = zeros(newPopObjJoined.n_rows,nEqConstraints);
	  for(int i=0; i < eqConstraints.n_rows; i++){
	    eqConstraints.row(i) = trans(this->evaluateEqConstraints(trans(newPopObjJoined.cols(0,nInputs-1).row(i)),eqConstraintPointers,eqConstraintParameterPointers));
	  }
	}

	//contraint handling operator
	d = this->constraintOperation(C,d,ineqConstraints,eqConstraints);

	//fitness assignment
	fit = this->fitnessAssignment(C,d.cols(0,nInputs-1),d.cols(nInputs,nInputs+nObj-1),true,ii);
	writeToFile("ddd.txt",fit);
	
      }

      iGen++;//update generation count
      this->setCurGen(iGen);
      
      stringstream iGenStream;
      iGenStream << iGen;
      genString = fileString;
      genString.append(iGenStream.str());
      string popString = genString;
      popString.append("_pop");
      
      //write current population to file ending in _pop 
      writeToFile(popString.c_str(),fit);     
      
      mat zrs = zeros(newPopObj.n_rows,2);
      mat largeArchive = join_cols(join_rows(newPopObj,zrs), archive); // merged archive of current best solutions
      
	    
      /////re-evaluation slows down optimization but is needed for dynamic penalty
      for(int i=0; i < largeArchive.n_rows && (largeArchive(i,0) < 1e10); i++){
	largeArchive.row(i).cols(nInputs,nInputs+nObj-1) = trans(this->evaluateObjectives(trans(largeArchive.cols(0,nInputs-1).row(i)),objFuncPointers,objFuncParameterPointers));
      }

      //if there are inequality constraints, evaluate them
      if(nIneqConstraints >0){
        ineqConstraints = zeros(largeArchive.n_rows,nIneqConstraints); 
	for(int i=0; i < ineqConstraints.n_rows; i++){
	  ineqConstraints.row(i) = trans(this->evaluateIneqConstraints(trans(largeArchive.cols(0,nInputs-1).row(i)),ineqConstraintPointers,ineqConstraintParameterPointers));
	}
      }
      //if there are equality constraints, evaluate them
      if(nEqConstraints >0){
	eqConstraints = zeros(largeArchive.n_rows,nEqConstraints);
	for(int i=0; i < eqConstraints.n_rows; i++){
	  eqConstraints.row(i) = trans(this->evaluateEqConstraints(trans(largeArchive.cols(0,nInputs-1).row(i)),eqConstraintPointers,eqConstraintParameterPointers));
	}
      }

      //constraint handling operator      
      largeArchive = this->constraintOperation(C,largeArchive,ineqConstraints,eqConstraints);
      
      //fitness assignment
      largeArchive = this->fitnessAssignment(C,largeArchive.cols(0,nInputs-1),largeArchive.cols(nInputs,nInputs+nObj-1));
      archive = largeArchive.rows(0,10*popSize-1);
      
      for(int i=0; i < archive.n_rows; i++){
	if(archive(i,nInputs+nObj) != 1)//only keep the rank 1 as best solutions so far
	  archive.row(i) = 1e20*ones(1,archive.n_cols);
      }
      
      int archiverows = 0;//remove empty rows in archive matrix
      for(archiverows = 0; archiverows < archive.n_rows; archiverows++){
	if(archive(archiverows,0) == 1e20){
	  break;
	}
      }
      ofstream arc;
      string arcString = genString;
      arcString.append("_pareto");
      //write the best solutions so far to file ending in _pareto      
      arc.open(arcString.c_str());
      arc << archive.rows(0,archiverows-1) ;
      arc.close();

      /////
      string NDHVstring = fileString;
      NDHVstring.append("_NDHV");
      ofstream nd;
      
      //convergence check for 1 objective
      if(nObj ==1){
	if(ii == startCheck)
	  prev = archive(0,nInputs+nObj-1);
	
	if(ii > startCheck){
	  if(ii>= startCheck + ndhv_tail){
	    if( checkObj.max() < EPS){
	      cout << "checkObj" << checkObj;
	      converged = 1;
	    }
	    else {
	      checkObj(((ii-startCheck)%ndhv_tail)) = abs(archive(0,nInputs+nObj-1)-prev)/abs(prev);
	      cout << (ii-startCheck)%ndhv_tail<<":"<<archive(0,nInputs+nObj-1) << ":" <<prev << endl;
	    }
	  }
	  else {
	    checkObj(((ii-startCheck)%ndhv_tail)) = abs(archive(0,nInputs+nObj-1)-prev)/abs(prev);
	    cout << (ii-startCheck)%ndhv_tail<<":"<<archive(0,nInputs+nObj-1) << ":" <<prev << endl;
	    
	  }
	  prev = archive(0,nInputs+nObj-1);
	}
      }
      
      if (nObj == 1){
	double monitor = archive(0,nInputs+nObj-1);
	nd.open("_1objvalues.dat",ios::app);
	nd << monitor << endl;
	nd.close();
      }
      else if(nObj> 1 ){
	//all other cases nObj > 1					
	double ndhv = this->NDHV(fit.cols(nInputs,nInputs+1),outputLowerLims,outputUpperLims);
	
	nd.open(NDHVstring.c_str(),ios::app);
	nd << ndhv << endl;
	nd.close();
	
	
	// old convergence criterion, to be removed, only use G->hasConverged()
	if(ii == startCheck)
	  prev = ndhv;
	if(ii > startCheck){
	  if(ii>= startCheck + ndhv_tail){
	    if( checkObj.max() < EPS)
	      converged = 1;
	    else {

	      checkObj(((ii-startCheck)%ndhv_tail)) = abs(ndhv-prev)/abs(prev);
	    }
	  }
	  else {
	    checkObj(((ii-startCheck)%ndhv_tail)) = abs(ndhv-prev)/abs(prev);
	  }
	  prev = ndhv;
	}
	
	// Convergence == true iff ndhv has stabilized and crowding dist has stabilized
	converged = (converged && this->hasConverged());
      }
      
      fitness = fit.rows(0,popSize-1);
      
      mat pareto(fitness.n_rows,fitness.n_cols);
      pareto  = fitness;
      
      //writeToFile("_results.dat", pareto);
      
      genEndTime = time(NULL);
      
      cout << "End generation " << ii << ", total " << difftime(genEndTime, genStartTime) << " seconds" << endl;
      
    }
    
    mat pareto(fitness.n_rows,fitness.n_cols);
    pareto  = fitness;
    //    writeToFile("_results.dat", pareto);
    
        
    endTime = time(NULL);
    cout << "Total run time: " << difftime(endTime, startTime) << endl;
    
}
