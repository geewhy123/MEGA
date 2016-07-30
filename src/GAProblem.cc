/****************************************
 *file:GAProblem.cc
 * contains most functions of the GAProblem class 
 * including specific operators of the GA
 * 
 *
 * 
 * Gary Yan, Peter Zhang
 * June 13,2012
 ****************************************/


#include <iostream>
#include <armadillo>
#include <stdlib.h>
#include "../include/GAProblem.hh"
#include <float.h>
#include <queue>
using namespace std;
using namespace arma;



/*--- Set and Get methods ----
 * Member variables are private, user can only access 
 * them through the public set/get methods
 */

void GAProblem::setnInputs(int nInputs){
	this->nInputs = nInputs;
}
int GAProblem::getnInputs(){
	return this->nInputs;
}
void GAProblem::setInputLowerLims(vec inputLowerLims){
	this->inputLowerLims = inputLowerLims;
}
vec GAProblem::getInputLowerLims(){
	return this->inputLowerLims;
}
void GAProblem::setInputUpperLims(vec inputUpperLims){
	this->inputUpperLims = inputUpperLims;
}
vec GAProblem::getInputUpperLims(){
	return this->inputUpperLims;
}

void GAProblem::setOutputLowerLims(vec outputLowerLims){
	this->outputLowerLims = outputLowerLims;
}
vec GAProblem::getOutputLowerLims(){
	return this->outputLowerLims;
}
void GAProblem::setOutputUpperLims(vec outputUpperLims){
	this->outputUpperLims =outputUpperLims;
}
vec GAProblem::getOutputUpperLims(){
	return this->outputUpperLims;
}
void GAProblem::setnObj(int nObj){
	this->nObj = nObj;
}
int GAProblem::getnObj(){
	return this->nObj;
}

double GAProblem::getDLim(){
  return this->d_lim;
}
void GAProblem::setDLim(double DLim){
  this->d_lim = DLim;
}
int GAProblem::getL(){
  return this->L;
}
void GAProblem::setL(int L){
  this->L = L;
}
vec GAProblem::getMaxCrowdingDist(){
  return this->maxCrowdingDist;
}
void GAProblem::setMaxCrowdingDist(vec maxCrowdingDist){
  this->maxCrowdingDist = maxCrowdingDist;
}

void GAProblem::setObjFuncPointers(double(**objFuncPointers)(vec x,void* P)){
  this->objFuncPointers = objFuncPointers;
}
double (**GAProblem::getObjFuncPointers())(vec x, void* P){                
  return this->objFuncPointers;
}                      
void GAProblem::setIneqConstraintPointers(double(**ineqConstraintPointers)(vec x,void* Pg)){
  this->ineqConstraintPointers = ineqConstraintPointers;
}
double (**GAProblem::getIneqConstraintPointers())(vec x, void* Pg){
  return this->ineqConstraintPointers;
}
void GAProblem::setEqConstraintPointers(double(**eqConstraintPointers)(vec x,void* Ph)){
  this->eqConstraintPointers = eqConstraintPointers;
}
double (**GAProblem::getEqConstraintPointers())(vec x, void* Ph){
  return this->eqConstraintPointers;
}
void GAProblem::setObjFuncParameterPointers(void** objFuncParameterPointers){
  this->objFuncParameterPointers = objFuncParameterPointers;
}
void ** GAProblem::getObjFuncParameterPointers(){
  return this->objFuncParameterPointers;
}
void GAProblem::setIneqConstraintParameterPointers(void** ineqConstraintParameterPointers){
  this->ineqConstraintParameterPointers = ineqConstraintParameterPointers;
}
void ** GAProblem::getIneqConstraintParameterPointers(){
  return this->ineqConstraintParameterPointers;
}
void GAProblem::setEqConstraintParameterPointers(void** eqConstraintParameterPointers){
  this->eqConstraintParameterPointers = eqConstraintParameterPointers;
}

void ** GAProblem::getEqConstraintParameterPointers(){
  return this->eqConstraintParameterPointers;
}

void GAProblem::setCurrentPop(mat pop){
  this->currentPop = pop;
}

mat GAProblem::getCurrentPop(){
  return this->currentPop;
}


void GAProblem::setConstraintType(string s){
  this->constraintType = s;
}
string GAProblem::getConstraintType(){
  return this->constraintType;
}


void GAProblem::setnIneqConstr(int n){
  this->nIneqConstr = n;
}
int GAProblem::getnIneqConstr(){
  return this->nIneqConstr;
}
void GAProblem::setnEqConstr(int n){
  this->nEqConstr = n;
}
int GAProblem::getnEqConstr(){
  return this->nEqConstr;
}
void GAProblem::setCurGen(int iGen){
  this->currentGen = iGen;
}
int GAProblem::getCurGen(){
  return this->currentGen;
}

//penalty coefficients is an array of (nIneqConstr + nEqConstr ) x nObj
// with each element representing how much the inequality or equality constraint
// penalized the objective
// rows relating to inequality constraints are placed first, followed by equality 
// constraint rows.
void GAProblem::setPenaltyCoefficients(mat P){
  this->penaltyCoefficients = P;
}
mat GAProblem::getPenaltyCoefficients(){
  return this->penaltyCoefficients;
}

// Convergence criterion:
// The standard deviation of d_l over past L runs is less than d_lim
// d_l is the max crowding distance at run l; L and d_lim are constants. 
// For example, L = 100, d_lim = 0.02 (normalized)
bool GAProblem::hasConverged()
{
  cout << "Convergence: delta_l" << stddev(maxCrowdingDist) << endl;
  return (mean(maxCrowdingDist) > 0.0) &&
    (mean(maxCrowdingDist) < 1.0) &&
    (stddev(maxCrowdingDist) /*/ mean(maxCrowdingDist)*/ < d_lim);
}



/* InitializePopulation function
 * Input: Instance of GAControl class
 * Output: matrix of size Popsize x nInputs of random input values
 * between lower and upper limit 
 * Random seed is specified as a member variable of GAProblem instance
 * use time(NULL) for random
 */
mat GAProblem::initializePopulation(GAControl C){
  // When checking whether different NHDV algorithms give the same values,
  // fix the seed so that the results are the same.
  //C.setInitialSeed(time(0));// clock time as the random seed
  //C.setInitialSeed(100);
	
  //	d_lim = d_lim_input;
  //	L = L_input;
  int L = this->getL();
  double d_lim = this->getDLim();
  int seed = C.getInitialSeed();
  this->maxCrowdingDist =ones<vec>(L);
  //  cout << maxCrowdingDist<<":::::"<< L;
	for (int i = 0; i < L; i ++)
		this->maxCrowdingDist(i) = (i % 2 == 0) ? 0 : 1; // creating an outrageous crowding distance standard deviation
	/*for (int i = 0; i < L; i ++)
	  cout << maxCrowdingDist(i) << " ";*/
	//cout << endl;
	
	srand(seed);
	int popSize = C.getPopSize();
	int nInputs = this->getnInputs();
	
	mat initPopInput=zeros(popSize,nInputs);
	
	vec inputi = zeros(popSize); 
	for(int i=0; i < nInputs; i++){
		double low = this->getInputLowerLims()[i];
		double upp = this->getInputUpperLims()[i];
		inputi = low+(upp - low) * randu<vec>(popSize); 
		initPopInput.col(i) = inputi;// randomize each input column		
	}
	
	return initPopInput;
}

struct objP;

/* evaluateObjectives function
 * Inputs: x - vector of input values
 *         *f - function pointer taking input values 
 *         objParams - pointer to function parameters
 * Output: objective value evaluated for the scalar function f
 */
vec GAProblem::evaluateObjectives(vec x, double(**f)(vec x,void* p),void** objParams){

  int nObj = this->getnObj();
  int iGen = this->getCurGen();	
  vec result = zeros(nObj);
    
  for(int i=0; i < nObj; i++){
    /*    if(this->getConstraintType().compare("PENALTY")==0){
      //      cout << "calling penalized objectives..." << endl;
      //cout <<"calling g ineq constr..." << endl;
      vec g = evaluateIneqConstraints(x,this->getIneqConstraintPointers(),this->getIneqConstraintParameterPointers());
      //cout << "calling h eq constr..." << endl;
      vec h = evaluateEqConstraints(x,this->getEqConstraintPointers(),this->getEqConstraintParameterPointers());
      //cout <<" calling Q...." <<endl;
      result(i) = Q[i](f[i](x,objParams[i]),g,h,iGen);

    }
    else*/
      result(i) =   f[i](x,objParams[i]);
      
  }
  return result;
}

/* evaluateIneqConstraints function                                                                                                                   * Inputs: x - vector of input values                                                                                                                 *         *g - function pointer taking input values                                                                                                  *         constrParams - pointer to function parameters               
 * Outputs: vector of constraint values evaluated for the scalar function g  (each individual)
 */
vec GAProblem::evaluateIneqConstraints(vec x, double(**g)(vec x, void* p),void** constrParams){
  int nConstr = this->getnIneqConstr();
  vec result = zeros(nConstr);
  for(int i=0; i < nConstr; i++)
    result(i) = g[i](x,constrParams[i]);
	return result;
}

/* evaluateEqConstraints function                                                                                                                   * Inputs: x - vector of input values                                                                                                                 *         *g - function pointer taking input values                                                                                                  *         constrParams - pointer to function parameters               
 * Outputs: vector of constraint values evaluated for the scalar function h  (each individual)
 */
vec GAProblem::evaluateEqConstraints(vec x, double(**h)(vec x, void* p),void** constrParams){
  int nConstr = this->getnEqConstr();
  vec result = zeros(nConstr);
  for(int i=0; i < nConstr; i++)
    result(i) = h[i](x,constrParams[i]);
  return result;
}

// Constraint operation
// called before fitness assignment
// Only implemented now for penalty function
// output:
// - matrix of population and penalized objective functions of size (popsize x (nInputs + nObj))
// inputs:
// - GAControl *C : pointer to GAControl instance with parameters
// - mat ineqConstraints : matrix of inequality constraint evaluations for each individual ( size (popSize x nIneqConstr))
// - mat eqConstraints: matrix of equality constraint evaluations for each individual ( size (popSize x nEqConstr))
mat GAProblem::constraintOperation(GAControl* C,mat popAndObj,mat ineqConstraints,mat eqConstraints){
  if(this->getConstraintType().substr(0,7).compare("PENALTY") ==0){
    return this->penalizeObjectives(C, popAndObj,ineqConstraints,eqConstraints); 
  }
}

//implements the penalty method for constraint handling (static or dynamic)
//the form of the static penalty value p is :
//p(x) = f(x) + c * h(x) ^2 * (i/N)^2, 
// is called when constraintType = "PENALTY-DYNAMIC"
//
//and static penalty is:
//p(x) = f(x) + c * max(g(x),0) ^2 
//is called otherwise when constraintType starts with "PENALTY"
//where i is the current generation, and N is the total number of generations 
// f(x) is the objective value, g(x) is the inequality constraint, h(x) is the equality constraint
//
// the quadratic dependence of the generation number for the dynamic penalty is take from :
// "Varying fitness functions in genetic algorithms: studying the rate of increase of the dynamic penalty terms"
// by Kazarlis and Petridis, 1998, cited in 
// "Theoretical and numerical constraint-handling techniques used with evolutionary algorithms: a survey of the state of the art"
// by Coello , 2001: 2.2 Dynamic Penalties: equation (17)
// 
// Output:
// -  matrix of population and penalized objective functions of size (popsize x (nInputs + nObj))
// Inputs:
// - GAControl *C : pointer to GAControl instance with parameters
// - mat ineqConstraints : matrix of inequality constraint evaluations for each individual ( size (popSize x nIneqConstr))
// - mat eqConstraints: matrix of equality constraint evaluations for each individual ( size (popSize x nEqConstr))

mat GAProblem::penalizeObjectives(GAControl* C,mat popAndObj,mat ineqConstraints,mat eqConstraints){
  
  mat objectives = zeros(popAndObj.n_rows,nObj);
  objectives = popAndObj.cols(nInputs,nInputs+nObj-1);
  
  mat penaltyCoefficients = this->getPenaltyCoefficients();
  double penalty = 0.0;
  for(int i=0; i < objectives.n_rows; i++){
    for(int j=0 ; j < nObj; j++){
      for(int k1=0; k1 < nIneqConstr ; k1++){

	// add ineq constraint contributions
	if((this->getConstraintType().size() > 7)&& (this->getConstraintType().substr(8,14).compare("DYNAMIC") ==0)){ // dynamic penalty
	  penalty = pow(max(ineqConstraints(i,k1),0.0),2)*penaltyCoefficients(k1,j)*pow(((double)this->getCurGen())/C->getnGens(),2);
	}
	else{//static penalty
	  penalty = pow(max(ineqConstraints(i,k1),0.0),2)*penaltyCoefficients(k1,j);
	}
	objectives(i,j) = objectives(i,j) + penalty;
      }
      for(int k2=0; k2< nEqConstr; k2++){
	// add eq constraint contributions
	if(this->getConstraintType().substr(8,14).compare("DYNAMIC") ==0)//dynamic penalty
	  penalty = pow(eqConstraints(i,k2),2)*penaltyCoefficients(k2+nIneqConstr,j)*pow(((double)this->getCurGen())/C->getnGens(),2);
	else{//static penalty
	  penalty = pow(eqConstraints(i,k2),2)*penaltyCoefficients(k2+nIneqConstr,j);
	}
	objectives(i,j) = objectives(i,j) + penalty;
	
      }
    }
  }

  return join_rows(popAndObj.cols(0,nInputs-1),objectives);
    
}


//helper function
mat GAProblem::sortByVec(vec sortBy,mat input, int asc){
	mat sorted(input.n_rows,input.n_cols);
	uvec index = sort_index(sortBy,asc);
	
	for(int i=0; i < input.n_rows; i++){
		sorted.row(i) = input.row(index(i));
	}
	
	return sorted;
}

/* fitnessAssignment function
 * 
 * - calls different variants of GA according to input variable fitnessAssignmentType of the GAControl class
 * - currently only implements NSGA-II
 * - need other types of GA (NSGAII, SPEA, etc) later
 *
 * Inputs: 
 * - GAControl * C : pointer to GAControl instance with parameters
 * - popInput : matrix of size popSize x nInputs of input values 
 * - objVals : matrix of size popSize x nObj of corresponding objective values
 * - isPureNewPop : boolean value that evaluates to true if the population is new
 * - genCount : current generation number
 * Outputs:
 * - matrix of size popSize x (nInputs + nObj + 2) , last 2 columns are rank and crowding distance values
 * 
 */
mat GAProblem::fitnessAssignment(GAControl* C, mat popInput,mat objVals, bool isPureNewPop, int genCount){
  string s = C->getFitnessAssignmentType();

  if(s.compare("NSGA-II") ==0){
return fitnessAssignmentNSGAII(popInput,objVals,isPureNewPop,genCount);
    
  }//else if...other variants SPEA, etc .. (to be implemented)
  else{
 return fitnessAssignmentNSGAII(popInput,objVals,isPureNewPop,genCount);
  }
}


/* fitnessAssignmentNSGAII function
 * - implementation of the NSGAII variant of GA
 * - follows algorithm by kangal
 * - follows implementation in MATLAB code 
 * Inputs: 
 * - GAControl * C : pointer to GAControl instance with parameters
 * - popInput : matrix of size popSize x nInputs of input values 
 * - objVals : matrix of size popSize x nObj of corresponding objective values
 * - isPureNewPop : boolean value that evaluates to true if the population is new
 * - genCount : current generation number
 * Outputs:
 * - matrix of size popSize x (nInputs + nObj + 2) , last 2 columns are rank and crowding distance values
 * 
 */
mat GAProblem::fitnessAssignmentNSGAII(mat popInput, mat objVals, bool isPureNewPop, int genCount){//, const string type){
  
  mat popAndObj = join_rows(popInput,objVals); //combine input and objective values into one matrix of size popSize x (nInputs + nObj)
  int nRows = popInput.n_rows;
  int nObj = objVals.n_cols;
  int nInputs = popInput.n_cols;

  //for one objective,only objective value determines rank.  Crowding distance irrelevant
  if(nObj == 1){
    cout << "1obj fitness" << endl;
    mat popObjRankCrowd = zeros<mat>(nRows,nInputs+nObj+2);
    popObjRankCrowd.cols(0,nInputs+nObj-1) = popAndObj;
    popObjRankCrowd = this->sortByVec(popObjRankCrowd.col(nInputs+nObj-1),popObjRankCrowd,0);
    for(int i=0; i < popObjRankCrowd.n_rows ; i++){
      popObjRankCrowd(i,nInputs+nObj) = i+1;
    }

    popObjRankCrowd.col(nInputs+nObj+1) = DBL_MAX*ones<vec>(nRows,1);
    
    return popObjRankCrowd;
    
  }
  
  //for nObj > 1

  vec individualN = zeros<vec>(nRows,1);
  mat individualP = zeros<mat>(nRows,nRows);
  mat front = zeros<mat>(nRows,nRows);  
  int domLess,domEqual,domMore,count=0;
  int fcount = 0, fx=0,p,q;
  mat rank = zeros<mat>(nRows,1);
  int nFront=1;
  popAndObj = join_rows(popAndObj,rank);

  //assign rank
  for (int i=0; i < nRows; i++){ //loop over each individual p
    individualN(i) = 0; 
    count = 0;
    for (int j=0; j < nRows; j++){
      domLess = 0;
      domEqual = 0;
      domMore = 0;
      for(int k=0; k < nObj; k++){
	
	if(popAndObj(i, nInputs +k) < popAndObj(j, nInputs + k))
	  domLess++;
	else if(popAndObj(i,nInputs +k)== popAndObj(j, nInputs + k))
	  domEqual++;
	else 
	  domMore++;
      }
      if (domLess == 0 && domEqual != nObj)
	individualN(i)++; //count number of other individuals that dominate p
      else if(domMore == 0 && domEqual != nObj){
	
	individualP(count,i) = j+1; //keeps a list of other individuals dominated by p
	count++;
      }
    }
    if(individualN(i) == 0){
      popAndObj(i,nObj + nInputs) = 1; 
      front(fx++,fcount) = i+1;
      //assign fronts
    }
  }
  
  
  for(int l = 1; l < nRows; l++){
    p=0;
    for(int m=0; m < nRows; m++){
      if( front(m,l-1) == 0) 
	break;
      for(int n = 0; n < nRows; n++){
	if(individualP(n,front(m,l-1)-1) ==0) 
	  break;
	individualN(individualP(n,front(m,l-1)-1)-1)--;
	
	if(individualN(individualP(n,front(m,l-1)-1)-1) == 0){
	  
	  front(p++,l) = individualP(n,front(m,l-1)-1);
	  nFront = l+1;
	}
      }
      
    }  
  }
  
  for(int ii = 1; ii < nRows ; ii++){
    if(front(0,ii) == 0) 
      break;
    for(int jj=0; jj < nRows ; jj++){
      if(front(jj,ii) == 0) 
	break;
      popAndObj(front(jj,ii)-1,nObj+nInputs)=ii+1;
    }
  }
  
  /*----------------
    /evaluate crowding distances
    /---------------------*/
  
  /*sort based on rank*/
  
  
  mat newSort = zeros(popAndObj.n_rows,popAndObj.n_cols);
  
  newSort = this->sortByVec(popAndObj.col(nInputs+nObj),popAndObj,0);
  
  
  mat result(0, newSort.n_cols+1);
  int indCount = 0;
  int end = 0;
  
  for (int i=0; i < nFront ; i++){//loop over number of fronts  
    uvec q1 = find(newSort.col(nInputs+nObj) ==i+1 );
    end = q1.n_rows;
    vec temp=zeros(q1.n_rows);
    mat newSub = newSort.rows(q1(0),q1(q1.n_rows-1));
    newSub.insert_cols(newSub.n_cols,1);
    
    for(int k=0; k < nObj ; k++){
      if(temp.n_rows == 1){
	
	newSub(0,nInputs+nObj+1) = DBL_MAX;
	indCount = q1.n_rows;
	break;
      }
      if(temp.n_rows == 2){
	// for 2 individuals in the same rank, set both crowding dist to inf
	newSub(0,nInputs+nObj+1) = DBL_MAX;
	newSub(1,nInputs+nObj+1) = DBL_MAX;
	
	indCount = q1.n_rows;
	break;
      }
      
      newSub = this->sortByVec(newSub.col(nInputs+k),newSub,0);
      
      newSub(0,nInputs+nObj+1) = DBL_MAX; // set individuals with extreme objective values to have crowding dist of inf
      newSub(newSub.n_rows-1,nInputs+nObj+1) = DBL_MAX;// set individuals with extreme objective values to have crowding dist of inf
      
      
      for(int j=1; j < newSub.n_rows-1; j++){
	//calculates L1 distance between adjacent individuals for each objective and calculates the sum	
	temp(j) = (newSub(j+1,nInputs+k)-newSub(j-1,nInputs+k))/(this->getOutputUpperLims()[k]-this->getOutputLowerLims()[k]);
	
	newSub(j,nInputs+nObj+1) += temp(j);
      }
      
    }
    
    //see convergence criterion in "A steady performance criterion for Pareto-based evolutionary algorithms" Rudenko and Schoenauer
    if (isPureNewPop && i == 0) // For convergence criterion, only check the crowding distance of rank 1 front
      {
	double maxCD = 0.0;
	for (int j = 1; j < newSub.n_rows - 1; j ++) // disregarding the head and tail
	  {
	    maxCD = (newSub(j, nInputs + nObj + 1) != DBL_MAX) && (newSub(j, nInputs + nObj + 1) > maxCD) ? 
	      newSub(j, nInputs + nObj + 1) : maxCD;
	  }
	
	//cout << "maxCD " << maxCD << endl;
	
	this->maxCrowdingDist(genCount % L) = maxCD;
	
	// testing
	//cout << "fitnessAssignment gen " << genCount << " rank " << i << endl;
	for (int j = 0; j < maxCrowdingDist.n_rows; j ++)
	  ;//cout << maxCrowdingDist(j) << " ";
	//cout << endl;
      }
    
    newSub = this->sortByVec(newSub.col(nInputs+nObj+1),newSub,1);
    result = join_cols(result,newSub);    
    indCount += temp.n_rows;
  }

  return result;
  
}

// parentSelection function
// - calls different variants of parentSelection operator according to input variable parentSelectionType of the GAControl class
// - currently only implements binary tournament (BINTOUR)
// - can implement other types (constrained BINTOUR,roulette, etc)later
// 
// Outputs
// - matrix of selected parents for the mating pool
// Inputs
// - GAControl * C - pointer to GAControl instance with parameters
// - mat popAndFitness - population, objective and fitness values of the current population, size (popSize x (nInputs + nObj + 2)
mat GAProblem::parentSelection(GAControl* C,mat popAndFitness){
  
    string s = C->getParentSelectionType();

    if(s.compare("BIN-TOUR")==0){
      return parentSelectionBINTOUR(C,popAndFitness);


    }//else if...                                                                                                                                    
    else{
      return parentSelectionBINTOUR(C,popAndFitness);
    }
}

//parentSelectionBINTOUR function
// - implements the binary tournament selection for selecting parents without replacement
// Outputs:
// - matrix of size nParents x (nInputs+ nObj) containing parents selected for recombination
// Inputs:
// - GAControl * C : pointer to GAControl instance containing parameters
// - mat popAndFitness: matrix of population, objective, and fitness values sorted based on rank, and 
//   in each rank based on crowding distance in decreasing order.  Individuals higher up in the matrix 
//   (with lower row index) are fitter than others. size (popSize x (nInputs + nObj +2))
mat GAProblem::parentSelectionBINTOUR(GAControl*C, mat popAndFitness){
	int nInputs = this->getnInputs();
	int nObj = this->getnObj();
	parentSelectionParams params = C->getParentSelectionParams();
	mat parents = zeros(params.nParents,popAndFitness.n_cols);
	
	//number of parents to be selected
	for(int i=0; i < params.nParents; i++){
		int size = popAndFitness.n_rows;
		//pick two rows randomly
		int index1=rand() % size ;
		int index2 = rand() % size ;
		
		//if repeated, select again
		while(index1 == index2)
			index2 = rand()% size ;
		
		mat first = popAndFitness.row(index1);
		mat second = popAndFitness.row(index2);
		//rewrite compare function here for constrained bin-tour		
		if(index1 < index2){
			parents.row(i) = first;
			popAndFitness.shed_row(index1);
		}
		else{
			parents.row(i) = second;
			popAndFitness.shed_row(index2);
		}
		
	}
	
	return parents.cols(0,nInputs+nObj-1);
}
//true if first better than second, false otherwise
//bool GAProblem::binTourCompare(){
//  if(this->getConstraintType().compare("") ==0);


//recombination function
// - calls different variants of recombination operator according to input variable recombinationType of the GAControl class
// - currently only implements simulated binary crossover (SBX) , referenced "Simulated Binary Crossover for Continuous Search Space" by Deb
// and Agrawal, 1995
// - can implement other types later
//
// Outputs:
// - matrix of offspring generated by parents
// Inputs:
// - GAControl* C : pointer to GAControl instance containing parameters
// - mat parents: pool of selected parents from parentSelection operator
mat GAProblem::recombination(GAControl* C, mat parents){
  string s = C->getRecombinationType();
  //if (strcmp(s, "SBX")==0){
  if(s.compare("SBX") ==0){
    return recombinationSBX(C,parents);

  }//else if...(other recombination methods later)                                                                                                                                       
  else{

    return recombinationSBX(C,parents);
  }
}


// recombinationSBX function
// implements the simulated binary crossover (SBX) recombination method
// Outputs:
// - matrix of offspring generated by parents
// Inputs:
// - GAControl* C : pointer to GAControl instance containing parameters
// - mat parents: pool of selected parents from parentSelection operator
// check: what happens if recombination probability test fails?
// - pick another 2
// - take parents as offspring
// - do not pick again
// currently approach 1 is taken, but should be validated , since it should be the same as setting recombination probability to 1?
mat GAProblem::recombinationSBX(GAControl* C, mat parents){
  recombinationParams params = C->getRecombinationParams();
	int nParents = parents.n_rows;
	
	int nOffspring = params.nOffspring;
	int nInputs = this->getnInputs();
	
	mat offspring = zeros(params.nOffspring+nOffspring%2,nInputs);
	
	int distIndex = C->getDistIndex();
	
	int nNewInd = 0; 
	vec lo = this->getInputLowerLims();
	vec hi = this->getInputUpperLims();
	mat u = zeros(1,nInputs);
	mat bq = zeros(1,nInputs);
	mat child1 = zeros(1,nInputs);
	mat child2 = zeros(1,nInputs);
	int curindex = 0;
	//pick 2 parents, generate 2 offspring
	for(int i=0; i < ceil(params.nOffspring/2.0); i++){
		if((double)rand()/RAND_MAX < params.recombProb){    
			
			int index1 = rand()%nParents;
			
			int index2 = rand()%nParents;
			while(index1 == index2)
				index2 = rand()%nParents;
			
			mat parent1 = parents.row(index1);
			mat parent2 = parents.row(index2);
			
			// worry about repeated selection?
			
			for(int j=0; j < nInputs; j++){
				u(j) = (double)rand()/RAND_MAX;
				
				if(u(j) <= 0.5)
					bq(j) = pow((2*u(j)),(1.0/(distIndex+1)));
				else bq(j) = pow(1.0/(2*(1-u(j))),(1.0/(distIndex+1)));
				
				child1(j) = 0.5*(((1+bq(j))*parent1(j)) + (1-bq(j))*parent2(j));
				child2(j) = 0.5*(((1-bq(j))*parent1(j)) + (1+bq(j))*parent2(j));
				
				if(child1(j) > hi(j))
					child1(j) = hi(j);
				else if(child1(j) < lo(j))
					child1(j) = lo(j);
				
				if(child2(j) > hi(j))
					child2(j) = hi(j);
				else if(child2(j) < lo(j))
					child2(j) = lo(j);
				
			}
			nNewInd = nNewInd+2;    
			
			offspring.row(curindex+i) = child1;
			offspring.row(curindex+i+1) = child2;
			curindex = curindex+1;
		}
		else i--;// pick another 2

	}
	
	//if nOffspring is odd, delete the last one to match the number
	if(ceil(nOffspring/2.0) != nOffspring/2.0){
		offspring.shed_row(offspring.n_rows-1);
		nNewInd = nNewInd-1;
	}

	return offspring;
}


// mutation function
// - calls different variants of mutation operator according to input variable mutationType of the GAControl class
// - currently only implements Polynomial mutation
// - can implement other types later
//
// Outputs:
// - matrix of mutated offspring 
// Inputs:
// - GAControl* C : pointer to GAControl instance containing parameters
// - mat offspring: matrix of previously generated offspring to be mutated
mat GAProblem::mutation(GAControl *C, mat offspring){
  string s = C->getMutationType();
  if(s.compare("POLYNOMIAL")==0){
    return mutationPOLY(C,offspring);

  }//else if...(other mutation methods)                                         
  else{
    return mutationPOLY(C,offspring);
  }
}


// mutationPOLY function
// performs polynomial mutation on each offspring
// Output :
// - matrix of mutated offspring: size equal to input offspring
// Inputs:
// - GAControl* C - pointer to GAControl instance containing parameters
// - mat offpspring - matrix of input offspring individuals , size (nOffspring x nInputs)
mat GAProblem::mutationPOLY(GAControl* C, mat offspring){
  mutationParams params = C->getMutationParams();	
	
  int popSize = offspring.n_rows;
  
  int nInputs = offspring.n_cols; 
  double prob = params.mutationProb;
  int distIndex = C->getDistIndex();
  vec lo = this->getInputLowerLims();
  vec hi = this->getInputUpperLims();
  int nMutations = 0;
    
  vec r = zeros(nInputs);
  vec delta = zeros(nInputs);
  mat newIndividual = zeros(1,nInputs);
  mat newOffspring = offspring;
  
  // loop over the population of solutions
  for (int i = 0; i<popSize; i++){
    // perform mutation with probability = params.prob
    if ((double)rand()/RAND_MAX < prob){
      // mutation is performed independently on each input variable
      for (int j = 0 ; j< nInputs; j++){
	r(j) = (double)rand()/RAND_MAX;
	
	if (r(j) < 0.5)
	  delta(j) = pow((2*r(j)),(1.0/(distIndex+1))) - 1;
	else
	  delta(j) = 1 - pow((2*(1 - r(j))),(1.0/(distIndex+1)));
	
	
	// Generate the corresponding child element
	newIndividual(0,j) = offspring(i,j) + delta(j);
	
	// Make sure that the generated element is within the decision space.
				
	if (newIndividual(0,j) > hi(j))
	  newIndividual(0,j) = hi(j);
	else if (newIndividual(0,j) < lo(j))
	  newIndividual(0,j) = lo(j);
	
      }
      nMutations++;
      
      newOffspring.row(i) = newIndividual.row(0);
      
    }
  }
  return newOffspring;
}




// ------------helper functions for old NDHV (currently unused)-------------

// helper function, generates full factorial needed for NDHV
mat fullfact(vec input ){
	int prod = 1;
	int count=1;
	int bcount = 1;
	for(int i=0; i < input.n_rows; i++)
		prod *= input(i);
	mat result = zeros(prod,input.n_rows);
	
	for(int j=0; j < result.n_cols;j++){
		count = pow(2.0,j);
		for(int k=0; k < result.n_rows; k++){
			if(k ==0)
				result(k,j) = 1;
			else if((k/count)%2==0){
				result(k,j) = 1;
			}
			else {
				result(k,j) = 2;
			}
		}
	}
		
	return result;
}

// helper function - sorts a matrix by values of a column
mat matSortCol(mat input){
	int n = input.n_cols;
	//vec index;
	mat sorted = zeros(input.n_rows,input.n_cols);
	for(int i=0; i < n; i++){
		uvec index = sort_index(input.col(i));
		for(int j=0; j < input.n_rows ; j++){
			sorted(j,i) = input(index(j),i);
		}
	}
	sorted.print();
	return sorted;
}

// helper function - imitates matlab functions - convert matrix to a cell
field<mat> mat2cell(mat vecSortedPts,int a, vec b){
	field<mat> result(1,b.n_rows);
	result(0) = vecSortedPts.cols(0,b(0)-1);
	int sum = b(0);
	for(int i=1; i < b.n_rows; i++){
		result(i) = vecSortedPts.cols(sum,sum+b(i)-1);
		sum += b(i);  
	}
	return result;
}

// helper function - 2D grid generation
field<mat> ndgrid2(field<mat> x){
	field<mat> result(1,x.n_cols);  
	int n = x.n_cols;
	
	
	if (n ==2){
		mat a(x(0).n_cols,x(1).n_cols);
		mat b(x(0).n_cols,x(1).n_cols);
		
		for(int i=0; i < a.n_cols; i++){
			a.col(i) = trans(x(0));
		}
		for(int j=0; j < b.n_rows; j++){
			b.row(j) = x(1);
		}  
		result(0) = a;
		result(1) = b;
	}
	return result; 
}
// helper function - 3D grid generation
field<cube> ndgrid3(field<mat> x)
{
	field<cube> result(1,x.n_cols);
	int n = x.n_cols;
	if(n==3){
	  
	  cube a(x(0).n_cols,x(1).n_cols,x(2).n_cols);
	  cube b(x(0).n_cols,x(1).n_cols,x(2).n_cols);
	  cube c(x(0).n_cols,x(1).n_cols,x(2).n_cols);
	  
	  for(int i=0; i < a.n_slices; i++){
	    for(int ii = 0; ii < a.n_cols; ii++)
	      a.slice(i).col(ii) = trans(x(0));
	  }
	  for(int j=0; j < b.n_slices; j++){
	    for(int jj= 0 ; jj < b.n_rows; jj++)
	      b.slice(j).row(jj) = x(1);
	  }
	  for(int k=0; k < c.n_slices; k++){
	    c.slice(k) = x(2)(k)*ones(c.n_rows,c.n_cols);
	  }
	  result(0) = a;
	  result(1) = b;
	  result(2) = c;
	  
	}
	return result;
	
}
// helper function - evaluates polynomial represented by vector of coefficients
double polyval(vec a,double n){
	double sum= 0;
	for(int i=0; i < a.n_elem; i++){
		sum += a(i)*pow((double) n, (int) a.n_elem-i-1);
	}
	return sum;
}


double product(mat a){
	double prod = 1.0;
	for(int i=0; i < a.n_cols; i++)
		prod*= a(0,i);
	return prod;
	
}
mat power(double x, mat p){
	mat result(1,p.n_cols);
	for(int i=0; i < p.n_cols ; i++){
		result(0,i) = pow(x,p(i));
    }
	return result;
}
mat findsub(mat A, urowvec b){
	mat result(1,b.n_cols);
	for(int i=0; i < b.n_cols; i++)
		result(0,i) = A(b(i));
	return result;
}
mat mod(mat A, int x){
	mat result(1,A.n_cols);
	for(int i=0; i < A.n_cols; i++)
		result(0,i) = (int)A(i) % x;
	return result;
}
//-----------end of old NDHV helper functions



/* NDHV 2obj function - sums rectangles for special case of 2 objectives
 * Inputs : points - matrix of size n x nObj of points on the Pareto Front 
 *			^^ are we considering only the Pareto or the entire population? Is Pareto the "best solutions from all generations"?
 *          refMin - vector of size nObj x 1 of lower limits on objective functions
 *          refMax - vector of size nObj x 1 of upper limits on objective functions
 * Output : NDHV - non dominated hypervolume value
 *
 */
double GAProblem::NDHV(mat points, vec refMin,vec refMax){
	
	time_t start_time, end_time;
	start_time = time (NULL);
		
	// Start of New Method for 2 Objectives
	double ndVol = 0.0;
	int m = points.n_cols;

	vec tlPt = zeros(2);
	vec brPt = zeros(2);
	
	tlPt(0) = refMin(0);
	tlPt(1) = refMax(1);
	brPt(0) = refMax(0);
	brPt(1) = refMin(1);
	
	points = join_cols(points, tlPt.t());
	points = join_cols(points, brPt.t());
	
	mat sortedPoints = sortByVec(points.col(0), points, 0);
	
	int numPoints = sortedPoints.n_rows;
	rowvec p1 = sortedPoints.row(0);
	rowvec p2 = sortedPoints.row(1);
	
	int i = 1;
	while (i < numPoints)
	{
		p2 = sortedPoints.row(i);
		if(p2(1) >= p1(1))
		{
			i ++;
		}
		else 
		{
			ndVol += p1(1) * (p2(0) - p1(0));
			p1 = sortedPoints.row(i);
			i ++;
		}
	}
	// End of new method for 2 objectives
	
	
	 
	// Start of Original Method
	/*
	int m = points.n_cols;
	
	vec dummy = 2*ones(m);
	rowvec refMin2 = refMin.t();
	rowvec refMax2 = refMax.t();

	mat aux = fullfact(dummy);
	aux = aux-ones(aux.n_rows,aux.n_cols);
	
	mat deltas = refMax2-refMin2;
	aux = aux.rows(1,aux.n_rows-2);

	mat auxPoints= zeros(aux.n_rows,aux.n_cols);

	for(int j=0; j < aux.n_rows ; j++)
		auxPoints.row(j) = refMin2 + deltas % aux.row(j);
		
	points = join_cols(points,auxPoints);
		
	int n = points.n_rows;
	m = points.n_cols;
	double ndVol = 0.0;
	
	mat sortedPoints = matSortCol(points);

	sortedPoints.reshape(1,m*n);
		
	field<mat> pointGrid(1,m);
	pointGrid(0) = sortedPoints;
	vec b = n*ones(m);
	field<mat> auxVar = mat2cell(sortedPoints,1,b);
    
	
	field<cube> pointC;
	field<mat> pointM;
	
	if(this->getnObj() == 3)
		pointC = ndgrid3(auxVar);
	else if (this->getnObj() == 2)
		pointM = ndgrid2(auxVar); 
	
	
	int nCells=  pow((double) n-1,m);
	int nNodes = pow((double) n,m);
	
	queue<int> theQ;  // queue of volumes to be added
	theQ.push(1);
	vec a = ones(m);
	int idxNextCorner = polyval(a,n);

	mat visited = ones(1,1);
	int i=0;
	int idxc1,idxc2;

	mat corner1(1,0);
	mat corner2(1,0);
	mat append = ones(1,1);
	uvec idx1(n);
	urowvec idx0 =  trans(idx1);

	mat p0(1,m);
	for(int k = 0; k < n; k++)
		idx0(k) = k;
	for(int k=0; k < m; k++)
		p0(0,k) = k+1;
	
	urowvec idx = idx0;
	mat p = p0;
	
	int nonDominated;
	double cellVol=0.0;
	mat idxNeighborCells(1,m);
	mat one = ones(1,p.n_cols);
	
    while(!theQ.empty()){
		i = theQ.front();
		
		
		theQ.pop();
		idxc1 = i-1;
		idxc2 = i+idxNextCorner-1;
		
		
		corner1.clear();
		corner2.clear();
		for(int j=0; j < m; j++){
			if(this->getnObj() == 3){
				corner1 = join_rows(corner1,(pointC(0,j))(idxc1)*append);
				corner2 = join_rows(corner2,(pointC(0,j))(idxc2)*append);
			}
			else if(this->getnObj() == 2){
				corner1 = join_rows(corner1,(pointM(0,j))(idxc1)*append);
				corner2 = join_rows(corner2,(pointM(0,j))(idxc2)*append);
			}
		}
		
		idx = idx0;
		
		for(int j=0; j < m; j++){
			uvec found = find(findsub(points.col(j),idx) < corner2(j));
			idx = trans(found);
		}
		
		if(idx.is_empty())
			nonDominated = 1;
		else nonDominated = 0;
		
		if(nonDominated){
			mat absolute = abs(corner2-corner1);
			cellVol = product(absolute);
			ndVol += cellVol;
			
			p = p0;      
			
			idxNeighborCells = i+power(n,p-one);
			
			
			idxNeighborCells = findsub(idxNeighborCells,trans(find(mod(idxNeighborCells,n)!=0)));
			
			idxNeighborCells = findsub(idxNeighborCells,trans(find(idxNeighborCells<=nNodes-idxNextCorner)));
			
			
			for(int j=0; j < idxNeighborCells.size(); j++){
				urowvec test = trans(find(idxNeighborCells(j) == visited));
				
				if(test.n_elem == 0){
					theQ.push(idxNeighborCells(j));
					
					visited.insert_cols(visited.n_elem,idxNeighborCells(j)*append);
					
				}
				
			}
			
		}
		
		
	}
	// End of original method
	 */
	
	end_time = time(NULL);
	cout << "NDHV time: " << difftime(end_time, start_time) << " seconds" << endl;
	
	return ndVol;
}


// showAll function 
//prints the variables of the GAProblem class to standard output
// Outputs: none
// Inputs: none
void GAProblem::showAll(){
  cout << "-----GAProblem Instance-----" << endl;
  cout << "pop:" << pop << endl;
  cout << "popObj:" << popObj << endl;
  cout << "nInputs:" << nInputs<< endl;
  cout << "nObj:" << nObj<< endl;
  cout <<"inputLowerLims:" << inputLowerLims<<endl;
  cout <<"inputUpperLims:" << inputUpperLims << endl;
  cout <<"outputLowerLims:"<<  outputLowerLims << endl;
  cout << "outputUpperLims:"<< outputUpperLims << endl;

  cout <<"d_lim" << d_lim << endl;
  cout << "L:" <<L; // see initializePopulation for description of converged and d_lim                                                    
  cout << "maxCrowdingDist:"<<maxCrowdingDist<<endl; // vector of L numbers from past L generations, each number is the max crowding dist in each geneneration.                                                                                                                                  


  for(int i=0; i < this->getnObj(); i++)
    cout <<"void **Ps:" <<objFuncParameterPointers[i] << ":" ;
cout << endl;

 for(int j=0; j< this->getnObj();j++)
   cout << "FPtrs:" << objFuncPointers[j]<<":"  ;
 cout << endl;

 for(int i=0; i < this->getnIneqConstr();i++)
   cout << "void**Pg" << ineqConstraintParameterPointers[i]<< ":";
 cout << endl;
 for(int i=0; i < this->getnIneqConstr();i++)
  cout << "GPtrs" << ineqConstraintPointers[i]<< ":";
 cout << endl;

 for(int i=0; i < this->getnEqConstr();i++)
   cout <<" void**Ph" << eqConstraintParameterPointers[i] << ":";
 cout << endl;
 for(int i=0; i < this->getnEqConstr();i++)
   cout << "HPtrs:" << eqConstraintPointers[i]<< ":";
 cout << endl;
 
 cout << "constraintType:" << this->getConstraintType()<< endl;
 cout << "nIneqConstr:"<< this->getnIneqConstr()<<endl;
 cout << "nEqConstr:" << this->getnEqConstr() <<endl;
 cout << "penaltyCoefficients:" << this->getPenaltyCoefficients() << endl;
 
}
