/****************************************
 *file:GAControl.cc
 * contains most functions and variables of the GAControl class 
 * 
 *
 * 
 * Gary Yan, Peter Zhang
 * June 13,2012
 ****************************************/

#include <iostream>
#include <armadillo>
#include <stdlib.h>
#include "../include/GAControl.hh"
#include "../include/parentSelectionParams.hh"
#include <float.h>


using namespace std;
using namespace arma;


// get and set methods - since member variables are declared
// private, these methods are used for accessing and writing 
// to them, so that they are not accidently overwritten 

void GAControl::setPopSize(int popSize){
  this->popSize = popSize;
}

int GAControl::getPopSize(){
  return this->popSize;
}
void GAControl::setInitialSeed(int seed){
  this->initialSeed = seed;
}

int GAControl::getInitialSeed(){
  return this->initialSeed;
}
parentSelectionParams GAControl::getParentSelectionParams(){
  return this->params;
}
void GAControl::setParentSelectionParams(parentSelectionParams param){
  this->params = param;
}
recombinationParams GAControl::getRecombinationParams(){
  return this->rparams;
}
void GAControl::setRecombinationParams(recombinationParams rparam){
  this->rparams = rparam;
}

mutationParams GAControl::getMutationParams(){
  return this->mparams;
}
void GAControl::setMutationParams(mutationParams mparam){
  this->mparams = mparam;
}

void GAControl::showAll(){

  cout <<"-----GAControl-----" << endl;
  cout <<"popSize:" << popSize << endl;
  cout << "initialSeed:"<< initialSeed << endl;
  cout << "pparams:"<< params.nParents<<endl;
  cout <<"rparams:"<< rparams.nOffspring<<":"<<rparams.recombProb<<endl;
  cout << "mparams:"<< mparams.mutationProb<<endl;
  cout << "Output:" << fileName << endl;
  cout << "distIndex:"<<distIndex << endl;
}

void GAControl::setOutputFile(string s){
  this-> fileName = s;
}
 string GAControl::getOutputFile(){
  return this->fileName;
}
void GAControl::setFitnessAssignmentType(string s){
  this->fitnessAssignmentType = s;
}
string GAControl::getFitnessAssignmentType(){
  return this->fitnessAssignmentType;
}
void GAControl::setRecombinationType(string s){
  this->recombinationType = s;
}
string GAControl::getRecombinationType(){
  return this->recombinationType;
}
void GAControl::setMutationType(string s){
  this->mutationType = s;
}
string GAControl::getMutationType(){
  return this->mutationType;
}
void GAControl::setParentSelectionType(string s){
  this->parentSelectionType = s;
}
string GAControl::getParentSelectionType(){
  return this->parentSelectionType;
}

double GAControl::getDistIndex(){
  return this->distIndex;
}
void GAControl::setDistIndex(double distIndex){
  this->distIndex = distIndex;
}

double GAControl::getnGens(){
  return this->nGens;
}
void GAControl::setnGens(int nGens){
  this->nGens = nGens;
}

GAControl::GAControl(){
  //  int rootn = (int)round(sqrt(G->getnInputs())+1.5);
  //cout << "nInputs:"<<G->getnInputs()<<"rootn:"<<rootn;
  //this->nGens = 2*rootn;
}
GAControl::~GAControl(){

}
