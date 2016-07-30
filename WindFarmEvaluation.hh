#ifndef WindFarmEvaluation_H
#define WindFarmEvaluation_H

#include <vector>
#include <armadillo>
using namespace std;
using namespace arma;

double powerCurve(double speed);
struct objP{
  double z; // turbine height
  double z0;
  int numDir;
  int numTurbines;
  double r; // turbine radius
  double drot; // turbine diameter
  vec windDirections;
  mat windSpeeds;
  mat windProb;
};
struct objP2{
  ;
};
struct objP3{
  ;
  };

//TODO: the turbine coordinates into 3 dimension vector, and remove the height variable
class WindFarmEvaluation
{
public:
	WindFarmEvaluation();

	// Calculation of the 
	//coordinates: all x, followed by all y
	//theta: wind direction, counterclockwise, pointing EAST is theta = 0
	//z: turbine height, can be a vector later
	//drot: rotor diameter
	//z0: roughness of terrain
	//powerCurve: pointer to the function that calculates power based on wind speed
	double calculatePower(vec coordinates,  // turbine locations
						  vec windDirections, 
						  mat windSpeeds,
						  mat windProb,
						  double z,
						  double drot,
						  double z0,
						  double (*powerCurve)(double));

	// coordinates: 2 dimensional coordinates of turbines, x coordinates, followed by y, followed by z
	// z: height of all turbines - can be put into the previous vector
	// receptorCoordinates: 3 dimensional coordinates of receptors, x coordinates, followed by y, followed by z
	// Returns a vector of sound levels at each receptor
	vec calculateSoundField(vec coordinates, 
							double z,
							vec receptorCoordinates,
							double G,
							double Gs,
							double Gr);
	//helper function for sound calculation
	vector <double> Aground(double xt1, double zt1, double xr1, double zr1, double G, double Gs, double Gr);

	//Helper functions for now, use a linear algebra library later
	double dotProduct(vector <double> x, vector <double> y);
	vector <double> vecSub(vector <double> x, vector <double> y);
	vector <double> normalize(vector <double> x);
	void printVector(vector <double> x);
};

#endif
