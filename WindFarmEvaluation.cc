#include "WindFarmEvaluation.hh"
#include <math.h>
#include <iostream>
#include <cmath>
#include <armadillo>
#define PI 3.1415926535897932 // less expensive than acos(-1)
using namespace std;
using namespace arma;

WindFarmEvaluation::WindFarmEvaluation()
{}

//TODO expanding windSpeed and theta to a more general case - for example 36 bins
//TODO combine the turbine height into the coordinates vector so that we can have variable turbine heights later

//TODO optimize the code
// method 1: reduce the time in half by checking i and j only once (still n^2 run time)
// method 2: calculate the "field of influece" for each turbine and remember; then use the polygon check for downstream/upstream (might be slow for the first evalution, but faster later)
double WindFarmEvaluation::calculatePower(vec coordinates, 
					  vec windDirections, 
					  mat windSpeeds,
					  mat windProb,
					  double z, // can later change to a vector with variable height
					  double drot, 
					  double z0,
					  double (*powerCurve)(double))
{
        //initialization
	double r = drot/2;
	int numTurbines = coordinates.n_elem / 2;
	double totalPower = 0.0;
	double upperbound = 0.0;
	double alpha = 0.5 / log(z/z0); // assume turbine height is constant for now

	
	for (unsigned w = 0; w < windDirections.n_elem; w ++)
	{
		double theta = windDirections(w);

		mat downstream = ones(numTurbines, numTurbines) * -1; //1 if i is in downstream of j, 0 otherwise, -1 if not initialized

		for (int i = 0; i < numTurbines; i ++)
		{
			for (int j = 0; j < numTurbines; j ++)
			{
				if (downstream(i,j) == -1) //not already initialized
				{	
					/*
					 \      /
					  \  * /<-- centre of turbine i
					   \__/ <-- turbine j
						\/  
					*/
					downstream(i,j) = 0;

					if (i != j)
					{
						double al = atan(alpha); //half angle of the wake spread

						vec pt1(2); // coordinates of the "left" end of turbine j
						vec pt2(2); // coordinates of the "right" end of turbine j
						pt1(0) = coordinates(j) + r * cos(theta + PI/2);
						pt2(0) = coordinates(j) + r * cos(theta - PI/2);
						pt1(1) = coordinates(j + numTurbines) + r * sin(theta + PI/2);
						pt2(1) = coordinates(j + numTurbines) + r * sin(theta - PI/2);

						vec r1(2); //vector pointing from pt1 to turbine i
						vec r2(2);
						r1(0) = coordinates(i) - pt1(0);
						r1(1) = coordinates(i + numTurbines) - pt1(1);
						r2(0) = coordinates(i) - pt2(0);
						r2(1) = coordinates(i + numTurbines) - pt2(1);
						
						vec r1_hat = r1 / norm(r1, 2); // normalized r1
						vec r2_hat = r2 / norm(r2, 2);

						double halfangle = (al + PI / 2) / 2;

						vec pt1_c(2); // centre line at pt1 between wake left boundary and turbine
						vec pt2_c(2);
						pt1_c(0) = cos(theta - PI/2 + halfangle);
						pt1_c(1) = sin(theta - PI/2 + halfangle);
						pt2_c(0) = cos(theta + PI/2 - halfangle);
						pt2_c(1) = sin(theta + PI/2 - halfangle);

						if (dot(pt1_c, r1_hat) >= cos(halfangle) && dot(pt2_c, r2_hat) >= cos(halfangle))
						{
							downstream(i, j) = 1;
							downstream(j, i) = 0; // two turbines cannot be in downstream of each other at the same time
						}
					}
				}
			}
		}


		vec effSpeed (numTurbines); //stores final effective wind speed at i
		vec defSumSqr = zeros(numTurbines); // sum of the sqr of deficits (0~1) from multiple wakes
		
		// Calculate the total power produced by turbine i
		// First, calculate the total wake effects on turbine i
		// Second, calculate the resulting wind speed and power output
		for (int i = 0; i < numTurbines; i ++)
		{
			for (int j = 0; j < numTurbines; j ++)
			{
				if (j != i && downstream(i,j) == 1) // if downstream and different turbines
				{
					double dd = (coordinates(i) - coordinates(j)) * cos(theta) + (coordinates(i + numTurbines) - coordinates(j + numTurbines)) * sin(theta);
					double deficit = (2.0 * 0.32679) * pow(r / (r + alpha * dd) , 2.0);
					defSumSqr(i) += pow(deficit, 2);
				}
			}
		}

		for (unsigned l = 0; l < windSpeeds.n_cols; l ++)
		{
			for (int i = 0; i < numTurbines; i ++)
			{
				effSpeed[i] = windSpeeds(w,l) * (1 - sqrt(defSumSqr(i)));
				totalPower += powerCurve(effSpeed(i)) * windProb(w,l);
				upperbound += powerCurve(windSpeeds(w,l)) * windProb(w,l);
			}
		}
	}

	double efficiency = totalPower / upperbound;

	//cout << "Total power: " << totalPower << " (" << efficiency * 100 << "%)" << endl;

	return  totalPower;
}

vec WindFarmEvaluation::calculateSoundField(vec coordinates, 
											double z,
											vec receptorCoordinates,
											double G,
											double Gs,
											double Gr)
{	
	//Sound Power of each point source
	double Lw = 105; // in dB (from UMass presentation)

	//Atmospheric attenuation coefficient for octave bands of noise at 20 Degree Celsius [dB/km]
	vector <double> alpha (8);
	alpha[0] = 0.1; // Nominal Frequency of 63Hz
	alpha[1] = 0.3; // Nominal Frequency of 125Hz
	alpha[2] = 1.1; // Nominal Frequency of 250Hz
	alpha[3] = 2.8; // Nominal Frequency of 500Hz
	alpha[4] = 5.0; // Nominal Frequency of 1000Hz
	alpha[5] = 9.0; // Nominal Frequency of 2000Hz
	alpha[6] = 22.9; // Nominal Frequency of 4000Hz
	alpha[7] = 76.6; // Nominal Frequency of 8000Hz

	// Standard A Weighting 
	vector <double> Af (8);
	Af[0] = -26.22; // Nominal Frequency of 63Hz
	Af[1] = -16.19; // Nominal Frequency of 125Hz
	Af[2] = -8.68; // Nominal Frequency of 250Hz
	Af[3] = -3.25; // Nominal Frequency of 500Hz
	Af[4] = 0; // Nominal Frequency of 1000Hz
	Af[5] = 1.2; // Nominal Frequency of 2000Hz
	Af[6] = 0.96; // Nominal Frequency of 4000Hz
	Af[7] = -1.14; // Nominal Frequency of 8000Hz

	// number of receptors
	unsigned numR = receptorCoordinates.n_elem / 3;
	// number of turbines
	unsigned numT = coordinates.n_elem / 2; // TODO change to 3D and change this

	// reference distance in calculating the Geometrical Divergence
	double d_o = 1; // in m

	// Equivalent continuous A-weighted downwind sound pressure level for receptors
	vec Lat (numR);

	for (unsigned i = 0; i < numR; i ++)
	{
		double OuterSummation = 0.0;
		for (unsigned j = 0; j < numT; j ++)
		{
			vec turbinePosition (3);
			vec receptorPosition (3);

			turbinePosition(0) = coordinates(j);
			turbinePosition(1) = coordinates(j + numT);
			turbinePosition(2) = z;
			
			receptorPosition(0) = receptorCoordinates(i);
			receptorPosition(1) = receptorCoordinates(i + numR);
			receptorPosition(2) = receptorCoordinates(i + numR * 2);

			//obtain distance between source and receptor: 3D
			double d = norm(turbinePosition - receptorPosition, 2);

			//attenuation due to geometrical divergence, in dB
			double Adiv = 20 * log10(d / d_o) + 11;

			//calculate Aatm for 8 nominal frequencies, attenuation due to atm absorption
			vector <double> Aatm(8);
			for (int k = 0; k < 8; k ++)
				Aatm[k] = alpha[k] * d / 1000;
			
			//double G = 0;
			//double Gs = 0;
			//double Gr = 0;

			vector <double> Agr = Aground(coordinates[j], z, receptorCoordinates[i], receptorCoordinates[i + numR * 2], G, Gs, Gr);

			//Total attenuation at 8 nominal frequencies
			vector <double> A (8);
			for (int k = 0; k < 8; k ++)
				A[k] = Adiv + Aatm[k] + Agr[k];

			// Calculate DOmega
			//double hm = 0.5 * (z + receptorCoordinates[i + 2 * numR]);
			//double dp = sqrt(pow(receptorCoordinates[i] - coordinates[j], 2) + pow(receptorCoordinates[i + numR] - coordinates[j + numT], 2));
			//double DOmega = 10 * log10(1 + ((pow(d,2) / (pow(dp,2) + pow(2 * hm, 2)))));
			double DOmega = 0.0;

			// Question: why "for k 1:1"?

			// Cauculate the equivalent continuous downwind octave-band sound pressure level at receiver location Lft_DW
			vector <double> Lft (8);
			for (int k = 0; k < 8; k ++)
				Lft[k] = Lw + DOmega - A[k];

			// Calculate the summation terms inside the Lat term
			// Let InnerSummation = sum the Lft and Af for 8 nominal frequencies of each wind turbine
			double InnerSummation = 0;
			for (int k = 0; k < 8; k ++)
			{
				double exponent = 0.1 * (Lft[k] + Af[k]);
				InnerSummation += pow(10, exponent);
			}

			OuterSummation += InnerSummation;
		}
		Lat(i) = 10 * log10(OuterSummation);
	}

	return Lat;
}

vector <double> WindFarmEvaluation::Aground(double xt1, double zt1, double xr1, double zr1, double G, double Gs, double Gr)
{
	vector <double> agr (8);

	double hs = zt1;
	double dp = abs(xt1 - xr1);
	double as = 1.5 + 3 * exp(-0.12 * (pow(hs - 5, 2))) * (1 - exp(-dp / 50)) + 5.7 * exp(-0.09 * (pow(hs, 2))) * (1 - exp((-2.8 * pow(10.0, -6.0)) * pow(dp, 2)));
	double bs = 1.5 + 8.6 * exp(-0.09 * pow(hs, 2)) * (1 - exp(-dp/50));
	double cs = 1.5 + 14 * exp(-0.46 * pow(hs, 2)) * (1 - exp(-dp/50));
	double ds = 1.5 + 5 * exp(-0.9 * pow(hs, 2)) * (1 - exp(-dp/50));

	vector <double> As (8);
	
	As[0] = -1.5;
	As[1] = -1.5 * Gs * as;
	As[2] = -1.5 + Gs * bs;
	As[3] = -1.5 + Gs * cs;
	As[4] = -1.5 + Gs * ds;
	As[5] = -1.5 * (1 - Gs);
	As[6] = -1.5 * (1 - Gs);
	As[7] = -1.5 * (1 - Gs);

	double hr = zr1;
	double ar = 1.5 + 3.0 * exp(-0.12 * pow((hr - 5), 2)) * (1.0 - exp(-dp / 50)) + 5.7 * exp((-0.09 * pow(hr, 2))) * (1 - exp((-2.8 * pow(10.0, -6)) * pow(dp, 2)));
	double br = 1.5 + 8.6 * exp(-0.09 * pow(hr, 2)) * (1 - exp(-dp / 50));
	double cr = 1.5 + 14.0 * exp(-0.46 * pow(hr, 2)) * (1 - exp(-dp / 50));
	double dr = 1.5 + 5.0 * exp(-0.9 * pow(hr, 2)) * (1 - exp(-dp / 50));

	vector <double> Ar (8);
	Ar[0] = -1.5;
	Ar[1] = -1.5 * Gr * ar;
	Ar[2] = -1.5 + Gr * br;
	Ar[3] = -1.5 + Gr * cr;
	Ar[4] = -1.5 + Gr * dr;
	Ar[5] = -1.5 * (1 - Gr);
	Ar[6] = -1.5 * (1 - Gr);
	Ar[7] = -1.5 * (1 - Gr);

	double q;

	if (dp <= 30 * (hs + hr))
		q = 0.0;
	else
		q = 1 - (30 * (hs + hr) / dp);

	vector <double> Am (8);
	Am[0] = -3 * (pow(q, 2));
	Am[1] = -3 * q * (1 - G);
	Am[2] = -3 * q * (1 - G);
	Am[3] = -3 * q * (1 - G);
	Am[4] = -3 * q * (1 - G);
	Am[5] = -3 * q * (1 - G);
	Am[6] = -3 * q * (1 - G);
	Am[7] = -3 * q * (1 - G);

	for (int k = 0; k < 8; k ++)
		agr[k] = As[k] + Ar[k] + Am[k];

	return agr;
}

double WindFarmEvaluation::dotProduct(vector <double> x, vector <double> y)
{
	double sum = 0.0;
	for (unsigned i = 0; i < x.size(); i ++)
		sum += x[i] * y[i];
	return sum;
}

vector <double> WindFarmEvaluation::vecSub(vector <double> x, vector <double> y)
{
	vector <double> ret (x.size());
	for (unsigned i = 0; i < ret.size(); i ++)
		ret[i] = x[i] - y[i];
	return ret;
}

vector <double> WindFarmEvaluation::normalize(vector <double> x)
{
	double card = 0.0;
	for (unsigned i = 0; i < x.size(); i ++)
		card += pow(x[i],2);
	card = sqrt(card);
	
	for (unsigned i = 0; i < x.size(); i ++)
		x[i] = x[i] / card;
	
	return x;
}

void WindFarmEvaluation::printVector(vector <double> x)
{
	for (unsigned i = 0; i < x.size(); i ++)
		cout << x[i] << " ";
	cout << endl;
}

