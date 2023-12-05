#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include "particle.h"
#include "nucleusColl.h"
#include <nlohmann/json.hpp>
using namespace std;
using json = nlohmann::json;

const double gamma = 0.57721566;

const double stepSize = 1e-4; // cm

/*************************** Soft Collision Functions ********************************/

double factorial(double num);

double Si(double y);

double Ci(double y);

double f1(double y, double v);

double f2(double y);

double bethe(Particle projectile, Material target);

double vavilovFunc(double thickness, double Q, Particle projectile, Material target);


/*************************************************************************************/

/************************ Angular Distribution Functions *****************************/

double meanSquareScatterAngle(Particle projectile, Material target);

double angularPDF(double alpha, double theta);

double sampleScatterAngle(double crossSection, Particle projectile, Material target);

/*************************************************************************************/

/*********************** Light Particle Interactions *********************************/
double sigmaFinder(json data, double energy);

double* lightParticleCSBasic(Particle projectile, string targetMaterial);
/*************************************************************************************/

int main(){

	double seed = time(NULL);

	srand(seed);


	return 0;
}

double factorial(int num){

	double value = 1.0;

	for (int i = 2; i <= num; i++){

		value = value * ((double) i);

	}

	return value;

}

double Si(double y){

	double value = 0;

	for(int i = 0; i < 11; i++){

		value += pow(-1, i) * pow(y, (2 * i + 1)) / ((2 * i + 1) * factorial((2 * i + 1)));
	}

	return value;

}

double Ci(double y){

	double value = 0;

	for(int i = 0; i < 11; i++){

		value += pow(-1, i) * pow(y, (2 * i)) / ((2 * i) * factorial((2 * i)));
	}

	value += (gamma + log(y));

	return value;

}

double f1(double y, double v){
	
	double value = pow(v, 2) * (log(y) - Ci(y)) - y * Si(y) - cos(y);

	return value;

}

double f2(double y){

	double value = y * (log(y) - Ci(y)) + pow(v, 2) * Si(y) + sin(y);

	return value;

}

double bethe(Particle projectile, Material target);

double vavilovFunc(double thickness, double Q, Particle projectile, Material target){

	double charge = (double) projectile.getZ();

	double velocity = vel(sqrt(pow(projectile.xMom, 2) + pow(projectile.yMom, 2) + pow(projectile.zMom, 2)), projectile.getMass()); // vel from nucleusColl.h

	double ksi = 0.15343 * pow(charge, 2) * target.ZEFF * target.density / pow(velocity, 2) / target.AEFF;

	double beta = bethe(projectile, target);

	double delE = Q;

	double k = ksi / delE;

	double lambda = (Q - beta * thickness) / delE - k * (1 + pow(velocity, 2) - gamma);

	double alpha = exp(k * (1 + gamma * pow(velocity, 2))) / (delE * 3.141592653);

	double count = (0.001 / 2) * exp(k * f1(0.0)) * cos(k * f2(0.0));

	count += (0.001 / 2) * exp(k * f1(100.0)) * cos(lambda * 100 + k * f2(100.0));

	for (double i = 0.001; i < 100.0; i += 0.001){

		count += 0.001 * exp(k * f1(i)) * cos(lambda * i + k * f2(i));

	}

	return (alpha * count);


}


double meanSquareScatterAngle(Particle projectile, Material target){

	double Z = target.ZEff;

	double velocity = vel(sqrt(pow(projectile.xMom, 2) + pow(projectile.yMom, 2) + pow(projectile.zMom, 2)), projectile.getMass()); // vel from nucleusColl.h
	double kineticEnergy = projectile.getMass() / sqrt(1 + pow(velocity, 2));

	double minTheta = 3.723e-3 * sqrt(Z, (1.0 / 3.0)) / sqrt(kineticEnergy * (kineticEnergy + 2 * projectile.getMass()));

	double value = 4 * pow(minTheta, 2) * log(183 / pow(Z, (1.0 / 3.0)));

	return value;


}

double angularPDF(double alpha, double theta){

	double value = (2 * theta / alpha) * exp(-1 * pow(theta, 2) / alpha);

	return value;

}

double sampleScatterAngle(double crossSection, Particle projectile, Material target){

	double alpha = crossSection * stepSize * meanSquareScatterAngle(projectile, target);

	double maxAngle = sqrt(alpha / 2.0);

	double eta1 = ((double) rand() / (RAND_MAX));

	double eta2 = ((double) rand() / (RAND_MAX));

	double x = 3.141592653 * eta1;

	while ((eta2 * angularPDF(alpha, maxAngle)) > angularPDF(alpha, x)){

		eta1 = ((double) rand() / (RAND_MAX));

		eta2 = ((double) rand() / (RAND_MAX));

		x = 3.141592653 * eta1;


	}

	return x;

}

double sigmaFinder(json data, double energy){
	double eFixed = energy * 1e+6; // eV
	int trueLine = 0;
	int lineBelow = 0;
	double sigma = 0.0;
	double m = 0.0;
	int datasetSize = data["datasets"][0]["pts"].size();
	if (eFixed > data["datasets"][0]["pts"][datasetSize-1]["E"]){
		sigma = data["datasets"][0]["pts"][datasetSize-1]["Sig"];
	}
	else{
		for (int i = 0; i < datasetSize; i++){
			if (data["datasets"][0]["pts"][i]["E"] >= eFixed){
				trueLine = i;
				lineBelow = i - 1;
				break;
			}
		}
	}

	if (data["datasets"][0]["pts"][trueLine]["E"] == eFixed){
		sigma = data["datasets"][0]["pts"][trueLine]["Sig"];
	}

	else{
		m = (data["datasets"][0]["pts"][trueLine]["Sig"] - data["datasets"][0]["pts"][lineBelow]["Sig"]) / (data["datasets"][0]["pts"][trueLine]["E"] - data["datasets"][0]["pts"][lineBelow]["E"]);

		sigma = m * (eFixed - data["datasets"][0]["pts"][lineBelow]["E"]) + data["datasets"][0]["pts"][lineBelow]["Sig"];


	}

	return sigma;
}

double* lightParticleCSBasic(Particle projectile, string targetMaterial){

	string particle;

	if (projectile.PID == "2112"){
		particle = "N";
	}

	else{
		particle = "P";
	}

	string fileNameInit = targetMaterial + "_" + particle;

	// init sigmaTOT
	
	double sigmaTOT;

	fileInput = "/ParticleCS/" + fileNameInit + "_TOT.json";

	ifstream f(fileInput);

	json j = json::parse(f);

	sigmaTOT = sigmaFinder(j, projectile.getEnergy());

	f.close();

	// init sigmaEL
	
	double sigmaEL;

	fileInput = "/ParticleCS/" + fileNameInit + "_EL.json";

	ifstream f(fileInput);

	j = json::parse(f);

	sigmaEL = sigmaFinder(j, projectile.getEnergy());

	f.close();

	// assemble array
	
	double sigmaNON = sigmaTOT - sigmaEL;

	double outputCS[3] = {sigmaTOT, sigmaEL, sigmaNON};

	return outputCS;
}
