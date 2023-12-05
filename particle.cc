#include <iostream>
#include <string>
#include <fstream>
#include "particle.h"
using namespace std;

double Particle::getMass(){

	double mass;

	if (PID.size() == 4){

		if (PID == "2212"){
			mass = 938.27;
		}

		else{
			mass = 939.57;
		}

	}
	
	else if (PID.size() == 10){

		string tempPID;
		string line;
		int curLine = 0;
		int searchLine = 0;

		ifstream fin;

		fin.open("iMass.txt");

		while(getline(fin, line)){
			curLine++;
			if(line.find(PID, 0) != string::npos){
				searchLine = curLine;
			}
		}

		fin.close();

		fin.open("iMass.txt");

		for(int i = 0; i < searchLine; i++){
			fin >> tempPID >> mass;
		}

		fin.close();
	}

	else{
		mass = 0;
	}


	return mass;
}

int Particle::getA(){
	
	int A;

	if (PID.size() == 10){
		int tens = ((int) PID[4]) - ((int) '0');
		int ones = ((int) PID[5]) - ((int) '0');

		A = tens * 10 + ones;
	}

	else if (PID.size() == 4){
		A = 1;
	}

	else{
		A = 0;
	}

	return A;

}

int Particle::getZ(){

	int Z;

	if (PID.size() == 10){
		int tens = ((int) PID[7]) - ((int) '0');
		int ones = ((int) PID[8]) - ((int) '0');
		Z = tens * 10 + ones;
	}

	else if (PID == "2212"){
		Z = 1;
	}

	else{
		Z = 0;
	}

	return Z;

}

double Particle::getEnergy(){

	double momRad = sqrt(pow(xMom, 2) + pow(yMom, 2) + pow(zMom, 2));

	double velocity = sqrt(pow(momRad, 2) / (pow(getMass(), 2) + pow(momRad, 2)));

	double correction = momRad / (getMass() * velocity);

	double T = correction * getMass();

	return T;

}

void Particle::readout(){
	cout << "PID: " << PID << endl;
	cout << "Mass: " << getMass() << " MeV/c^2" << endl;
	cout << "Position: (" << xCord << ", " << yCord << ", " << zCord << ") fm" << endl;
	cout << "Momentum: (" << xMom << ", " << yMom << ", " << zMom << ") MeV/c" << endl;
}
