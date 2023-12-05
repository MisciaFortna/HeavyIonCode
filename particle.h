#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <string>
#include <fstream>
using namespace std;

class Particle{
	public:
		string PID = "0";
		double weight = 1.0;
		double xCord = 0.0;
		double yCord = 0.0;
		double zCord = 0.0;

		double xMom = 0.0;
		double yMom = 0.0;
		double zMom = 0.0;



		double getMass();

		int getA();

		int getZ();

		double getEnergy();

		void readout();

};

#endif
