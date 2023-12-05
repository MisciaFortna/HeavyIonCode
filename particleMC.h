//particleMC.h
#include <string>
#include <vector>
#include "particle.h"
using namespace std;

class Material {
	public:
		vector<string> elements("H");
		vector<double> elementRatios(1.0);
		double density = 0.0;			// g/cc
		double meanExcitationEnergy = 0.0;	// eV
		double atomicMass = 0.0;		// g/mol
		double ZEff = 0.0;
		double AEff = 0.0;

};

class Target {
	public:
		Material composition;
		double HLX = 0.0;	// cm
		double HLY = 0.0;	// cm
		double HLZ = 0.0;	// cm

		double xCord = 0.0; // cm
		double yCord = 0.0; // cm
		double zCord = 0.0; // cm

		int inTarget(Particle p);

};

class Source {
	private:
		string PID = "0";	// 10LZZZAAAI for heavy ions
		double mass = 0.0;	// MeV/c^2
		int Z = 0;		// Atomic Number
		int A = 0;		// Atomic Mass
		int N = 0;		// Number of Neutrons

	public:

		double posX = 0.0;	// cm Describes the average initial x position for spawn sampling
		double posY = 0.0;	// cm Describes the average initial y position for spawn sampling
		double posZ = 0.0;	// cm
		double xCutoff = 1.0;	// cm Describes half of x range cutoff for spawn sampling
		double yCutoff = 1.0;	// cm Describes half of y range cutoff for spawn sampling
		double angCutoff = 0.2;	// deg Describes half of polar angle cutoff for spawn sampling
		double energy = 0.0;	// MeV

};


