//Nucleus Coll Code
#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

class Particle{
	public:
		string PID = "0";
		double mass = 0.0;
		double xCord = 0.0;
		double yCord = 0.0;
		double zCord = 0.0;

		double xMom = 0.0;
		double yMom = 0.0;
		double zMom = 0.0;
};


// all in fm so 1e-15

double* particleSpawn(string PID, double avgR, double avgP);

double crossSection(double energy);

int finder(int el, vector<int> listSet);

double vel(double p, double m); // given in units of c

void coll(Particle &nuc1, Particle &nuc2);

void dispLog(int counter, Particle sys[], int time);

double distance(Particle p1, Particle p2);

double dHalf(Particle p1, Particle p2);

double dubU(Particle p1, Particle p2);

double tripU(Particle p1, Particle p2, Particle p3);

double U(int count, Particle set[]);

void readout(Particle p);

int elChecker(vector<int> v, int x);

int R(Particle p1, Particle p2);

int P(Particle p1, Particle p2);

int main(){

	srand(time(NULL));
	
	double bMax = 3.19848732601 - 0.8;

	double impactParam = bMax * ((double) rand()) / (RAND_MAX);
	double stepSize = 3.32236936788 / 2.5;
	int i;

	double mass = 12 * 931.49410242; // mass of carbon atom

	double energy = 200 * 12; // KE of carbon atom

	double labP = mass * sqrt(pow((energy / mass), 2) + (2 * energy / mass));

	double comPProjectile = -1 * labP / 2;
	double comPTarget = labP / 2;

	double comR = 0.0;

	double* tempResults;

	Particle nuc1[12];

	Particle nuc2[12];

	for (i = 0; i <= 5; i++){
		nuc1[i].PID = "2212";
		nuc2[i].PID = "2212";
		nuc1[i].mass = 939.57; // MeV/c^2
		nuc2[i].mass = 939.57; // MeV/c^2
	}

	for (i = 6; i <= 11; i++){
		nuc1[i].PID = "2112";
		nuc2[i].PID = "2112";
		nuc1[i].mass = 938.28; // MeV/c^2
		nuc2[i].mass = 938.28; // MeV/c^2
	}


	for (i = 0; i <= 11; i++){
		tempResults = particleSpawn(nuc1[i].PID, comR, comPProjectile);

		nuc1[i].xCord = tempResults[0];
		nuc1[i].yCord = tempResults[1];
		nuc1[i].zCord = tempResults[2];
		nuc1[i].xMom = tempResults[3];
		nuc1[i].yMom = tempResults[4];
		nuc1[i].zMom = tempResults[5];
	}

	for (i = 0; i <= 11; i++){
		tempResults = particleSpawn(nuc2[i].PID, comR, comPTarget);

		nuc2[i].xCord = tempResults[0];
		nuc2[i].yCord = tempResults[1];
		nuc2[i].zCord = tempResults[2];
		nuc2[i].xMom = tempResults[3];
		nuc2[i].yMom = tempResults[4];
		nuc2[i].zMom = tempResults[5];
	}

	// change to COM Frame
	for (i = 0; i <= 11; i++){
		nuc1[i].yCord += (impactParam / 2);
		nuc1[i].zCord -= (5 * stepSize / 2);

		nuc2[i].yCord -= (impactParam / 2);
		nuc2[i].zCord += (5 * stepSize / 2);
	}

	// Collision Work
	// Init Field/Coll System
	Particle sys[24];

	int j = 0;

	for (i = 0; i <= 11; i++){
		sys[i] = nuc1[i];
		j = i + 12;
		sys[j] = nuc2[i];
	}

	int count = sizeof(sys) / sizeof(Particle);

	// step 1: Mean Field Movement
	//
	
	Particle tempSet1[24];
	Particle tempSet2[24];

	double h = 1e-6;
	double pXDot = 0.0;
	double pYDot = 0.0;
	double pZDot = 0.0;


	for (int t = 0; t < 6; t++){

	for (i = 0; i < count; i++){

		pXDot = 0.0;
		pYDot = 0.0;
		pZDot = 0.0;
/*
		for (i = 0; i < count; i++){
			readout(sys[i]);
		}**/

		copy(begin(sys), end(sys), begin(tempSet1));
		copy(begin(sys), end(sys), begin(tempSet2));


		tempSet1[i].xCord += h;
		pXDot = -1 * (U(count, tempSet1) - U(count, sys)) / h;
		//cout << pXDot << endl;
		tempSet1[i].xCord -= h;
		//copy(begin(sys), end(sys), begin(tempSet1));

		tempSet1[i].yCord += h;
		pYDot = -1 * (U(count, tempSet1) - U(count, sys)) / h;
		tempSet1[i].yCord -= h;
		//copy(begin(sys), end(sys), begin(tempSet1));

		tempSet1[i].zCord += h;
		pZDot = -1 * (U(count, tempSet1) - U(count, sys)) / h;
		tempSet1[i].zCord -= h;
		//copy(begin(sys), end(sys), begin(tempSet1));

		tempSet2[i].xMom += stepSize * pXDot;
		tempSet2[i].yMom += stepSize * pYDot;
		tempSet2[i].zMom += stepSize * pZDot;

	}


	for (i = 0; i < (count); i++){

		sys[i].xCord += stepSize * sys[i].xMom / sqrt(pow(sys[i].mass, 2) + pow(sys[i].xMom, 2));
		sys[i].xMom = tempSet2[i].xMom;

		sys[i].yCord += stepSize * sys[i].yMom / sqrt(pow(sys[i].mass, 2) + pow(sys[i].yMom, 2));
		sys[i].yMom = tempSet2[i].yMom;

		sys[i].zCord += stepSize * sys[i].zMom / sqrt(pow(sys[i].mass, 2) + pow(sys[i].zMom, 2));
		sys[i].zMom = tempSet2[i].zMom;
	}
/*
	cout << "\nValues Post Field" << endl;
	for (i = 0; i < count; i++){
		cout << "Particle " << i + 1 << ": Pos -> (" << sys[i].xCord << ", " << sys[i].yCord << ", " << sys[i].zCord << ")" << endl;
		cout << "Particle " << i + 1 << ": Mom -> (" << sys[i].xMom << ", " << sys[i].yMom << ", " << sys[i].zMom << ")" << endl;
	}**/

	// step 2: Collision Test
	vector<int> listOfColl{0};

	for (i = 1; i < count; i++){
		if (distance(sys[0], sys[i]) <= dHalf(sys[0], sys[i])){
			coll(sys[0], sys[i]);
			listOfColl.insert(listOfColl.end(), 0);
			listOfColl.insert(listOfColl.end(), i);
			break;
		}
	}

	if (listOfColl.size() == 3)
		listOfColl.erase(listOfColl.begin());


	for (i = 1; i < count; i++){
		for (int j = i + 1; j < count; j++){
			if ((1 - finder(i, listOfColl)) && (1 - finder(j, listOfColl)) && (distance(sys[i], sys[j]) <= dHalf(sys[i], sys[j]))){
				coll(sys[i], sys[j]);
				listOfColl.insert(listOfColl.end(), i);
				listOfColl.insert(listOfColl.end(), j);
				break;
			}
		}
	}
/*
	cout << "\nList of Collisions" << endl;
	for (i = 0; i < listOfColl.size(); i++){
		cout << listOfColl[i] << endl;
	}**/

	for (i = 0; i < count; i++){
		if(abs(sys[i].xMom) <= 1e-5){
			sys[i].xMom = 0;
		}
		
		if(abs(sys[i].yMom) <= 1e-5){
			sys[i].yMom = 0;
		}

		if(abs(sys[i].zMom) <= 1e-5){
			sys[i].zMom = 0;
		}
	}

//	for (int i: listOfColl)
//		cout << i << endl;

	//dispLog(count, sys, t);

	}

	// Frag Identification
	vector<vector<int>> fragSet; // sets of fragments (based on elements of sys)
	vector<int> pairList; // elements of sys
	vector<int> tempList;

	// Identification of initial frag (sys[0])
	pairList.insert(pairList.begin(), 0);

	int tempCounter = 1;
	int tempCounter2 = 1;
	int tempPCounter = 0;
	int tempNCounter = 0;

	for(i = 1; i < count; i++){
		if (R(sys[0], sys[i]) && P(sys[0], sys[i])){
			pairList.insert(pairList.end(), i);
		}
	}

	for (i = 0; i < pairList.size(); i++){
		if(sys[pairList[i]].PID == "2112"){
			tempNCounter++;
		}
		else{
			tempPCounter++;
		}
	}

	if (((tempNCounter == 0) && (tempNCounter != 0)) || ((tempNCounter != 0) && (tempPCounter == 0))){
		for (i = 0; i < pairList.size(); i++){
			tempList.insert(tempList.begin(), pairList[i]);
			fragSet.insert(fragSet.end(), tempList);
			tempList.clear();
		}
	}

	else{
		fragSet.insert(fragSet.begin(), pairList);
	}

	pairList.clear();
	tempPCounter = 0;
	tempNCounter = 0;

	for (i = 1; i < (count - 2); i++){
		for (j = 0; j < fragSet.size(); j++){
			if (elChecker(fragSet[j], i)){
				tempCounter = 0;
				break;
			}
		}
		if (tempCounter){
			pairList.insert(pairList.begin(), i);
			for (j = i + 1; j < count; j++){

				for (int k = 0; k < fragSet.size(); k++){
					if (elChecker(fragSet[k], j)){
						tempCounter2 = 0;
						break;
					}
				}


				if (tempCounter2){
					if (R(sys[i], sys[j]) && P(sys[i], sys[j])){
						pairList.insert(pairList.end(), j);
					}
				}

				tempCounter2 = 1;
			}

			for (i = 0; i < pairList.size(); i++){
				if(sys[pairList[i]].PID == "2112"){
					tempNCounter++;
				}
				else{
					tempPCounter++;
				}
			}

			if (((tempNCounter == 0) && (tempNCounter != 0)) || ((tempNCounter != 0) && (tempPCounter == 0))){
				for (i = 0; i < pairList.size(); i++){
					tempList.insert(tempList.begin(), pairList[i]);
					fragSet.insert(fragSet.end(), tempList);
					tempList.clear();
				}
			}
			else{
				fragSet.insert(fragSet.end(), pairList);
			}

			pairList.clear();
			tempPCounter = 0;
			tempNCounter = 0;
		}
		tempCounter = 1;
	}

	pairList.clear();
	tempCounter = 1;
	tempCounter2 = 1;

	for (i = 0; i < fragSet.size(); i++){
		if (elChecker(fragSet[i], count - 2)){
			tempCounter = 0;
			break;
		}
		if (elChecker(fragSet[i], count - 1)){
			tempCounter2 = 0;
			break;
		}
	}


	if (tempCounter && tempCounter2){
		if (R(sys[count - 2], sys[count - 1]) && P(sys[count - 2], sys[count - 1])){
			pairList.insert(pairList.begin(), count - 2);
			pairList.insert(pairList.end(), count - 1);

			for (i = 0; i < pairList.size(); i++){
				if(sys[pairList[i]].PID == "2112"){
					tempNCounter++;
				}
				else{
					tempPCounter++;
				}
			}

			if (((tempNCounter == 0) && (tempNCounter != 0)) || ((tempNCounter != 0) && (tempPCounter == 0))){
				for (i = 0; i < pairList.size(); i++){
					tempList.insert(tempList.begin(), pairList[i]);
					fragSet.insert(fragSet.end(), tempList);
					tempList.clear();
				}
			}
			else{
				fragSet.insert(fragSet.end(), pairList);
			}

			tempNCounter = 0;
			tempPCounter = 0;
			pairList.clear();
		}

		else{
			pairList.insert(pairList.begin(), count - 2);
			fragSet.insert(fragSet.end(), pairList);
			pairList.clear();

			pairList.insert(pairList.end(), count - 1);
			fragSet.insert(fragSet.end(), pairList);
			pairList.clear();
		}
	}

	else if (tempCounter){
		pairList.insert(pairList.begin(), count - 2);
		fragSet.insert(fragSet.end(), pairList);
		pairList.clear();
	}

	else if (tempCounter2){
		pairList.insert(pairList.end(), count - 1);
		fragSet.insert(fragSet.end(), pairList);
		pairList.clear();
	}

	int setSize = fragSet.size();
	vector<Particle> fragList(setSize);

	int nCounter = 0;
	int pCounter = 0;
	int Z = pCounter;
	int A = pCounter + nCounter;
	string pidVal = "1000";

	for (i = 0; i < setSize; i++){
		if (fragSet[i].size() == 1){
			fragList[i].PID = sys[fragSet[i][0]].PID;
			fragList[i].mass = sys[fragSet[i][0]].mass;
		}
		else{
			nCounter = 0;
			pCounter = 0;
			for (j = 0; j < fragSet[i].size(); j++){
				if (sys[fragSet[i][j]].PID == "2112"){
					nCounter++;
				}
				else{
					pCounter++;
				}
			}

			Z = pCounter;
			A = nCounter + pCounter;
			
			pidVal = "1000";

			if (Z < 10){
				pidVal += "0";
				pidVal += to_string(Z);
			}

			else{
				pidVal += to_string(Z);
			}

			if (A < 10){
				pidVal += "00";
				pidVal += to_string(A);
				pidVal += "0";
			}

			else{
				pidVal += "0";
				pidVal += to_string(A);
				pidVal += "0";
			}

			fragList[i].PID = pidVal;

		}
	}

	// averaged positions and momenta
	double totMass;
	double invRelMass;
	double xSum;
	double ySum;
	double zSum;
	double pXSum;
	double pYSum;
	double pZSum;
	int fragIndex;

	for (i = 0; i < setSize; i++){
		totMass = 0.0;
		invRelMass = 0.0;
		xSum = 0.0;
		ySum = 0.0;
		zSum = 0.0;
		pXSum = 0.0;
		pYSum = 0.0;
		pZSum = 0.0;

		for (j = 0; j < fragSet[i].size(); j++){
			fragIndex = fragSet[i][j];
			totMass += sys[fragIndex].mass;
			invRelMass += pow(sys[fragIndex].mass, -1);
			xSum += (sys[fragIndex].mass * sys[fragIndex].xCord);
			ySum += (sys[fragIndex].mass * sys[fragIndex].yCord);
			zSum += (sys[fragIndex].mass * sys[fragIndex].zCord);

			pXSum += sys[fragIndex].xMom;
			pYSum += sys[fragIndex].yMom;
			pZSum += sys[fragIndex].zMom;
		}

		fragList[i].xCord = xSum / totMass;
		fragList[i].yCord = ySum / totMass;
		fragList[i].zCord = zSum / totMass;

		fragList[i].xMom = pow(invRelMass, -1) * pXSum / totMass;
		fragList[i].yMom = pow(invRelMass, -1) * pYSum / totMass;
		fragList[i].zMom = pow(invRelMass, -1) * pZSum / totMass;
	}

	for (i = 0; i < fragList.size(); i++){
		cout << "Fragment " << i + 1 << ": " << fragList[i].PID << endl;
		cout << "\tPosition: (" << fragList[i].xCord << ", " << fragList[i].yCord << ", " << fragList[i].zCord << ")" << endl;
		cout << "\tMomentum: (" << fragList[i].xMom << ", " << fragList[i].yMom << ", " << fragList[i].zMom << ")" << endl;
		cout << endl;
	}


	return 0;
}


double* particleSpawn(string PID, double avgR, double avgP){

	static double values[6];

	double width = 1.7;// fm

	if (PID == "2112")
		width = 3.32236936788 - 0.8;
	else
		width = 3.19848732601 - 0.8;

	double hBar = 1.05471817e-34; //Js
	double pi = 3.141592653;
	double sigmaR = width; // fm
	double sigmaP = hBar/(width * 1e-15); //Ns

	double xCord = 0.0;
	double yCord = 0.0;
	double zCord = 0.0;

	double xMom = 0.0;
	double yMom = 0.0;
	double zMom = 0.0;

	double sampR = 0.0;
	double sampP = 0.0;

	double eta1 = 0.0;
	double eta2 = 0.0;

	// init sampR and sampP
	
	eta1 = ((double) rand()) / (RAND_MAX);
	eta2 = ((double) rand()) / (RAND_MAX);

//	sampR = sqrt(-2 * log(eta2)) * cos(pow((pi * hBar),3) * eta1) * sigmaR + avgR;
//	sampP = (sqrt(-2 * log(eta2)) * sin(pow((pi * hBar),3) * eta1) * sigmaP) / 5.34428565e-22; // MeV/c
	sampR = sqrt(-2 * log(eta2)) * cos(2 * pi * eta1) * sigmaR + avgR;
	sampP = (sqrt(-2 * log(eta2)) * sin(2 * pi * eta1) * sigmaP) / 5.34428565e-22; // MeV/c

	// init coords
	
	zCord = sampR * ((double) rand()) / (RAND_MAX);
	yCord = (sqrt(pow(sampR,2) - pow(zCord,2))) * ((double) rand()) / (RAND_MAX);
	xCord = sqrt(pow(sampR,2) - pow(zCord,2) - pow(yCord,2));

	// init momenta
	
	//xMom = sampP * ((double) rand()) / (RAND_MAX);
	//yMom = (sqrt(pow(sampP,2) - pow(xMom,2))) * ((double) rand()) / (RAND_MAX);
	//zMom = sqrt(pow(sampP,2) - pow(xMom,2) - pow(yMom,2));
	
	zMom = (2 * sampP) * ((double) rand()) / (RAND_MAX) + (avgP - sampP);
	double diff = zMom - avgP;
	yMom = sqrt(pow(sampP, 2) - pow(diff, 2)) * ((double) rand()) / (RAND_MAX);
	zMom = sqrt(pow(sampP, 2) - pow(diff, 2) - pow(yMom, 2));
		

	values[0] = xCord;
	values[1] = yCord;
	values[2] = zCord;
	values[3] = xMom;
	values[4] = yMom;
	values[5] = zMom;

	return values;

}

double crossSection(double energy){

	double cs = (1 + (5 / energy)) * (40 + 109 * cos(0.199 * sqrt(energy)) * exp(-0.451 * pow((energy - 25), 0.258)));

	return cs;

}

int finder(int el, vector<int> listSet){

	int temp = 0;

	for(int i = 0; i <= listSet.size(); i++){

		if (el == listSet[i]){
			temp++;
		}
	}

	if (temp > 0){
		return 1;
	}

	else{
		return 0;
	}
}

double vel(double p, double m){ // vel in units of c

	double v = p / sqrt(pow(m,2) + pow(p,2));
	return v;

}

void coll(Particle &nuc1, Particle &nuc2){
	double mEff = 1 / ((1 / nuc1.mass) + (1 / nuc2.mass));

	double r1 = sqrt(pow(nuc1.xCord, 2) + pow(nuc1.yCord, 2) + pow(nuc1.zCord, 2));
	double r2 = sqrt(pow(nuc2.xCord, 2) + pow(nuc2.yCord, 2) + pow(nuc2.zCord, 2));
	double denom = r2 - r1;

	vector<double> n{0.0, 0.0, 0.0};
	n[0] = (nuc2.xCord - nuc1.xCord) / denom;
	n[1] = (nuc2.yCord - nuc1.yCord) / denom;
	n[2] = (nuc2.zCord - nuc1.zCord) / denom;

	vector<double> v1{vel(nuc1.xMom, nuc1.mass), vel(nuc1.yMom, nuc1.mass), vel(nuc1.zMom, nuc1.mass)};
	vector<double> v2{vel(nuc2.xMom, nuc2.mass), vel(nuc2.yMom, nuc2.mass), vel(nuc2.zMom, nuc2.mass)};

	double u12 = 0.0;

	for(int i = 0; i < 3; i++){
		u12 += n[i] * (v1[i] - v2[i]);
	}

	double restitution = 1.0;

	double J = (1.0 + restitution) * mEff * u12;

	for(int i = 0; i < 3; i++){
		v1[i] -= (J / mEff * n[i]);
		v2[i] += (J / mEff * n[i]);
		if (v1[i] >= 1.0 || v1[i] <= -1.0){
			v1[i] = 0.8 * ((double) rand()) / (RAND_MAX);
		}
		if (v2[i] >= 1.0 || v2[i] <= -1.0){
			v2[i] = 0.8 * ((double) rand()) / (RAND_MAX);
		}
	}

	nuc1.xMom = nuc1.mass * pow(v1[0], 2) / sqrt(1 - pow(v1[0], 2));
	nuc1.yMom = nuc1.mass * pow(v1[1], 2) / sqrt(1 - pow(v1[1], 2));
	nuc1.zMom = nuc1.mass * pow(v1[2], 2) / sqrt(1 - pow(v1[2], 2));

	nuc2.xMom = nuc2.mass * pow(v2[0], 2) / sqrt(1 - pow(v2[0], 2));
	nuc2.yMom = nuc2.mass * pow(v2[1], 2) / sqrt(1 - pow(v2[1], 2));
	nuc2.zMom = nuc2.mass * pow(v2[2], 2) / sqrt(1 - pow(v2[2], 2));
}

void dispLog(int counter, Particle sys[], int time){

	for (int i = 0; i < counter; i++){

		//cout << "\n Particles at time " << time + 1 << endl;

		cout << i + 1 << ", " << sys[i].xCord << ", " << sys[i].yCord << ", " << sys[i].zCord << ", " << sys[i].xMom << ", " << sys[i].yMom << ", " << sys[i].zMom << endl;

	}


}


double distance(Particle p1, Particle p2){

	double r = sqrt(pow((p2.xCord - p1.xCord), 2) + pow((p2.yCord - p1.yCord), 2) + pow((p2.zCord - p1.zCord), 2));

	return r;
}

double dHalf(Particle p1, Particle p2){

	double mEff = 1 / ((1 / p1.mass) + (1 / p2.mass));

	vector<double> v1{vel(p1.xMom, p1.mass), vel(p1.yMom, p1.mass), vel(p1.zMom, p1.mass)};

	vector<double> v2{vel(p2.xMom, p2.mass), vel(p2.yMom, p2.mass), vel(p2.zMom, p2.mass)};

	vector<double> relV{(v2[0] - v1[0]), (v2[1] - v1[1]), (v2[2] - v1[2])};

	double vScal = sqrt(pow(relV[0], 2) + pow(relV[1], 2) + pow(relV[2], 2));

	double T = 0.5 * mEff * pow(vScal, 2);

	double p = sqrt(T * (T + 2 * mEff));

	double rootS = sqrt(2 * mEff * (mEff + sqrt(pow(p, 2) + pow(mEff, 2))));

	return (sqrt(crossSection(rootS) / 31.41592653));

}

double dubU(Particle p1, Particle p2){

	double sigmaR = 1.7;

	double PI = 3.141592653;

	vector<double> rio{p1.xCord, p1.yCord, p1.zCord};
	vector<double> rjo{p2.xCord, p2.yCord, p2.zCord};

	double square = -2 * (rio[0] * rjo[0] + rio[1] * rjo[1] + rio[2] * rjo[2]);

	for (int i = 0; i <= rio.size(); i++){
		square += pow(rio[i], 2);
		square += pow(rjo[i], 2);
	}

	double t1 = -2.9273566254e+21;

	double value = t1 / (pow((4 * PI * pow(sigmaR, 2)), 1.5)) * exp(-1 / (4 * pow(sigmaR, 2)) * square) * 6.24e-18;

	return value; // MeV
}

double tripU(Particle p1, Particle p2, Particle p3){

	double sigmaR = 1.7;

	double PI = 3.141592653;

	vector<double> rio{p1.xCord, p1.yCord, p1.zCord};
	vector<double> rjo{p2.xCord, p2.yCord, p2.zCord};
	vector<double> rko{p3.xCord, p3.yCord, p3.zCord};

	double squareLeft = -2 * (rio[0] * rjo[0] + rio[1] * rjo[1] + rio[2] * rjo[2]);

	for (int i = 0; i <= rio.size(); i++){
		squareLeft += pow(rio[i], 2);
		squareLeft += pow(rjo[i], 2);
	}

	double squareRight = -2 * (rio[0] * rko[0] + rio[1] * rko[1] + rio[2] * rko[2]);

	for (int i = 0; i <= rio.size(); i++){
		squareRight += pow(rio[i], 2);
		squareRight += pow(rko[i], 2);
	}

	double denom = pow((2 * PI * pow(sigmaR, 2)), 1.5) * pow(3, 1.5);

	double t2 = -3.6591957818e+21;

	double value = t2 / denom * exp(-1 / (4 * pow(sigmaR, 2)) * (squareLeft + squareRight)) * 6.24e-18;

	return value; // MeV
}

double U(int count, Particle set[]){

	double totalU = 0.0;
	double totalUij = 0.0;
	double totalUijk = 0.0;

	for (int i = 0; i < count; i++){
		for (int j = i + 1; j < count; j++){
			totalUij += dubU(set[i], set[j]);
			for (int k = j + 1; k < count; k++){
				totalUijk += tripU(set[i], set[j], set[k]);
			}
		}
	}

	totalU = totalUij + totalUijk;

	return totalU;
}

void readout(Particle p){
	cout << "PID: " << p.PID << endl;
	cout << "Mass: " << p.mass << " MeV/c^2" << endl;
	cout << "Position: (" << p.xCord << ", " << p.yCord << ", " << p.zCord << ") fm" << endl;
	cout << "Momentum: (" << p.xMom << ", " << p.yMom << ", " << p.zMom << ") MeV/c" << endl;
}

int elChecker(vector<int> v, int x){
	if(find(v.begin(), v.end(), x) != v.end()){
		return 1;
	}
	else{
		return 0;
	}
}

int R(Particle p1, Particle p2){
	double r0 = 3.5; // fm
	double xDiff = p1.xCord - p2.xCord;
	double yDiff = p1.yCord - p2.yCord;
	double zDiff = p1.zCord - p2.zCord;

	double r = sqrt(pow(xDiff, 2) + pow(yDiff, 2) + pow(zDiff, 2));

	if (r < r0){
		return 1;
	}
	else{
		return 0;
	}
}

int P(Particle p1, Particle p2){
	double p0 = 250; // MeV/c
	double xDiff = p1.xMom - p2.xMom;
	double yDiff = p1.yMom - p2.yMom;
	double zDiff = p1.zMom - p2.zMom;

	double p = sqrt(pow(xDiff, 2) + pow(yDiff, 2) + pow(zDiff, 2));

	if (p < p0){
		return 1;
	}
	else{
		return 0;
	}
}
