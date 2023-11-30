#ifndef NUCLEUSCOLL_H
#define NUCLEUSCOLL_H
#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include "particle.h"
using namespace std;

vector<Particle> nuColl(Particle projectile, string targetPID);

double* particleSpawn(Particle nucleus, string PID, double avgR, double avgP);

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

int elChecker(vector<int> v, int x);

int R(Particle p1, Particle p2);

int P(Particle p1, Particle p2);

double framePoly(double x, double c[5]);

double Rp(Particle p);

double Rn(Particle p);



#endif
