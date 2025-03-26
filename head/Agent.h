#pragma once

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}


#include <stdio.h>
#include <stdlib.h>
#include "MatrixCPU.h"

class Agent
{
	int _id;
	float _pLim[2];
	float _cost[2];
	int _nVoisin;
	MatrixCPU* _omega = nullptr;
	float _Lb;
	float _Ub;
	MatrixCPU* _x0 = nullptr;
	int _nAgent;
	int _type;



public:
	static const int AGENT_CONSUMER = 1;
	static const int AGENT_GENERATOR = 2;
	static const int AGENT_PROSUMER = 3;

	Agent();
	Agent(int id, float pLim1, float pLim2, float cost1, float cost2, int nVoisin, MatrixCPU* connec, int nAgent,int type);
	Agent(const Agent& a);
	Agent& operator= (const Agent& a);

	void setAgent(int id, float pLim1, float pLim2, float cost1, float cost2, int nVoisin, MatrixCPU* connec, int nAgent,int type);
	void updateP0(float P0, float dP);
	void generateOmega(MatrixCPU* omega, MatrixCPU* connect, int id, int nAgent, int nVoisin);

	MatrixCPU getVoisin();
	int getId();
	float getPmin();
	float getPmax();
	float getA();
	float getB();
	int getNVoisin();
	float getLb();
	float getUb();
	int getType();


	void display() const;

	~Agent();
};

