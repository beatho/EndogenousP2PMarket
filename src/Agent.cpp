#include "../head/Agent.h"


Agent::Agent()
{
	_id = 0;
	_pLim[0] = 0;
	_pLim[1] = 0;
	_cost[0] = 0;
	_cost[1] = 0;
	_nVoisin = 0;
	_nAgent = 0;
	_omega = new MatrixCPU(1, 1);
	_x0 = new MatrixCPU(1, 1);
	_type = 0; 
	_Ub = 0;
	_Lb = 0;
}

Agent::Agent(int id, float pLim1, float pLim2, float cost1, float cost2, int nVoisin, MatrixCPU* connec, int nAgent,int type)
{
	_id = id;
	_pLim[0] = pLim1;
	_pLim[1] = pLim2;
	_cost[0] = cost1;
	_cost[1] = cost2;
	_nVoisin = nVoisin;
	_nAgent = nAgent;
	
	if (nVoisin > 0) {
		_omega = new MatrixCPU(nVoisin, 1);
		generateOmega(_omega, connec, id, nAgent, nVoisin);
		_x0 = new MatrixCPU(1, nVoisin);
	}
	
	_type = type;
	switch (type) {
	case AGENT_CONSUMER:
		_Ub = 0;
		_Lb = pLim1;
		break;
	case AGENT_GENERATOR:
		_Ub = pLim2;
		_Lb = 0;
		break;
	case AGENT_PROSUMER:
		_Ub = pLim2;
		_Lb = pLim1;
		break;
	default:
		_Ub = 0;
		_Lb = 0;
		return;
	}
}

Agent::Agent(const Agent& a)
{
	_id = a._id;
	_pLim[0] = a._pLim[0];
	_pLim[1] = a._pLim[1];
	_cost[0] = a._cost[0];
	_cost[1] = a._cost[1];
	_nVoisin = a._nVoisin;
	_nAgent = a._nAgent;
	if (_nVoisin > 0) {
		_omega = new MatrixCPU(*(a._omega));
		_x0 = new MatrixCPU(*(a._x0));
	}
	_type = a._type;
	_Ub = a._Ub;
	_Lb = a._Lb;
}


Agent& Agent::operator=(const Agent& a)
{
	_id = a._id;
	_pLim[0] = a._pLim[0];
	_pLim[1] = a._pLim[1];
	_cost[0] = a._cost[0];
	_cost[1] = a._cost[1];
	_nVoisin = a._nVoisin;
	_nAgent = a._nAgent;
	DELETEB(_omega);
	DELETEB(_x0);
	if (_nVoisin > 0) {
		_omega = new MatrixCPU(*(a._omega));
		_x0 = new MatrixCPU(*(a._x0));
	}
	_type = a._type;
	_Ub = a._Ub;
	_Lb = a._Lb;
	return *this;
}


void Agent::setAgent(int id, float pLim1, float pLim2, float cost1, float cost2, int nVoisin, MatrixCPU* connec, int nAgent,int type)
{
	_id = id;
	_pLim[0] = pLim1;
	_pLim[1] = pLim2;
	_cost[0] = cost1;
	_cost[1] = cost2;
	_nVoisin = nVoisin;
	_nAgent = nAgent;
	DELETEB(_omega);
	DELETEB(_x0);
	
	if (nVoisin > 0) {
		_omega = new MatrixCPU(nVoisin, 1);
		generateOmega(_omega, connec, id, nAgent, nVoisin);
		_x0 = new MatrixCPU(1, nVoisin);
	}
	
	
	
	_type = type;
	switch (type) {
	case AGENT_CONSUMER:
		_Ub = 0;
		_Lb = pLim1;
		break;
	case AGENT_GENERATOR:
		_Ub = pLim2;
		_Lb = 0;
		break;
	case AGENT_PROSUMER:
		_Ub = pLim2;
		_Lb = pLim1;
		break;
	default:
		return;
	}
}

void Agent::updateP0(float P0, float dP)
{

	_pLim[0] = -(1 + dP) * P0;
	_pLim[1] = -(1 - dP) * P0;
	_cost[0] = 1;
	_cost[1] = P0 * _cost[0];
	
	switch (_type) {
	case AGENT_CONSUMER:
		_Ub = 0;
		_Lb = _pLim[0];
		break;
	case AGENT_GENERATOR:
		_Ub = _pLim[1];
		_Lb = 0;
		break;
	case AGENT_PROSUMER:
		_Ub = _pLim[1];
		_Lb = _pLim[0];
		break;
	default:
		return;
	}
}

void Agent::generateOmega(MatrixCPU* omega, MatrixCPU* connect, int id, int nAgent, int nVoisin)
{
	int k = 0;
	for (int n = 0; n < nAgent; n++) 
	{
		if (connect->get(id, n) > 0) {
			omega->set(k, 0, n);
			k++;
		}
	}
}



MatrixCPU Agent::getVoisin()
{
	return *_omega;
}

int Agent::getId()
{
	return _id;
}

float Agent::getPmin()
{
	return _pLim[0];
}

float Agent::getPmax()
{
	return _pLim[1];
}

float Agent::getA()
{
	return _cost[0];
}

float Agent::getB()
{
	return _cost[1];
}

int Agent::getNVoisin()
{
	return _nVoisin;
}

float Agent::getLb()
{
	return _Lb;
}

float Agent::getUb()
{
	return _Ub;
}

int Agent::getType()
{
	return _type;
}

void Agent::display() const
{
	std::cout << "Agent " << _id << "of type " << _type << std::endl;
	std::cout << "Power : [ " << _pLim[0] << " , " << _pLim[1] << " ]" << std::endl;
	std::cout << "Cost function : [ " << _cost[0] << " , " << _cost[1] << " ]" << std::endl;

	if (_nVoisin > 0) {
		std::cout << "Peers : " << std::endl;
		_omega->display();
	}
	else {
		std::cout << "No peers" << std::endl;
	}
	printf("end agent \n"); // seem to have at least one printf somewhere to use OSQP, so it is here...
}


Agent::~Agent() {
#ifdef DEBUG_DESTRUCTOR
	std::cout << "destruction agent" << std::endl;
#endif 
	DELETEB(_omega);
	DELETEB(_x0);
}