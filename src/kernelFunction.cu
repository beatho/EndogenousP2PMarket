
#include "../head/kernelFunction.cuh"


// OPF ADMM

__global__ void updateQ(float* Q, float* X, float* MU, float _rho, int sizeOPF) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;


	for (int j = index; j < sizeOPF; j += step) {
		Q[j] = -(MU[j] + _rho * X[j]); //Q[i].set(j, 0, -(Mu[i].get(j, 0) + _rho * X[i].get(j, 0)));
	}
}


__global__ void updateMUGPU(float* Mu, float* Y, float* X, float rho, int sizeOPF) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;


	for (int j = index; j < sizeOPF; j += step) {
		Mu[j] = Mu[j] + rho * (X[j] - Y[j]);
	}
}

__global__ void removeLossAgent(float* _nAgentByBus, float* CoresAgentBusBegin) {

	int thIdx = threadIdx.x;
	if (thIdx == 0) {
		_nAgentByBus[0] = _nAgentByBus[0] - 1;
		CoresAgentBusBegin[0] = 1;
	}
}


__global__ void initVoltageBound(float* VlimReal, float* Vlim, float* constraintLo, float* constraintUp, float* nChild, int nBus) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int b = index; b < nBus; b += step) {
		float nb = sqrtf((nChild[b] + 1) / 2);
		float ub = constraintUp[b + nBus];
		float lb = constraintLo[b + nBus];

		VlimReal[b] = lb;
		VlimReal[b + nBus] = ub;

		Vlim[b] = lb * lb * nb;
		Vlim[b + nBus] = ub * ub * nb;
	}
}
__global__ void divideMultiplyByNagentByBus(float* Apt1, float* Apt2, float* PnTilde, float* PnTmin, float* PnTmax, float* nAgentByBus, float rhol, int nBus) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int b = index; b < nBus; b += step) {
		int nb = nAgentByBus[b];
		Apt1[b] = nb * rhol;
		Apt1[b + nBus] = nb * rhol;

		Apt2[b] = nb * nb * rhol;
		Apt2[b + nBus] = nb * nb * rhol;

		if (nb > 0) {
			PnTilde[b] = PnTilde[b] / nb;
			PnTilde[b + nBus] = PnTilde[b + nBus] / nb;

			PnTmin[b] = PnTmin[b] / nb;
			PnTmin[b + nBus] = PnTmin[b + nBus] / nb;

			PnTmax[b] = PnTmax[b] / nb;
			PnTmax[b + nBus] = PnTmax[b + nBus] / nb;
		}
	}
}


__global__ void initDFSPQ(float* X, float* Pb, float* nChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, int nBus) {


	extern __shared__ int globalMemory[];
	int* ChildsSh = globalMemory;
	bool* hasfinished = (bool*)(&(ChildsSh[nBus]));
	/*__shared__ int ChildsSh[2];
	__shared__ bool hasfinished[2];*/


	__shared__ bool notfinished;

	int bus = threadIdx.x; // one block and _blocksize = nBus !!!!

	if (bus == 0) {
		notfinished = true;
	}
	__syncthreads();

	if (bus < nBus) {
		hasfinished[bus] = false;
		int indice = indiceBusBegin[bus];
		int indiceChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
		int nb = nChild[bus];
		bool mustCompute = (nb == 0);
		ChildsSh[bus] = (bus < (nBus - 1)) ? Childs[bus] : 0;
		__syncthreads();
		while (notfinished) {
			__syncthreads();
			if (mustCompute) { // divergent mais on n'y peut rien
				float p = Pb[bus];
				float q = Pb[bus + nBus];
				for (int i = 0; i < nb; i++) {
					int c = ChildsSh[indiceChild + i];
					int indiceBusChild = indiceBusBegin[c];
					p += X[indiceBusChild];
					q += X[indiceBusChild + 1];
				}
				X[indice] = (bus > 0) * p;
				X[indice + 1] = (bus > 0) * q;
				float Si = p * p + q * q;
				X[indice + 2] = (bus > 0) * Si / X[indice + 3];
				hasfinished[bus] = true;
				if (bus == 0) {
					notfinished = false;
				}
				//notfinished = (bus != 0); // tous essaie d'ecrire la m�me chose sauf quand 0 sera tout seul �tant le premier noeud, le seul sans anc�tre
			}
			__syncthreads();
			// trouver qui doit tourner � la prochaine boucle
			mustCompute = !(hasfinished[bus]);

			for (int i = 0; i < nb; i++) {
				int c = ChildsSh[indiceChild + i];

				mustCompute = (mustCompute && hasfinished[c]); // il suffit qu'un enfant n'a pas fini pour que cela soit false
			}
			__syncthreads();
		}
	}

}
__global__ void initDFSPQ(float* X, float* nChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, int nBus) {


	extern __shared__ int globalMemory[];
	int* ChildsSh = globalMemory;
	bool* hasfinished = (bool*)(&(ChildsSh[nBus]));
	/*__shared__ int ChildsSh[2];
	__shared__ bool hasfinished[2];*/


	__shared__ bool notfinished;

	int bus = threadIdx.x; // one block and _blocksize = nBus !!!!

	if (bus == 0) {
		notfinished = true;
	}
	__syncthreads();

	if (bus < nBus) {
		hasfinished[bus] = false;
		int indice = indiceBusBegin[bus];
		int indiceChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
		int nb = nChild[bus];
		bool mustCompute = (nb == 0);
		ChildsSh[bus] = (bus < (nBus - 1)) ? Childs[bus] : 0;
		while (notfinished) {
			__syncthreads();
			if (mustCompute) { // divergent mais on n'y peut rien
				float p = X[indice + 4];
				float q = X[indice + 5];
				for (int i = 0; i < nb; i++) {
					int c = ChildsSh[indiceChild + i];
					int indiceBusChild = indiceBusBegin[c];
					p += X[indiceBusChild];
					q += X[indiceBusChild + 1];
				}
				X[indice] = (bus > 0) * p;
				X[indice + 1] = (bus > 0) * q;
				float Si = p * p + q * q;
				X[indice + 2] = (bus > 0) * Si / X[indice + 3];
				hasfinished[bus] = true;
				if (bus == 0) {
					notfinished = false;
				}
				//notfinished = (bus != 0); // tous essaie d'ecrire la m�me chose sauf quand 0 sera tout seul �tant le premier noeud, le seul sans anc�tre
			}
			__syncthreads();
			// trouver qui doit tourner � la prochaine boucle
			mustCompute = !(hasfinished[bus]);

			for (int i = 0; i < nb; i++) {
				int c = ChildsSh[indiceChild + i];
				mustCompute = (mustCompute && hasfinished[c]); // il suffit qu'un enfant n'a pas fini pour que cela soit false
			}
			__syncthreads();
		}
	}

}


__global__ void initPosAgent(float* PosAgent, float* nAgentByBus, float* CoresAgentBusBegin, float* CoresAgentBus) {

	int bus = blockIdx.x; // un bloc par bus
	int nB = nAgentByBus[bus];
	int index = threadIdx.x;
	int step = blockDim.x;
	int begin = CoresAgentBusBegin[bus];
	for (int i = index; i < nB; i += step) {
		int n = CoresAgentBus[begin + i];
		PosAgent[n] = i; // pas du tout coalescent
	}


}


__global__ void initPQAgentV(float* X, float* indiceBusBegin, float* CoresAgentBus, float* nAgentByBus, float* beginBusAgent, float* Pn, int nAgent) {
	int bus = blockIdx.x;
	int thI = threadIdx.x;
	int step = blockDim.x;
	int begin = beginBusAgent[bus];
	int nB = nAgentByBus[bus];
	int fin = begin + nB;
	int indiceBus = indiceBusBegin[bus];

	if (thI == 0) {
		X[indiceBus + 3] = 1; // vi
		X[indiceBus + 4 + 2 * nB] = 1; // vai
	}

	for (int i = begin + thI; i < fin; i += step) { // ecriture coalecente mais pas lecture
		int agent = CoresAgentBus[i];

		X[indiceBus + 4 + thI] = Pn[agent]; 
		X[indiceBus + 4 + nB + thI] = Pn[agent + nAgent]; 

	}
}
__global__ void initPQAgent(float* X, float* indiceBusBegin, float* CoresAgentBus, float* nAgentByBus, float* beginBusAgent, float* Pn, int nAgent) {
	int bus = blockIdx.x;
	int thI = threadIdx.x;
	int step = blockDim.x;
	int begin = beginBusAgent[bus];
	int nB = nAgentByBus[bus];
	int fin = begin + nB;
	int indiceBus = indiceBusBegin[bus];


	for (int i = begin + thI; i < fin; i += step) { // ecriture coalecente mais pas lecture
		int agent = CoresAgentBus[i];

		X[indiceBus + 4 + thI] = Pn[agent]; // pi = sum(pn)
		X[indiceBus + 4 + nB + thI] = Pn[agent + nAgent]; // qi = sum(pn)

	}
}

__global__ void initPQV(float* X, float* indiceBusBegin, float* nAgentByBus, float* PnTilde, int nBus) {
	int bus = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = bus; i < nBus; i += step) { // pas coalescent
		int indiceBus = indiceBusBegin[i];
		X[indiceBus + 3] = 1; // vi

		// pi & qi
		int Nb = nAgentByBus[i];
		X[indiceBus + 4] = Nb * PnTilde[i]; // pi = sum(pn)
		X[indiceBus + 5] = Nb * PnTilde[i + nBus]; // qi = sum(pn)

	}
}
__global__ void initPQ(float* X, float* indiceBusBegin, float* nAgentByBus, float* PnTilde, int nBus) {
	int bus = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = bus; i < nBus; i += step) { // pas coalescent
		int indiceBus = indiceBusBegin[i];

		// pi & qi
		int Nb = nAgentByBus[i];
		X[indiceBus + 4] = Nb * PnTilde[i]; // pi = sum(pn)
		X[indiceBus + 5] = Nb * PnTilde[i + nBus]; // qi = sum(pn)

	}
}


__global__ void defineSizeBig(float* sizeOPFADMMbig, float* nChild, float* CoresBusBegin, float* sizeOPFADMM, float* CoresBusBeginBig, float* nAgentByBus) {

	int bus = blockIdx.x; // un bloc par bus
	int thIdx = threadIdx.x;
	int step = blockDim.x;

	int nC = nChild[bus];
	int sizeOPF = nC * 3 + 5 + 2 * nAgentByBus[bus];
	int debut = CoresBusBegin[bus];

	if (thIdx == 0) {
		sizeOPFADMM[bus] = sizeOPF;
	}

	for (int i = thIdx; i < sizeOPF; i += step) {
		sizeOPFADMMbig[i + debut] = sizeOPF;
		CoresBusBeginBig[i + debut] = debut;
	}

}

__global__ void defineSizeBig(float* sizeOPFADMMbig, float* nChild, float* CoresBusBegin, float* sizeOPFADMM, float* CoresBusBeginBig, float* nAgentByBus, int lossType, int nBus, int nAgent) {

	int bus = blockIdx.x; // un bloc par bus
	int thIdx = threadIdx.x;
	int step = blockDim.x;

	if (bus < nBus) {
		int nC = nChild[bus];
		int sizeOPF = nC * 3 + 5 + 2 * nAgentByBus[bus];
		int debut = CoresBusBegin[bus];

		if (thIdx == 0) {
			sizeOPFADMM[bus] = sizeOPF;
		}

		for (int i = thIdx; i < sizeOPF; i += step) {
			sizeOPFADMMbig[i + debut] = sizeOPF;
			CoresBusBeginBig[i + debut] = debut;
		}
	}
	else { // loss
		if (lossType) { // Current
			int sizeOPF = nBus + 2;
			int debut = CoresBusBegin[bus];
			if (thIdx == 0) {
				sizeOPFADMM[bus] = sizeOPF;
			}

			for (int i = thIdx; i < sizeOPF; i += step) {
				sizeOPFADMMbig[i + debut] = sizeOPF;
				CoresBusBeginBig[i + debut] = debut;
			}
		}
		else { // POWER
			int sizeOPF = 2 * nAgent;
			int debut = CoresBusBegin[bus];
			if (thIdx == 0) {
				sizeOPFADMM[bus] = sizeOPF;
			}

			for (int i = thIdx; i < sizeOPF; i += step) {
				sizeOPFADMMbig[i + debut] = sizeOPF;
				CoresBusBeginBig[i + debut] = debut;
			}
		}

	}
	

}
__global__ void defineSizeBig(float* sizeOPFADMMbig, float* nChild, float* CoresBusBegin, float* sizeOPFADMM, float* CoresBusBeginBig) {

	int bus = blockIdx.x; // un bloc par bus
	int thIdx = threadIdx.x;
	int step = blockDim.x;

	int nC = nChild[bus];
	int sizeOPF = nC * 3 + 7;
	int debut = CoresBusBegin[bus];

	if (thIdx == 0) {
		sizeOPFADMM[bus] = sizeOPF;
	}

	for (int i = thIdx; i < sizeOPF; i += step) {
		sizeOPFADMMbig[i + debut] = sizeOPF;
		CoresBusBeginBig[i + debut] = debut;
	}

}


__global__ void updateXOPFADMM(float* X, float* Chat, float* Vbound, float* nAgentByBus, float* nChild, float* indiceBusBegin, float* CoresChatBegin,
	float* CoresAgentBusBegin, float* CoresAgentBus, float* Cost1, float* Cost2, float* Pmin, float* Pmax, float rho, int nBus, int nAgent, bool Lagrange) {

	int bus = blockIdx.x;
	int index = threadIdx.x;
	int step = blockDim.x;
	int beginChat = CoresChatBegin[bus];
	double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, x3min, x3max, gamma, k2; // double peut �tre necessaire
	double c1122; // c3 : voltage -> indice + 3, c4 : current -> indice + 2;
	double coefPoly2[2];
	double root2[4];
	double root3[4];
	double root4[4];
	double coefPoly3[3];
	int typeSol = 0;
	int BestRoot = 0;
	double bestGamma = -1;
	double p = 0;

	int nRoot = 0;


	int begining = indiceBusBegin[bus];
	int nC = nChild[bus];

	bool goodSol = false;
	k2 = sqrt(2.0 / (nC + 1));
	if (index == 0)
	{
		/*if (bus == 0) {
			goodSol = true;
			x1 = 0;
			x2 = 0;
			x4 = 0;
			x3 = 1 / k2;
			gamma = 0;
		}
		else {*/

		c1 = -2 * Chat[beginChat];
		c2 = -2 * Chat[beginChat + 1];
		c4 = -2 * Chat[beginChat + 2];
		c3 = -2 * Chat[beginChat + 3] / k2;

		c1122 = c1 * c1 + c2 * c2;
		x3min = Vbound[bus];
		x3max = Vbound[bus + nBus];

		// case without constraint

		x1 = -c1 / 2;
		x2 = -c2 / 2;
		x3 = -c3 / 2;
		x4 = -c4 / 2;

		lambdaUp = 0;
		lambdaLo = 0;
		if (x3 < x3min) {
			x3 = x3min;
			lambdaLo = (2 * x3 + c3);
		}
		else if (x3 > x3max) {
			x3 = x3max;
			lambdaUp = -(2 * x3 + c3);
		}
		gamma = k2 * x4 - (x1 * x1 + x2 * x2) / x3; // ce n'est pas vraiment gamma, doit être positif

		if (gamma >= 0) {
			// the solution is good !
			goodSol = true;
		}
		else {
			if (c1122 == 0) {
				x4 = 0;
				goodSol = true;
			}
			if (gamma > bestGamma) {
				typeSol = 1;
				bestGamma = gamma;
			}
		}
		//}
		if (!goodSol) {
			x3 = x3max;

			coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
			coefPoly2[1] = 1 / x3;
			coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
			coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);
			nRoot = resolveRealPolynome3without2termGPU(root2, coefPoly2[0], coefPoly2[1]);

			for (int n = 0; n < nRoot; n++) {
				p = root2[n];

				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaUp = -(2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

				if (gamma >= 0 && lambdaUp >= 0) {
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && lambdaUp > bestGamma) {
					typeSol = 2;
					bestGamma = min(gamma, lambdaUp);
					BestRoot = n;
				}

			}
		}
			// case x3 = x3min lambdaUp = 0
		if (!goodSol) {
			x3 = x3min;

			coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
			coefPoly2[1] = 1 / x3;
			coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
			coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);
			nRoot = resolveRealPolynome3without2termGPU(root3, coefPoly2[0], coefPoly2[1]);

			for (int n = 0; n < nRoot; n++) {
				p = root3[n];
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaLo = (2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

				if (gamma >= 0 && lambdaLo >= 0) {
					// the solution is good !
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && lambdaLo > bestGamma) {
					typeSol = 3;
					bestGamma = min(gamma, lambdaLo);
					BestRoot = n;
				}
			}
		}
			// case xmin<x3<xmax lambdaLo = 0 lambdaUp = 0
		if (!goodSol) {

			coefPoly3[0] = c1122 / k2 * (2 * c3 / k2 - c4);
			coefPoly3[1] = (c3 - 2 * c4 / k2);
			coefPoly3[2] = -1;
			coefPoly3[0] = coefPoly3[0] * k2 * k2 / (c1122 * c1122);
			coefPoly3[1] = coefPoly3[1] * k2 * k2 / (c1122 * c1122);
			coefPoly3[2] = coefPoly3[2] * k2 * k2 / (c1122 * c1122);

			nRoot = resvolveRealPolynome4without2termGPU(root4, coefPoly3[0], coefPoly3[1], coefPoly3[2], Lagrange);

			for (int n = 0; n < nRoot; n++) {
				p = root4[n];
				x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;

				if (gamma >= 0 && x3 <= x3max && x3 >= x3min) {
					// the solution is good !
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma) {
					typeSol = 4;
					bestGamma = min(min(gamma, (x3max - x3)), (x3 - x3min));
					BestRoot = n;
				}
			}
		}
		if (!goodSol) {

			if (typeSol == 1) {
				// case without constraint
				x1 = -c1 / 2;
				x2 = -c2 / 2;
				x3 = -c3 / 2;
				x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
				x4 = -c4 / 2; // ou  (x1 * x1 + x2 * x2) / (k2 * x3)
			}
			else {
				if (typeSol == 2) {
					x3 = x3max;
					p = root2[BestRoot];
				}
				else if (typeSol == 3) {
					x3 = x3min;
					p = root3[BestRoot];
				}
				else if (typeSol == 4) {
					p = root4[BestRoot];
					x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
					x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
				}
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
			}
		}

		X[begining] = x1;
		X[begining + 1] = x2;
		X[begining + 2] = x4;
		X[begining + 3] = x3 * k2;

	}
	int nb = nAgentByBus[bus];
	int beginAgent = CoresAgentBusBegin[bus];
	for (int i = index; i < nb; i += step) {
		int n = CoresAgentBus[i + beginAgent];
		float ub = Pmax[n];
		float lb = Pmin[n];
		float pn = (rho * Chat[beginChat + 4 + i] - Cost2[n]) / (Cost1[n] + rho);
		pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;


		ub = Pmax[n + nAgent];
		lb = Pmin[n + nAgent];
		float qn = (rho * Chat[beginChat + 4 + nb + i] - Cost2[n + nAgent]) / (Cost1[n + nAgent] + rho);
		qn = (ub - qn) * (qn > ub) + (lb - qn) * (qn < lb) + qn;

		// pn & qn
		X[begining + 4 + i] = pn;
		X[begining + 4 + nb + i] = qn;
	}


	// X =  {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}	

}

__global__ void updateXOPFADMM(float* X, float* Chat, float* Vbound, float* PnMoy, float* nAgentByBus, float* nChild, float* indiceBusBegin, int nBus, bool Lagrange) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, x3min, x3max, gamma, k2; // double peut �tre necessaire
	double c1122; // c3 : voltage -> indice + 3, c4 : current -> indice + 2;
	double coefPoly2[2];
	double root2[4];
	double root3[4];
	double root4[4];
	double coefPoly3[3];
	int typeSol = 0;
	int BestRoot = 0;
	double bestGamma = -1;
	double p = 0;

	int nRoot = 0;

	for (int b = index + 1; b < nBus; b += step) {
		int begining = indiceBusBegin[b];
		int nC = nChild[b];

		bool goodSol = false;
		k2 = sqrt(2.0 / (nC + 1));
		/*if (b == 0) { // slack bus divergent, mais on n'y peut rien
			goodSol = true;
			x1 = 0;
			x2 = 0;
			x4 = 0;
			x3 = 1 / k2;
			gamma = 0;
		}*/
		//else {

		c1 = -2 * Chat[b];
		c2 = -2 * Chat[b + nBus];
		c4 = -2 * Chat[b + 2 * nBus];
		c3 = -2 * Chat[b + 3 * nBus] / k2;

		c1122 = c1 * c1 + c2 * c2;
		x3min = Vbound[b];
		x3max = Vbound[b + nBus];

		// case without constraint

		x1 = -c1 / 2;
		x2 = -c2 / 2;
		x3 = -c3 / 2;
		x4 = -c4 / 2;

		lambdaUp = 0;
		lambdaLo = 0;
		if (x3 < x3min) {
			x3 = x3min;
			lambdaLo = (2 * x3 + c3);
		}
		else if (x3 > x3max) {
			x3 = x3max;
			lambdaUp = -(2 * x3 + c3);
		}
		gamma = k2 * x4 - (x1 * x1 + x2 * x2) / x3; // ce n'est pas vraiment gamma, doit être positif

		if (gamma >= 0) {
			// the solution is good !
			goodSol = true;
		}
		else {
			if (c1122 == 0) {
				x4 = 0;
				goodSol = true;
			}
			if (gamma > bestGamma) {
				typeSol = 1;
				bestGamma = gamma;
			}
		}
		//}
		if (!goodSol) {
			x3 = x3max;

			coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
			coefPoly2[1] = 1 / x3;
			coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
			coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);
			nRoot = resolveRealPolynome3without2termGPU(root2, coefPoly2[0], coefPoly2[1]);

			for (int n = 0; n < nRoot; n++) {
				p = root2[n];

				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaUp = -(2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

				if (gamma >= 0 && lambdaUp >= 0) {
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && lambdaUp > bestGamma) {
					typeSol = 2;
					bestGamma = min(gamma, lambdaUp);
					BestRoot = n;
				}

			}
		}
		// case x3 = x3min lambdaUp = 0
		if (!goodSol) {
			x3 = x3min;

			coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
			coefPoly2[1] = 1 / x3;
			coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
			coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);
			nRoot = resolveRealPolynome3without2termGPU(root2, coefPoly2[0], coefPoly2[1]);

			for (int n = 0; n < nRoot; n++) {
				p = root3[n];
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaLo = (2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

				if (gamma >= 0 && lambdaLo >= 0) {
					// the solution is good !
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && lambdaLo > bestGamma) {
					typeSol = 3;
					bestGamma = min(gamma, lambdaLo);
					BestRoot = n;
				}
			}
		}
		// case xmin<x3<xmax lambdaLo = 0 lambdaUp = 0
		if (!goodSol) {

			coefPoly3[0] = c1122 / k2 * (2 * c3 / k2 - c4);
			coefPoly3[1] = (c3 - 2 * c4 / k2);
			coefPoly3[2] = -1;
			coefPoly3[0] = coefPoly3[0] * k2 * k2 / (c1122 * c1122);
			coefPoly3[1] = coefPoly3[1] * k2 * k2 / (c1122 * c1122);
			coefPoly3[2] = coefPoly3[2] * k2 * k2 / (c1122 * c1122);

			nRoot = resvolveRealPolynome4without2termGPU(root4, coefPoly3[0], coefPoly3[1], coefPoly3[2], Lagrange);

			for (int n = 0; n < nRoot; n++) {
				p = root4[n];
				x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;

				if (gamma >= 0 && x3 <= x3max && x3 >= x3min) {
					// the solution is good !
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma) {
					typeSol = 4;
					bestGamma = min(min(gamma, (x3max - x3)), (x3 - x3min));
					BestRoot = n;
				}
			}
		}


		if (!goodSol) {

			if (typeSol == 1) {
				// case without constraint
				x1 = -c1 / 2;
				x2 = -c2 / 2;
				x3 = -c3 / 2;
				x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
				x4 = -c4 / 2; // ou  (x1 * x1 + x2 * x2) / (k2 * x3)
			}
			else {
				if (typeSol == 2) {
					x3 = x3max;
					p = root2[BestRoot];
				}
				else if (typeSol == 3) {
					x3 = x3min;
					p = root3[BestRoot];
				}
				else if (typeSol == 4) {
					p = root4[BestRoot];
					x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
					x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
				}
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
			}
		}
		// X =  {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}	
		X[begining] = x1;
		X[begining + 1] = x2;
		X[begining + 2] = x4;
		X[begining + 3] = x3 * k2;

	}

	for (int b = index; b < nBus; b += step) {
		int begining = indiceBusBegin[b];
		X[begining + 4] = PnMoy[b] * nAgentByBus[b];
		X[begining + 5] = PnMoy[b + nBus] * nAgentByBus[b];
	}
}


__global__ void updateXEndoMarket(float* X, float* Chat, float* Vbound, float* nChild, float* CoresChatBegin, float* indiceBusBegin, int nBus) {
	
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	
	double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, x3min, x3max, gamma, k2; // double peut �tre necessaire
	double c1122; // c3 : voltage -> indice + 3, c4 : current -> indice + 2;
	double coefPoly2[2];
	double root2[4];
	double root3[4];
	double root4[4];
	double coefPoly3[3];
	int typeSol = 0;
	int BestRoot = 0;
	double bestGamma = -1;
	double p = 0;

	int nRoot = 0;


	bool goodSol = false;
	
	for (int bus = index + 1; bus < nBus; bus += step) {
		int beginChat = CoresChatBegin[bus];
		int begining = indiceBusBegin[bus];
		int nC = nChild[bus];
		k2 = sqrt(2.0 / (nC + 1));

		c1 = -2 * Chat[beginChat];
		c2 = -2 * Chat[beginChat + 1];
		c4 = -2 * Chat[beginChat + 2];
		c3 = -2 * Chat[beginChat + 3] / k2;

		c1122 = c1 * c1 + c2 * c2;
		x3min = Vbound[bus];
		x3max = Vbound[bus + nBus];

		// case without constraint

		x1 = -c1 / 2;
		x2 = -c2 / 2;
		x3 = -c3 / 2;
		x4 = -c4 / 2;

		lambdaUp = 0;
		lambdaLo = 0;
		if (x3 < x3min) {
			x3 = x3min;
			lambdaLo = (2 * x3 + c3);
		}
		else if (x3 > x3max) {
			x3 = x3max;
			lambdaUp = -(2 * x3 + c3);
		}

		gamma = k2 * x4 - (x1 * x1 + x2 * x2) / x3; // ce n'est pas vraiment gamma, doit être positif

		if (gamma >= 0) {
			// the solution is good !
			goodSol = true;
		}
		else {
			if (c1122 == 0) {
				x4 = 0;
				goodSol = true;
			}
			else if (gamma > bestGamma) {
				typeSol = 1;
				bestGamma = gamma;
			}
		}

		if (!goodSol) {
			x3 = x3max;

			coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
			coefPoly2[1] = 1 / x3;
			coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
			coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);
			nRoot = resolveRealPolynome3without2termGPU(root2, coefPoly2[0], coefPoly2[1]);

			for (int n = 0; n < nRoot; n++) {
				p = root2[n];

				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaUp = -(2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

				if (gamma >= 0 && lambdaUp >= 0) {
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && lambdaUp > bestGamma) {
					typeSol = 2;
					bestGamma = min(gamma, lambdaLo);
					BestRoot = n;
				}

			}
			// case x3 = x3min lambdaUp = 0
		}
		if (!goodSol) {
			x3 = x3min;

			coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
			coefPoly2[1] = 1 / x3;
			coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
			coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);
			nRoot = resolveRealPolynome3without2termGPU(root3, coefPoly2[0], coefPoly2[1]);

			for (int n = 0; n < nRoot; n++) {
				p = root3[n];
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaLo = (2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

				if (gamma >= 0 && lambdaLo >= 0) {
					// the solution is good !
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && lambdaLo > bestGamma) {
					typeSol = 3;
					bestGamma = min(gamma, lambdaLo);
					BestRoot = n;
				}
			}
		}
		// case xmin<x3<xmax lambdaLo = 0 lambdaUp = 0
		if (!goodSol) {
			coefPoly3[0] = c1122 / k2 * (2 * c3 / k2 - c4);
			coefPoly3[1] = (c3 - 2 * c4 / k2);
			coefPoly3[2] = -1;
			coefPoly3[0] = coefPoly3[0] * k2 * k2 / (c1122 * c1122);
			coefPoly3[1] = coefPoly3[1] * k2 * k2 / (c1122 * c1122);
			coefPoly3[2] = coefPoly3[2] * k2 * k2 / (c1122 * c1122);

			nRoot = resvolveRealPolynome4without2termGPU(root4, coefPoly3[0], coefPoly3[1], coefPoly3[2]);

			for (int n = 0; n < nRoot; n++) {
				p = root4[n];
				x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;

				if (gamma >= 0 && x3 <= x3max && x3 >= x3min) {
					// the solution is good !
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma) {
					typeSol = 4;
					bestGamma = min(min(gamma, (x3max - x3)), (x3 - x3min));
					BestRoot = n;
				}
			}
		}

		if (!goodSol) {

			if (typeSol == 1) {
				// case without constraint
				x1 = -c1 / 2;
				x2 = -c2 / 2;
				x3 = -c3 / 2;
				x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
				x4 = -c4 / 2; // ou  (x1 * x1 + x2 * x2) / (k2 * x3)
			}
			else {
				if (typeSol == 2) {
					x3 = x3max;
					p = root2[BestRoot];
				}
				else if (typeSol == 3) {
					x3 = x3min;
					p = root3[BestRoot];
				}
				else if (typeSol == 4) {
					p = root4[BestRoot];
					x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
					x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
				}
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
			}
		}

		X[begining] = x1;
		X[begining + 1] = x2;
		X[begining + 2] = x4;
		X[begining + 3] = x3 * k2;
		


	}
	

}

__global__ void updateXOPFADMMCons(float* X, float* Pn, float* Chat, float* Vbound, float* nAgentByBus, float* nChild, float* indiceBusBegin, float* CoresChatBegin,
	float* CoresAgentBusBegin, float* CoresAgentBus, float* Cost1, float* Cost2, float* Pmin, float* Pmax, float rho, int losstype, int nBus, int nAgent, bool Lagrange) {

	int bus = blockIdx.x + 1; // pas faire le bus 0
	int index = threadIdx.x;
	int step = blockDim.x;
	
	double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, x3min, x3max, gamma, k2; // double peut �tre necessaire
	double c1122; // c3 : voltage -> indice + 3, c4 : current -> indice + 2;
	double coefPoly2[2];
	double root2[4];
	double root3[4];
	double root4[4];
	double coefPoly3[3];
	int typeSol = 0;
	int BestRoot = 0;
	double bestGamma = -1;
	double p = 0;

	int nRoot = 0;


	int begining = indiceBusBegin[bus];
	int nC = nChild[bus];
	int beginChat = CoresChatBegin[bus];

	bool goodSol = false;
	k2 = sqrt(2.0 / (nC + 1));
	if (bus < nBus) {
		if (index == 0)
		{
			c1 = -2 * Chat[beginChat];
			c2 = -2 * Chat[beginChat + 1];
			c4 = -2 * Chat[beginChat + 2];
			c3 = -2 * Chat[beginChat + 3] / k2;

			c1122 = c1 * c1 + c2 * c2;
			x3min = Vbound[bus];
			x3max = Vbound[bus + nBus];

			// case without constraint

			x1 = -c1 / 2;
			x2 = -c2 / 2;
			x3 = -c3 / 2;
			x4 = -c4 / 2;

			lambdaUp = 0;
			lambdaLo = 0;
			if (x3 < x3min) {
				x3 = x3min;
				lambdaLo = (2 * x3 + c3);
			}
			else if (x3 > x3max) {
				x3 = x3max;
				lambdaUp = -(2 * x3 + c3);
			}
			gamma = k2 * x4 - (x1 * x1 + x2 * x2) / x3; // ce n'est pas vraiment gamma, doit être positif

			if (gamma >= 0) {
				// the solution is good !
				goodSol = true;
			}
			else {
				if (c1122 == 0) {
					x4 = 0;
					goodSol = true;
				}
				if (gamma > bestGamma) {
					typeSol = 1;
					bestGamma = gamma;
				}
			}
			//}
			if (!goodSol) {
				x3 = x3max;

				coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
				coefPoly2[1] = 1 / x3;
				coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
				coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);
				nRoot = resolveRealPolynome3without2termGPU(root2, coefPoly2[0], coefPoly2[1]);

				for (int n = 0; n < nRoot; n++) {
					p = root2[n];

					x1 = p * c1 * x3;
					x2 = p * c2 * x3;
					x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
					gamma = (2 * x4 + c4) / k2;
					lambdaUp = -(2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

					if (gamma >= 0 && lambdaUp >= 0) {
						goodSol = true;
						break;
					}
					if (gamma > bestGamma && lambdaUp > bestGamma) {
						typeSol = 2;
						bestGamma = min(gamma, lambdaUp);
						BestRoot = n;
					}

				}
			}
			// case x3 = x3min lambdaUp = 0
			if (!goodSol) {
				x3 = x3min;

				coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
				coefPoly2[1] = 1 / x3;
				coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
				coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);
				nRoot = resolveRealPolynome3without2termGPU(root3, coefPoly2[0], coefPoly2[1]);

				for (int n = 0; n < nRoot; n++) {
					p = root3[n];
					x1 = p * c1 * x3;
					x2 = p * c2 * x3;
					x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
					gamma = (2 * x4 + c4) / k2;
					lambdaLo = (2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

					if (gamma >= 0 && lambdaLo >= 0) {
						// the solution is good !
						goodSol = true;
						break;
					}
					if (gamma > bestGamma && lambdaLo > bestGamma) {
						typeSol = 3;
						bestGamma = min(gamma, lambdaLo);
						BestRoot = n;
					}
				}
			}
			// case xmin<x3<xmax lambdaLo = 0 lambdaUp = 0
			if (!goodSol) {

				coefPoly3[0] = c1122 / k2 * (2 * c3 / k2 - c4);
				coefPoly3[1] = (c3 - 2 * c4 / k2);
				coefPoly3[2] = -1;
				coefPoly3[0] = coefPoly3[0] * k2 * k2 / (c1122 * c1122);
				coefPoly3[1] = coefPoly3[1] * k2 * k2 / (c1122 * c1122);
				coefPoly3[2] = coefPoly3[2] * k2 * k2 / (c1122 * c1122);

				nRoot = resvolveRealPolynome4without2termGPU(root4, coefPoly3[0], coefPoly3[1], coefPoly3[2], Lagrange);

				for (int n = 0; n < nRoot; n++) {
					p = root4[n];
					x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
					x1 = p * c1 * x3;
					x2 = p * c2 * x3;
					x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
					gamma = (2 * x4 + c4) / k2;

					if (gamma >= 0 && x3 <= x3max && x3 >= x3min) {
						// the solution is good !
						goodSol = true;
						break;
					}
					if (gamma > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma) {
						typeSol = 4;
						bestGamma = min(min(gamma, (x3max - x3)), (x3 - x3min));
						BestRoot = n;
					}
				}
			}
			if (!goodSol) {

				if (typeSol == 1) {
					// case without constraint
					x1 = -c1 / 2;
					x2 = -c2 / 2;
					x3 = -c3 / 2;
					x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
					x4 = -c4 / 2; // ou  (x1 * x1 + x2 * x2) / (k2 * x3)
				}
				else {
					if (typeSol == 2) {
						x3 = x3max;
						p = root2[BestRoot];
					}
					else if (typeSol == 3) {
						x3 = x3min;
						p = root3[BestRoot];
					}
					else if (typeSol == 4) {
						p = root4[BestRoot];
						x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
						x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
					}
					x1 = p * c1 * x3;
					x2 = p * c2 * x3;
					x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				}
			}

			X[begining] = x1;
			X[begining + 1] = x2;
			X[begining + 2] = x4;
			X[begining + 3] = x3 * k2;

		}
		int nb = nAgentByBus[bus];
		int beginAgent = CoresAgentBusBegin[bus];
		for (int i = index; i < nb; i += step) {
			int n = CoresAgentBus[i + beginAgent];
			float ub = Pmax[n];
			float lb = Pmin[n];
			float pn = (rho * Chat[beginChat + 4 + i] - Cost2[n]) / (Cost1[n] + rho);
			pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;


			ub = Pmax[n + nAgent];
			lb = Pmin[n + nAgent];
			float qn = (rho * Chat[beginChat + 4 + nb + i] - Cost2[n + nAgent]) / (Cost1[n + nAgent] + rho);
			qn = (ub - qn) * (qn > ub) + (lb - qn) * (qn < lb) + qn;

			// pn & qn
			X[begining + 4 + i] = pn;
			X[begining + 4 + nb + i] = qn;
			Pn[n] = pn;
			Pn[n + nAgent] = qn;
		}
	}
	else { 
		if (index == 0) { 
			// bus des pertes
			float pn = ( rho * Chat[beginChat] - Cost2[0]) / (Cost1[0] + rho);
			float qn = ( rho * Chat[beginChat + 1] - Cost2[nAgent]) / (Cost1[nAgent] + rho);
			int offset = (losstype == 0) * nAgent + (losstype == 1) * 1;
			X[begining] = pn;
			X[begining + offset] = qn;
		}
		// puissance sur le bus 0
		int nb = nAgentByBus[0];
		int beginAgent = CoresAgentBusBegin[0];
		int begin0 = indiceBusBegin[0];
		int beginChat = CoresChatBegin[0];

		for (int i = index; i < nb; i += step) {
			int n = CoresAgentBus[i + beginAgent];
			float ub = Pmax[n];
			float lb = Pmin[n];
			float pn = (rho * Chat[beginChat + 4 + i] - Cost2[n]) / (Cost1[n] + rho);
			pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;


			ub = Pmax[n + nAgent];
			lb = Pmin[n + nAgent];
			float qn = (rho * Chat[beginChat + 4 + nb + i] - Cost2[n + nAgent]) / (Cost1[n + nAgent] + rho);
			qn = (ub - qn) * (qn > ub) + (lb - qn) * (qn < lb) + qn;

			// pn & qn
			X[begin0 + 4 + i] = pn;
			X[begin0 + 4 + nb + i] = qn;
			Pn[n] = pn;
			Pn[n + nAgent] = qn;
		}



	}
	


	// X =  {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}	

}


__global__ void updateXPnOPFADMMCons(float* X, float* Pn, float* Chat, float* nAgentByBus, float* indiceBusBegin, float* CoresChatBegin,
	float* CoresAgentBusBegin, float* CoresAgentBus, float* Cost1, float* Cost2, float* Pmin, float* Pmax, float rho, int losstype, int nBus, int nAgent) {

	int bus = blockIdx.x;
	int index = threadIdx.x;
	int step = blockDim.x;

	int begining = indiceBusBegin[bus];
	int beginChat = CoresChatBegin[bus];
		
	if (bus < nBus) {
		int nb = nAgentByBus[bus];
		int beginAgent = CoresAgentBusBegin[bus];
		for (int i = index; i < nb; i += step) {
			int n = CoresAgentBus[i + beginAgent];
			float ub = Pmax[n];
			float lb = Pmin[n];
			float pn = (rho * Chat[beginChat + 4 + i] - Cost2[n]) / (Cost1[n] + rho);
			pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;


			ub = Pmax[n + nAgent];
			lb = Pmin[n + nAgent];
			float qn = (rho * Chat[beginChat + 4 + nb + i] - Cost2[n + nAgent]) / (Cost1[n + nAgent] + rho);
			qn = (ub - qn) * (qn > ub) + (lb - qn) * (qn < lb) + qn;

			// pn & qn
			X[begining + 4 + i] = pn;
			X[begining + 4 + nb + i] = qn;
			Pn[n] = pn;
			Pn[n + nAgent] = qn;
		}
	}
	else {
		if (index == 0) {
			// bus des pertes
			float pn = (rho * Chat[beginChat] - Cost2[0]) / (Cost1[0] + rho);
			float qn = (rho * Chat[beginChat + 1] - Cost2[nAgent]) / (Cost1[nAgent] + rho);
			int offset = (losstype == 0) * nAgent + (losstype == 1) * 1;
			X[begining] = pn;
			X[begining + offset] = qn;
		}
	
	}

	// X =  {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}	

}



__global__ void updateXPn(float* X, float* Pn, float* P, float* nVoisin, float* indiceBusBegin, float* nAgentByBus, float* CoresAgentBusBegin, float* CoresAgentBus, int lossType, int nAgent, int nBus) {
	
	int bus = blockIdx.x; // un bloc par bus
	int index = threadIdx.x;
	int begin = indiceBusBegin[bus];
	int step = blockDim.x;
	
	if (bus == nBus) { // bus des pertes
		if (index == 0) {
			int indice = (lossType == 0) * (nAgent - 1) + 1;
			float pn = P[0] * nVoisin[0];
			float qn = P[nAgent] * nVoisin[nAgent];
			Pn[0] = pn;
			Pn[nAgent] = qn;
			X[begin] = pn;
			X[begin + indice] = qn;
		}
	}
	else {
		
		int Nb = nAgentByBus[bus];
		int beginAgent = CoresAgentBusBegin[bus];
		for (int In = index; In < Nb; In += step) {
			int n = CoresAgentBus[In + beginAgent];

			float pn = P[n] * nVoisin[n];
			float qn = P[n + nAgent] * nVoisin[n + nAgent];

			Pn[n] = pn;
			Pn[nAgent + n] = qn;
			// pn & qn
			X[begin + 4 + In] = pn;
			X[begin + 4 + Nb + In] = qn;
		}

	}



}

__global__ void communicateX(float* X, float* nChild, float* Ancestor, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* nAgentByBus, int nBus) {
	int bus = blockIdx.x;
	int index = threadIdx.x;
	int step = blockDim.x;

	int indice = indiceBusBegin[bus];
	int indiceChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
	int nb = nChild[bus];
	int nAgent = nAgentByBus[bus];

	if (index == 0) { // Vai
		int Ai = Ancestor[bus];
		int indiceAi = bus > 0 ? indiceBusBegin[Ai] : 0;
		X[indice + 4 + 2 * nAgent] = bus > 0 ? X[indiceAi + 3] : 1;
	}


	for (int voisin = index; voisin < nb; voisin += step) { //  coalescent en ecriture mais pas en lecture 
		int c = Childs[indiceChild + voisin];
		int indiceBusChild = indiceBusBegin[c];
		X[indice + 5 + 2 * nAgent + voisin] = X[indiceBusChild];
		X[indice + 5 + 2 * nAgent + nb + voisin] = X[indiceBusChild + 1];
		X[indice + 5 + 2 * nAgent + 2 * nb + voisin] = X[indiceBusChild + 2];
	}

}
__global__ void communicateX(float* X, float* nChild, float* Ancestor, float* Childs, float* indiceBusBegin, float* indiceChildBegin, int nBus) {
	int bus = blockIdx.x;
	int index = threadIdx.x;
	int step = blockDim.x;

	int indice = indiceBusBegin[bus];
	int indiceChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
	int nb = nChild[bus];

	if (index == 0) { // Vai
		int Ai = Ancestor[bus];
		int indiceAi = bus > 0 ? indiceBusBegin[Ai] : 0;
		X[indice + 6] = bus > 0 ? X[indiceAi + 3] : 1;
	}


	for (int voisin = index; voisin < nb; voisin += step) { // pas coalescent et am�liorable
		int c = Childs[indiceChild + voisin];
		int indiceBusChild = indiceBusBegin[c];
		X[indice + 7 + voisin] = X[indiceBusChild]; // un truc du style X[indice + 7 + voisin], X[indice + 7 + nb + voisin] et X[indice + 7 + 2*nb + voisin] serait coalescent
		X[indice + 7 + nb + voisin] = X[indiceBusChild + 1];
		X[indice + 7 + 2 * nb + voisin] = X[indiceBusChild + 2];
	}
	/*
	for (int i = 0; i < _nBus; i++) {

			if (i > 0) {
				int Ai = Ancestor.get(i, 0);
				X[i].set(6, 0, X[Ai].get(2, 0));
			}

			int m = nChild.get(i, 0);
			for (int j = 0; j < m; j++) {
				int c = Childs[i].get(j, 0);
				X[i].set(7 + 3 * j, 0, X[c].get(0, 0));
				X[i].set(8 + 3 * j, 0, X[c].get(1, 0));
				X[i].set(9 + 3 * j, 0, X[c].get(3, 0));

			}
		}

	*/


}

__global__ void communicateX(float* X, float* nChild, float* Ancestor, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* nAgentByBus, float* CoresBusAgent, float* PosAgent, int Losstype, int nBus, int nAgent) {
	int bus = blockIdx.x;
	int index = threadIdx.x;
	int step = blockDim.x;

	int indice = indiceBusBegin[bus];
	if (bus < nBus) {
		int indiceChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
		int nb = nChild[bus];
		int nAgent = nAgentByBus[bus];

		if (index == 0 && bus > 0) { // Vai
			int Ai = Ancestor[bus];
			int indiceAi = indiceBusBegin[Ai];
			X[indice + 4 + 2 * nAgent] = X[indiceAi + 3];
		}


		for (int voisin = index; voisin < nb; voisin += step) { //  coalescent en ecriture mais pas en lecture 
			int c = Childs[indiceChild + voisin];
			int indiceBusChild = indiceBusBegin[c];
			X[indice + 5 + 2 * nAgent + voisin] = X[indiceBusChild];
			X[indice + 5 + 2 * nAgent + nb + voisin] = X[indiceBusChild + 1];
			X[indice + 5 + 2 * nAgent + 2 * nb + voisin] = X[indiceBusChild + 2];
		}
	}
	else {
		if (Losstype == 1) {// Current
			for (int bus2 = index; bus2 < nBus; bus2 += step) {
				int indicebus = indiceBusBegin[bus2];
				X[indice + 2 + bus2] = X[indicebus + 2];
			}
		}
		else { // Puissance
			for (int n = 1 + index; n < nAgent; n += step) {
				int bus2 = CoresBusAgent[n];
				int In = PosAgent[n];
				int nAgentBus = nAgentByBus[bus2];

				X[indice + n] = X[bus + 4 + In];
				X[indice + n + nAgent] = X[bus + 4 + nAgentBus + In];
			}/**/
		}
	}
}


__global__ void setPnFromX(float* Pn, float* X, float* indiceBusBegin, float* CoresAgentBus, float* nAgentByBus, float* beginBusAgent, int nAgent) {
	int bus = blockIdx.x;
	int thI = threadIdx.x;
	int step = blockDim.x;
	int begin = beginBusAgent[bus];
	int nB = nAgentByBus[bus];
	int fin = begin + nB;
	int indiceBus = indiceBusBegin[bus];


	for (int i = begin + thI; i < fin; i += step) { // ecriture coalecente mais pas lecture
		int agent = CoresAgentBus[i];

		Pn[agent] = X[indiceBus + 4 + thI];
		Pn[agent + nAgent] = X[indiceBus + 4 + nB + thI];

	}
}





// FIN OPF ADMM


// Debut PF

__global__ void calcWinterCar(float* Pinter, float* Qinter, float* VoltageCar, float* Glin, float* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B) {


	int index = threadIdx.x;
	int step = blockDim.x;
	int i = blockIdx.x;
	extern __shared__ float shE[];
	int begin = CoresBusLin[i];
	int end = begin + nLines[blockIdx.x];
	int B2 = 2 * B;

	for (int n = index; n < B2; n += step)
	{
		shE[n] = VoltageCar[n];
	}
	__syncthreads();

	for (int l = begin + index; l < end; l += step) {
		int k = CoresVoiLin[l];
		float g = Glin[l];
		float b = Blin[l];

		float a = g * shE[k] - b * shE[k + B];
		float c = b * shE[k] + g * shE[k + B];

		Pinter[l] = shE[i] * a + shE[i + B] * c;
		Qinter[l] = shE[i + B] * a - shE[i] * c;

	}

}


__global__ void initECar(float* VoltageRealIm, float* E, int B) {


	int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
	int size = gridDim.x * blockDim.x;
	for (int i = thIdx; i < B; i += size) {
		float V0 = E[i + B];
		float theta0 = E[i];
		VoltageRealIm[i] = V0 * cos(theta0);
		VoltageRealIm[i + B] = V0 * sin(theta0);
	}
}


__global__ void initECar(float* VoltageRealImD, float v0, float w0, int B) {


	int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
	int size = gridDim.x * blockDim.x;
	for (int i = thIdx; i < B; i += size) {
		
		VoltageRealImD[i] = v0;
		VoltageRealImD[i + B] = w0;
	}
}

__global__ void calcEGPU(float* E, float* VoltageRealIm, int B) {

	int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
	int size = gridDim.x * blockDim.x;

	for (int i = thIdx; i < B; i += size) {
		float Rev = VoltageRealIm[i];
		float Imv = VoltageRealIm[i + B];


		E[i + B] = sqrt(Rev * Rev + Imv * Imv);
		E[i] = atan2(Imv, Rev);
	}

}


__global__ void calculYGPU(float* Y, float* E, float* Voltage, float* Blin2, float* Glin2, float* CoresLineBus, int B, int L) {

	int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
	int size = gridDim.x * blockDim.x;

	for (int b = thIdx; b < 2 * B; b += size) {
		Y[b] = E[b];
	}
	for (int l = thIdx; l < L; l += size) {
		int busFrom = CoresLineBus[l];
		int busTo = CoresLineBus[L + l];

		float ei = Voltage[busFrom];
		float fi = Voltage[busFrom + B];
		float ej = Voltage[busTo];
		float fj = Voltage[busTo + B];
		float Blin = Blin2[l];
		float Glin = Glin2[l];
		float Pij = (ei * ei + fi * fi - ei * ej - fi * fj) * Glin + (ei * fj - ej * fi) * Blin;
		Y[2 * B + l] = Pij;
	}

}


__global__ void calculYGPU(float* Y, float* E, float* Blin2, float* Glin2, float* CoresLineBus, int B, int L) {

	int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
	int size = gridDim.x * blockDim.x;

	for (int b = thIdx; b < 2 * B; b += size) {
		Y[b] = E[b];
	}
	for (int l = thIdx; l < L; l += size) {
		int busFrom = CoresLineBus[l];
		int busTo = CoresLineBus[L + l];
		float vi = E[busFrom + B];
		float thetai = E[busFrom];
		float vj = E[busTo + B];
		float thetaj = E[busTo];

		float ei = vi * cos(thetai);
		float fi = vi * sin(thetai);
		float ej = vj * cos(thetaj);
		float fj = vj * sin(thetaj);
		float Blin = Blin2[l];
		float Glin = Glin2[l];
		float Pij = (ei * ei + fi * fi - ei * ej - fi * fj) * Glin + (ei * fj - ej * fi) * Blin;
		
		Y[2 * B + l] = Pij;
	}

}


__global__ void calculPhiGPU(float* Phi, float* E, float* Blin2, float* Glin2, float* CoresLineBus, int B, int L) {

	int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
	int size = gridDim.x * blockDim.x;

	for (int l = thIdx; l < L; l += size) {
		int busFrom = CoresLineBus[l];
		int busTo = CoresLineBus[L + l];
		float vi = E[busFrom + B];
		float thetai = E[busFrom];
		float vj = E[busTo + B];
		float thetaj = E[busTo];

		float ei = vi * cos(thetai);
		float fi = vi * sin(thetai);
		float ej = vj * cos(thetaj);
		float fj = vj * sin(thetaj);
		float Blin = Blin2[l];
		float Glin = Glin2[l];
		float Pij = (ei * ei + fi * fi - ei * ej - fi * fj) * Glin + (ei * fj - ej * fi) * Blin;
		Phi[2 * B + l] = Pij;
	}

}



__global__ void setY(float* Y, float* E, float* Phi, int B, int L) {
	
	int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
	int size = gridDim.x * blockDim.x;

	for (int b = thIdx; b < 2 * B; b += size) {
		Y[b] = E[b];
	}
	for (int l = thIdx; l < L; l += size) {
		Y[2 * B + l] = Phi[l];
	}
}
// Fin PF
/*

double b = coef[0];
double d = coef[1];
double e = coef[2];
int nRoot = 0;

if (b * b * b + 8 * d == 0) {
	//if (abs(b * b * b + 8 * d) < 0.00000001) {
	// passage de p^4 + b p^3 + d p + e -> a p^4 + b p^2 + c = 0
	double B = -3 * b * b / 8;
	double C = -3 * b * b * b * b / 256 - b * d / 4 + e;

	double delta = B * B - 4 * C;
	//std::cout << "Delta " << delta;
	if (delta == 0) {
		double z = -B / 2;
		nRoot = 2;
		//std::cout << " z " << z << std::endl;
		root[0] = sqrt(z);
		root[1] = -sqrt(z);
		return nRoot;
	}
	else if (delta > 0) {
		double z1 = (-B + sqrt(delta)) / 2;
		double z2 = (-B - sqrt(delta)) / 2;
		//std::cout << " z1 " << z1 << " z2 " << z2 << std::endl;
		if (z1 >= 0) {
			root[0] = sqrt(z1);
			root[1] = -sqrt(z1);
			nRoot = 2;
		} if (z2 >= 0) {
			root[nRoot] = sqrt(z2);
			root[nRoot + 1] = -sqrt(z2);
			nRoot += 2;
		}
		return nRoot;
	}
	else { // delta < 0
		//std::cout << "pas de racines réelle !!!! rip, on tente le pas bicarré" << std::endl;
	}
}

// for the lambda polynome
double coef2[2];
double rootlambda[3];
coef2[0] = (2 * b * d - 8 * e) / 8;
coef2[1] = -(b * b * e + d * d) / 8;
int nRootlambda = resolveRealPolynome3without2term(rootlambda, coef2);




for (int i = 0; i < nRootlambda; i++) {
	double lambda0 = rootlambda[i];
	//std::cout << "poly3 " << coef2[0] * lambda0 + coef2[1] + lambda0 * lambda0 * lambda0 << std::endl;
	double mu1 = 2 * lambda0 + (b * b) / 4;
	if (mu1 > 0) {
		double mu0 = sqrt(mu1);
		double DeltaP = -2 * lambda0 + 2 * (d - b * lambda0) / mu0 + b * mu0 + b * b / 2;
		double DeltaM = -2 * lambda0 - 2 * (d - b * lambda0) / mu0 - b * mu0 + b * b / 2;
		if (DeltaP >= 0) {
			root[nRoot] = (-mu0 + sqrt(DeltaP)) / 2 - b / 4;
			root[nRoot + 1] = (-mu0 - sqrt(DeltaP)) / 2 - b / 4;
			nRoot = nRoot + 2;
			//std::cout << "  Dp   ";
		}
		if (DeltaM >= 0) {
			root[nRoot] = (mu0 + sqrt(DeltaM)) / 2 - b / 4;
			root[nRoot + 1] = (mu0 - sqrt(DeltaM)) / 2 - b / 4;
			nRoot = nRoot + 2;
			//std::cout << "  DM   ";
		}
		if (nRoot > 0) {
			//std::cout << "poly4 " << coef[0] << " " <<  coef[1] << " " << coef[2] << std::endl;
			return nRoot;
		}
	}

}
double lambda0 = rootlambda[0];
double mu0 = sqrt(2 * lambda0 + (b * b) / 4);
double DeltaP = -2 * lambda0 + 2 * (d - b * lambda0) / mu0 + b * mu0 + b * b / 2;
double DeltaM = -2 * lambda0 - 2 * (d - b * lambda0) / mu0 - b * mu0 + b * b / 2;
std::cout << "no real root found " << abs(b * b * b + 8 * d) << " " << lambda0 << " " << mu0 << " " << DeltaP << " " << DeltaM << std::endl;



*/







