
#include "../head/TestKernel.cuh"


// --------------------------------- Q part -------------------------------------------------------

float testCalculQpart(int method) {

	std::string fileName = "TempsQpart"+std::to_string(method) + ".csv";
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	float time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 6;
	const int nNLine = 6;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000 };//, 10000

	int nLine[nNLine] = { 10, 100, 500, 1000, 5000, 10000 }; // 10000 };
	int blockSize = 256;
	float values[nSimu];
	
	MatrixCPU temps(nNAgent*nNLine, nSimu, 0);
	MatrixCPU nAgentMat(1, nNAgent);
	MatrixCPU nLineMat(1, nNLine);
	for (int j = 0; j < nSimu; j++) {
		values[j] = (float)(rand()) / rand();
	}
	int indice = 0;
	for (int i = 0; i < nNAgent; i++) {
		nAgentMat.set(0, i, nAgent[i]);
		for (int j = 0; j < nNLine; j++) {
			nLineMat.set(0, j, nLine[j]);
			std::cout << "iteration (" << i << ", " << j << ") nAgent " << nAgent[i] << " nline " << nLine[j] << std::endl;
			for (int simu = 0; simu < nSimu; simu++) {
				MatrixGPU alpha(nLine[j], nAgent[i], values[j], 1);
				MatrixGPU Qpart(nLine[j], nAgent[i], 0, 1);
				MatrixGPU alphaTrans(nAgent[i], nLine[j], values[j], 1);
				MatrixGPU QpartTrans(nAgent[i], nLine[j], 0, 1);
				//alpha.setRand(values[j]);
				//alphaTrans.setTrans(&alpha);
				cudaDeviceSynchronize();
				time = 0;
				int numBlocks;
				for (int repet = 0; repet < nRepet; repet++) {
					MatrixGPU alphaCopy(alpha);
					MatrixGPU QpartCopy(Qpart);
					MatrixGPU alphaTransCopy(alphaTrans);
					MatrixGPU QpartTransCopy(QpartTrans);
					int N = nAgent[i];
					int L = nLine[j];
					cudaDeviceSynchronize();
					switch (method)
					{
					case 0:
						numBlocks = L;
						a = std::chrono::high_resolution_clock::now();
						calculQpartLineBloc <<<numBlocks, blockSize, N * sizeof(float) >> > (QpartCopy._matrixGPU, alphaCopy._matrixGPU, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 1:
						numBlocks = N;
						a = std::chrono::high_resolution_clock::now();
						calculQpartAgentBloc <<<numBlocks, blockSize >> > (QpartCopy._matrixGPU, alphaCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 2:
						numBlocks = L;
						a = std::chrono::high_resolution_clock::now();
						calculQpartLineBlocTrans << <numBlocks, blockSize, N * sizeof(float) >> > (QpartTransCopy._matrixGPU, alphaTransCopy._matrixGPU, N, L);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 3: 
						numBlocks = N;
						a = std::chrono::high_resolution_clock::now();
						calculQpartAgentBlocTrans << <numBlocks, blockSize >> > (QpartTransCopy._matrixGPU, alphaTransCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 4:
						numBlocks = L;
						a = std::chrono::high_resolution_clock::now();
						calculQpartLineBlocReverse << <numBlocks, blockSize, N * sizeof(float) >> > (QpartCopy._matrixGPU, alphaCopy._matrixGPU, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 5:
						numBlocks = L;
						a = std::chrono::high_resolution_clock::now();
						calculQpartLineBlocReverseTrans << <numBlocks, blockSize, N * sizeof(float) >> > (QpartTransCopy._matrixGPU, alphaTransCopy._matrixGPU, N, L);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 6:
						numBlocks = L;
						a = std::chrono::high_resolution_clock::now();
						calculQpartLineBlocReverseBis << <numBlocks, blockSize, N * sizeof(float) >> > (QpartCopy._matrixGPU, alphaCopy._matrixGPU, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 7:
						numBlocks = L;
						a = std::chrono::high_resolution_clock::now();
						calculQpartLineBlocReverseBisTrans << <numBlocks, blockSize, N * sizeof(float) >> > (QpartTransCopy._matrixGPU, alphaTransCopy._matrixGPU, N, L);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					default:
						return 0;
						break;
					}
					time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
				}
				temps.set(indice, simu, (float)time / nRepet);
			}
			indice++;
		}
	}
	nAgentMat.saveCSV(fileName, mode);
	nLineMat.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	return temps.sum();
}

__global__ void calculQpartLineBloc(float* Qpart, float* alpha, const int N)
{
	int index = threadIdx.x;
	int step = blockDim.x;
	int l = blockIdx.x;
	extern __shared__ float shAlpha[];

	for (int n = index; n < N; n += step)
	{
		shAlpha[n] = alpha[l * N + n];
	}
	__syncthreads();

	for (int n = index; n < N; n += step)
	{
		float s = 0;
		for (int j = n + 1; j < N; j++) {
			s += shAlpha[j]; // c'est moche cet accès de mémoire partagée
		}
		Qpart[l*N + n] = s;
	}
}

__global__ void calculQpartAgentBloc(float* Qpart, float* alpha, const int L, const int N) {
	int index = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	

	for (int l = index; l < L; l += step)
	{
		float s = 0;
		for (int j = n + 1; j < N; j++) {
			s += alpha[l*N+j]; 
		}
		Qpart[l*N + n] = s;
	}
}

__global__ void calculQpartLineBlocTrans(float* Qpart, float* alpha, const int N, const int nLine)
{
	int index = threadIdx.x;
	int step = blockDim.x;
	int l = blockIdx.x;
	extern __shared__ float shAlpha[];

	for (int n = index; n < N; n += step)
	{
		shAlpha[n] = alpha[ n*nLine + l]; // moche
	}
	__syncthreads();

	for (int n = index; n < N; n += step)
	{
		float s = 0;
		for (int j = n + 1; j < N; j++) {
			s += shAlpha[j]; // c'est moche cet accès de mémoire partagée
		}
		Qpart[ n*nLine +l] = s; // moche
	}
}

__global__ void calculQpartLineBlocReverseTrans(float* Qpart, float* alpha, const int N, const int nLine)
{
	int index = threadIdx.x;
	int step = blockDim.x;
	int l = blockIdx.x;
	extern __shared__ float shAlpha[];

	for (int n = index; n < N; n += step)
	{
		shAlpha[n] = alpha[n * nLine + l]; // moche
	}
	__syncthreads();
	float s_pre = 0;
	int n_pre = N - 1;
	for (int n = (N - index - 1); n >= 0; n -= step)
	{
		float s = 0;
		for (int j = n_pre; j > n; j--) {
			s += shAlpha[j]; // c'est moche cet accès de mémoire partagée
		}
		s = s + s_pre;
		Qpart[n * nLine + l] = s;
		s_pre = s;
		n_pre = n;
	}
}

__global__ void calculQpartAgentBlocTrans(float* Qpart, float* alpha, const int L, const int N) {
	int index = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	

	for (int l = index; l < L; l += step)
	{
		float s = 0;
		for (int j = n + 1; j < N; j++) {
			s += alpha[j * L + l]; // nombre de boucle dépend du bloc pas du thread, accès coalescent
		}
		Qpart[ n * L + l] = s;
	}
}

__global__ void calculQpartLineBlocReverse(float* Qpart, float* alpha, const int N) // est ce que cela marche ???????
{
	int index = threadIdx.x;
	int step = blockDim.x;
	int l = blockIdx.x;
	extern __shared__ float shAlpha[];

	for (int n = index; n < N; n += step)
	{
		shAlpha[n] = alpha[l * N + n];
	}
	__syncthreads();

	float s_pre = 0;
	int n_pre = N-1;
	for (int n = (N-index-1); n >=0; n -= step)
	{
		float s = 0;
		for (int j = n_pre; j > n; j--) {
			s += shAlpha[j]; // c'est moche cet accès de mémoire partagée
		}
		s = s + s_pre;
		Qpart[l * N + n] = s; 
		s_pre = s;
		n_pre = n;
	}
}

__global__ void calculQpartLineBlocReverseBis(float* Qpart, float* alpha, const int N) // est ce que cela marche ???????
{
	int index = threadIdx.x;
	int step = blockDim.x;
	int l = blockIdx.x;
	extern __shared__ float shAlpha[];

	for (int n = index; n < N; n += step)
	{
		shAlpha[n] = alpha[l * N + n];
	}
	__syncthreads();

	for (int n = (N - index - 1); n >= 0; n -= step)
	{
		float s = 0;
		for (int j = N-1; j > n; j--) {
			s += shAlpha[j]; // c'est moche cet accès de mémoire partagée
		}
		Qpart[l * N + n] = s;
	}
}

__global__ void calculQpartLineBlocReverseBisTrans(float* Qpart, float* alpha, const int N, const int nLine) // est ce que cela marche ???????
{
	int index = threadIdx.x;
	int step = blockDim.x;
	int l = blockIdx.x;
	extern __shared__ float shAlpha[];

	for (int n = index; n < N; n += step)
	{
		shAlpha[n] = alpha[l * N + n];
	}
	__syncthreads();

	for (int n = (N - index - 1); n >= 0; n -= step)
	{
		float s = 0;
		for (int j = N - 1; j > n; j--) {
			s += shAlpha[j]; // c'est moche cet accès de mémoire partagée
		}
		Qpart[n * nLine + l] = s;
	}
}

// --------------------------------- alpha -------------------------------------------------------

float testCalculAlpha(int method)
{
	std::string fileName = "TempsAlpha" + std::to_string(method) + ".csv";
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 6;
	const int nNLine = 6;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000 };
	int nLine[nNLine] = { 10, 100, 500, 1000, 5000, 10000 };
	int blockSize = 256;
	float values[nSimu];
	float values2[nSimu];

	MatrixCPU temps(nNAgent * nNLine, nSimu, 0);
	MatrixCPU nAgentMat(1, nNAgent);
	MatrixCPU nLineMat(1, nNLine);
	for (int j = 0; j < nSimu; j++) {
		values[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
	}
	int indice = 0;
	for (int i = 0; i < nNAgent; i++) {
		nAgentMat.set(0, i, nAgent[i]);
		for (int j = 0; j < nNLine; j++) {
			nLineMat.set(0, j, nLine[j]);
			//std::cout << "iteration (" << i << ", " << j << ") nAgent " << nAgent[i] << " nline " << nLine[j] << std::endl;
			for (int simu = 0; simu < nSimu; simu++) {
				MatrixGPU G(nLine[j], nAgent[i], values[j], 1);
				MatrixGPU Pn(nAgent[i], 1, values2[j], 1);
				MatrixGPU GTrans(nAgent[i], nLine[j], values[j], 1);
				MatrixGPU alpha(nLine[j], nAgent[i], 0, 1);
				MatrixGPU alphaTrans(nAgent[i], nLine[j], 0, 1);
				//G.setRand(values[j]);
				//GTrans.setTrans(&G);
				//Pn.setRand(values2[j]);
				cudaDeviceSynchronize();
				time = 0;
				int numBlocks;
				
				for (int repet = 0; repet < nRepet; repet++) {
					MatrixGPU alphaCopy(alpha);
					MatrixGPU PnCopy(Pn);
					MatrixGPU GCopy(G);
					MatrixGPU alphaTransCopy(alphaTrans);
					MatrixGPU GTransCopy(GTrans);
					int N = nAgent[i];
					int L = nLine[j];
					const int nThread = 16;
					const int bx = (N + nThread - 1) / nThread;
					const int by = (L + nThread - 1) / nThread;
					dim3 dimBlock(nThread, nThread);
					dim3 gridBlock(bx, by);
					dim3 gridBlockTrans(by, bx);
				
					cudaDeviceSynchronize();
					switch (method)
					{
					case 0:
						numBlocks = N;
						a = std::chrono::high_resolution_clock::now();
						updateAlphaSh << <numBlocks, blockSize >> > (alphaCopy._matrixGPU, GCopy._matrixGPU, PnCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 1:
						a = std::chrono::high_resolution_clock::now();
						updateAlpha2D << <gridBlock, dimBlock >> > (alphaCopy._matrixGPU, GCopy._matrixGPU, PnCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 2:
						numBlocks = (N * L + blockSize - 1) / blockSize;
						a = std::chrono::high_resolution_clock::now();
						updateAlpha1D << <numBlocks, blockSize >> > (alphaCopy._matrixGPU, GCopy._matrixGPU, PnCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 3: 
						numBlocks = N;
						a = std::chrono::high_resolution_clock::now();
						updateAlphaShTrans << <numBlocks, blockSize >> > (alphaTransCopy._matrixGPU, GTransCopy._matrixGPU, PnCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 4: 
						a = std::chrono::high_resolution_clock::now();
						updateAlpha2DTrans << <gridBlockTrans, dimBlock >> > (alphaTransCopy._matrixGPU, GTransCopy._matrixGPU, PnCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 5: 
						numBlocks = (N * L + blockSize - 1) / blockSize;
						a = std::chrono::high_resolution_clock::now();
						updateAlpha1DTrans << <numBlocks, blockSize >> > (alphaTransCopy._matrixGPU, GTransCopy._matrixGPU, PnCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
						
					default:
						return 0;
						break;
					}
					time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
				}
				temps.set(indice, simu, (float)time / nRepet);
			}
			indice++;
		}
	}
	nAgentMat.saveCSV(fileName, mode);
	nLineMat.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	return temps.sum();
}


__global__ void updateAlphaSh(float* alpha, float* G, float* Pn, const int nLine, const int nAgent)
{
	// un bloc par agent
	int index = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	__shared__ float shPn;
	if (index == 0) {
		shPn = Pn[n];
	}
	__syncthreads();
	for (int l = index; l < nLine; l += step)
	{
		alpha[l * nAgent + n] = G[l * nAgent + n] * shPn;
	}
}

__global__ void updateAlpha2D(float* alpha, float* G, float* Pn, const int nLine, const int nAgent)
{
	int indexX = threadIdx.x + blockIdx.x*blockDim.x;
	int stepX = blockDim.x*gridDim.x;
	int indexY = threadIdx.y + blockIdx.y * blockDim.y;
	int stepY = blockDim.y * gridDim.y;
	
	for (int n = indexX; n < nAgent; n += stepX)
	{
		float PnLocal = Pn[n];
		for (int l = indexY; l < nLine; l += stepY)
		{
			alpha[l * nAgent + n] = G[l * nAgent + n] * PnLocal;
		}
	}
}

__global__ void updateAlpha1D(float* alpha, float* G, float* Pn, const int nLine, const int nAgent)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int step = blockDim.x * gridDim.x;
	int N = nAgent * nLine;

	for (int i = index; i < N; i += step)
	{
		int k = i % nAgent;
		alpha[i] = G[i] * Pn[k];
	}
}

__global__ void updateAlphaShTrans(float* alpha, float* G, float* Pn, const int nLine, const int nAgent)
{
	// un bloc par agent
	int index = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	__shared__ float shPn;
	if (index == 0) {
		shPn = Pn[n];
	}
	__syncthreads();
	for (int l = index; l < nLine; l += step)
	{
		alpha[ n * nLine + l ] = G[n * nLine + l] * shPn;
	}
}

__global__ void updateAlpha2DTrans(float* alpha, float* G, float* Pn, const int nLine, const int nAgent)
{
	// alpha et G en (n,l)
	int indexX = threadIdx.x + blockIdx.x * blockDim.x;
	int stepX = blockDim.x * gridDim.x;
	int indexY = threadIdx.y + blockIdx.y * blockDim.y;
	int stepY = blockDim.y * gridDim.y;

	for (int n = indexY; n < nAgent; n += stepY)
	{
		float PnLocal = Pn[n];
		for (int l = indexX; l < nLine; l += stepX)
		{
			alpha[n*nLine +l] = G[n * nLine + l] * PnLocal;
		}
	}
}

__global__ void updateAlpha1DTrans(float* alpha, float* G, float* Pn, const int nLine, const int nAgent) {

	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int step = blockDim.x * gridDim.x;
	int N = nAgent * nLine;

	for (int i = index; i < N; i += step)
	{
		int k = i / nAgent;
		alpha[i] = G[i] * Pn[k];
	}

}

// --------------------------------- Cpa -------------------------------------------------------

float testCalculCpa(int method) {
	std::string fileName = "TempsCpa" + std::to_string(method) + ".csv";
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 6;
	const int nNLine = 6;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000 };
	int nLine[nNLine] = { 10, 100, 500, 1000, 5000, 10000 };
	int blockSize = 256;
	float values[nSimu];
	float values2[nSimu];

	MatrixCPU temps(nNAgent * nNLine, nSimu, 0);
	MatrixCPU nAgentMat(1, nNAgent);
	MatrixCPU nLineMat(1, nNLine);
	for (int j = 0; j < nSimu; j++) {
		values[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
	}
	int indice = 0;
	for (int i = 0; i < nNAgent; i++) {
		nAgentMat.set(0, i, nAgent[i]);
		for (int j = 0; j < nNLine; j++) {
			nLineMat.set(0, j, nLine[j]);
			//std::cout << "iteration (" << i << ", " << j << ") nAgent " << nAgent[i] << " nline " << nLine[j] << std::endl;
			for (int simu = 0; simu < nSimu; simu++) {
				MatrixGPU Cp2(nAgent[i], 1, 0, 1);
				MatrixGPU tempL1(nLine[j], 1, values[j], 1);
				MatrixGPU G(nLine[j], nAgent[i], values2[j], 1);
				MatrixGPU GTrans(nAgent[i], nLine[j], values2[j], 1);
				//tempL1.setRand(values[j]);
				//G.setRand(values2[j]);
				//GTrans.setTrans(&G);
				cudaDeviceSynchronize();
				time = 0;
				int numBlocks;
				for (int repet = 0; repet < nRepet; repet++) {
					MatrixGPU Cp2Copy(Cp2);
					MatrixGPU tempL1Copy(tempL1);
					MatrixGPU GCopy(G);
					MatrixGPU GTransCopy(GTrans);
					int N = nAgent[i];
					int L = nLine[j];
					cudaDeviceSynchronize();
					switch (method)
					{
					case 0:
						numBlocks = N;
						a = std::chrono::high_resolution_clock::now();
						updateCp2aTest<256> << <numBlocks, blockSize >> > (Cp2Copy._matrixGPU, tempL1Copy._matrixGPU, GCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 1:
						numBlocks = N;
						a = std::chrono::high_resolution_clock::now();
						updateCp2aTestTrans<256> << <numBlocks, blockSize >> > (Cp2Copy._matrixGPU, tempL1Copy._matrixGPU, GTransCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					default:
						return 0;
						break;
					}
					time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
				}
				temps.set(indice, simu, (float)time / nRepet);
			}
			indice++;
		}
	}
	nAgentMat.saveCSV(fileName, mode);
	nLineMat.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	return temps.sum();
}

template <unsigned int blockSize>
__global__ void updateCp2aTest(float* Cp2, float* diffKappa, float* G, const int nLine, const int nAgent) {
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {

		float t = G[j * nAgent + n] * diffKappa[j];
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduceTest<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		Cp2[n] = shArr[0];
	}
}

template <unsigned int blockSize>
__global__ void updateCp2aTestTrans(float* Cp2, float* diffKappa, float* G, const int nLine, const int nAgent) {
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {

		float t = G[ n * nLine + j] * diffKappa[j];
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduceTest<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		Cp2[n] = shArr[0];
	}
}

// --------------------------------- Cpb -------------------------------------------------------

float testCalculCpb(int method)
{
	std::string fileName = "TempsCpb" + std::to_string(method) + ".csv";
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 6;
	const int nNLine = 6;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000 };
	int nLine[nNLine] = { 10, 100, 500, 1000, 5000, 10000 };
	int blockSize = 256;
	float values[nSimu];
	float values2[nSimu];

	MatrixCPU temps(nNAgent * nNLine, nSimu, 0);
	MatrixCPU nAgentMat(1, nNAgent);
	MatrixCPU nLineMat(1, nNLine);
	for (int j = 0; j < nSimu; j++) {
		values[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
	}
	int indice = 0;
	for (int i = 0; i < nNAgent; i++) {
		nAgentMat.set(0, i, nAgent[i]);
		for (int j = 0; j < nNLine; j++) {
			nLineMat.set(0, j, nLine[j]);
			int N = nAgent[i];
			int L = nLine[j];
			//std::cout << "iteration (" << i << ", " << j << ") nAgent " << nAgent[i] << " nline " << nLine[j] << std::endl;
			for (int simu = 0; simu < nSimu; simu++) {
				MatrixGPU tempN1(N, 1, 0, 1);
				MatrixGPU Qpart(L, N, values[simu], 1);
				MatrixGPU QpartTrans(N, L, values[simu], 1);
				MatrixGPU G(L, N, values2[simu], 1);
				MatrixGPU GTrans(N, L, values2[simu], 1);
				//Qpart.setRand(values[simu]); en faisant de l'aléatoire, la mesure "plante" au bout d'un moment...
				//G.setRand(values2[simu]);
				//GTrans.setTrans(&G);
				//QpartTrans.setTrans(&Qpart);
				cudaDeviceSynchronize();
				time = 0;
				int numBlocks;
				for (int repet = 0; repet < nRepet; repet++) {
					MatrixGPU tempN1Copy(tempN1);
					MatrixGPU QpartCopy(Qpart);
					MatrixGPU QpartTransCopy(QpartTrans);
					MatrixGPU GCopy(G);
					MatrixGPU GTransCopy(GTrans);
					
					cudaDeviceSynchronize();
					switch (method)
					{
					case 0:
						numBlocks = N;
						a = std::chrono::high_resolution_clock::now();
						updateCp2bTest<256> << <numBlocks, blockSize >> > (tempN1Copy._matrixGPU, GCopy._matrixGPU, QpartCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 1:
						numBlocks = N;
						a = std::chrono::high_resolution_clock::now();
						updateCp2bTestTrans<256> << <numBlocks, blockSize >> > (tempN1Copy._matrixGPU, GTransCopy._matrixGPU, QpartTransCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					default:
						return 0;
						break;
					}
					time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
				}
				temps.set(indice, simu, (float)time / nRepet);
			}
			indice++;
		}
	}
	nAgentMat.saveCSV(fileName, mode);
	nLineMat.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	return temps.sum();
}
template <unsigned int blockSize>
__global__ void updateCp2bTest(float* tempN1, float* G, float* Qpart, const int nLine, const int nAgent)
{
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {

		float t = G[j * nAgent + n] * Qpart[j * nAgent + n];
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduceTest<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		tempN1[n] = 2 * shArr[0];
	}

}

template <unsigned int blockSize>
__global__ void updateCp2bTestTrans(float* tempN1, float* G, float* Qpart, const int nLine, const int nAgent)
{
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {

		float t = G[n * nLine + j] * Qpart[n * nLine + j];
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduceTest<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		tempN1[n] = 2 * shArr[0];
	}

}

// --------------------------------- Q tot -------------------------------------------------------

float testCalculQtot(int method)
{
	std::string fileName = "TempsQtot" + std::to_string(method) + ".csv";
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 6;
	const int nNLine = 6;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000 };
	int nLine[nNLine] = { 10, 100, 500, 1000, 5000, 10000 };
	int blockSize = 256;
	float values[nSimu];
	float values2[nSimu];

	MatrixCPU temps(nNAgent * nNLine, nSimu, 0);
	MatrixCPU nAgentMat(1, nNAgent);
	MatrixCPU nLineMat(1, nNLine);
	for (int j = 0; j < nSimu; j++) {
		values[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
	}
	int indice = 0;
	for (int i = 0; i < nNAgent; i++) {
		nAgentMat.set(0, i, nAgent[i]);
		for (int j = 0; j < nNLine; j++) {
			nLineMat.set(0, j, nLine[j]);
			//std::cout << "iteration (" << i << ", " << j << ") nAgent " << nAgent[i] << " nline " << nLine[j] << std::endl;
			for (int simu = 0; simu < nSimu; simu++) {
				MatrixGPU Qtot(nLine[j], 1, 0, 1);
				MatrixGPU Qpart(nLine[j], nAgent[i], values[j], 1);
				MatrixGPU QpartTrans(nAgent[i], nLine[j], values[j], 1);
				MatrixGPU alpha(nLine[j], nAgent[i], values2[j], 1);
				MatrixGPU alphaTrans(nAgent[i], nLine[j], values2[j], 1);
				//Qpart.setRand(values[j]);
				//alpha.setRand(values2[j]);
				//QpartTrans.setTrans(&Qpart);
				//alphaTrans.setTrans(&alpha);
				
				cudaDeviceSynchronize();
				time = 0;
				int numBlocks;
				for (int repet = 0; repet < nRepet; repet++) {
					MatrixGPU QtotCopy(Qtot);
					MatrixGPU QpartCopy(Qpart);
					MatrixGPU QpartTransCopy(QpartTrans);
					MatrixGPU alphaCopy(alpha);
					MatrixGPU alphaTransCopy(alphaTrans);
					int N = nAgent[i];
					int L = nLine[j];
					cudaDeviceSynchronize();
					switch (method)
					{
					case 0:
						a = std::chrono::high_resolution_clock::now();
						QtotCopy.sum(&alphaCopy);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 1:
						numBlocks = ( L + blockSize - 1) / blockSize;
						a = std::chrono::high_resolution_clock::now();
						updateQtotTest << <numBlocks, blockSize >> > (QtotCopy._matrixGPU, QpartCopy._matrixGPU, alphaCopy._matrixGPU, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 2:
						numBlocks = (L + blockSize - 1) / blockSize;
						a = std::chrono::high_resolution_clock::now();
						updateQtotTestTrans << <numBlocks, blockSize >> > (QtotCopy._matrixGPU, QpartTransCopy._matrixGPU, alphaTransCopy._matrixGPU, L);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					default:
						return 0;
						break;
					}
					time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
				}
				temps.set(indice, simu, (float)time / nRepet);
			}
			indice++;
		}
	}
	nAgentMat.saveCSV(fileName, mode);
	nLineMat.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	return temps.sum();
}

__global__ void updateQtotTest(float* Qtot, float* Qpart, float* alpha, const int nLine, const int nAgent) {

	
	int thIdx = threadIdx.x + blockIdx.x * blockDim.x;
	int step = blockDim.x * gridDim.x;
	
	for (int l = thIdx; l < nLine; l += step) {

		Qtot[l] = Qpart[l * nAgent] + alpha[l * nAgent];
	}
}

__global__ void updateQtotTestTrans(float* Qtot, float* Qpart, float* alpha, const int nLine) {

	
	int thIdx = threadIdx.x + blockIdx.x * blockDim.x;
	int step = blockDim.x * gridDim.x;

	for (int l = thIdx; l < nLine; l += step) {

		Qtot[l] = Qpart[l] + alpha[l];
	}
}

template <unsigned int blockSize>
__device__ void warpReduceTest(volatile float* sdata, unsigned int tid) {
	if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
	if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
	if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
	if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
	if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
	if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
}

// --------------------------------- Cp -------------------------------------------------------
float testCalculCp(int method) {
	std::string fileName = "TempsCp" + std::to_string(method) + ".csv";
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 6;
	const int nNLine = 6;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000 };
	int nLine[nNLine] = { 10, 100, 500, 1000, 5000, 10000 };
	int blockSize = 256;
	float values[nSimu];
	float values2[nSimu];
	float values3[nSimu];
	float values4[nSimu];

	MatrixCPU temps(nNAgent * nNLine, nSimu, 0);
	MatrixCPU nAgentMat(1, nNAgent);
	MatrixCPU nLineMat(1, nNLine);
	for (int j = 0; j < nSimu; j++) {
		values[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
		values3[j] = (float)(rand()) / rand();
		values4[j] = (float)(rand()) / rand();
	}
	int indice = 0;
	for (int i = 0; i < nNAgent; i++) {
		nAgentMat.set(0, i, nAgent[i]);
		for (int j = 0; j < nNLine; j++) {
			nLineMat.set(0, j, nLine[j]);
			std::cout << "iteration (" << i << ", " << j << ") nAgent " << nAgent[i] << " nline " << nLine[j] << std::endl;
			for (int simu = 0; simu < nSimu; simu++) {
				MatrixGPU Cp2(nAgent[i], 1, 0, 1);
				MatrixGPU tempL1(nLine[j], 1, values[j], 1);
				MatrixGPU G(nLine[j], nAgent[i], values2[j], 1);
				MatrixGPU GTrans(nAgent[i], nLine[j], values2[j], 1);
				MatrixGPU Qpart(nLine[j], nAgent[i], values3[j], 1);
				MatrixGPU QpartTrans(nAgent[i], nLine[j], values3[j], 1);
				MatrixGPU nVoisin(nAgent[i], 1, nAgent[i], 1);
				float rho1 = values4[j];
				//tempL1.setRand(values[j]);
				//G.setRand(values2[j]);
				//GTrans.setTrans(&G);
				cudaDeviceSynchronize();
				time = 0;
				int numBlocks;
				for (int repet = 0; repet < nRepet; repet++) {
					MatrixGPU Cp2Copy(Cp2);
					MatrixGPU tempL1Copy(tempL1);
					MatrixGPU GCopy(G);
					MatrixGPU GTransCopy(GTrans);
					MatrixGPU QpartCopy(Qpart);
					MatrixGPU QpartTransCopy(QpartTrans);
					MatrixGPU nVoisinCopy(nVoisin);
					int N = nAgent[i];
					int L = nLine[j];
					cudaDeviceSynchronize();
					switch (method)
					{
					case 0:
						numBlocks = N;
						a = std::chrono::high_resolution_clock::now();
						updateCp2Test<256> << <numBlocks, blockSize >> > (Cp2Copy._matrixGPU, tempL1Copy._matrixGPU, GCopy._matrixGPU, QpartCopy._matrixGPU, nVoisinCopy._matrixGPU, rho1, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 1:
						numBlocks = N;
						a = std::chrono::high_resolution_clock::now();
						updateCp2TestTrans<256> << <numBlocks, blockSize >> > (Cp2Copy._matrixGPU, tempL1Copy._matrixGPU, GTransCopy._matrixGPU, QpartTransCopy._matrixGPU, nVoisinCopy._matrixGPU, rho1, L, N);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					default:
						return 0;
						break;
					}
					time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
				}
				temps.set(indice, simu, (float)time / nRepet);
			}
			indice++;
		}
	}
	nAgentMat.saveCSV(fileName, mode);
	nLineMat.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	return temps.sum();
}

template <unsigned int blockSize>
__global__ void updateCp2Test(float* Cp2, float* diffKappa, float* G, float* Qpart, float* nVoisin, float rho1, const int nLine, const int nAgent) {
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {

		float t = G[j * nAgent + n] * (diffKappa[j] + 2 * Qpart[j * nAgent + n]);
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduceTest<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		Cp2[n] = rho1 * nVoisin[n] * shArr[0];
	}
}

template <unsigned int blockSize>
__global__ void updateCp2TestTrans(float* Cp2, float* diffKappa, float* G, float* Qpart,float* nVoisin, float rho1, const int nLine, const int nAgent) {
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {
		float Gloc = G[n * nLine + j];
		float dKloc = diffKappa[j];
		float Q = Qpart[n * nLine + j];
		float t = Gloc * (dKloc + 2 * Q);
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduceTest<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		Cp2[n] = rho1 * nVoisin[n] * shArr[0];
	}
}

// --------------------------------- RexX -------------------------------------------------------
float testCalculResX(int method)
{
	std::string fileName = "TempsResX" + std::to_string(method) + ".csv";
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 1;
	const int nNLine = 8;
	const int nSimu = 100;
	const int nRepet = 10;
	int nLine[nNLine] = { 10, 100, 500, 1000, 5000, 10000, 100000, 1000000 };
	int blockSize = 256;
	float values[nSimu];
	float values2[nSimu];
	float values3[nSimu];
	float values4[nSimu];

	MatrixCPU temps(nNAgent * nNLine, nSimu, 0);
	MatrixCPU nAgentMat(1, nNAgent);
	MatrixCPU nLineMat(1, nNLine);
	for (int j = 0; j < nSimu; j++) {
		values[j] =  (float)(rand()) / rand() * (0.5 - rand() % 2);
		values2[j] = (float)(rand()) / rand() * (0.5 - rand() % 2);
		values3[j] = (float)(rand()) / rand() * (0.5 - rand() % 2);
		values4[j] = (float)(rand()) / rand() * (0.5 - rand() % 2);
	}
	int indice = 0;
	for (int i = 0; i < nNAgent; i++) {
		for (int j = 0; j < nNLine; j++) {
			nLineMat.set(0, j, nLine[j]);
			int numBlocks = ceil((nLine[j] + blockSize - 1) / blockSize);
			std::cout << "iteration (" << i << ", " << j << ") "  << " nline " << nLine[j] << std::endl;
			for (int simu = 0; simu < nSimu; simu++) {
				MatrixGPU res(nLine[j], 1, 0, 1);
				MatrixGPU tempL2(nLine[j], 1, 0, 1);
				MatrixGPU kappa1(nLine[j], 1, values[j], 1);
				MatrixGPU kappa2(nLine[j], 1, values2[j], 1);
				MatrixGPU kappaPre1(nLine[j], 1, values3[j], 1);
				MatrixGPU kappaPre2(nLine[j], 1, values4[j], 1);
				cudaDeviceSynchronize();
				time = 0;
				for (int repet = 0; repet < nRepet; repet++) {
					MatrixGPU resCopy(res);
					MatrixGPU tempL2Copy(tempL2);
					MatrixGPU kappa1Copy(kappa1);
					MatrixGPU kappa2Copy(kappa2);
					MatrixGPU kappaPre1Copy(kappaPre1);
					MatrixGPU kappaPre2Copy(kappaPre2);
					int L = nLine[j];
					cudaDeviceSynchronize();
					switch (method)
					{
					case 0:
						a = std::chrono::high_resolution_clock::now();
						resCopy.set(&kappa1Copy);
						tempL2.set(&kappa2Copy);
						kappaPre1Copy.projectNeg();
						kappaPre2Copy.projectNeg();
						resCopy.projectNeg();
						tempL2.projectNeg();
						resCopy.subtract(&kappaPre1Copy);
						tempL2.subtract(&kappaPre2Copy);
						resCopy.multiplyT(&resCopy);
						tempL2.multiplyT(&tempL2);
						resCopy.add(&tempL2);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 1:
						a = std::chrono::high_resolution_clock::now();
						updateResXTest << <numBlocks, blockSize >> > (resCopy._matrixGPU, kappa1Copy._matrixGPU, kappa2Copy._matrixGPU, kappaPre1Copy._matrixGPU, kappaPre2Copy._matrixGPU, L);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					default:
						return 0;
						break;
					}
					time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
				}
				temps.set(indice, simu, (float)time / nRepet);
			}
			indice++;
		}
	}
	nAgentMat.saveCSV(fileName, mode);
	nLineMat.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	return temps.sum();
}

__global__ void updateResXTest(float* res, float* Kappa1, float* Kappa2, float* KappaPre1, float* KappaPre2, const int nLine) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < nLine; l += step)
	{
		float k1 = (Kappa1[l] < 0) * Kappa1[l];
		float k2 = (Kappa2[l] < 0) * Kappa2[l];
		float k1pre = (KappaPre1[l] < 0) * KappaPre1[l];
		float k2pre = (KappaPre2[l] < 0) * KappaPre2[l];

		k1 -= k1pre;
		k2 -= k2pre;

		res[l] = k1 * k1 + k2 * k2;	
	}
}

// --------------------------------- LAMBDABt1 -------------------------------------------------------

float testCalculLAMBDABt1(int method)
{
	std::string fileName = "TempsLAMBDA" + std::to_string(method) + ".csv";
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 6;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000 };
	int ntrade[nNAgent];
	int blockSize = 256;
	float values[nSimu];
	float values2[nSimu];
	float rhos[nSimu];
	MatrixCPU temps(nNAgent, nSimu, 0);
	MatrixCPU nAgentMat(1, nNAgent, 0);
	MatrixCPU nTradeMat(1, nNAgent, 0);
	for (int j = 0; j < nSimu; j++) {
		values[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
		rhos[j] = (float)(rand() % 100) / rand();
	}
	for (int i = 0; i < nNAgent; i++) {
		ntrade[i] = nAgent[i] * nAgent[i] / 2;
		int numBlocks = ceil((ntrade[i] + blockSize - 1) / blockSize);
		nAgentMat.set(0, i, nAgent[i]);
		nTradeMat.set(0, i, ntrade[i]);
		MatrixGPU Bt1(ntrade[i], 1, 0, 1);
		MatrixGPU CoresLinAgent(ntrade[i], 1);
		MatrixGPU CoresLinVoisin(ntrade[i], 1);
		MatrixGPU CoresMatLin(nAgent[i], nAgent[i], -1);
		MatrixGPU CoresLinTrans(ntrade[i], 1);
		int indice = 0;
		int idVoisin = 0;

		for (int idAgent = 0; idAgent < nAgent[i]; idAgent++) {
			int Nvoisinmax = nAgent[i] / 2;
			if (idAgent < Nvoisinmax) {
				idVoisin = 0;
			}
			else {
				idVoisin = nAgent[i] / 2;
			}
			for (int voisin = idVoisin; voisin < Nvoisinmax; voisin++) {
				CoresLinAgent.set(indice, 0, idAgent);
				CoresLinVoisin.set(indice, 0, voisin);
				CoresMatLin.set(idAgent, voisin, indice);
				indice = indice + 1;
			}
		}
		for (int lin = 0; lin < ntrade[i]; lin++) {
			int i = CoresLinAgent.get(lin, 0);
			int j = CoresLinVoisin.get(lin, 0);
			int k = CoresMatLin.get(j, i);
			CoresLinTrans.set(lin, 0, k);
		}
		CoresLinAgent.transferGPU();
		CoresLinVoisin.transferGPU();
		CoresMatLin.transferGPU();
		CoresLinTrans.transferGPU();
		//std::cout << "iteration (" << i << ", " << j << ") nAgent " << nAgent[i] << " nline " << nLine[j] << std::endl;
		for (int simu = 0; simu < nSimu; simu++) {
			MatrixGPU LAMBDALin(ntrade[i], 1, values[simu], 1);
			MatrixGPU trade(ntrade[i], 1, values2[simu], 1);
			float rho = rhos[simu];
			cudaDeviceSynchronize();
			time = 0;
			for (int repet = 0; repet < nRepet; repet++) {
				MatrixGPU LAMBDALinCopy(LAMBDALin);
				MatrixGPU tradeCopy(trade);
				MatrixGPU CoresLinTransCopy(CoresLinTrans);
				MatrixGPU Bt1Copy(Bt1);
				int M = ntrade[i];
				cudaDeviceSynchronize();
				switch (method)
				{
				case 0:
					a = std::chrono::high_resolution_clock::now();
					updateLAMBDAGPUTest << <numBlocks, blockSize >> > (LAMBDALinCopy._matrixGPU, tradeCopy._matrixGPU, rho, CoresLinTransCopy._matrixGPU, M);
					updateBt1GPUTest << <numBlocks, blockSize >> > (Bt1Copy._matrixGPU, tradeCopy._matrixGPU, rho, LAMBDALinCopy._matrixGPU, CoresLinTransCopy._matrixGPU, M);
					cudaDeviceSynchronize();
					b = std::chrono::high_resolution_clock::now();
					break;
				case 1:
					a = std::chrono::high_resolution_clock::now();
					updateLAMBDABt1GPUTest << <numBlocks, blockSize >> > (Bt1Copy._matrixGPU, LAMBDALinCopy._matrixGPU, tradeCopy._matrixGPU, rho, CoresLinTransCopy._matrixGPU, M);
					cudaDeviceSynchronize();
					b = std::chrono::high_resolution_clock::now();
					break;
				default:
					return 0;
					break;
				}
				time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
			}
			temps.set(i, simu, (float)time / nRepet);
		}
	}
	nAgentMat.saveCSV(fileName, mode);
	nTradeMat.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	return temps.sum();
}


__global__ void updateLAMBDAGPUTest(float* LAMBDALin, float* tradeLin, float rho, float* CoresLinTrans, int const N)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < N; l += step)
	{
		float m = LAMBDALin[l];
		int k = CoresLinTrans[l];
		LAMBDALin[l] = m + 0.5 * rho * (tradeLin[l] + tradeLin[k]);
	}
}
__global__ void updateBt1GPUTest(float* Bt1, float* tradeLin, float rho, float* LAMBDA, float* CoresLinTrans, int const N)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < N; l += step)
	{
		int k = CoresLinTrans[l];
		Bt1[l] = 0.5 * (tradeLin[l] - tradeLin[k]) - LAMBDA[l] / rho;
	}

}

__global__ void updateLAMBDABt1GPUTest(float* Bt1, float* LAMBDA, float* tradeLin, float rho, float* CoresLinTrans, int const N) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < N; l += step)
	{
		int k = CoresLinTrans[l];
		float m = LAMBDA[l] + 0.5 * rho * (tradeLin[l] + tradeLin[k]);
		Bt1[l] = 0.5 * (tradeLin[l] - tradeLin[k]) - m / rho;
		LAMBDA[l] = m;
	}


}

// --------------------------------- Chat -------------------------------------------------------
float testCalculChat(int method, int blockSize, int repartition) {
	std::string fileName = "TempsCalculChat" + std::to_string(method) + "_" + std::to_string(blockSize) + "_" + std::to_string(repartition) + ".csv";
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNBus = 8;
	const int nSimu = 100;
	const int nRepet = 10;
	int nBus[nNBus] = { 10, 100, 500, 1000, 5000, 10000, 50000, 100000 };

	float values[nSimu];
	float values2[nSimu];
	
	float rhos[nSimu];
	MatrixCPU temps(nNBus, nSimu, 0);
	MatrixCPU nBusMat(1, nNBus, 0);
	for (int j = 0; j < nSimu; j++) {
		values[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
		rhos[j] = (float)(rand() % 100) / rand();
	}
	for (int i = 0; i < nNBus; i++) {
		std::cout << "iteration (" << i << ") nBus " << nBus[i] << std::endl;
		nBusMat.set(0, i, nBus[i]);
		int numBlocksB = ceil((nBus[i] + blockSize - 1) / blockSize);
		int sizeOPFTotal = 10 * nBus[i] - 3;
		MatrixGPU ChatCopy(6, nBus[i], 0, 1);
		MatrixGPU Bpt2Copy(2 * nBus[i], 1, 0, 1);
		MatrixGPU nAgentByBus(nBus[i], 1, 2, 1);
		MatrixGPU nChild(nBus[i], 1);
		MatrixGPU Ancestor(nBus[i], 1);
		MatrixGPU PosChild(nBus[i], 1);
		MatrixGPU Childs(nBus[i], 1);
		MatrixGPU indiceBusBegin(nBus[i], 1);
		MatrixGPU indiceChildBegin(nBus[i], 1);
		int debut;
		switch (repartition)
		{
		case 0: // tous sur 0
			
			nChild.set(0, 0, nBus[i] - 1);
			Ancestor.set(0, 0, -1);
			debut = 3 * (nBus[i] - 1) + 7;
			for (int b = 1; b < nBus[i]; b++) {
				PosChild.set(b, 0, b - 1);
				Childs.set(b, 0, b + 1);
				indiceChildBegin.set(b, 0, nBus[i]);
				indiceBusBegin.set(b, 0, debut);
				debut += 7;
			}
			break;
		case 1:
			// ligne
			debut = 0;
			for (int b = 0; b < nBus[i]; b++) {
				nChild.set(b, 0, 1);
				Ancestor.set(b, 0, b - 1);
				//PosChild.set(b, 0, 0);
				Childs.set(b, 0, b + 1);
				indiceBusBegin.set(b, 0, debut);
				indiceChildBegin.set(b, 0, b);
				debut += 10;
			}
			break;
		default:
			// n enfant par bus
			int nStep = log(nBus[i]) / log(repartition);
			int BusBegin = 0;
			debut = 0;
			Ancestor.set(0, 0, -1);
			for (int n = 0; n < nStep - 1; n++) {
				int nBusStep = pow(repartition, n);
				for (int ancestor = 0; ancestor < nBusStep; ancestor++) {
					int idAncestor = BusBegin + ancestor; // à changer ?
					nChild.set(idAncestor, 0, repartition);
					indiceBusBegin.set(idAncestor, 0, debut);
					indiceChildBegin.set(idAncestor, 0, idAncestor * repartition);
					debut += 3 * repartition + 7;
					for (int b = 0; b < repartition; b++) {
						int idBus = BusBegin + nBusStep + ancestor * repartition + b; // à changer ?
						Ancestor.set(idBus, 0, idAncestor);
						PosChild.set(idBus, 0, b);
					}
				}
				BusBegin += pow(repartition, n);
			}
			

			int nBusRestant = nBus[i] - (pow(repartition, nStep) - 1);
			int nBAncestor = nBusRestant / repartition;
			std::cout << BusBegin << " " << nBusRestant << " " << nBAncestor << std::endl;
			for (int ancestor = 0; ancestor < nBAncestor; ancestor++) {
				int idAncestor = BusBegin + ancestor; // à changer ?
				indiceBusBegin.set(idAncestor, 0, debut);
				nChild.set(idAncestor, 0, repartition);
				indiceChildBegin.set(idAncestor, 0, idAncestor * repartition);
				debut += 3 * repartition + 7;
				for (int b = 0; b < repartition; b++) {
					int idBus = pow(repartition, nStep) - 1 + ancestor * repartition + b; // à changer ?
					Ancestor.set(idBus, 0, idAncestor);
					PosChild.set(idBus, 0, b);
				}
			}
			
			for (int b = BusBegin + nBAncestor; b < nBus[i]; b++) { // le reste n'a pas d'enfant
				indiceBusBegin.set(b, 0, debut);
				debut += 7;
			}

			nBusRestant = nBusRestant % repartition;

			for (int b = 0; b < nBusRestant; b++) {
				int idBus = nBus[i] - nBusRestant + b; // à changer ?
				Ancestor.set(idBus, 0, 0);
				PosChild.set(idBus, 0, b + repartition);
			}

			for (int b = 0; b < nBus[i]; b++) {
				Childs.set(b, 0, b + 1);
			}
			break;
		}
		nChild.transferGPU();
		Ancestor.transferGPU();
		PosChild.transferGPU();
		Childs.transferGPU();
		indiceBusBegin.transferGPU();
		indiceChildBegin.transferGPU();
		

		for (int simu = 0; simu < nSimu; simu++) {
			MatrixGPU Y(sizeOPFTotal, 1, values[simu], 1);
			MatrixGPU Mu(sizeOPFTotal, 1, values2[simu], 1);
			float _rho = rhos[simu];
			cudaDeviceSynchronize();
			time = 0;
			for (int repet = 0; repet < nRepet; repet++) {
				MatrixGPU YCopy(Y);
				MatrixGPU MuCopy(Mu);
				MatrixGPU nChildCopy(nChild);
				MatrixGPU AncestorCopy(Ancestor);
				MatrixGPU PosChildCopy(PosChild);
				MatrixGPU ChildsCopy(Childs);
				MatrixGPU indiceBusBeginCopy(indiceBusBegin);
				MatrixGPU indiceChildBeginCopy(indiceChildBegin);
				
				int B = nBus[i];
				int numBlock = B;
				
				cudaDeviceSynchronize();
				switch (method)
				{
				case 0:
					a = std::chrono::high_resolution_clock::now();
					
					switch (blockSize) {
					case 512:
						updateChatGPUTest<512> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, _rho, B);
						break;
					case 256:
						updateChatGPUTest<256> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, _rho, B);
						break;
					case 128:
						updateChatGPUTest<128> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, _rho, B);
						break;
					case 64:
						updateChatGPUTest< 64> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, _rho, B);
						break;
					case 32:
						updateChatGPUTest< 32> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, _rho, B);
						break;
					case 16:
						updateChatGPUTest< 16> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, _rho, B);
						break;
					case  8:
						updateChatGPUTest<  8> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, _rho, B);
						break;
					case  4:
						updateChatGPUTest<  4> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, _rho, B);
						break;
					case  2:
						updateChatGPUTest<  2> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, _rho, B);
						break;
					case  1:
						updateChatGPUTest<  1> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, _rho, B);
						break;
					}

					updateBpt2Test << < numBlocksB, blockSize >> > (Bpt2Copy._matrixGPU, ChatCopy._matrixGPU, nAgentByBus._matrixGPU, B);

					
					cudaDeviceSynchronize();
					b = std::chrono::high_resolution_clock::now();
					break;
				case 1:
					a = std::chrono::high_resolution_clock::now();
					switch (blockSize) {
					case 512:
						updateChatBpt2Test<512> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, Bpt2Copy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, nAgentByBus._matrixGPU, _rho, B);
						break;
					case 256:
						updateChatBpt2Test<256> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, Bpt2Copy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, nAgentByBus._matrixGPU, _rho, B);
						break;
					case 128:
						updateChatBpt2Test<128> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, Bpt2Copy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, nAgentByBus._matrixGPU, _rho, B);
						break;
					case 64:
						updateChatBpt2Test< 64> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, Bpt2Copy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, nAgentByBus._matrixGPU, _rho, B);
						break;
					case 32:
						updateChatBpt2Test< 32> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, Bpt2Copy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, nAgentByBus._matrixGPU, _rho, B);
						break;
					case 16:
						updateChatBpt2Test< 16> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, Bpt2Copy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, nAgentByBus._matrixGPU, _rho, B);
						break;
					case  8:
						updateChatBpt2Test<  8> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, Bpt2Copy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, nAgentByBus._matrixGPU, _rho, B);
						break;
					case  4:
						updateChatBpt2Test<  4> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, Bpt2Copy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, nAgentByBus._matrixGPU, _rho, B);
						break;
					case  2:
						updateChatBpt2Test<  2> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, Bpt2Copy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, nAgentByBus._matrixGPU, _rho, B);
						break;
					case  1:
						updateChatBpt2Test<  1> << <numBlock, blockSize >> > (ChatCopy._matrixGPU, Bpt2Copy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, nAgentByBus._matrixGPU, _rho, B);
						break;
					}
					cudaDeviceSynchronize();
					b = std::chrono::high_resolution_clock::now();
					break;
				case 2:
					a = std::chrono::high_resolution_clock::now();
					updateChatBpt2OneDoAll <<< numBlocksB, blockSize >> > (ChatCopy._matrixGPU, Bpt2Copy._matrixGPU, YCopy._matrixGPU, MuCopy._matrixGPU, nChildCopy._matrixGPU, AncestorCopy._matrixGPU, PosChildCopy._matrixGPU, ChildsCopy._matrixGPU, indiceBusBeginCopy._matrixGPU, indiceChildBeginCopy._matrixGPU, nAgentByBus._matrixGPU, _rho, B);
					cudaDeviceSynchronize();
					b = std::chrono::high_resolution_clock::now();
					break;
				default:
					return 0;
					break;
				}
				time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
			}
			temps.set(i, simu, (float)time / nRepet);
		}
	}
	nBusMat.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	return temps.sum();
}


template <unsigned int blockSize>
__global__ void updateChatGPUTest(float* Chat,  float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float _rho, int nBus) {

	int bus = blockIdx.x;
	int index = threadIdx.x;
	int step = blockDim.x;

	__shared__ float shArr[blockSize]; // c'est grand pour pas grand chose...


	int indice = indiceBusBegin[bus];
	int indiceChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
	int nb = nChild[bus];
	int Ai = Ancestor[bus];
	int c = posChild[bus];
	float var = 0;

	if (index < 6) {
		//float Phat, Qhat, lhat, phat, qhat;
		var = Y[indice + index] / 2 - MU[indice + index] / (2 * _rho);
		if (bus > 0) {
			if (index < 3) {
				int nAi = nChild[Ai];
				int indiceAncBus = indiceBusBegin[Ai] + 7 + nAi * index + c;
				//var = indiceAncBus;
				var += Y[indiceAncBus] / 2 - MU[indiceAncBus] / (2 * _rho);
			}
		}
	}
	float vhat = 0;
	float muhat = 0;
	for (int i = index; i < nb; i += step) {
		int c = Childs[indiceChild + i];
		int indiceBusChild = indiceBusBegin[c];
		muhat += MU[indiceBusChild + 6]; // pas du tout coalescent
		vhat += Y[indiceBusChild + 6]; // pas du tout coalescent
	}
	shArr[index] = vhat / (nb + 1) - muhat / (_rho * (nb + 1));
	__syncthreads();
	for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
		if (index < size) {
			shArr[index] += shArr[index + size];
		}
		__syncthreads();
	}

	if (index < 6) {
		if (index == 3) {
			var = shArr[0] + Y[indice + 3] / (nb + 1) - MU[indice + 3] / (_rho * (nb + 1)); //shArr[0];
		}
		Chat[index * nBus + bus] = var; // pas coalescent mais bon perdu pour perdu
	}
}
__global__ void updateBpt2Test(float* Bpt2, float* Chat, float* nAgentByBus, int nBus) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;


	for (int b = index; b < nBus; b += step) {
		int nA = nAgentByBus[b];

		Bpt2[b] = nA > 0 ? 2 * Chat[b + 4 * nBus] / nA : 0; // �criture coalescente et lecture coalescente
		Bpt2[b + nBus] = nA > 0 ? 2 * Chat[b + 5 * nBus] / nA : 0;

	}

}

template <unsigned int blockSize>
__global__ void updateChatBpt2Test(float* Chat, float* Bpt2, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* nAgentByBus, float _rho, int nBus) {
	int bus = blockIdx.x;
	int index = threadIdx.x;
	int step = blockDim.x;

	__shared__ float shArr[blockSize]; // c'est grand pour pas grand chose...


	int indice = indiceBusBegin[bus];
	int indiceChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
	int nb = nChild[bus];
	int Ai = Ancestor[bus];
	int c = posChild[bus];
	int nA = nAgentByBus[bus];
	float var = 0;

	if (index < 6) {
		//float Phat, Qhat, lhat, phat, qhat;
		var = Y[indice + index] / 2 - MU[indice + index] / (2 * _rho);
		if (bus > 0) {
			if (index < 3) {
				int nAi = nChild[Ai];
				int indiceAncBus = indiceBusBegin[Ai] + 7 + nAi * index + c;
				//var = indiceAncBus;
				var += Y[indiceAncBus] / 2 - MU[indiceAncBus] / (2 * _rho);
			}
		}
	}
	float vhat = 0;
	float muhat = 0;
	for (int i = index; i < nb; i += step) {
		int c = Childs[indiceChild + i];
		int indiceBusChild = indiceBusBegin[c];
		muhat += MU[indiceBusChild + 6]; // pas du tout coalescent
		vhat += Y[indiceBusChild + 6]; // pas du tout coalescent
	}
	shArr[index] = vhat / (nb + 1) - muhat / (_rho * (nb + 1));
	__syncthreads();
	for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
		if (index < size) {
			shArr[index] += shArr[index + size];
		}
		__syncthreads();
	}

	if (index < 6) {
		if (index == 3) {
			var = shArr[0] + Y[indice + 3] / (nb + 1) - MU[indice + 3] / (_rho * (nb + 1)); //shArr[0];
		}
		Chat[index * nBus + bus] = var; // pas coalescent mais bon perdu pour perdu
		if (index == 4) {
			Bpt2[bus] = nA > 0 ? 2 * var / nA : 0; // �criture coalescente et lecture coalescente
		} 
		if (index == 5) {
			Bpt2[bus + nBus] = nA > 0 ? 2 * var / nA : 0;
		}
	}
}


__global__ void updateChatBpt2OneDoAll(float* Chat, float* Bpt2, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* nAgentByBus, float _rho, int nBus) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;


	for (int bus = index; bus < nBus; bus += step) {
		int indice = indiceBusBegin[bus];
		int indiceChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
		int nb = nChild[bus];
		int Ai = Ancestor[bus];
		int c = posChild[bus];
		int nA = nAgentByBus[bus];
		float Phat, Qhat, lhat, phat, qhat, vhat, muhat;
		
		Phat = Y[indice] / 2     - MU[indice] / (2 * _rho);
		Qhat = Y[indice + 1] / 2 - MU[indice + 1] / (2 * _rho);
		lhat = Y[indice + 2] / 2 - MU[indice + 2] / (2 * _rho);
		phat = Y[indice + 4] / 2 - MU[indice + 4] / (2 * _rho);
		qhat = Y[indice + 5] / 2 - MU[indice + 5] / (2 * _rho);
		vhat = 0;
		muhat = 0;
		if (bus > 0) {
			int nAi = nChild[Ai];
			int indiceAncBus = indiceBusBegin[Ai] + 7 + c;
			Phat += Y[indiceAncBus] / 2 - MU[indiceAncBus] / (2 * _rho);
			Qhat += Y[indiceAncBus + nAi] / 2 - MU[indiceAncBus + nAi] / (2 * _rho);
			lhat += Y[indiceAncBus + 2 * nAi] / 2 - MU[indiceAncBus + 2 * nAi] / (2 * _rho);
		}
		for (int i = 0; i < nb; i++) {
			int c = Childs[indiceChild + i];
			int indiceBusChild = indiceBusBegin[c];
			muhat += MU[indiceBusChild + 6]; // pas du tout coalescent
			vhat += Y[indiceBusChild + 6]; // pas du tout coalescent
		}
		vhat = (vhat + Y[indice + 3]) / (nb + 1) - (muhat + MU[indice + 3]) / (_rho * (nb + 1));
		
		// ecriture coalescente 
		Chat[bus] = Phat;
		Chat[1 * nBus + bus] = Qhat;
		Chat[2 * nBus + bus] = lhat;
		Chat[3 * nBus + bus] = vhat;
		Chat[4 * nBus + bus] = phat;
		Chat[5 * nBus + bus] = qhat;
		Bpt2[bus] = nA > 0 ? 2 * phat / nA : 0; 
		Bpt2[bus + nBus] = nA > 0 ? 2 * qhat / nA : 0;
	}

}


// --------------------------------- probleme local -------------------------------------------------------

float testCalculPnShared(int method, int blockSize, int repartition) {
	std::string fileName = "TempsCalculPnShared" + std::to_string(method) + "_" + std::to_string(blockSize) + "_" + std::to_string(repartition) + ".csv";
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNBus = 7;
	const int nNAgent = 8;
	const int nSimu = 100;
	const int nRepet = 10;
	int nBus[nNBus] = { 20, 100, 500, 1000, 5000, 10000, 50000 };
	float nAgent[nNAgent] = { 0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10 };
	float epsL = -1; 
	int nIterL = 1000; // pour pas que le nombre d'itération varie en fonction de l'aléatoire, mais induit un biais, un bus avec moins d'agents devrait normalement 
	// nécessiter moins d'itérations. Donc pas généralisable au cas où on cherche la convergence...

	float values[nSimu];
	float values2[nSimu];
	float values3[nSimu];
	float rhos[nSimu];
	MatrixCPU temps(nNBus * nNAgent, nSimu, 0);
	MatrixCPU nBusMat(1, nNBus, 0);
	MatrixCPU nAgentMat(1, nNAgent);
	for (int j = 0; j < nSimu; j++) {
		values[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
		values3[j] = (float)(rand()) / rand();
		rhos[j] = (float)(rand() % 100) / rand();
	}
	for (int i = 0; i < nNBus; i++) {
		nBusMat.set(0, i, nBus[i]);
		for (int j = 0; j < nNAgent; j++) {
			int Agent = nAgent[j] * nBus[i];
			std::cout << "iteration (" << i << "," << j << ") nBus " << nBus[i] << " agent " << Agent << std::endl;
			nAgentMat.set(0, j, nAgent[j]);
			
			MatrixGPU nAgentByBus(nBus[i], 1);
			MatrixGPU CoresSoloBusAgent(nBus[i], 1, -1);
			MatrixGPU CoresAgentBus(Agent, 1);
			MatrixGPU CoresAgentBusBegin(nBus[i], 1);

		
			

			int nAgentByBusQ = Agent / nBus[i];
			int nAgentByBusR = Agent % nBus[i];

			int debut = 0;
			switch (repartition)
			{
			case 0: // tous sur 0
				nAgentByBus.set(0, 0, Agent);
				
				for (int b = 1; b < nBus[i]; b++) {
					CoresAgentBusBegin.set(b, 0, Agent);
				}

				for (int n = 0; n < Agent; n++) {
					CoresAgentBus.set(n, 0, Agent - 1 - n);
				}
			
				break;
			case 1: // equilibré
				std::cout << nAgentByBusQ << " " << nAgentByBusR << std::endl;
				for (int b = 0; b < nBus[i]; b++) {
					int nA = nAgentByBusQ + 1 * (nAgentByBusR > b);
					nAgentByBus.set(b, 0, nA);
					CoresAgentBusBegin.set(b, 0, debut);
					
					if (nA == 1) {
						CoresSoloBusAgent.set(b, 0, debut);
					}
					debut += nA;

					
				}
				for (int n = 0; n < Agent; n++) {
					CoresAgentBus.set(n, 0, Agent - 1 - n);
				}
				break;
			default: // random
			
				break;
			}
			CoresSoloBusAgent.transferGPU();
			CoresAgentBus.transferGPU();
			CoresAgentBusBegin.transferGPU();
			nAgentByBus.transferGPU();

			for (int simu = 0; simu < nSimu; simu++) {
				
				MatrixGPU Ap2(2 * Agent, 1, 0.1, 1);
				MatrixGPU Apt2(2 * Agent, 1, values2[simu], 1);
				MatrixGPU Cp(2 * Agent, 1, values3[simu], 1);
				MatrixGPU Pmin(2 * Agent, 1, values[simu], 1);
				MatrixGPU Pmax(2 * Agent, 1, values[simu] + 1, 1);
				MatrixGPU Apt1(2 * Agent, 1, values3[simu] - values2[simu] , 1);
				MatrixGPU Bpt2(2 * Agent, 1, values3[simu] + values2[simu], 1);
				float _rhol = rhos[simu];
				
				time = 0;
				for (int repet = 0; repet < nRepet; repet++) {
					MatrixGPU Ap2Copy(Ap2);
					MatrixGPU CpCopy(Cp);
					MatrixGPU PminCopy(Pmin);
					MatrixGPU PmaxCopy(Pmax);
					MatrixGPU Apt1Copy(Apt1);
					MatrixGPU Apt2Copy(Apt2);
					MatrixGPU Bpt2Copy(Bpt2);

					MatrixGPU nAgentByBusCopy(nAgentByBus);
					MatrixGPU CoresSoloBusAgentCopy(CoresSoloBusAgent);
					MatrixGPU CoresAgentBusCopy(CoresAgentBus);
					MatrixGPU CoresAgentBusBeginCopy(CoresAgentBusBegin);

					MatrixGPU PnCopy(2 * Agent, 1, 0, 1);
					MatrixGPU PnPreCopy(2 * Agent, 1, 0, 1);
					MatrixGPU PnMoyCopy(2 * nBus[i], 1, 0, 1);
					MatrixGPU PnTildeCopy(2 * nBus[i], 1, 0, 1);
					MatrixGPU MuLCopy(2 * nBus[i], 1, 0, 1);


					int B = nBus[i];
					int numBlocks = B;

					cudaDeviceSynchronize();
					switch (method)
					{
					case 0:
						a = std::chrono::high_resolution_clock::now();
					
						switch (blockSize) {
						case 512:
							updatePnPGPUSharedResidualTest<512> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case 256:
							updatePnPGPUSharedResidualTest<256> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);							break;
						case 128:
							updatePnPGPUSharedResidualTest<128> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);						break;
						case 64:
							updatePnPGPUSharedResidualTest< 64> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case 32:
							updatePnPGPUSharedResidualTest< 32> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case 16:
							updatePnPGPUSharedResidualTest< 16> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case  8:
							updatePnPGPUSharedResidualTest<  8> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case  4:
							updatePnPGPUSharedResidualTest<  4> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case  2:
							updatePnPGPUSharedResidualTest<  2> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case  1:
							updatePnPGPUSharedResidualTest<  1> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						}


						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 1:
						a = std::chrono::high_resolution_clock::now();
						switch (blockSize) {
						case 512:
							updatePnPGPUSharedResidualSameThreadTest<512> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case 256:
							updatePnPGPUSharedResidualSameThreadTest<256> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);							break;
						case 128:
							updatePnPGPUSharedResidualSameThreadTest<128> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);						break;
						case 64:
							updatePnPGPUSharedResidualSameThreadTest< 64> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case 32:
							updatePnPGPUSharedResidualSameThreadTest< 32> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case 16:
							updatePnPGPUSharedResidualSameThreadTest< 16> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case  8:
							updatePnPGPUSharedResidualSameThreadTest<  8> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case  4:
							updatePnPGPUSharedResidualSameThreadTest<  4> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case  2:
							updatePnPGPUSharedResidualSameThreadTest<  2> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						case  1:
							updatePnPGPUSharedResidualSameThreadTest<  1> << <numBlocks, blockSize >> > (PnCopy._matrixGPU, PnPreCopy._matrixGPU, PnMoyCopy._matrixGPU, PnTildeCopy._matrixGPU, MuLCopy._matrixGPU, nAgentByBusCopy._matrixGPU, _rhol, Ap2Copy._matrixGPU, CpCopy._matrixGPU,
								PminCopy._matrixGPU, PmaxCopy._matrixGPU, Apt1Copy._matrixGPU, Apt2Copy._matrixGPU, Bpt2Copy._matrixGPU, CoresSoloBusAgentCopy._matrixGPU, CoresAgentBusCopy._matrixGPU, CoresAgentBusBeginCopy._matrixGPU, epsL, nIterL, Agent, B);
							break;
						}
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					default:
						return 0;
						break;
					}
					time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
				}
				temps.set(i * nNAgent + j, simu, (float)time / nRepet);
			}
		}		
	}
	nBusMat.saveCSV(fileName, mode);
	nAgentMat.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	return temps.sum();
}



template <unsigned int blockSize>
__global__ void updatePnPGPUSharedResidualTest(float* Pn, float* PnPre, float* PnMoy, float* PnTilde, float* MUL, float* nAgentByBus, float _rhol, float* Ap2, float* Cp, float* Pmin,
	float* Pmax, float* Apt1, float* Apt2, float* Bpt2, float* CoresSoloBusAgent, float* CoresBusAgent, float* CoresBusAgentBegin, float eps, int nIterLMax, int nAgent, int nBus) {

	//Definition de toutes les variables locales
	int i = blockIdx.x; // c'est aussi l'identifiant du bus !
	unsigned int thIdx = threadIdx.x;

	// ne change pas

	float Ap2local;
	float Ap12local;
	float Cplocal;
	float Pminlocal;
	float Pmaxlocal;
	float Pnlocal; // change
	float Pnprelocal; // change

	float bpt, MULOCAL, moy, p;
	float m, r, ub, lb, t;
	// le changement doit �tre partag� par tous les threads du bloc

	__shared__ float MuShared[2];
	__shared__ float PnMoyShared[2];
	__shared__ float PnTildeShared[2];
	__shared__ bool mustContinue;

	// constant et commun � tous les thread d'un bloc
	__shared__ float Apt1Shared;
	__shared__ float Apt2Shared;
	__shared__ float Apt12Shared;
	__shared__ float Bpt2Shared[2];
	__shared__ int nAgentShared;
	__shared__ float at1Shared;

	__shared__ float shArrP[blockSize / 2 + 1];
	__shared__ float shArrQ[blockSize / 2 + 1];

	if (thIdx < (blockSize / 2 + 1)) {
		shArrP[thIdx] = 0;
		shArrQ[thIdx] = 0;
	}
	if (thIdx < 1) {
		Apt1Shared = Apt1[i]; // rho_l *Ni, m�me pour les 2
		Apt2Shared = Apt2[i]; // rho * Ni^2, m�me pour les 2
		Apt12Shared = Apt1Shared + Apt2Shared; // m�me pour les 2
		nAgentShared = nAgentByBus[i];
		at1Shared = _rhol;
		mustContinue = false;
	}


	if (thIdx < 2) {
		Bpt2Shared[thIdx] = Bpt2[i + nBus * thIdx];
		MuShared[thIdx] = MUL[i + nBus * thIdx];
		PnMoyShared[thIdx] = PnMoy[i + nBus * thIdx];
		PnTildeShared[thIdx] = PnTilde[i + nBus * thIdx];
	}
	__syncthreads();

	int iter = 0;
	if (nAgentShared > 0) { // sinon il n'y a rien � faire
		int indicePorQ = thIdx / nAgentShared; // 0 or 1
		const int CoresAgentLinLocal = CoresBusAgentBegin[i];
		const int j = CoresAgentLinLocal + thIdx;
		double res = 0;
		if (nAgentShared == 1) { // cas trivial s'il n'y a qu'un agent, la divergence est entre les blocs donc c'est ok
			if (thIdx < 2) {
				int agent = CoresSoloBusAgent[i];
				// Cplocal et Ap12local, Pmaxlocal, Pminlocal � definir
				Cplocal = Cp[agent + thIdx * nAgent];
				Ap2local = Ap2[agent + thIdx * nAgent];
				ub = Pmax[agent + thIdx * nAgent];
				lb = Pmin[agent + thIdx * nAgent];
				r = (Apt2Shared * Bpt2Shared[thIdx] - Cplocal) / (Apt2Shared + Ap2local); //pn = (_rho * Bpt2.get(b, 0) - Cost2.get(n, 0)) / ((_rho + Cost1.get(n, 0)));
				t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
				Pnlocal = t;
				Pnprelocal = t;
				PnMoyShared[thIdx] = t;
				PnTildeShared[thIdx] = t;// PnMoyShared[thIdx];
			}
		}
		else {
			Pnlocal = 0;
			if (thIdx < 2* nAgentShared)
			{
				int indice = j < nAgentShared ? j : j - nAgentShared;
				int agent = CoresBusAgent[indice] + indicePorQ * nAgent;
				// P & Q
				Ap2local = Ap2[agent];
				Ap12local = Ap2local + _rhol;
				Cplocal = Cp[agent];
				Pminlocal = Pmin[agent];
				Pmaxlocal = Pmax[agent];
				Pnlocal = Pn[agent];
			}

			//Calcul des it�rations

			for (iter = 0; iter < nIterLMax; iter++) {
				__syncthreads();
				if (thIdx < 2 * nAgentShared) {
					MULOCAL = MuShared[indicePorQ]; // il y a 2 valeurs c'est tr�s chiant, ce n'est pas broadcast
					moy = PnMoyShared[indicePorQ]; // avec de la chance cela fait 2 broadcasts
					p = PnTildeShared[indicePorQ];
					// P & Q
					Pnprelocal = Pnlocal;
					m = Pnlocal - moy + p - MULOCAL; // Pn.get(i, 0) - PnMoy.get(bus, 0) + PnTilde.get(bus, 0) - MuL.get(bus, 0);
					r = (m * at1Shared - Cplocal) / Ap12local; // pn = (Bp1.get(n, 0) * _rhol - Cost2.get(n, 0)) / Ap12.get(n, 0);
					ub = Pmaxlocal;
					lb = Pminlocal;
					t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
					Pnlocal = t;
					res = (double)t - Pnprelocal;
					res = res * res;
					if (res > eps) {
						mustContinue = true; // pas de race condition, car l'ordre n'importe pas,
					}
				}
				if (thIdx < nAgentShared) {
					shArrP[thIdx] = Pnlocal;
				}
				else if (thIdx < 2*nAgentShared)
				{
					shArrQ[thIdx - nAgentShared] = Pnlocal;
				}

				__syncthreads();
				if (blockSize >= 512) {
					if (thIdx < 128) {
						shArrP[thIdx] += shArrP[thIdx + 128];
					}
					else if (thIdx < 256)
					{
						shArrQ[thIdx - 128] += shArrQ[thIdx];
					}
					__syncthreads();
				}
				if (blockSize >= 256) {
					if (thIdx < 64) {
						shArrP[thIdx] += shArrP[thIdx + 64];
					}
					else if (thIdx < 128)
					{
						shArrQ[thIdx - 64] += shArrQ[thIdx];
					}
					__syncthreads();
				}
				if (blockSize >= 128) {
					if (thIdx < 32) {
						warpReduce<blockSize>(shArrP, thIdx);
					}
					else if (thIdx < 64)
					{
						warpReduce<blockSize>(shArrQ, thIdx - 32);
					}

				}
				__syncthreads();

				if (thIdx < 2) {
					moy = (shArrP[0] + thIdx * (shArrQ[0] - shArrP[0])) / nAgentShared;
					PnMoyShared[thIdx] = moy;
					bpt = moy + MuShared[thIdx]; //Bpt1.set(b, 0, MuL.get(b, 0) + PnMoy.get(b, 0));
					p = (Apt1Shared * bpt + Apt2Shared * Bpt2Shared[thIdx]) / Apt12Shared; //pn = (Bpt1.get(b, 0) * Apt1.get(b, 0) + Bpt2.get(b, 0) * Apt2.get(b, 0)) / Apt12.get(b, 0);
					PnTildeShared[thIdx] = p;
					res = p - moy;
					res = res * res;
					if (res > eps) {
						mustContinue = true;
					}
					MuShared[thIdx] = MuShared[thIdx] + moy - p; // mu = MuL.get(b, 0) + PnMoy.get(b, 0) - PnTilde.get(b, 0);
				}
				__syncthreads();
				if (!mustContinue) {
					break;
				}
				else {
					__syncthreads();
					if (thIdx == 0) {
						mustContinue = false;
					}
				}
			}
		}
		//Ecriture des it�rations
		__syncthreads();

		if (thIdx < 2*nAgentShared)
		{
			int indice = j - indicePorQ * nAgentShared;
			int agent = CoresBusAgent[indice] + indicePorQ * nAgent;

			Pn[agent] = Pnlocal;
			PnPre[agent] = Pnprelocal;

		}
		if (thIdx < 2) {
			PnMoy[blockIdx.x + thIdx * nBus] = PnMoyShared[thIdx];// TMoyShared;
			PnTilde[blockIdx.x + thIdx * nBus] = PnTildeShared[thIdx];// PShared;
			MUL[blockIdx.x + thIdx * nBus] = MuShared[thIdx];// MuShared;
		}
	}
}



template <unsigned int blockSize>
__global__ void updatePnPGPUSharedResidualSameThreadTest(float* Pn, float* PnPre, float* PnMoy, float* PnTilde, float* MUL, float* nAgentByBus, float _rhol, float* Ap2, float* Cp, float* Pmin,
	float* Pmax, float* Apt1, float* Apt2, float* Bpt2, float* CoresSoloBusAgent, float* CoresBusAgent, float* CoresBusAgentBegin, float eps, int nIterLMax, int nAgent, int nBus) {

	//Definition de toutes les variables locales
	int i = blockIdx.x; // c'est aussi l'identifiant du bus !
	unsigned int thIdx = threadIdx.x;

	// ne change pas

	float Ap2local[2];
	float Ap12local[2];
	float Cplocal[2];
	float Pminlocal[2];
	float Pmaxlocal[2];
	float Pnlocal[2]; // change
	float Pnprelocal[2]; // change

	float bpt, MULOCAL, moy, p;
	float m, r, ub, lb, t;
	// le changement doit �tre partag� par tous les threads du bloc

	__shared__ float MuShared[2];
	__shared__ float PnMoyShared[2];
	__shared__ float PnTildeShared[2];
	__shared__ bool mustContinue;

	// constant et commun � tous les thread d'un bloc
	__shared__ float Apt1Shared;
	__shared__ float Apt2Shared;
	__shared__ float Apt12Shared;
	__shared__ float Bpt2Shared[2];
	__shared__ int nAgentShared;
	__shared__ float at1Shared;

	__shared__ float shArrP[blockSize / 2 + 1];
	__shared__ float shArrQ[blockSize / 2 + 1];

	if (thIdx < (blockSize / 2 + 1)) {
		shArrP[thIdx] = 0;
		shArrQ[thIdx] = 0;
	}
	if (thIdx < 1) {
		Apt1Shared = Apt1[i]; // rho_l *Ni, m�me pour les 2
		Apt2Shared = Apt2[i]; // rho * Ni^2, m�me pour les 2
		Apt12Shared = Apt1Shared + Apt2Shared; // m�me pour les 2
		nAgentShared = nAgentByBus[i];
		at1Shared = _rhol;
		mustContinue = false;
	}


	if (thIdx < 2) {
		Bpt2Shared[thIdx] = Bpt2[i + nBus * thIdx];
		MuShared[thIdx] = MUL[i + nBus * thIdx];
		PnMoyShared[thIdx] = PnMoy[i + nBus * thIdx];
		PnTildeShared[thIdx] = PnTilde[i + nBus * thIdx];
	}
	__syncthreads();

	int iter = 0;
	if (nAgentShared > 0) { // sinon il n'y a rien � faire
		const int CoresAgentLinLocal = CoresBusAgentBegin[i];
		const int j = CoresAgentLinLocal + thIdx;
		double res = 0;
		if (nAgentShared == 1) { // cas trivial s'il n'y a qu'un agent, la divergence est entre les blocs donc c'est ok
			if (thIdx == 0) {
				int agent = CoresSoloBusAgent[i];
				// Cplocal et Ap12local, Pmaxlocal, Pminlocal � definir
				Cplocal[0] = Cp[agent];
				Ap2local[0] = Ap2[agent];
				ub = Pmax[agent];
				lb = Pmin[agent];
				r = (Apt2Shared * Bpt2Shared[0] - Cplocal[0]) / (Apt2Shared + Ap2local[0]); //pn = (_rho * Bpt2.get(b, 0) - Cost2.get(n, 0)) / ((_rho + Cost1.get(n, 0)));
				t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
				Pnlocal[0] = t;
				Pnprelocal[0] = t;
				PnMoyShared[0] = t;
				PnTildeShared[0] = t;// PnMoyShared[thIdx];

				// Q 
				Cplocal[1] = Cp[agent + nAgent];
				Ap2local[1] = Ap2[agent + nAgent];
				ub = Pmax[agent + nAgent];
				lb = Pmin[agent + nAgent];
				r = (Apt2Shared * Bpt2Shared[1] - Cplocal[1]) / (Apt2Shared + Ap2local[1]); //pn = (_rho * Bpt2.get(b, 0) - Cost2.get(n, 0)) / ((_rho + Cost1.get(n, 0)));
				t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
				Pnlocal[1] = t;
				Pnprelocal[1] = t;
				PnMoyShared[1] = t;
				PnTildeShared[1] = t;// PnMoyShared[thIdx];
			}
		}
		else {
			Pnlocal[0] = 0;
			Pnlocal[1] = 0;
			if (thIdx < nAgentShared)
			{
				int agent = CoresBusAgent[j];
				// P & Q
				Ap2local[0] = Ap2[agent];
				Ap12local[0] = Ap2local[0] + _rhol;
				Cplocal[0] = Cp[agent];
				Pminlocal[0] = Pmin[agent];
				Pmaxlocal[0] = Pmax[agent];
				Pnlocal[0] = Pn[agent];

				Ap2local[1] = Ap2[agent + nAgent];
				Ap12local[1] = Ap2local[1] + _rhol;
				Cplocal[1] = Cp[agent + nAgent];
				Pminlocal[1] = Pmin[agent + nAgent];
				Pmaxlocal[1] = Pmax[agent + nAgent];
				Pnlocal[1] = Pn[agent + nAgent];
			}

			//Calcul des it�rations

			for (iter = 0; iter < nIterLMax; iter++) {
				__syncthreads();
				if (thIdx < nAgentShared) {
					// P
					MULOCAL = MuShared[0]; 
					moy = PnMoyShared[0]; 
					p = PnTildeShared[0];
					
					Pnprelocal[0] = Pnlocal[0];
					m = Pnlocal[0] - moy + p - MULOCAL; // Pn.get(i, 0) - PnMoy.get(bus, 0) + PnTilde.get(bus, 0) - MuL.get(bus, 0);
					r = (m * at1Shared - Cplocal[0]) / Ap12local[0]; // pn = (Bp1.get(n, 0) * _rhol - Cost2.get(n, 0)) / Ap12.get(n, 0);
					ub = Pmaxlocal[0];
					lb = Pminlocal[0];
					t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
					Pnlocal[0] = t;
					res = (double)t - Pnprelocal[0];
					res = res * res;
					if (res > eps) {
						mustContinue = true; // pas de race condition, car l'ordre n'importe pas,
					}
					// Q
					MULOCAL = MuShared[1];
					moy = PnMoyShared[1];
					p = PnTildeShared[1];

					Pnprelocal[1] = Pnlocal[1];
					m = Pnlocal[1] - moy + p - MULOCAL; // Pn.get(i, 0) - PnMoy.get(bus, 0) + PnTilde.get(bus, 0) - MuL.get(bus, 0);
					r = (m * at1Shared - Cplocal[1]) / Ap12local[1]; // pn = (Bp1.get(n, 0) * _rhol - Cost2.get(n, 0)) / Ap12.get(n, 0);
					ub = Pmaxlocal[1];
					lb = Pminlocal[1];
					t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
					Pnlocal[1] = t;
					res = (double)t - Pnprelocal[1];
					res = res * res;
					if (res > eps) {
						mustContinue = true; // pas de race condition, car l'ordre n'importe pas,
					}

				}
				
				shArrP[thIdx] = Pnlocal[0];
				shArrQ[thIdx] = Pnlocal[1];
				
				
				__syncthreads();
				if (blockSize >= 512) {
					if (thIdx < 256) {
						shArrP[thIdx] += shArrP[thIdx + 256];
						shArrQ[thIdx] += shArrQ[thIdx + 256];
					}
					__syncthreads();
				}
				if (blockSize >= 256) {
					if (thIdx < 128) {
						shArrP[thIdx] += shArrP[thIdx + 128];
						shArrQ[thIdx] += shArrQ[thIdx + 128];
					}
					__syncthreads();
				}
				if (blockSize >= 128) {
					if (thIdx < 64) {
						shArrP[thIdx] += shArrP[thIdx + 64];
						shArrQ[thIdx] += shArrQ[thIdx + 64];
					}
					__syncthreads();
				}
				if (blockSize >= 64) {
					if (thIdx < 32) {
						warpReduce<blockSize>(shArrP, thIdx);
						warpReduce<blockSize>(shArrQ, thIdx);
					}
				}
				__syncthreads();

				if (thIdx == 0) {
					// P
					moy = shArrP[0] / nAgentShared;
					PnMoyShared[0] = moy;
					bpt = moy + MuShared[0]; //Bpt1.set(b, 0, MuL.get(b, 0) + PnMoy.get(b, 0));
					p = (Apt1Shared * bpt + Apt2Shared * Bpt2Shared[0]) / Apt12Shared; //pn = (Bpt1.get(b, 0) * Apt1.get(b, 0) + Bpt2.get(b, 0) * Apt2.get(b, 0)) / Apt12.get(b, 0);
					PnTildeShared[0] = p;
					res = p - moy;
					res = res * res;
					if (res > eps) {
						mustContinue = true;
					}
					MuShared[0] = MuShared[0] + moy - p; // mu = MuL.get(b, 0) + PnMoy.get(b, 0) - PnTilde.get(b, 0);
					// Q
					moy = shArrQ[0] / nAgentShared;
					PnMoyShared[1] = moy;
					bpt = moy + MuShared[1]; //Bpt1.set(b, 0, MuL.get(b, 0) + PnMoy.get(b, 0));
					p = (Apt1Shared * bpt + Apt2Shared * Bpt2Shared[1]) / Apt12Shared; //pn = (Bpt1.get(b, 0) * Apt1.get(b, 0) + Bpt2.get(b, 0) * Apt2.get(b, 0)) / Apt12.get(b, 0);
					PnTildeShared[1] = p;
					res = p - moy;
					res = res * res;
					if (res > eps) {
						mustContinue = true;
					}
					MuShared[1] = MuShared[1] + moy - p; // mu = MuL.get(b, 0) + PnMoy.get(b, 0) - PnTilde.get(b, 0);

				}
				__syncthreads();
				if (!mustContinue) {
					break;
				}
				else {
					__syncthreads();
					if (thIdx == 0) {
						mustContinue = false;
					}
				}
			}
		}
		//Ecriture des it�rations
		__syncthreads();

		if (thIdx < nAgentShared)
		{
			int agent = CoresBusAgent[j];

			Pn[agent] = Pnlocal[0];
			PnPre[agent] = Pnprelocal[0];
			Pn[agent + nAgent] = Pnlocal[1];
			PnPre[agent + nAgent] = Pnprelocal[1];

		}
		if (thIdx == 0) {
			PnMoy[blockIdx.x] = PnMoyShared[0];// TMoyShared;
			PnTilde[blockIdx.x] = PnTildeShared[0];// PShared;
			MUL[blockIdx.x] = MuShared[0];// MuShared;
			PnMoy[blockIdx.x + nBus] = PnMoyShared[1];// TMoyShared;
			PnTilde[blockIdx.x + nBus] = PnTildeShared[1];// PShared;
			MUL[blockIdx.x + nBus] = MuShared[1];// MuShared;
		}
	}
}


// --------------------------------- Voltage GS -------------------------------------------------------


float testCalculVGS(int method) {
	std::string fileName = "TempsVGS" + std::to_string(method) + ".csv";
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNBus = 3;
	const int nNLine = 5;
	const int nSimu = 100;
	const int nRepet = 10;
	int nBus[nNBus] = { 10, 100, 500};
	// nombre ligne : radial, 1/2 milieu,  milieu , 3/2 milieu  full connecté
	int nLineMax = 70000; // pouréviter les badalloc
	int blockSize = 256;
	float values[nSimu];
	float values2[nSimu];

	MatrixCPU temps(nNBus * nNLine, nSimu, 0);
	MatrixCPU nBusMat(1, nNBus);
	MatrixCPU nLineMat(1, nNLine * nNBus);
	for (int j = 0; j < nSimu; j++) {
		values[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
	}
	int indice = 0;
	for (int i = 0; i < nNBus; i++) {
		nBusMat.set(0, i, nBus[i]);
		for (int j = 0; j < nNLine; j++) {
			int nLine = 0;
			int nLineFully = nBus[i] * (nBus[i] - 1) / 2;
			nLineFully = MYMIN(nLineFully, nLineMax);
			switch (j)
			{
			case 0:
				nLine = nBus[i] - 1;
				break;
			case 1:
				nLine = 0.25 * nLineFully;
				break;
			case 2:
				nLine =  0.5 * nLineFully;
				break;
			case 3:
				nLine = 0.75 * nLineFully;
				break;
			case 4:
				nLine = nLineFully;
				break;
			default:
				throw std::invalid_argument("wrong number of differents number of line");
			}
			nLineMat.set(0, indice, nLine);
			std::cout << indice << " iteration (" << i << ", " << j << ") nBus " << nBus[i] << " nline " << nLine << std::endl;
			for (int simu = 0; simu < nSimu; simu++) {
				int B = nBus[i];
				MatrixGPU VoltageRealIm(2 * B, 1, values[j] + values2[j], 1);
				MatrixGPU RMGgrid(B + 2 * nLine, 1, values[j], 1);
				MatrixGPU RPGgrid(B + 2 * nLine, 1, values2[j], 1);
				MatrixGPU W0(2 * B, 1, values[j] - values2[j], 1);
				MatrixGPU Rgrid(B + 2 * nLine, 1, values2[j] * 2, 1);
				MatrixGPU Xgrid(B + 2 * nLine, 1, values[j] *2, 1);
				MatrixGPU CoresTrans(B + 2 * nLine, 1, 0, 1);
				

				time = 0;
				int numBlocks;
				for (int repet = 0; repet < nRepet; repet++) {
					StudyCase cas;
					if (j) {
						cas.genGridHTB(B, nLine, 1, 10, 1);
					}
					else {
						cas.genGridBT(B, B, B, 10, 1);
					}
					
					MatrixGPU VoltageRealImCopy(VoltageRealIm);
					MatrixGPU RMGgridCopy(RMGgrid);
					MatrixGPU RPGgridCopy(RPGgrid);
					MatrixGPU W0Copy(W0);
					MatrixGPU RgridCopy(Rgrid);
					MatrixGPU XgridCopy(Xgrid);
			
					MatrixGPU CoresVoiLin(cas.getCoresVoiLin());
					MatrixGPU CoresBusLin(cas.getCoresBusLin());
					MatrixGPU nLines(cas.getNLines());
					
					if (method) {
						CoresTrans.transferCPU();
						int* decompte = new int[B];
						for (int i = 0; i < B; i++) {
							decompte[i] = 0;
						}


						for (int i = 0; i < B; i++) {
							int begin = CoresBusLin.get(i, 0);
							for (int l = begin + 1; l < (begin + nLines.get(i, 0)); l++) { // l = (j, i)
								int j = CoresVoiLin.get(l, 0);
								CoresTrans.set(l, 0, CoresBusLin.get(j, 0) + 1 + decompte[j]);
								decompte[j]++;
							}
						}


						DELETEA(decompte);

						
						CoresTrans.transferGPU();
					}
					CoresVoiLin.transferGPU();
					CoresBusLin.transferGPU();
					nLines.transferGPU();
					

					
					int L = nLine;
					int BL2 = B + 2 * L;
					int B2 = 2 * B;
					cudaDeviceSynchronize();
					switch (method)
					{
					case 0:
						numBlocks = 1;
						a = std::chrono::high_resolution_clock::now();
						calculVoltOneStepTest<256> << <1, blockSize, B2 * sizeof(float) >> > (VoltageRealImCopy._matrixGPU, W0Copy._matrixGPU, RgridCopy._matrixGPU, XgridCopy._matrixGPU, RMGgridCopy._matrixGPU, RPGgridCopy._matrixGPU, CoresVoiLin._matrixGPU, CoresBusLin._matrixGPU, nLines._matrixGPU, B);
						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 1:
						numBlocks = B;
						a = std::chrono::high_resolution_clock::now();
						calculVoltStep1Test<256> << <numBlocks, blockSize, B2 * sizeof(float) >> > (VoltageRealImCopy._matrixGPU, W0Copy._matrixGPU, RgridCopy._matrixGPU, XgridCopy._matrixGPU, RMGgridCopy._matrixGPU, RPGgridCopy._matrixGPU, CoresVoiLin._matrixGPU, CoresBusLin._matrixGPU, nLines._matrixGPU, B);
						calculVoltStep2bisTest << <1, blockSize, 2 * (BL2 + B) * sizeof(float) >> > (VoltageRealImCopy._matrixGPU, RMGgridCopy._matrixGPU, RPGgridCopy._matrixGPU, CoresVoiLin._matrixGPU, CoresBusLin._matrixGPU, nLines._matrixGPU, CoresTrans._matrixGPU, B, BL2);

						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					case 2:
						numBlocks = B;
						a = std::chrono::high_resolution_clock::now();
						calculVoltStep1Test<256> << <numBlocks, blockSize, B2 * sizeof(float) >> > (VoltageRealImCopy._matrixGPU, W0Copy._matrixGPU, RgridCopy._matrixGPU, XgridCopy._matrixGPU, RMGgridCopy._matrixGPU, RPGgridCopy._matrixGPU, CoresVoiLin._matrixGPU, CoresBusLin._matrixGPU, nLines._matrixGPU, B);
						calculVoltStep2Test << <1, blockSize, 2 * (BL2 + B) * sizeof(float) >> > (VoltageRealImCopy._matrixGPU, RMGgridCopy._matrixGPU, RPGgrid._matrixGPU, CoresVoiLin._matrixGPU, CoresBusLin._matrixGPU, nLines._matrixGPU, CoresTrans._matrixGPU, B);

						cudaDeviceSynchronize();
						b = std::chrono::high_resolution_clock::now();
						break;
					default:
						return 0;
						break;
					}
					time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
				}
				temps.set(indice, simu, (float)time / nRepet);
			}
			indice++;
		}
	}
	nBusMat.saveCSV(fileName, mode);
	nLineMat.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	return temps.sum();
}

template <unsigned int blockSize>
__global__  void calculVoltOneStepTest(float* VoltageRealIm, float* W0, float* Rgrid, float* Xgrid, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B) {

	int thIdx = threadIdx.x;
	//int i = blockIdx.x; un seul bloc
	int size = blockDim.x;

	extern __shared__ float Voltage[];
	__shared__ float shArr[blockSize];
	__shared__ float shArr2[blockSize];
	for (int k = thIdx; k < 2 * B; k += size) {
		Voltage[k] = VoltageRealIm[k];
	}
	__syncthreads();

	for (int i = 1; i < B; i++) { // calcul de i = iter
		float sum = 0;
		float sum2 = 0;
		int begin = CoresBusLin[i];
		int end = begin + nLines[i];

		for (int l = begin + thIdx + 1; l < end; l += size) {
			int k = CoresVoiLin[l];
			sum  += (RMGgrid[l] * Voltage[k] - RPGgrid[l] * Voltage[k + B]);
			sum2 += (RPGgrid[l] * Voltage[k] + RMGgrid[l] * Voltage[k + B]);
		}
		shArr[thIdx] = sum;
		shArr2[thIdx] = sum2;
		__syncthreads();

		if (blockSize >= 512) {
			if (thIdx < 256) {
				shArr[thIdx] += shArr[thIdx + 256];
				shArr2[thIdx] += shArr2[thIdx + 256];
			}
			__syncthreads();
		}
		if (blockSize >= 256) {
			if (thIdx < 128) {
				shArr[thIdx] += shArr[thIdx + 128];
				shArr2[thIdx] += shArr2[thIdx + 128];
			}
			__syncthreads();
		}
		if (blockSize >= 128) {
			if (thIdx < 64) {
				shArr[thIdx] += shArr[thIdx + 64];
				shArr[thIdx] += shArr2[thIdx + 64];
			} __syncthreads();
		}
		if (thIdx < 32) {
			warpReduce<blockSize>(shArr, thIdx);
			warpReduce<blockSize>(shArr2, thIdx);
		}

		if (thIdx == 0) {
			float vi = Voltage[i];
			float wi = Voltage[i + B];
			float r = Rgrid[i];
			float x = Xgrid[i];
			float W0_local = W0[i];
			float W0B_local = W0[i + B];
			                                                    
			float norm = vi * vi + wi * wi;
			float c = (W0_local * vi + W0B_local * wi) / norm;
			float d = (W0_local * wi - W0B_local * vi) / norm;



			VoltageRealIm[i] = -shArr[0] + c * r - d * x;
			VoltageRealIm[i + B] = -shArr2[0] + d * r + c * x;
		}
		__syncthreads();
	}
	for (int k = thIdx; k < 2 * B; k += size) {
		VoltageRealIm[k] = Voltage[k];
	}
	__syncthreads();


}


__global__ void calculVoltStep2Test(float* VoltageRealIm, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B) {

	/* int thIdx = threadIdx.x;
	 int i = blockIdx.x;
	 int size = blockDim.x;
	 int begin = CoresBusLin[i]; // k = CoresBusLin[iter]; !!!
	 int end = begin + nLines[i]; // k + nLines[iter]
	 // il ne faut pas Ypj/Ypp mais bien Yjp/Yjj, donc il faut savoir quel voisin est p...


	if(i>iter){
		 for (int voisin = thIdx + begin + 1; voisin < end; voisin += size) {
			 int p = CoresVoiLin[voisin];

			 if (p == iter ) { // pour trouver quel indice est p, c'est plutôt nul

				 float db1 = RMGgrid[voisin] * VoltageRealIm[iter] - RPGgrid[voisin] * VoltageRealIm[iter + B];
				 float db2 = RPGgrid[voisin] * VoltageRealIm[iter] + RMGgrid[voisin] * VoltageRealIm[iter + B];

				 VoltageRealIm[i] = VoltageRealIm[i] - db1;
				 VoltageRealIm[i + B] = VoltageRealIm[i + B] - db2;
			 }
		 }
	}*/

	int thIdx = threadIdx.x;
	//int i = blockIdx.x; un seul bloc
	int size = blockDim.x;

	extern __shared__ float Voltage[];

	for (int k = thIdx; k < 2 * B; k += size) {
		Voltage[k] = VoltageRealIm[k];
	}
	__syncthreads();


	for (int iter = 0; iter < B - 1; iter++) {
		int begin = CoresBusLin[iter]; // k = CoresBusLin[iter]; !!!
		int end = begin + nLines[iter]; // k + nLines[iter]
		float ei = Voltage[iter];
		float fi = Voltage[iter + B];

		for (int l = thIdx + begin + 1; l < end; l += size) { // voisin
			int j = CoresVoiLin[l];

			if (j > iter) {
				int lTrans = CoresTrans[l]; // accès pas du tout coalescent !!!

				float ri = RMGgrid[lTrans]; // accès pas du tout coalescent mais c'est sur la mémoire partagé
				float li = RPGgrid[lTrans];

				float db1 = ri * ei - li * fi;
				float db2 = li * ei + ri * fi;

				Voltage[j] = Voltage[j] - db1;
				Voltage[j + B] = Voltage[j + B] - db2;
			}
		}
		__syncthreads();
	}

	for (int k = thIdx; k < 2 * B; k += size) {
		VoltageRealIm[k] = Voltage[k];
	}
	__syncthreads();


}

__global__ void calculVoltStep2bisTest(float* VoltageRealIm, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B, int BL2) {


	int thIdx = threadIdx.x;
	//int i = blockIdx.x; un seul bloc
	int size = blockDim.x;


	extern __shared__ float RI[];

	float* Voltage = &RI[2 * BL2];
	for (int l = thIdx; l < BL2; l += size) { // coalecent
		RI[l] = RMGgrid[l];
		RI[l + BL2] = RPGgrid[l];
	}
	for (int k = thIdx; k < 2 * B; k += size) {
		Voltage[k] = VoltageRealIm[k];
	}
	__syncthreads();


	for (int iter = 0; iter < B - 1; iter++) {
		int begin = CoresBusLin[iter]; // k = CoresBusLin[iter]; !!!
		int end = begin + nLines[iter]; // k + nLines[iter]

		float ei = Voltage[iter];
		float fi = Voltage[iter + B];

		for (int l = thIdx + begin + 1; l < end; l += size) { // voisin
			int j = CoresVoiLin[l];

			if (j > iter) {
				int lTrans = CoresTrans[l];

				float ri = RI[lTrans]; // accès pas du tout coalescent mais c'est sur la mémoire partagé
				float li = RI[lTrans + BL2];


				float db1 = ri * ei - li * fi;
				float db2 = li * ei + ri * fi;

				Voltage[j] = Voltage[j] - db1;
				Voltage[j + B] = Voltage[j + B] - db2;
			}
		}
		__syncthreads();
	}
	for (int k = thIdx; k < 2 * B; k += size) {
		VoltageRealIm[k] = Voltage[k];
	}
	__syncthreads();

}



template <unsigned int blockSize>
__global__ void calculVoltStep1Test(float* VoltageRealIm, float* W0, float* Rgrid, float* Xgrid, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B) {

	__shared__ float shArr[blockSize];
	__shared__ float shArr2[blockSize];
	extern __shared__ float shE[];
	int thIdx = threadIdx.x;
	int i = blockIdx.x;
	int step = blockSize;
	int begin = CoresBusLin[i];
	int end = begin + nLines[i];
	int B2 = 2 * B;

	if (i != 0) {
		for (int n = thIdx; n < B2; n += step)
		{
			shE[n] = VoltageRealIm[n];
		}
		__syncthreads();
		float sum = 0;
		float sum2 = 0;
		for (int l = begin + thIdx + 1; l < end; l += step) {
			int k = CoresVoiLin[l];
			if (k > i) {
				sum -= (RMGgrid[l] * shE[k] - RPGgrid[l] * shE[k + B]);
				sum2 -= (RPGgrid[l] * shE[k] + RMGgrid[l] * shE[k + B]);
			}
		}


		shArr[thIdx] = sum;
		shArr2[thIdx] = sum2;
		__syncthreads();

		if (blockSize >= 512) {
			if (thIdx < 256) {
				shArr[thIdx] += shArr[thIdx + 256];
				shArr2[thIdx] += shArr2[thIdx + 256];
			}
			__syncthreads();
		}
		if (blockSize >= 256) {
			if (thIdx < 128) {
				shArr[thIdx] += shArr[thIdx + 128];
				shArr2[thIdx] += shArr2[thIdx + 128];
			}
			__syncthreads();
		}
		if (blockSize >= 128) {
			if (thIdx < 64) {
				shArr[thIdx] += shArr[thIdx + 64];
				shArr[thIdx] += shArr2[thIdx + 64];
			} __syncthreads();
		}
		if (thIdx < 32) {
			warpReduce<blockSize>(shArr, thIdx);
			warpReduce<blockSize>(shArr2, thIdx);
		}
		if (thIdx == 0) {
			float vi = shE[i];
			float wi = shE[i + B];
			float r = Rgrid[i];
			float x = Xgrid[i];
			float W0_local = W0[i];
			float W0B_local = W0[i + B];

			float norm = vi * vi + wi * wi;
			float c = (W0_local * vi + W0B_local * wi) / norm;
			float d = (W0_local * wi - W0B_local * vi) / norm;



			VoltageRealIm[i] = shArr[0] + c * r - d * x;
			VoltageRealIm[i + B] = shArr2[0] + d * r + c * x;
		}
	}
}
