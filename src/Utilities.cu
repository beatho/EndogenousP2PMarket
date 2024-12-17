

#include "../head/Utilities.cuh"
#define PI 3.14159265359
#define ITERNEWTON 50
#define EPSNEWTON 0.00000001
#define F3(a, b, c, x) (x*x*x + a *x*x + b*x + c)
#define F4(a, b, c, d, x) (x*x*x*x + a *x*x*x + b*x*x + c*x + d)

#define FPRIM3(a, b, x) (3*x*x + 2*a*x + b)
#define FPRIM4(a, b, c, x) (4*x*x*x + a*3*x*x + b*2*x + c)

#define FSECON3(a, x) (6*x + 2*a)
#define FSECON4(a, b, x) (12*x*x + a*6*x + b*2)

#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)

template <typename T>
void check(T err, const char* const func, const char* const file,
	const int line)
{
	if (err != cudaSuccess)
	{
		std::cerr << "CUDA Runtime Error at: " << file << ":" << line
			<< std::endl;
		std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
		// We don't exit when we encounter CUDA errors in this example.
		// std::exit(EXIT_FAILURE);
	}
}


void checkLast(const char* const file, const int line)
{
	cudaDeviceSynchronize();
	cudaError_t err{ cudaGetLastError() };
	if (err != cudaSuccess)
	{
		std::cerr << "CUDA Runtime Error at: " << file << ":" << line
			<< std::endl;
		std::cerr << cudaGetErrorString(err) << std::endl;
		// We   exit when we encounter CUDA errors in this example.
		std::exit(EXIT_FAILURE);
	}
}






int resolveRealPolynome3without2term(double* root, double* coef) {
	/*
	* return : the number of real root for the polynome x^3 + px + q = 0
	* root : is a array of size 3
	* coeff : is a array of size 2 (p and q of x^3 + px + q = 0)
	* if complexe root, root[1] is the real part and root[2] the imaginary part 
	* must add -b/3a if bx^2 is not null at the begining and coefPolynome3From4to2coef is used
	*/
	double p = coef[0];
	double q = coef[1];
	double Delta = -4 * p*p*p - 27 *q*q;
	if (Delta == 0) {
		root[0] = -3 * q / (2 * p);
		root[1] = -2 * root[0];
		root[2] = root[0];
		return 2;
	}
	else if (Delta < 0) {
		double z0 = cbrt((-q + sqrt(-Delta / 27.0)) / 2.0) + cbrt((-q - sqrt(-Delta / 27.0)) / 2.0);
		root[0] = z0; // b2 = z0

		double c2 = p + z0 * z0;
		double delta2 = z0 * z0 - 4 * c2; // négatif normalement
		root[1] = -z0 / 2; // partie réelle de la racine double
		root[2] = sqrt(-delta2) / 2;
		return 1;
	}
	else {
		double r = (3 * q * sqrt(3)) / (2 * p * sqrt(-p)) ;
		for (int k = 0; k < 3; k++) {
			root[k] = 2.0 * sqrt(-p / 3.0) * cos((acos(r) + 2.0 * k * PI) / 3.0	);
		}

		return 3;
	}
}

int resolveRealPolynome3Newton(double* root, double* coef, double init)
{
	/*
	* return : the number of real root for the polynome x^3 + s x^2 + px + q = 0
	* root : is a array of size 3
	* coeff : is a array of size 3 (s and p and q of x^3 + s x^2 + px + q = 0 )
	*/
	
	double x_i = 0;
	int nRoot = 0;
	int i = 0;
	double x_pre = init;
	if (coef[2] != 0) {
		// solve Newton
				
		double eps = 2 * EPSNEWTON;

		while (i<ITERNEWTON && eps>EPSNEWTON)
		{
			//std::cout << i << " : " << x_i << std::endl;
			x_i = x_pre - (F3(coef[0], coef[1], coef[2], x_pre)) / (FPRIM3(coef[0], coef[1], x_pre));
			eps = (x_i - x_pre) * (x_i - x_pre);
			x_pre = x_i;
			i++;
		}	
		//std::cout << i << " " << eps << " " << x_i << std::endl;
		
	}
	if (i == ITERNEWTON) {
		i = 0; // on réessaie avec une autre init
		double b = coef[0];
		double c = coef[1];
		double d = coef[2];

		double p = (3 * c - b * b) / (3);
		double q = (2 * b * b * b - 9 * b * c + 27 * d) / (27);

		if (q == 0) {
			x_i = 0 - b / 3;
		}
		else if (F3(0, p, q, init) > 0) {
			x_pre = findAntpoly3Neg(p, q);
			if (F3(0, p, q, x_pre) >= 0) {
				std::cout << "polynome " << p << " " << q << std::endl;
				std::cout << "problème on the solution neg " << x_pre << " " << F3(0, p, q, x_pre) << std::endl;
			}
		}
		else {
			x_pre = findAntpoly3Pos(p, q);
			if (F3(0, p, q, x_pre) <= 0) {
				std::cout << "polynome " << p << " " << q << std::endl;
				std::cout << "problème on the solution pos " << x_pre << " " << F3(0, p, q, x_pre) << std::endl;
			}
		}
		double eps = 2 * EPSNEWTON;

		while (i<ITERNEWTON && eps>EPSNEWTON)
		{
			//std::cout << i << " : " << x_i << std::endl;
			x_i = x_pre - (F3(0, p, q, x_pre)) / (FPRIM3(0, p, x_pre));
			eps = (x_i - x_pre) * (x_i - x_pre);
			x_pre = x_i;
			i++;
		}
		x_i = x_i - b / 3;

		if (i == ITERNEWTON) { // racine non trouvé ...
			std::cout << "***************** Cela n'a pas marché **********" << std::endl;
			double coefTemp[2];

			coefTemp[0] = p;
			coefTemp[1] = q;
			nRoot = resolveRealPolynome3without2term(root, coefTemp);
			std::cout << "analytics method used, nRoot = " << nRoot << std::endl;
			for (int k = 0; k < nRoot; k++) {
				root[k] = root[k] - b / 3;
			}

			return nRoot;
		}
	}
	
	root[nRoot] = x_i;
	nRoot++;
	
	// second degré
	// x^2 + b x + c = 0 tel que  x^3 + s x^2 + px + q = (x-x_i)(x^2 + b x + c)
	double B = coef[0] + x_i;
	double C = coef[1] + x_i * B;


	double delta = B * B - 4 * C;
	
	if (delta == 0) {
		double z = -B / 2;
		root[nRoot] = z;
		nRoot++;
		//root[1] = z;
		//nRoot++;
		//std::cout << " z " << z << std::endl;
		return nRoot;
	}
	else if (delta > 0) {
		double z1 = (-B + sqrt(delta)) / 2;
		double z2 = (-B - sqrt(delta)) / 2;
		root[nRoot] = z1;
		nRoot++;
		root[nRoot] = z2;
		nRoot++;
		return nRoot;
	}
	else { // delta < 0
		//std::cout << "pas d'autres racines réelle !!!! " << std::endl;
	}
return nRoot;
	
	
	
}


int resolveRealPolynome3Laguerre(double* root, double* coef, double init)
{
	/*
	* return : the number of real root for the polynome x^3 + s x^2 + px + q = 0
	* root : is a array of size 3
	* coeff : is a array of size 3 (s and p and q of x^3 + s x^2 + px + q = 0 )
	*/

	double x_i = 0;
	int nRoot = 0;
	int n = 3; // degré du polynome
	int i = 0;
	if (coef[2] != 0) {
		// solve Laguerre
		double x_pre = init;

		
		double eps = 2 * EPSNEWTON;

		while (i<ITERNEWTON && eps>EPSNEWTON)
		{
			std::cout << i << " : " << x_i << std::endl;
			double p = F3(coef[0], coef[1], coef[2], x_pre);
			double p2 = FPRIM3(coef[0], coef[1], x_pre);
			double p3 = FSECON3(coef[0], x_pre);
			double S1 = p / p2;
			double S2 = p3 / p - S1 * S1;

			if (S1 > 0) {
				x_i = x_pre - n / (S1 + sqrt((1 - n) * (n * S2 + S1 * S1)));
			}
			else {
				x_i = x_pre - n / (S1 - sqrt((1 - n) * (n * S2 + S1 * S1)));
			}

			eps = (x_i - x_pre) * (x_i - x_pre);
			x_pre = x_i;
			i++;
		}
		std::cout << i << " " << eps << " " << x_i << std::endl;
	}
	if (i == ITERNEWTON) { // racine non trouvé ...
		double coefTemp[2];
		double b = coef[0];
		double c = coef[1];
		double d = coef[2];

		double p = (3 * c - b * b) / (3);
		double q = (2 * b * b * b - 9 * b * c + 27 * d) / (27);
		coefTemp[0] = p;
		coefTemp[1] = q;
		nRoot = resolveRealPolynome3without2term(root, coefTemp);
		std::cout << "analytics method used, nRoot = " << nRoot << std::endl;
		for (int k = 0; k < nRoot; k++) {
			root[k] = root[k] - b / 3;
		}

		return nRoot;
	} 

	root[nRoot] = x_i;
	nRoot++;

	// second degré
	// x^2 + b x + c = 0 tel que  x^3 + s x^2 + px + q = (x-x_i)(x^2 + b x + c)
	double B = coef[0] + x_i;
	double C = coef[1] + x_i * B;


	double delta = B * B - 4 * C;

	if (delta == 0) {
		double z = -B / 2;
		root[nRoot] = z;
		nRoot++;
		//root[1] = z;
		//nRoot++;
		//std::cout << " z " << z << std::endl;
		return nRoot;
	}
	else if (delta > 0) {
		double z1 = (-B + sqrt(delta)) / 2;
		double z2 = (-B - sqrt(delta)) / 2;
		root[nRoot] = z1;
		nRoot++;
		root[nRoot] = z2;
		nRoot++;
		return nRoot;
	}
	else { // delta < 0
		//std::cout << "pas d'autres racines réelle !!!! " << std::endl;
	}


	return nRoot;
}


int resolveRealPolynome3Halley(double* root, double* coef, double init)
{
	/*
	* return : the number of real root for the polynome x^3 + s x^2 + px + q = 0
	* root : is a array of size 3
	* coeff : is a array of size 3 (s and p and q of x^3 + s x^2 + px + q = 0 )
	*/

	double x_i = 0;
	int nRoot = 0;
	int n = 3; // degré du polynome
	int i = 0;
	double x_pre = init;
	
	if (coef[2] != 0) {
		// solve Laguerre
		
		double eps = 2 * EPSNEWTON;

		while (i<ITERNEWTON && eps>EPSNEWTON)
		{
			//std::cout << i << " : " << x_i << std::endl;
			double p = F3(coef[0], coef[1], coef[2], x_pre);
			double p2 = FPRIM3(coef[0], coef[1], x_pre);
			double p3 = FSECON3(coef[0], x_pre);
			
			x_i = x_pre - (2 * p * p2) / (2 * p2 * p2 - p * p3);

			eps = (x_i - x_pre) * (x_i - x_pre);
			x_pre = x_i;
			i++;
		}
	}
	if (i == ITERNEWTON) {
		i = 0; // on réessaie avec une autre init
		double b = coef[0];
		double c = coef[1];
		double d = coef[2];

		double p = (3 * c - b * b) / (3);
		double q = (2 * b * b * b - 9 * b * c + 27 * d) / (27);

		if (q == 0) {
			x_i = 0 - b / 3;
		}
		else if (F3(0, p, q, init) > 0) {
			x_pre = findAntpoly3Neg(p, q);
			if (F3(0, p, q, x_pre) >= 0) {
				std::cout << "polynome " << p << " " << q << std::endl;
				std::cout << "problème on the solution neg " << x_pre << " " << F3(0, p, q, x_pre) << std::endl;
			}
		}
		else {
			x_pre = findAntpoly3Pos(p, q);
			if (F3(0, p, q, x_pre) <= 0) {
				std::cout << "polynome " << p << " " << q << std::endl;
				std::cout << "problème on the solution pos " << x_pre << " " << F3(0, p, q, x_pre) << std::endl;
			}
		}
		double eps = 2 * EPSNEWTON;

		while (i<ITERNEWTON && eps>EPSNEWTON)
		{
			//std::cout << i << " : " << x_i << std::endl;
			double p1 = F3(0, p, q, x_pre);
			double p2 = FPRIM3(0, p, x_pre);
			double p3 = FSECON3(0, x_pre);

			x_i = x_pre - (2 * p1 * p2) / (2 * p2 * p2 - p1 * p3);

			eps = (x_i - x_pre) * (x_i - x_pre);
			x_pre = x_i;
			i++;
		}
		x_i = x_i - b / 3;

		if (i == ITERNEWTON) { // racine non trouvé ...
			std::cout << "***************** Cela n'a pas marché **********" << std::endl;
			double coefTemp[2];
			
			coefTemp[0] = p;
			coefTemp[1] = q;
			nRoot = resolveRealPolynome3without2term(root, coefTemp);
			std::cout << "analytics method used, nRoot = " << nRoot << std::endl;
			for (int k = 0; k < nRoot; k++) {
				root[k] = root[k] - b / 3;
			}

			return nRoot;
		}
	}
	
	
	
	

	root[nRoot] = x_i;
	nRoot++;

	// second degré
	// x^2 + b x + c = 0 tel que  x^3 + s x^2 + px + q = (x-x_i)(x^2 + b x + c)
	double B = coef[0] + x_i;
	double C = coef[1] + x_i * B;


	double delta = B * B - 4 * C;

	if (delta == 0) {
		double z = -B / 2;
		root[nRoot] = z;
		nRoot++;
		//root[1] = z;
		//nRoot++;
		//std::cout << " z " << z << std::endl;
		return nRoot;
	}
	else if (delta > 0) {
		double z1 = (-B + sqrt(delta)) / 2;
		double z2 = (-B - sqrt(delta)) / 2;
		root[nRoot] = z1;
		nRoot++;
		root[nRoot] = z2;
		nRoot++;
		return nRoot;
	}
	else { // delta < 0
		//std::cout << "pas d'autres racines réelle !!!! " << std::endl;
	}


	return nRoot;
}


int resvolveRealPolynome4without2term(double* root, double* coef)
{
	/*
	* return : the number of real root for the polynome x^4 + bx^3 + dx + e = 0
	* root : is a array of size 4
	* coeff : is a array of size 3 (b,d,e)
	*/
	double b = coef[0];
	double d = coef[1];
	double e = coef[2];
	int nRoot = 0;

	if (b * b * b + 8 * d == 0) {
		//if (abs(b * b * b + 8 * d) < 0.00000001) {
		// passage de p^4 + b p^3 + d p + e -> a p^4 + b p^2 + c = 0
		double B = -3 * b * b / 8;
		double C = -3* b*b*b*b/256 - b*d/4 + e;

		double delta = B * B - 4 * C;
		//std::cout << "Delta " << delta;
		if (delta  == 0) {
			double z = - B / 2;
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
	//std::cout << "no real root found " << abs(b * b * b + 8 * d) << " " << lambda0 << " " << mu0 << " " << DeltaP << " " << DeltaM << std::endl;
	
	return nRoot;
}

void coefPolynome3From4to2coef(double* coef4, double* coef2)
{

	double a = coef4[0];

	if (a == 0) {
		throw std::invalid_argument("must be a thrid degree polynome, a!=0 ");
	}

	double b = coef4[1];
	double c = coef4[2];
	double d = coef4[3];

	coef2[0] = (3 * a * c - b * b) / (3 * a * a);
	coef2[1] = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);

}


int resolveRealPolynome4Newton(double* root, double* coef, double init)
{
	/*
	* return : the number of real root for the polynome x^4 + bx^3 + c x^2  + dx + e = 0
	* root : is a array of size 4
	* coeff : is a array of size 4 (b, c, d,e)
	*/
	double coef3[3];
	double x_i = 0;
	
	if (coef[3] != 0) {
		// solve Newton
		double x_pre = init;

		int i = 0;
		double eps = 2 * EPSNEWTON;

		while (i<ITERNEWTON && eps>EPSNEWTON)
		{
			//std::cout << i << " : " << x_i << std::endl;
			x_i = x_pre - (F4(coef[0], coef[1], coef[2], coef[3], x_pre)) / (FPRIM4(coef[0], coef[1], coef[2], x_pre));
			eps = (x_i - x_pre) * (x_i - x_pre);
			x_pre = x_i;
			i++;
		}
	}
	
	// troisième degré
	// x^2 + b x + c = 0 tel que  x^3 + s x^2 + px + q = (x-x_i)(x^2 + b x + c)
	coef3[0] = coef[0] + x_i;
	coef3[1] = coef[1] + x_i * coef3[0];
	coef3[2] = coef[2] + x_i * coef3[1];

	//std::cout << "Nouveau Polynome " << coef3[0] << " " << coef3[1] << " " << coef3[2] << std::endl;

	int nRoot = resolveRealPolynome3Newton(root, coef3, 0);
	root[nRoot] = x_i;
	nRoot++;

	return nRoot;
}


int resolveRealPolynome4Halley(double* root, double* coef, double init)
{
	/*
	* return : the number of real root for the polynome x^4 + bx^3 + c x^2  + dx + e = 0
	* root : is a array of size 4
	* coeff : is a array of size 4 (b, c, d,e)
	*/

	double x_i = 0;

	int n = 3; // degré du polynome
	if (coef[2] != 0) {
		 
		double x_pre = init;

		int i = 0;
		double eps = 2 * EPSNEWTON;

		while (i<ITERNEWTON && eps>EPSNEWTON)
		{
			double p = F4(coef[0], coef[1], coef[2], coef[3], x_pre);
			double p2 = FPRIM4(coef[0], coef[1], coef[2], x_pre);
			double p3 = FSECON4(coef[0], coef[1], x_pre);

			x_i = x_pre - (2 * p * p2) / (2 * p2 * p2 - p * p3);

			eps = (x_i - x_pre) * (x_i - x_pre);
			x_pre = x_i;
			i++;
		}
	}


	// troisième degré
	double coef3[3];


	double b = coef[0] + x_i;
	double c = coef[1] + x_i * b;
	double d = coef[2] + x_i * c;

	/*double p = (3  * c - b * b) / (3);
	//double q = (2 * b * b * b - 9  * b * c + 27  * d) / (27);
	//coef3[0] = 0;
	//coef3[1] = p;
	//coef3[2] = q;
	std::cout << "Nouveau Polynome " << coef3[0] << " " << coef3[1] << " " << coef3[2] << std::endl;
	std::cout << "Nouveau Polynome ou " << b << " " << c << " " << d << std::endl;*/
	coef3[0] = b;
	coef3[1] = c;
	coef3[2] = d;

	int nRoot = resolveRealPolynome3Halley(root, coef3, 0);
	/*for (int i = 0; i < nRoot; i++) {
		root[i] = root[i] - b / 3;
	}*/
	root[nRoot] = x_i;
	nRoot++;
	

	return nRoot;
}

int resvolveRealPolynome4without2termLagrange(double* root, double* coef) {
	/*
	* return : the number of real root for the polynome x^4 + bx^3 + dx + e = 0
	* root : is a array of size 4
	* coeff : is a array of size 3 (b,d,e)
	*/
	double b = coef[0];
	double d = coef[1];
	double e = coef[2];
	int nRoot = 0;

	// il faut passer de b d e à p q t, c'est le coef devant z^3 qui doit etre nul pas celui devant z^2 !!!

	double p = -3.0 * b * b / (8.0);
	double q = d + (b * b * b) / 8.0;
	double t = e - b * d / 4.0 - 3 * b * b * b * b / 256;

	//t = e/a - b*d/(4*a^2) + c*b^2/(16*a^3) - 3*b^4/(256 * a^4);


	int signe = 1 - 2 * (q > 0);

	/*double b2 = 2 * p;
	double c2 = (p * p - 4 * t);
	double d2 = -q * q;

	std::cout << b2 << " " << c2 << " " << d2 << std::endl;*/

	double coef2[2];
	double rootlambda[3];
	//coef2[0] = -b2*b2 / 3.0 + c2;
	//coef2[1] = (b2 / 27.0) * (2 * b2 *b2 - 9 * c2) + d2;
	coef2[0] = -4.0 * p * p / 3.0 + p * p - 4.0 * t;
	coef2[1] =  2.0 * p / 27.0 * (36.0 * t - p *p) - q * q;


	//std::cout << b2 /27.0 << " " << (2 * b2 * b2 - 9 * c2)<< " " << d2 << std::endl;
	//std::cout << coef2[0] << " " << coef2[1] << std::endl;

	int nRootlambda = resolveRealPolynome3without2term(rootlambda, coef2);
	// il y a eu un changement de variable donc il faut décaller les racines
	double offset = 2.0 * p / 3.0; // b/3a du poly du 2nd ordre

	
	if (nRootlambda == 1) { // une réelle et 2 compl
		rootlambda[0] = signe * sqrt(rootlambda[0] - offset); // racine réelle
		rootlambda[1] = rootlambda[1] - offset; // Partie relle de la racine 


		double terme = sqrt((rootlambda[1] + sqrt(rootlambda[1] * rootlambda[1] + rootlambda[2] * rootlambda[2])) / 2.0);
		
		root[0] = 0.5 * (rootlambda[0] + 2 * terme);
		root[1] = 0.5 * (rootlambda[0] - 2 * terme);

		for (int k = 0; k < 2; k++)
		{
			root[k] = root[k] - b / 4;
		}

		return 2;

	}
	else { // 3 réelles
		for (int i = 0; i < 3; i++) {
			if (rootlambda[i] - offset < 0) {
				return 0; // que des racines négatives -> pas de racines réelles
			}
			else {
				rootlambda[i] = signe * sqrt(rootlambda[i] - offset);
				
			}
		}
		
		root[0] = 0.5 * ( rootlambda[0] + rootlambda[1] + rootlambda[2]);
		root[1] = 0.5 * ( rootlambda[0] - rootlambda[1] - rootlambda[2]);
		root[2] = 0.5 * (-rootlambda[0] + rootlambda[1] - rootlambda[2]);
		root[3] = 0.5 * (-rootlambda[0] - rootlambda[1] + rootlambda[2]);

		for (int k = 0; k < 4; k++)
		{
			root[k] = root[k] - b / 4;
		}

		return 4;

	}
}


int resolveRealPolynome3without2termEigen(double* root, double* coef) {
	/*
	* return : the number of real root for the polynome x^3 + px + q = 0
	* root : is a array of size 3
	* coeff : is a array of size 2 (p and q of x^3 + px + q = 0)
	*/
	int nRoot = 0;

	Eigen::Vector4d coeff(coef[1], coef[0], 0 , 1); //double coef4[4] = { 1, -5, 3, 1 };
	Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
	solver.compute(coeff);
	const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType& r = solver.roots();

	//std::cout << r << std::endl;

	for (int k = 0; k < 3; k++) {
		if (r(k).imag() == 0) { // racine réelle
			root[nRoot] = r(k).real();
			nRoot++;
		}
	}

	return nRoot;
}

int resvolveRealPolynome4without2termEigen(double* root, double* coef)
{
	/*
	* return : the number of real root for the polynome x^4 + bx^3 + dx + e = 0
	* root : is a array of size 4
	* coeff : is a array of size 3 (b,d,e)
	*/

	int nRoot = 0;
	Eigen::VectorXd coeff(5);
	coeff(0) = coef[2];
	coeff(1) = coef[1];
	coeff(2) = 0;
	coeff(3) = coef[0];
	coeff(4) = 1;

	Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
	solver.compute(coeff);
	const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType& r = solver.roots();

	for (int k = 0; k < 4; k++) {
		if (r(k).imag() == 0) { // racine réelle
			root[nRoot] = r(k).real();
			nRoot++;
		}
	}

	return nRoot;

}


int resvolveRealPolynome4without2term(double* root, double* coef, bool Lagrange) {
	if (Lagrange) {
		return resvolveRealPolynome4without2termLagrange(root, coef);
	}
	else {
		return resvolveRealPolynome4without2term(root, coef);
	}
}


__host__ __device__ double findAntpoly3Neg(double p, double q) {
	/* 
		find one x where f(x) <0
		with f(x) = x^3 + px + q
	*/

	if (q > 0) {
		if (p > 0) {
			return -(q / p);
		}
		else {
			return -1.26*MAX(MAX(1 , -  p),  q);
		}
	}
	else {
		return 0;
	}


}
__host__ __device__ double findAntpoly3Pos(double p, double q) {

	/*
		find one x where f(x) >0
		with f(x) = x^3 + px + q
	*/

	if (q < 0) {
		if (p > 0) {
			return -(q / p);
		}
		else {
			return 1.26*MAX(MAX(1, -p), -q);
		}
	}
	else {
		return 0;
	}

}


__device__ int resolveRealPolynome3without2termGPU(double* root, double p, double q) {

	double Delta = -4.0 * p * p * p - 27.0 * q * q;
	if (Delta == 0) {
		root[0] = -3.0 * q / (2.0 * p);
		root[1] = -2.0 * root[0];
		return 2;
	}
	else if (Delta < 0) {
		double z0 =  cbrt((-q + sqrt(-Delta / 27.0)) / 2.0) + cbrt((-q - sqrt(-Delta / 27.0)) / 2.0);
		root[0] = z0; // b2 = z0
		
		double c2 = p + z0*z0;
		double delta2 = z0 * z0 - 4 * c2; // négatif normalement
		root[1] = -z0 / 2; // partie réelle de la racine double
		root[2] = sqrt(-delta2) / 2;
		return 1;
	}
	else {

		for (int k = 0; k < 3; k++) {
			double r = (3.0 * q * sqrt(3.0)) / (2.0 * p * sqrt(-p));
			r = -1.0 * (r < -1.0) + 1.0 * (r > 1.0) + r * (r<1.0 && r>-1.0);
			root[k] = 2.0 * sqrt(-p / 3.0) * cos((acos(r) + 2.0 * k * PI) / 3.0);
		}

		return 3;
	}
}

/**/
__device__ int resvolveRealPolynome4without2termGPU(double* root, double b, double d, double e)
{
	/*
	* return : the number of real root for the polynome x^4 + bx^3 + dx + e = 0
	* root : is a array of size 4
	* coeff : is a array of size 3 (b,d,e)
	*/
	int nRoot = 0;

	if (b * b * b + 8.0 * d == 0) {

		// passage de p^4 + b p^3 + d p + e -> a p^4 + b p^2 + c = 0
		double B = -3.0 * b * b / 8.0;
		double C = -3.0 * b * b * b * b / 256.0 - b * d / 4.0 + e;

		double delta = B * B - 4.0 * C;
		//std::cout << "Delta " << delta;
		if (delta == 0) {
			double z = -B / 2.0;
			nRoot = 2;
			//std::cout << " z " << z << std::endl;
			root[0] = sqrt(z);
			root[1] = -sqrt(z);
			return nRoot;
		}
		else if (delta > 0) {
			double z1 = (-B + sqrt(delta)) / 2.0;
			double z2 = (-B - sqrt(delta)) / 2.0;
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
	}


	double rootlambda[3];
	double coef2_0 = (2.0 * b * d - 8.0 * e) / 8.0;
	double coef2_1 = -(b * b * e + d * d) / 8.0;
	int nRootlambda = resolveRealPolynome3without2termGPU(rootlambda, coef2_0, coef2_1);

	for (int i = 0; i < nRootlambda; i++) {
		double lambda0 = rootlambda[i];
		double mu1 = 2.0 * lambda0 + (b * b) / 4.0;
		if (mu1 > 0) {
			double mu0 = sqrt(mu1);
			double DeltaP = -2.0 * lambda0 + 2.0 * (d - b * lambda0) / mu0 + b * mu0 + b * b / 2.0;
			double DeltaM = -2.0 * lambda0 - 2.0 * (d - b * lambda0) / mu0 - b * mu0 + b * b / 2.0;
			if (DeltaP >= 0) {
				root[nRoot] = (-mu0 + sqrt(DeltaP)) / 2.0 - b / 4.0;
				root[nRoot + 1] = (-mu0 - sqrt(DeltaP)) / 2.0 - b / 4.0;
				nRoot = nRoot + 2;
			}
			if (DeltaM >= 0) {
				root[nRoot] = (mu0 + sqrt(DeltaM)) / 2.0 - b / 4.0;
				root[nRoot + 1] = (mu0 - sqrt(DeltaM)) / 2.0 - b / 4.0;
				nRoot = nRoot + 2;
			}
			if (nRoot > 0) {
				return nRoot;
			}
		}
	}
	return nRoot;
}

__device__ int resvolveRealPolynome4without2termGPULagrange(double* root, double b, double d, double e) {

	int nRoot = 0;

	// il faut passer de b d e à p q t, c'est le coef devant z^3 qui doit etre nul pas celui devant z^2 !!!

	double p = -3.0 * b * b / (8.0);
	double q = d + (b * b * b) / 8.0;
	double t = e - b * d / 4.0 - 3 * b * b * b * b / 256;

	//t = e/a - b*d/(4*a^2) + c*b^2/(16*a^3) - 3*b^4/(256 * a^4);


	int signe = 1 - 2 * (q > 0);

	double rootlambda[3];
	double coef2_0 = -4.0 * p * p / 3.0 + p * p - 4.0 * t;
	double coef2_1 = 2.0 * p / 27.0 * (36.0 * t - p * p) - q * q;


	
	int nRootlambda = resolveRealPolynome3without2termGPU(rootlambda, coef2_0, coef2_1);
	// il y a eu un changement de variable donc il faut décaller les racines
	double offset = 2.0 * p / 3.0; // b/3a du poly du 2nd ordre


	if (nRootlambda == 1) { // une réelle et 2 compl
		rootlambda[0] = signe * sqrt(rootlambda[0] - offset); // racine réelle
		rootlambda[1] = rootlambda[1] - offset; // Partie relle de la racine 


		double terme = sqrt((rootlambda[1] + sqrt(rootlambda[1] * rootlambda[1] + rootlambda[2] * rootlambda[2])) / 2.0);
	
		root[0] = 0.5 * (rootlambda[0] + 2 * terme);
		root[1] = 0.5 * (rootlambda[0] - 2 * terme);

		for (int k = 0; k < 2; k++)
		{
			root[k] = root[k] - b / 4;
		}

		return 2;

	}
	else { // 3 réelles
		for (int i = 0; i < 3; i++) {
			if (rootlambda[i] - offset < 0) {
				return 0; // que des racines négatives -> pas de racines réelles
			}
			else {
				rootlambda[i] = signe * sqrt(rootlambda[i] - offset);

			}
		}
		root[0] = 0.5 * (rootlambda[0] + rootlambda[1] + rootlambda[2]);
		root[1] = 0.5 * (rootlambda[0] - rootlambda[1] - rootlambda[2]);
		root[2] = 0.5 * (-rootlambda[0] + rootlambda[1] - rootlambda[2]);
		root[3] = 0.5 * (-rootlambda[0] - rootlambda[1] + rootlambda[2]);

		for (int k = 0; k < 4; k++)
		{
			root[k] = root[k] - b / 4;
		}

		return 4;
	}
}

__device__ int resvolveRealPolynome4without2termGPU(double* root, double b, double d, double e, bool Lagrange) {
	if (Lagrange) {
		return resvolveRealPolynome4without2termGPULagrange(root, b, d, e);
	}
	else {
		return resvolveRealPolynome4without2termGPU(root, b, d, e);
	}
}


__device__ void coefPolynome3From4to2coefGPU(double* coef4, double* coef2) {
	double a = coef4[0];

	double b = coef4[1];
	double c = coef4[2];
	double d = coef4[3];

	coef2[0] = (3.0 * a * c - b * b) / (3.0 * a * a);
	coef2[1] = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d) / (27.0 * a * a * a);
}


__global__ void resolveSeveralRealPolynome3termGPU(double* nRoot, double* roots, double* coefs, int nPoly) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = index; i < nPoly; i += step) {
		double coefsLocal2[2];
		double rootsLocal[3];
		double coefsLocal4[4];
		int nRootLocal = 0;
		for (int j = 0; j < 4; j++) {
			coefsLocal4[j] = coefs[j * nPoly + i];
		}
		if (coefsLocal4[0] == 0) { // polynone of degre 2
			if (coefsLocal4[1] == 0) { // polynom of degre 1
				if (coefsLocal4[2] != 0) { // no const
					nRootLocal = 1;
					rootsLocal[0] = -coefsLocal4[3] / coefsLocal4[2];
				}
				nRootLocal = 0;
			}
			else {
				// la flemme WIP
			}
		}
		else {
			
			
			coefPolynome3From4to2coefGPU(coefsLocal4, coefsLocal2);
			
			

			nRootLocal = resolveRealPolynome3without2termGPU(rootsLocal, coefsLocal2[0], coefsLocal2[1]);
			for (int k = 0; k < nRootLocal; k++) {
				rootsLocal[k] += -coefsLocal4[1] / (3 * coefsLocal4[0]);
			}
		}
		for (int j = 0; j < nRootLocal; j++) {
			roots[j * nPoly + i] = rootsLocal[j];
		}
		nRoot[i] = nRootLocal;

	}
}

__global__ void resolveSeveralRealPolynome4WO2termGPU(double* nRoot, double* roots, double* coefs, int nPoly) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = index; i < nPoly; i += step) {
		double rootsLocal[4];
		double coefsLocal3[3];
		int nRootLocal = 0;
		for (int j = 0; j < 4; j++) {
			coefsLocal3[j] = coefs[j * nPoly + i];
		}
		
		nRootLocal = resvolveRealPolynome4without2termGPU(rootsLocal, coefsLocal3[0], coefsLocal3[2], coefsLocal3[3]);
		
		
		for (int j = 0; j < nRootLocal; j++) {
			roots[j * nPoly + i] = rootsLocal[j];
		}
		nRoot[i] = nRootLocal;

	}
}


__global__ void resolveSeveralRealPolynome4WO2termGPULagrange(double* nRoot, double* roots, double* coefs, int nPoly) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = index; i < nPoly; i += step) {
		double rootsLocal[4];
		double coefsLocal3[3];
		int nRootLocal = 0;
		for (int j = 0; j < 4; j++) {
			coefsLocal3[j] = coefs[j * nPoly + i];
		}

		nRootLocal = resvolveRealPolynome4without2termGPULagrange(rootsLocal, coefsLocal3[0], coefsLocal3[2], coefsLocal3[3]);


		for (int j = 0; j < nRootLocal; j++) {
			roots[j * nPoly + i] = rootsLocal[j];
		}
		nRoot[i] = nRootLocal;

	}
}



__device__ int resolveRealPolynome3GPU(double* root, double b, double c, double d) {

	double x_i = 0;
	int nRoot = 0;
	int i = 0;
	double x_pre = 0;

	if (d != 0) {
		// solve Laguerre

		double eps = 2 * EPSNEWTON;

		while (i<ITERNEWTON && eps>EPSNEWTON)
		{
			//std::cout << i << " : " << x_i << std::endl;
			double p = F3(b, c, d, x_pre);
			double p2 = FPRIM3(b, c, x_pre);
			double p3 = FSECON3(b, x_pre);

			x_i = x_pre - (2 * p * p2) / (2 * p2 * p2 - p * p3);

			eps = (x_i - x_pre) * (x_i - x_pre);
			x_pre = x_i;
			i++;
		}
	}
	if (i == ITERNEWTON) {
		i = 0; // on réessaie avec une autre init

		double p = (3 * c - b * b) / (3);
		double q = (2 * b * b * b - 9 * b * c + 27 * d) / (27);

		if (F3(0, p, q, 0) > 0) {
			x_pre = findAntpoly3Neg(p, q);
		}
		else {
			x_pre = findAntpoly3Pos(p, q);
		}
		double eps = 2 * EPSNEWTON;

		while (i<ITERNEWTON && eps>EPSNEWTON)
		{
			//std::cout << i << " : " << x_i << std::endl;
			double p1 = F3(0, p, q, x_pre);
			double p2 = FPRIM3(0, p, x_pre);
			double p3 = FSECON3(0, x_pre);

			x_i = x_pre - (2 * p1 * p2) / (2 * p2 * p2 - p1 * p3);

			eps = (x_i - x_pre) * (x_i - x_pre);
			x_pre = x_i;
			i++;
		}
		x_i = x_i - b / 3;

		if (i == ITERNEWTON) { // racine non trouvé ...
			double coefTemp[2];

			coefTemp[0] = p;
			coefTemp[1] = q;
			nRoot = resolveRealPolynome3without2termGPU(root, p, q);
			for (int k = 0; k < nRoot; k++) {
				root[k] = root[k] - b / 3;
			}

			return nRoot;
		}
	}
	root[nRoot] = x_i;
	nRoot++;

	// second degré
	// x^2 + b x + c = 0 tel que  x^3 + s x^2 + px + q = (x-x_i)(x^2 + b x + c)
	double B = b + x_i;
	double C = c + x_i * B;


	double delta = B * B - 4 * C;

	if (delta == 0) {
		double z = -B / 2;
		root[nRoot] = z;
		nRoot++;
		return nRoot;
	}
	else if (delta > 0) {
		double z1 = (-B + sqrt(delta)) / 2;
		double z2 = (-B - sqrt(delta)) / 2;
		root[nRoot] = z1;
		nRoot++;
		root[nRoot] = z2;
		nRoot++;
		return nRoot;
	}
	
	return nRoot;

}/**/
__device__ int resvolveRealPolynome4GPU(double* root, double b, double c, double d, double e) {
	double x_i = 0;

	 
	if (e != 0) {
		 
		double x_pre = 0;

		int i = 0;
		double eps = 2 * EPSNEWTON;

		while (i<ITERNEWTON && eps>EPSNEWTON)
		{
			double p = F4(b, c, d, e, x_pre);
			double p2 = FPRIM4(b, c, d, x_pre);
			double p3 = FSECON4(b, c, x_pre);

			x_i = x_pre - (2 * p * p2) / (2 * p2 * p2 - p * p3);

			eps = (x_i - x_pre) * (x_i - x_pre);
			x_pre = x_i;
			i++;
		}
	}


	// troisième degré
	 

	double B = b + x_i;
	double C = c + x_i * B;
	double D = d + x_i * C;

	/*double p = (3  * c - b * b) / (3);
	//double q = (2 * b * b * b - 9  * b * c + 27  * d) / (27);
	//coef3[0] = 0;
	//coef3[1] = p;
	//coef3[2] = q;
	std::cout << "Nouveau Polynome " << coef3[0] << " " << coef3[1] << " " << coef3[2] << std::endl;
	std::cout << "Nouveau Polynome ou " << b << " " << c << " " << d << std::endl;*/
 

	int nRoot = resolveRealPolynome3GPU(root, B, C, D);
	/*for (int i = 0; i < nRoot; i++) {
		root[i] = root[i] - b / 3;
	}*/
	root[nRoot] = x_i;
	nRoot++;


	return nRoot;

}



__global__ void resolveSeveralRealPolynome3GPU(double* nRoot, double* roots, double* coefs, int nPoly) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = index; i < nPoly; i += step) {
		double coefsLocal3[3];
		double rootsLocal[3];
		double coefsLocal4[4];
		int nRootLocal = 0;
		for (int j = 0; j < 4; j++) {
			coefsLocal4[j] = coefs[j * nPoly + i];
		}
		coefsLocal3[0] = coefsLocal4[1] / coefsLocal4[0];
		coefsLocal3[1] = coefsLocal4[2] / coefsLocal4[0];
		coefsLocal3[2] = coefsLocal4[3] / coefsLocal4[0];

		nRootLocal = resolveRealPolynome3GPU(rootsLocal, coefsLocal3[0], coefsLocal3[1], coefsLocal3[2]);
		

		for (int j = 0; j < nRootLocal; j++) {
			roots[j * nPoly + i] = rootsLocal[j];
		}
		nRoot[i] = nRootLocal;
	}

}


__global__ void resolveSeveralRealPolynome4GPU(double* nRoot, double* roots, double* coefs, int nPoly) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = index; i < nPoly; i += step) {
		double rootsLocal[4];
		double coefsLocal4[4];
		int nRootLocal = 0;
		for (int j = 0; j < 4; j++) {
			coefsLocal4[j] = coefs[j * nPoly + i];
		}

		nRootLocal = resvolveRealPolynome4GPU(rootsLocal, coefsLocal4[0], coefsLocal4[1], coefsLocal4[2], coefsLocal4[3]);


		for (int j = 0; j < nRootLocal; j++) {
			roots[j * nPoly + i] = rootsLocal[j];
		}
		nRoot[i] = nRootLocal;

	}
}




