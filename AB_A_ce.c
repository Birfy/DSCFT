/////////AB2C4 dendritic triblock copolymer //////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include </export/home/cdy/fftw3/include/fftw3.h>
#include </data2/opt/include/gsl/gsl_blas.h>
#include </data2/opt/include/gsl/gsl_linalg.h>
#define MaxIT 50000 //Maximum iteration steps

#define Pi 3.141592653589

#define N_hist 50
/* Parameters used in Anderson convergence */
#define Del(k, i, n) del[(i) + NxNyNz * (n) + N_hist * NxNyNz * (k)]
#define Outs(k, i, n) outs[(i) + NxNyNz * (n) + N_hist * NxNyNz * (k)]
#define U(n, m) up[(m - 1) + (N_rec - 1) * (n - 1)]
#define V(n) vp[n - 1]
#define A(n) ap[n - 1]

double freeE(double *wA, double *wB, double *phA, double *phB, double *eta);
double getConc(double *phlA, double *phlB, double phs0, double *wA, double *wB);
void sovDifFft(double *g, double *w, double *qInt, double z, int ns, int sign, double epk);
void write_ph(double *phA, double *phB, double *wA, double *wB);
double error_cal(double *waDiffs, double *wbDiffs, double *wAs, double *wBs);
void update_flds_hist(double *waDiff, double *wbDiff, double *wAnew, double *wBnew, double *del, double *outs);
void Anderson_mixing(double *del, double *outs, int N_rec, double *wA, double *wB);

void init(int in, double *wA, double *wB);
void initW_Random(double *wA, double *wB);
void initW_BCC(double *wA, double *wB);

int ZDIMM, NsA1, NsB, NsA2;
double lx, ly, lz, ds0;
double hAB, fA1, fB, fA2;
double epA, epB;
int Nx, Ny, Nz, NxNyNz, Nzh1, NxNyNz1;
double *kxyz, dx, dy, dz;
int andmix, andmixsteps;
double wopt, wcmp, errMax;
double fAinit, fBinit;
int aismatrix;
double q1, q2;

char FEname[30], phname[30];

int main(int argc, char **argv)
{
	double *wA, *wB, *eta, *phA, *phB;
	int i, j, k, in, iseed = -3; //local_x_starti;
	double lylx, lzlx;
	long ijk;

	//MPI_Status status
	FILE *fp;
	time_t ts;
	iseed = time(&ts);
	srand48(iseed);

	//////////put in para ///////////////////////
	fp = fopen("para", "r");
	fscanf(fp, "%d", &in);
	fscanf(fp, "%d, %d", &andmix, &andmixsteps);
	fscanf(fp, "%lf, %lf", &wopt, &wcmp);
	fscanf(fp, "%lf", &errMax);
	fscanf(fp, "%lf, %lf", &lylx, &lzlx);
	fscanf(fp, "%d", &aismatrix);
	fscanf(fp, "%lf", &hAB);
	fscanf(fp, "%lf, %lf", &fA1, &fB);
	fscanf(fp, "%lf, %lf, %lf", &lx, &ly, &lz);
	fscanf(fp, "%d, %d, %d", &Nx, &Ny, &Nz);
	fscanf(fp, "%s", FEname); //output file name for parameters;
	fscanf(fp, "%s", phname); //output file name for configuration;
	fscanf(fp, "%lf", &ds0);
	fscanf(fp, "%lf, %lf", &epA, &epB);
	fclose(fp);

	NxNyNz = Nx * Ny * Nz;

	Nzh1 = Nz / 2 + 1;
	NxNyNz1 = Nx * Ny * Nzh1;
	double kx[Nx], ky[Ny], kz[Nz];

	wA = (double *)malloc(sizeof(double) * NxNyNz);
	wB = (double *)malloc(sizeof(double) * NxNyNz);
	phA = (double *)malloc(sizeof(double) * NxNyNz);
	phB = (double *)malloc(sizeof(double) * NxNyNz);
	eta = (double *)malloc(sizeof(double) * NxNyNz);
	kxyz = (double *)malloc(sizeof(double) * NxNyNz);

	if (lylx != 0)
		ly = lx * sqrt(lylx);
	if (lzlx != 0)
		lz = lx * sqrt(lzlx);
	dx = lx / Nx;
	dy = ly / Ny;
	dz = lz / Nz;

	fA2 = 1.0 - fA1 - fB;

	printf("hAB = %.3lf\n", hAB);
	printf("fA1 = %.3lf, fB = %.3lf, fA2 = %.3lf\n", fA1, fB, fA2);
	printf("dx = %.3lf, dy = %.3lf, dz = %.3lf\n", dx, dy, dz);

	NsA1 = ((int)(fA1 / ds0 + 1.0e-6));
	NsB = ((int)(fB / ds0 + 1.0e-6));
	NsA2 = ((int)(fA2 / ds0 + 1.0e-6));

	fp = fopen(FEname, "w");
	fprintf(fp, "Nx = %d, Ny = %d, Nz = %d\n", Nx, Ny, Nz);
	fprintf(fp, "hAB = %lf\n", hAB);
	fprintf(fp, "fA1 = %lf, fB = %lf, fA2 = %lf\n", fA1, fB, fA2);
	fprintf(fp, "NsA1 = %d, NsB = %d, NsA2 = %d\n", NsA1, NsB, NsA2);
	fprintf(fp, "dx = %.6lf, dy = %.6lf, dz = %.6lf\n", dx, dy, dz);
	fclose(fp);

	for (i = 0; i <= Nx / 2 - 1; i++)
		kx[i] = 2 * Pi * i * 1.0 / Nx / dx;
	for (i = Nx / 2; i < Nx; i++)
		kx[i] = 2 * Pi * (i - Nx) * 1.0 / dx / Nx;
	for (i = 0; i < Nx; i++)
		kx[i] *= kx[i];

	for (i = 0; i <= Ny / 2 - 1; i++)
		ky[i] = 2 * Pi * i * 1.0 / Ny / dy;
	for (i = Ny / 2; i < Ny; i++)
		ky[i] = 2 * Pi * (i - Ny) * 1.0 / dy / Ny;
	for (i = 0; i < Ny; i++)
		ky[i] *= ky[i];

	for (i = 0; i <= Nz / 2 - 1; i++)
		kz[i] = 2 * Pi * i * 1.0 / Nz / dz;
	for (i = Nz / 2; i < Nz; i++)
		kz[i] = 2 * Pi * (i - Nz) * 1.0 / dz / Nz;
	for (i = 0; i < Nz; i++)
		kz[i] *= kz[i];

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				ijk = (long)((i * Ny + j) * Nz + k);
				kxyz[ijk] = kx[i] + ky[j] + kz[k];
			}

	/***************Initialize wA, wB******************/
	if (aismatrix == 0)
	{
		fAinit = fA1 + fA2;
		fBinit = fB;
		init(in, wA, wB);
	}
	else if (aismatrix == 1)
	{
		fAinit = fB;
		fBinit = fA1 + fA2;
		init(in, wB, wA);
	}

	freeE(wA, wB, phA, phB, eta);

	free(wA);
	free(wB);
	free(phA);
	free(phB);
	free(eta);
	free(kxyz);

	return 1;
}

void init(int in, double *wA, double *wB)
{
	FILE *fp;
	long ijk;
	int i, j, k;
	double e1, e2, e3, e4;

	if (in == 0)
		initW_Random(wA, wB);
	else if (in == 233)
		initW_BCC(wA, wB);
	else if (in == 1)
	{
		fp = fopen("in.d", "r");
		for (i = 0; i < Nx; i++)
			for (j = 0; j < Ny; j++)
				for (k = 0; k < Nz; k++)
				{
					fscanf(fp, "%lf %lf %lf %lf", &e1, &e2, &e3, &e4);
					ijk = (long)((i * Ny + j) * Nz + k);
					wA[ijk] = e3;
					wB[ijk] = e4;
				}
		fclose(fp);
	}
	else if (in == 2)
	{
		fp = fopen("in.d", "r");
		for (i = 0; i < Nx; i++)
			for (j = 0; j < Ny; j++)
				for (k = 0; k < Nz; k++)
				{
					fscanf(fp, "%lf %lf %lf %lf", &e1, &e2, &e3, &e4);

					ijk = (long)((i * Ny + j) * Nz + k);
					wA[ijk] = hAB * e2;
					wB[ijk] = hAB * e1;
				}
		fclose(fp);
	}
}

void initW_Random(double *wA, double *wB)
{
	int i, j, k;
	long ijk;

	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				ijk = (long)((i * Ny + j) * Nz + k);

				wA[ijk] = hAB * fB + 0.10 * (drand48() - 0.5);
				wB[ijk] = hAB * (fA1 + fA2) + 0.10 * (drand48() - 0.5);
			}
	}
}

void initW_BCC(double *wA, double *wB)
{
	int i, j, k, nc, tag;
	long ijk;
	double xij, yij, zij;
	double xc[9], yc[9], zc[9];
	double xi, yj, zk, rij, r;
	double phat, phbt;
	FILE *fp;
	fp = fopen("init_BCC.dat", "w");

	r = pow((fAinit / (fAinit + fBinit) * lx * ly * lz / (4 * 4.0 * Pi / 3.0)), 1.0 / 3);

	xc[0] = 0.0;
	yc[0] = 0.0;
	zc[0] = 0.0;
	xc[1] = 0.0;
	yc[1] = 0.0;
	zc[1] = lz;
	xc[2] = 0.0;
	yc[2] = ly;
	zc[2] = 0.0;
	xc[3] = 0.0;
	yc[3] = ly;
	zc[3] = lz;
	xc[4] = lx;
	yc[4] = 0.0;
	zc[4] = lz;
	xc[5] = lx;
	yc[5] = ly;
	zc[5] = 0.0;
	xc[6] = lx;
	yc[6] = ly;
	zc[6] = lz;
	xc[7] = lx;
	yc[7] = 0.0;
	zc[7] = 0.0;
	xc[8] = lx / 2;
	yc[8] = ly / 2;
	zc[8] = lz / 2;

	for (i = 0; i < Nx; i++)
	{
		xi = i * dx;
		for (j = 0; j < Ny; j++)
		{
			yj = j * dy;
			for (k = 0; k < Nz; k++)
			{
				zk = k * dz;
				tag = 0;
				for (nc = 0; nc < 9; nc++)
				{
					xij = xi - xc[nc];
					yij = yj - yc[nc];
					zij = zk - zc[nc];

					rij = xij * xij + yij * yij + zij * zij;
					rij = sqrt(rij);
					if (rij < r)
					{
						tag = 1;
					}
				}
				phat = 0.0;
				phbt = 1.0;

				if (tag == 1)
				{
					phat = 1.0;
					phbt = 0.0;
				}
				ijk = (long)((i * Ny + j) * Nz + k);
				wA[ijk] = hAB * phbt + 0.040 * (drand48() - 0.5);
				wB[ijk] = hAB * phat + 0.040 * (drand48() - 0.5);
				if (aismatrix == 0)
					fprintf(fp, "%lf %lf %lf %lf\n", phat, phbt, wA[ijk], wB[ijk]);
				else if (aismatrix == 1)
					fprintf(fp, "%lf %lf %lf %lf\n", phbt, phat, wB[ijk], wA[ijk]);
			}
		}
	}
	fclose(fp);
}

//********************Output configuration******************************

void write_ph(double *phA, double *phB, double *wA, double *wB)
{
	int i, j, k;
	long ijk;
	FILE *fp = fopen(phname, "w");
	//	FILE *fpp=fopen("result.dat","w");
	//	fprintf(fp,"Nx=%d, Ny=%d, Nz=%d\n",Nx,Ny,Nz);
	//	fprintf(fp,"dx=%lf, dy=%lf, dz=%lf\n",dx,dy,dz);

	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny; j++)
		{
			for (k = 0; k < Nz; k++)
			{
				ijk = (long)((i * Ny + j) * Nz + k);
				fprintf(fp, "%lf %lf %lf %lf\n",
						phA[ijk], phB[ijk], wA[ijk], wB[ijk]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

//*************************************main loop****************************************

double freeE(double *wA, double *wB, double *phA, double *phB, double *eta)
{
	int i, j, k, iter, maxIter;
	long ijk;
	double freeEnergy, freeOld, qC;
	double freeW, freeAB, freeS, freeDiff;
	double Sm1, Sm2, beta, psum, fpsum, *psuC;
	double *waDiff, *wbDiff, inCompMax;
	double *del, *outs, *wAnew, *wBnew, err;
	int N_rec;
	FILE *fp;

	psuC = (double *)malloc(sizeof(double) * NxNyNz);
	waDiff = (double *)malloc(sizeof(double) * NxNyNz);
	wbDiff = (double *)malloc(sizeof(double) * NxNyNz);
	wAnew = (double *)malloc(sizeof(double) * NxNyNz);
	wBnew = (double *)malloc(sizeof(double) * NxNyNz);
	del = (double *)malloc(sizeof(double) * N_hist * 2 * NxNyNz);
	outs = (double *)malloc(sizeof(double) * N_hist * 2 * NxNyNz);

	Sm1 = 5.0e-6;
	Sm2 = 0.1e-7;
	maxIter = MaxIT;
	beta = 1.0;

	iter = 0;

	freeEnergy = 0.0;

	do
	{
		iter = iter + 1;

		getConc(phA, phB, 1.0, wA, wB);

		freeW = 0.0;
		freeAB = 0.0;
		freeS = 0.0;

		inCompMax = 0.0;

		for (ijk = 0; ijk < NxNyNz; ijk++)
		{
			eta[ijk] = (wA[ijk] + wB[ijk] - hAB) / 2.0;
			psum = 1.0 - phA[ijk] - phB[ijk];
			psuC[ijk] = psum;
			fpsum = fabs(psum);
			if (fpsum > inCompMax)
				inCompMax = fpsum;
			wAnew[ijk] = hAB * phB[ijk] + eta[ijk];
			wBnew[ijk] = hAB * phA[ijk] + eta[ijk];
			waDiff[ijk] = wAnew[ijk] - wA[ijk];
			wbDiff[ijk] = wBnew[ijk] - wB[ijk];
			waDiff[ijk] -= wcmp * psum;
			wbDiff[ijk] -= wcmp * psum;

			freeAB += (hAB * phA[ijk] * phB[ijk]);
			freeW -= (wA[ijk] * phA[ijk] + wB[ijk] * phB[ijk] + eta[ijk] * psum);
		}

		freeAB /= NxNyNz;
		freeW /= NxNyNz;

		freeS = - log(q1) - log(q2);

		freeOld = freeEnergy;
		freeEnergy = freeAB + freeW + freeS;

		//judge the error
		err = error_cal(waDiff, wbDiff, wA, wB);

		//update the history fields, and zero is new fields
		update_flds_hist(waDiff, wbDiff, wAnew, wBnew, del, outs);

		//if achieved some level, anderson-mixing, else simple-mixing
		if ((err > errMax || iter < 3 || andmix == 0) || (andmix == 2 && iter < andmixsteps))
		{
			for (ijk = 0; ijk < NxNyNz; ijk++)
			{
				wA[ijk] += wopt * waDiff[ijk];
				wB[ijk] += wopt * wbDiff[ijk];
			}
		}
		else
		{
			printf("iter  %4d  err  %.4f /***** enter Anderson mixing *****/\n", iter, err);
			N_rec = (iter - 1) < N_hist ? (iter - 1) : N_hist;
			Anderson_mixing(del, outs, N_rec, wA, wB);
		}

		//**** print out the free energy and error results ****

		if (iter == 1 || iter % 20 == 0 || iter >= maxIter)
		{
			if (iter == 1)
				fp = fopen("printout.txt", "w");
			else
				fp = fopen("printout.txt", "a");
			fprintf(fp, "%d\n", iter);
			fprintf(fp, "%10.8e, %10.8e, %10.8e, %10.8e, %e\n", freeEnergy, freeAB, freeW, freeS, inCompMax);
			fclose(fp);
		}
		printf(" %5d : %.8e, %.8e, %.8e\n", iter, freeEnergy, inCompMax, err);
		freeDiff = fabs(freeEnergy - freeOld);

		if (iter == 1 || iter % 100 == 0)
			write_ph(phA, phB, wA, wB);
	} while (iter < maxIter && (inCompMax > Sm1 || freeDiff > Sm2));

	fp = fopen("printout.txt", "a");
	fprintf(fp, "%d\n", iter);
	fprintf(fp, "%10.8e, %10.8e, %10.8e, %10.8e, %e\n", freeEnergy, freeAB, freeW, freeS, inCompMax);
	fclose(fp);
	//fclose(fp1);
	write_ph(phA, phB, wA, wB);

	free(psuC);
	free(waDiff);
	free(wbDiff);
	free(wAnew);
	free(wBnew);
	free(del);
	free(outs);

	return freeDiff;
}

double getConc(phlA, phlB, phs0, wA, wB) 
double *phlA, *phlB, phs0;
double *wA, *wB;
{
	int i, j, k, iz;
	long ijk, ijkiz;
	double *qA1, *qcA1, *qB, *qcB, *qA2, *qcA2;
	double ffl1, ffl2, *qInt, qtmp;

	qA1 = (double *)malloc(sizeof(double) * NxNyNz * (NsA1 + 1));
	qcA1 = (double *)malloc(sizeof(double) * NxNyNz * (NsA1 + 1));
	qB = (double *)malloc(sizeof(double) * NxNyNz * (NsB + 1));
	qcB = (double *)malloc(sizeof(double) * NxNyNz * (NsB + 1));
	qA2 = (double *)malloc(sizeof(double) * NxNyNz * (NsA2 + 1));
	qcA2 = (double *)malloc(sizeof(double) * NxNyNz * (NsA2 + 1));
	qInt = (double *)malloc(sizeof(double) * NxNyNz);

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		qInt[ijk] = 1.0;
	}

	sovDifFft(qA1, wA, qInt, fA1, NsA1, 1, epA);
	sovDifFft(qcB, wB, qInt, fB, NsB, -1, epB);
	sovDifFft(qA2, wA, qInt, fA2, NsA2, 1, epA);
	sovDifFft(qcA2, wA, qInt, fA2, NsA2, -1, epA);

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		qInt[ijk] = qA1[ijk * (NsA1 + 1) + NsA1];
	}

	sovDifFft(qB, wB, qInt, fB, NsB, 1, epB);

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		qInt[ijk] = qcB[ijk * (NsB + 1)];
	}

	sovDifFft(qcA1, wA, qInt, fA1, NsA1, -1, epA);

	q1 = 0.0;
	q2 = 0.0;

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		q1 += qcA1[ijk * (NsA1 + 1)];
		q2 += qcA2[ijk * (NsA2 + 1)];
	}

	q1 /= NxNyNz;
	q2 /= NxNyNz;

	ffl1 = phs0 / q1 * ds0;
	ffl2 = phs0 / q2 * ds0;

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		phlA[ijk] = 0.0;
		phlB[ijk] = 0.0;

		ZDIMM = NsA1 + 1;
		for (iz = 0; iz <= NsA1; iz++)
		{
			ijkiz = ijk * ZDIMM + iz;
			if (iz == 0 || iz == NsA1)
				phlA[ijk] += (0.50 * qA1[ijkiz] * qcA1[ijkiz] * ffl1);
			else
				phlA[ijk] += (qA1[ijkiz] * qcA1[ijkiz] * ffl1);
		}

		ZDIMM = NsB + 1;
		for (iz = 0; iz <= NsB; iz++)
		{
			ijkiz = ijk * ZDIMM + iz;
			if (iz == 0 || iz == NsB)
				phlB[ijk] += (0.50 * qB[ijkiz] * qcB[ijkiz] * ffl1);
			else
				phlB[ijk] += (qB[ijkiz] * qcB[ijkiz] * ffl1);
		}

		ZDIMM = NsA2 + 1;
		for (iz = 0; iz <= NsA2; iz++)
		{
			ijkiz = ijk * ZDIMM + iz;
			if (iz == 0 || iz == NsA2)
				phlA[ijk] += (0.50 * qA2[ijkiz] * qcA2[ijkiz] * ffl2);
			else
				phlA[ijk] += (qA2[ijkiz] * qcA2[ijkiz] * ffl2);
		}
	}
	free(qA1);
	free(qA2);
	free(qcA1);
	free(qcA2);
	free(qB);
	free(qcB);
	free(qInt);
}

void sovDifFft(double *g, double *w, double *qInt, double z, int ns, int sign, double epk)
{
	int i, j, k, iz;
	unsigned long ijk, ijkr;
	double dzc, *wdz;
	double *kxyzdz, dzc2;
	double *in;
	fftw_complex *out;
	fftw_plan p_forward, p_backward;

	wdz = (double *)malloc(sizeof(double) * NxNyNz);
	kxyzdz = (double *)malloc(sizeof(double) * NxNyNz);
	in = (double *)malloc(sizeof(double) * NxNyNz);

	out = (fftw_complex *)malloc(sizeof(fftw_complex) * NxNyNz1);
	dzc = z / ns;
	dzc2 = 0.50 * dzc;
	ZDIMM = ns + 1;
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				ijk = (i * Ny + j) * Nz + k;
				kxyzdz[ijk] = exp(-dzc * kxyz[ijk] * epk);
				wdz[ijk] = exp(-w[ijk] * dzc2);
			}
	p_forward = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, in, out, FFTW_ESTIMATE);
	p_backward = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, out, in, FFTW_ESTIMATE);
	if (sign == 1)
	{
		for (ijk = 0; ijk < NxNyNz; ijk++)
		{
			g[ijk * ZDIMM] = qInt[ijk];
		}

		for (iz = 1; iz <= ns; iz++)
		{
			for (ijk = 0; ijk < NxNyNz; ijk++)
			{
				in[ijk] = g[ijk * ZDIMM + iz - 1] * wdz[ijk];
			}

			fftw_execute(p_forward);

			for (i = 0; i < Nx; i++)
				for (j = 0; j < Ny; j++)
					for (k = 0; k < Nzh1; k++)
					{
						ijk = (i * Ny + j) * Nzh1 + k;
						ijkr = (i * Ny + j) * Nz + k;
						out[ijk][0] *= kxyzdz[ijkr]; //out[].re or .im for fftw2
						out[ijk][1] *= kxyzdz[ijkr]; //out[][0] or [1] for fftw3
					}

			fftw_execute(p_backward);

			for (ijk = 0; ijk < NxNyNz; ijk++)
			{
				g[ijk * ZDIMM + iz] = in[ijk] * wdz[ijk] / NxNyNz;
			}
		}
	}
	else
	{

		for (ijk = 0; ijk < NxNyNz; ijk++)
		{
			g[ijk * ZDIMM + ns] = qInt[ijk];
		}

		for (iz = ns - 1; iz >= 0; iz--)
		{
			for (ijk = 0; ijk < NxNyNz; ijk++)
			{
				in[ijk] = g[ijk * ZDIMM + iz + 1] * wdz[ijk];
			}

			fftw_execute(p_forward);

			for (i = 0; i < Nx; i++)
				for (j = 0; j < Ny; j++)
					for (k = 0; k < Nzh1; k++)
					{
						ijk = (i * Ny + j) * Nzh1 + k;
						ijkr = (i * Ny + j) * Nz + k;
						out[ijk][0] *= kxyzdz[ijkr];
						out[ijk][1] *= kxyzdz[ijkr];
					}

			fftw_execute(p_backward);

			for (ijk = 0; ijk < NxNyNz; ijk++)
			{
				//ijk=(i*Ny+j)*Nz2+k;
				//ijkr=(i*Ny+j)*Nz+k;
				g[ijk * ZDIMM + iz] = in[ijk] * wdz[ijk] / NxNyNz;
			}
		}
	}
	fftw_destroy_plan(p_forward);
	fftw_destroy_plan(p_backward);
	free(wdz);
	free(kxyzdz);
	free(in);
	free(out);
}

double error_cal(double *waDiffs, double *wbDiffs, double *wAs, double *wBs)
{
	double err_dif, err_w, err;
	int ijk;

	err = 0.0;
	err_dif = 0.0;
	err_w = 0.0;
	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		err_dif += pow(waDiffs[ijk], 2) + pow(wbDiffs[ijk], 2);
		err_w += pow(wAs[ijk], 2) + pow(wBs[ijk], 2);
	}
	err = err_dif / err_w;
	err = sqrt(err);

	return err;
}

void update_flds_hist(double *waDiff, double *wbDiff, double *wAnew, double *wBnew, double *del, double *outs)
{
	int ijk, j;

	for (j = N_hist - 1; j > 0; j--)
	{
		for (ijk = 0; ijk < NxNyNz; ijk++)
		{
			Del(0, ijk, j) = Del(0, ijk, j - 1);
			Del(1, ijk, j) = Del(1, ijk, j - 1);
			Outs(0, ijk, j) = Outs(0, ijk, j - 1);
			Outs(1, ijk, j) = Outs(1, ijk, j - 1);
		}
	}

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		Del(0, ijk, 0) = waDiff[ijk];
		Del(1, ijk, 0) = wbDiff[ijk];
		Outs(0, ijk, 0) = wAnew[ijk];
		Outs(1, ijk, 0) = wBnew[ijk];
	}
}

/*********************************************************************/
/*
  Anderson mixing [O(Nx)]

  CHECKED
*/

void Anderson_mixing(double *del, double *outs, int N_rec, double *wA, double *wB)
{
	int i, k, ijk;
	int n, m;
	double *up, *vp, *ap;
	int s;

	gsl_matrix_view uGnu;
	gsl_vector_view vGnu, aGnu;
	gsl_permutation *p;

	up = (double *)malloc(sizeof(double) * (N_rec - 1) * (N_rec - 1));
	vp = (double *)malloc(sizeof(double) * (N_rec - 1));
	ap = (double *)malloc(sizeof(double) * (N_rec - 1));

	/* 
	     Calculate the U-matrix and the V-vector 
     
		Follow Shuang, and add the A and B components together.
  	*/

	for (n = 1; n < N_rec; n++)
	{
		V(n) = 0.0;

		for (ijk = 0; ijk < NxNyNz; ijk++)
		{
			V(n) += (Del(0, ijk, 0) - Del(0, ijk, n)) * Del(0, ijk, 0);
			V(n) += (Del(1, ijk, 0) - Del(1, ijk, n)) * Del(1, ijk, 0);
		}

		for (m = n; m < N_rec; m++)
		{
			U(n, m) = 0.0;
			for (ijk = 0; ijk < NxNyNz; ijk++)
			{
				U(n, m) += (Del(0, ijk, 0) - Del(0, ijk, n)) * (Del(0, ijk, 0) - Del(0, ijk, m));
				U(n, m) += (Del(1, ijk, 0) - Del(1, ijk, n)) * (Del(1, ijk, 0) - Del(1, ijk, m));
			}
			U(m, n) = U(n, m);
		}
	}

	/* Calculate A - uses GNU LU decomposition for U A = V */

	uGnu = gsl_matrix_view_array(up, N_rec - 1, N_rec - 1);
	vGnu = gsl_vector_view_array(vp, N_rec - 1);
	aGnu = gsl_vector_view_array(ap, N_rec - 1);

	p = gsl_permutation_alloc(N_rec - 1);

	gsl_linalg_LU_decomp(&uGnu.matrix, p, &s);

	gsl_linalg_LU_solve(&uGnu.matrix, p, &vGnu.vector, &aGnu.vector);

	gsl_permutation_free(p);

	/* Update omega */

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		wA[ijk] = Outs(0, ijk, 0);
		wB[ijk] = Outs(1, ijk, 0);

		for (n = 1; n < N_rec; n++)
		{
			wA[ijk] += A(n) * (Outs(0, ijk, n) - Outs(0, ijk, 0));
			wB[ijk] += A(n) * (Outs(1, ijk, n) - Outs(1, ijk, 0));
		}
	}

	free(ap);
	free(vp);
	free(up);
}
