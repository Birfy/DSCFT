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

double freeE(double *wA, double *wB, double *phA, double *phB, double *eta, double *we);
double getConc(double *phlA, double *phlB, double phs0, double *wA, double *wB, double *we);
void sovDifFft(double *g, double *w, double *qInt, double z, int ns, int sign);
void write_ph(double *phA, double *phB, double *wA, double *wB);

void initW_Random(double *wA, double *wB);
void initW_w(double *we);

int ZDIMM, NsA, NsB;
double lx, ly, ds0;
double hAB, fA, fB;
int Nx, Ny, NxNy, Nyh1, NxNy1;
double *kxyz, dx, dy;
double wopt, wcmp;

char FEname[30], phname[30];

int main(int argc, char **argv)
{
	double *wA, *wB, *eta, *phA, *phB, *we;
	int i, j, k, iseed = -3; //local_x_starti;
	double lylx;
	long ijk;

	//MPI_Status status
	FILE *fp;
	time_t ts;
	iseed = time(&ts);
	srand48(iseed);

	//////////put in para ///////////////////////
	fp = fopen("para", "r");
	fscanf(fp, "%lf, %lf", &wopt, &wcmp);
	fscanf(fp, "%lf", &lylx);
	fscanf(fp, "%lf", &hAB);
	fscanf(fp, "%lf", &fA);
	fscanf(fp, "%lf, %lf", &lx, &ly);
	fscanf(fp, "%d, %d", &Nx, &Ny);
	fscanf(fp, "%s", FEname); //output file name for parameters;
	fscanf(fp, "%s", phname); //output file name for configuration;
	fscanf(fp, "%lf", &ds0);
	fclose(fp);

	NxNy = Nx * Ny;

	Nyh1 = Ny / 2 + 1;
	NxNy1 = Nx * Nyh1;
	double kx[Nx], ky[Ny];

	wA = (double *)malloc(sizeof(double) * NxNy);
	wB = (double *)malloc(sizeof(double) * NxNy);
    we = (double *)malloc(sizeof(double) * NxNy);
	phA = (double *)malloc(sizeof(double) * NxNy);
	phB = (double *)malloc(sizeof(double) * NxNy);
	eta = (double *)malloc(sizeof(double) * NxNy);
	kxyz = (double *)malloc(sizeof(double) * NxNy);

	if (lylx != 0)
		ly = lx * sqrt(lylx);
	dx = lx / Nx;
	dy = ly / Ny;
	fB = 1.0 - fA;

	printf("hAB = %.3lf\n", hAB);
	printf("fA = %.3lf, fB = %.3lf\n", fA, fB);
	printf("dx = %.3lf, dy = %.3lf\n", dx, dy);

	NsA = ((int)(fA / ds0 + 1.0e-6));
	NsB = ((int)(fB / ds0 + 1.0e-6));

	fp = fopen(FEname, "w");
	fprintf(fp, "Nx = %d, Ny = %d\n", Nx, Ny);
	fprintf(fp, "hAB = %lf\n", hAB);
	fprintf(fp, "fA = %lf, fB = %lf\n", fA, fB);
	fprintf(fp, "NsA = %d, NsB = %d\n", NsA, NsB);
	fprintf(fp, "dx = %.6lf, dy = %.6lf\n", dx, dy);
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

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			{
				ijk = (long)(i * Ny + j);
				kxyz[ijk] = kx[i] + ky[j];
			}

	/***************Initialize wA, wB******************/
	initW_Random(wA, wB);

    initW_w(we);

	freeE(wA, wB, phA, phB, eta, we);

    free(we);
	free(wA);
	free(wB);
	free(phA);
	free(phB);
	free(eta);
	free(kxyz);

	return 1;
}

void initW_Random(double *wA, double *wB)
{
	int i, j, k;
	long ijk;

	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny; j++)
			{
				ijk = (long)(i * Ny + j);

				wA[ijk] = hAB * fB + 0.10 * (drand48() - 0.5);
				wB[ijk] = hAB * fA + 0.10 * (drand48() - 0.5);
			}
	}
}

void initW_w(double *we)
{
    int i, j;
    long ijk;
    int sf[256];
    
    for (i = 0; i < 16; i++)
	{
		for (j = 0; j < 16; j++)
		{
			if (i % 2 == 0 || j % 2 == 0)
				sf[16 * i + j] = 1;
			else
				sf[16 * i + j] = 0;
		}
	}

    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
        {
            ijk = (long)(i * Ny + j);
            
            we[ijk] = 20.0 * sf[(i / (Nx / 16)) * 16 + j / (Ny / 16)];
        }
    }
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
				ijk = (long)(i * Ny + j);
				fprintf(fp, "%lf %lf %lf %lf\n",
						phA[ijk], phB[ijk], wA[ijk], wB[ijk]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

//*************************************main loop****************************************

double freeE(double *wA, double *wB, double *phA, double *phB, double *eta, double *we)
{
	int i, j, k, iter, maxIter;
	long ijk;
	double freeEnergy, freeOld, qC, freeContact;
	double freeW, freeAB, freeS, freeDiff;
	double Sm1, Sm2, beta, psum, fpsum;
	double *waDiff, *wbDiff, inCompMax;
	double *wAnew, *wBnew;
	FILE *fp;

	waDiff = (double *)malloc(sizeof(double) * NxNy);
	wbDiff = (double *)malloc(sizeof(double) * NxNy);
	wAnew = (double *)malloc(sizeof(double) * NxNy);
	wBnew = (double *)malloc(sizeof(double) * NxNy);

	Sm1 = 5.0e-6;
	Sm2 = 0.1e-7;
	maxIter = MaxIT;
	beta = 1.0;

	iter = 0;

	freeEnergy = 0.0;

	do
	{
		iter = iter + 1;

		qC = getConc(phA, phB, 1.0, wA, wB, we);

		freeW = 0.0;
		freeAB = 0.0;
		freeS = 0.0;
        freeContact = 0.0;

		inCompMax = 0.0;

		for (ijk = 0; ijk < NxNy; ijk++)
		{
			eta[ijk] = (wA[ijk] + wB[ijk] - hAB) / 2.0;
			psum = 1.0 - phA[ijk] - phB[ijk];
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
            freeContact = we[ijk] * (phB[ijk] - phA[ijk]);
		}

		freeAB /= NxNy;
		freeW /= NxNy;
        freeContact /= NxNy;

		freeS = -log(qC);

		freeOld = freeEnergy;
		freeEnergy = freeAB + freeW + freeS + freeContact;
		
		for (ijk = 0; ijk < NxNy; ijk++)
		{
			wA[ijk] += wopt * waDiff[ijk];
			wB[ijk] += wopt * wbDiff[ijk];
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
		
		printf(" %5d : %.8e, %.8e\n", iter, freeEnergy, inCompMax);
		
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

	free(waDiff);
	free(wbDiff);
	free(wAnew);
	free(wBnew);

	return freeDiff;
}

double getConc(double *phlA, double *phlB, double phs0, double *wA, double *wB, double *we) 
{
	int i, j, k, iz;
	long ijk, ijkiz;
	double *qA, *qcA, *qB, *qcB;
	double ql, ffl, *qInt, qtmp;

	qA = (double *)malloc(sizeof(double) * NxNy * (NsA + 1));
	qcA = (double *)malloc(sizeof(double) * NxNy * (NsA + 1));
	qB = (double *)malloc(sizeof(double) * NxNy * (NsB + 1));
	qcB = (double *)malloc(sizeof(double) * NxNy * (NsB + 1));
	qInt = (double *)malloc(sizeof(double) * NxNy);

    for (ijk = 0; ijk < NxNy; ijk++)
    {
        wA[ijk] -= we[ijk];
        wB[ijk] += we[ijk];
    }    

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = 1.0;
	}

	sovDifFft(qA, wA, qInt, fA, NsA, 1);
	sovDifFft(qcB, wB, qInt, fB, NsB, -1);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = qA[ijk * (NsA + 1) + NsA];
	}

	sovDifFft(qB, wB, qInt, fB, NsB, 1);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = qcB[ijk * (NsB + 1)];
	}

	sovDifFft(qcA, wA, qInt, fA, NsA, -1);

	ql = 0.0;
	for (ijk = 0; ijk < NxNy; ijk++)
	{
		ql += qcA[ijk * (NsA + 1)];
	}

	ql /= NxNy;

	ffl = phs0 / ql * ds0;

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		phlA[ijk] = 0.0;
		phlB[ijk] = 0.0;

		ZDIMM = NsA + 1;
		for (iz = 0; iz <= NsA; iz++)
		{
			ijkiz = ijk * ZDIMM + iz;
			if (iz == 0 || iz == NsA)
				phlA[ijk] += (0.50 * qA[ijkiz] * qcA[ijkiz] * ffl);
			else
				phlA[ijk] += (qA[ijkiz] * qcA[ijkiz] * ffl);
		}

		ZDIMM = NsB + 1;
		for (iz = 0; iz <= NsB; iz++)
		{
			ijkiz = ijk * ZDIMM + iz;
			if (iz == 0 || iz == NsB)
				phlB[ijk] += (0.50 * qB[ijkiz] * qcB[ijkiz] * ffl);
			else
				phlB[ijk] += (qB[ijkiz] * qcB[ijkiz] * ffl);
		}
	}
	free(qA);
	free(qcA);
	free(qB);
	free(qcB);
	free(qInt);
	free(wA);
	free(wB);

	return ql;
}

void sovDifFft(double *g, double *w, double *qInt, double z, int ns, int sign)
{
	int i, j, k, iz;
	unsigned long ijk, ijkr;
	double dzc, *wdz;
	double *kxyzdz, dzc2;
	double *in;
	fftw_complex *out;
	fftw_plan p_forward, p_backward;

	wdz = (double *)malloc(sizeof(double) * NxNy);
	kxyzdz = (double *)malloc(sizeof(double) * NxNy);
	in = (double *)malloc(sizeof(double) * NxNy);

	out = (fftw_complex *)malloc(sizeof(fftw_complex) * NxNy1);
	dzc = z / ns;
	dzc2 = 0.50 * dzc;
	ZDIMM = ns + 1;
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			{
				ijk = i * Ny + j;
				kxyzdz[ijk] = exp(-dzc * kxyz[ijk]);
				wdz[ijk] = exp(-w[ijk] * dzc2);
			}
	p_forward = fftw_plan_dft_r2c_2d(Nx, Ny, in, out, FFTW_ESTIMATE);
	p_backward = fftw_plan_dft_c2r_2d(Nx, Ny, out, in, FFTW_ESTIMATE);
	if (sign == 1)
	{
		for (ijk = 0; ijk < NxNy; ijk++)
		{
			g[ijk * ZDIMM] = qInt[ijk];
		}

		for (iz = 1; iz <= ns; iz++)
		{
			for (ijk = 0; ijk < NxNy; ijk++)
			{
				in[ijk] = g[ijk * ZDIMM + iz - 1] * wdz[ijk];
			}

			fftw_execute(p_forward);

			for (i = 0; i < Nx; i++)
				for (j = 0; j < Nyh1; j++)
					{
						ijk = i * Nyh1 + j;
						ijkr = i * Ny + j;
						out[ijk][0] *= kxyzdz[ijkr]; //out[].re or .im for fftw2
						out[ijk][1] *= kxyzdz[ijkr]; //out[][0] or [1] for fftw3
					}

			fftw_execute(p_backward);

			for (ijk = 0; ijk < NxNy; ijk++)
			{
				g[ijk * ZDIMM + iz] = in[ijk] * wdz[ijk] / NxNy;
			}
		}
	}
	else
	{

		for (ijk = 0; ijk < NxNy; ijk++)
		{
			g[ijk * ZDIMM + ns] = qInt[ijk];
		}

		for (iz = ns - 1; iz >= 0; iz--)
		{
			for (ijk = 0; ijk < NxNy; ijk++)
			{
				in[ijk] = g[ijk * ZDIMM + iz + 1] * wdz[ijk];
			}

			fftw_execute(p_forward);

			for (i = 0; i < Nx; i++)
				for (j = 0; j < Nyh1; j++)
					{
						ijk = i * Nyh1 + j;
						ijkr = i * Ny + j;
						out[ijk][0] *= kxyzdz[ijkr];
						out[ijk][1] *= kxyzdz[ijkr];
					}

			fftw_execute(p_backward);

			for (ijk = 0; ijk < NxNy; ijk++)
			{
				//ijk=(i*Ny+j)*Nz2+k;
				//ijkr=(i*Ny+j)*Nz+k;
				g[ijk * ZDIMM + iz] = in[ijk] * wdz[ijk] / NxNy;
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
