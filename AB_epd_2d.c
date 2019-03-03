#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include </export/home/cdy/fftw3/include/fftw3.h>
#include </data2/opt/include/gsl/gsl_blas.h>
#include </data2/opt/include/gsl/gsl_linalg.h>

#define Pi 3.141592653589

void freeE(double *wm, double *wp, double *phA, double *phB);
double getConc(double *phlA, double *phlB, double phs0, double *wm, double *wp);
void sovDifFft(double *g, double *w, double *qInt, double z, int ns, int sign, double epK);
void write_ph(int t, double *phA, double *phB, double *wm, double *wp);
void update_wm(double *phA, double *phB, double *wm);
void fftw(double *in, fftw_complex *in_k, int sign);
void wmtowp(double *wm, double *wp, double *phA, double *phB);

void init(int in, double *wm, double *wp);
void initW_Random(double *wm, double *wp);

int ZDIMM, NsA, NsB;
double lx, ly, ds0;
double hAB, fA, fB;
double epA, epB;
int Nx, Ny, NxNy, Nyh1, NxNy1;
double *kxyz, dx, dy;
double fAinit, fBinit;
int aismatrix;
double dt, rdt;
double q;
double *gfact;

int main(int argc, char **argv)
{
	double *wm, *wp, *phA, *phB;
	int i, j, k, in, iseed = -3;
	double lylx;
	long ijk;
	double ksq;

	FILE *fp;
	time_t ts;
	iseed = time(&ts);
	srand48(iseed);

	fp = fopen("para", "r");
	fscanf(fp, "%d", &in);
	fscanf(fp, "%lf", &lylx);
	fscanf(fp, "%d", &aismatrix);
	fscanf(fp, "%lf", &hAB);
	fscanf(fp, "%lf", &fA);
	fscanf(fp, "%lf, %lf", &lx, &ly);
	fscanf(fp, "%d, %d", &Nx, &Ny);
	fscanf(fp, "%lf", &ds0);
	fscanf(fp, "%lf, %lf", &epA, &epB);
	fscanf(fp, "%lf, %lf", &dt, &rdt);
	fclose(fp);

	NxNy = Nx * Ny;
	Nyh1 = Ny / 2 + 1;
	NxNy1 = Nx * Nyh1;
	double kx[Nx], ky[Ny];

	wm = (double *)malloc(sizeof(double) * NxNy);
	wp = (double *)malloc(sizeof(double) * NxNy);
	phA = (double *)malloc(sizeof(double) * NxNy);
	phB = (double *)malloc(sizeof(double) * NxNy);
	kxyz = (double *)malloc(sizeof(double) * NxNy);
	gfact = (double *)malloc(sizeof(double) * NxNy);

	if (lylx != 0)
		ly = lx * sqrt(lylx);

	dx = lx / Nx;
	dy = ly / Ny;

	fB = 1 - fA;

	printf("hAB = %.3lf\n", hAB);
	printf("fA = %.3lf, fB = %.3lf\n", fA, fB);
	printf("dx = %.3lf, dy = %.3lf\n", dx, dy);

	NsA = ((int)(fA / ds0 + 1.0e-6));
	NsB = ((int)(fB / ds0 + 1.0e-6));

	fp = fopen("fet.dat", "w");
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
			ksq = kx[i] + ky[j];

			if (ijk == 0)
				gfact[ijk] = 1.0;
			else
			{
				gfact[ijk] = (fA * ksq + exp(-ksq * fA) - 1.0);
				gfact[ijk] += (1.0 - exp(-ksq * fA)) * (1.0 - exp(-ksq * fB));
				gfact[ijk] += (fB * ksq + exp(-ksq * fB) - 1.0);
				gfact[ijk] *= (2.0 / ksq / ksq);
			}
			gfact[ijk] = dt / (1.0 + dt * gfact[ijk]);
		}

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
		{
			ijk = (long)(i * Ny + j);
			kxyz[ijk] = kx[i] + ky[j];
		}

	if (aismatrix == 0)
	{
		fAinit = fA;
		fBinit = fB;
		init(in, wm, wp);
	}
	else if (aismatrix == 1)
	{
		fAinit = fB;
		fBinit = fA;
		init(in, wm, wp);
	}

	freeE(wm, wp, phA, phB);

	free(wm);
	free(wp);
	free(phA);
	free(phB);
	free(kxyz);

	return 1;
}

void init(int in, double *wm, double *wp)
{
	FILE *fp;
	long ijk;
	int i, j, k;
	double e1, e2, e3, e4;

	if (in == 0)
		initW_Random(wm, wp);
	else if (in == 1)
	{
		fp = fopen("in.d", "r");
		for (i = 0; i < Nx; i++)
			for (j = 0; j < Ny; j++)
			{
				fscanf(fp, "%lf %lf %lf %lf", &e1, &e2, &e3, &e4);
				ijk = (long)(i * Ny + j);
				wm[ijk] = (e4 - e3) + hAB * (e2 - e1) + 0.1 * (drand48() - 0.5);
				wp[ijk] = 0.1 * (drand48() - 0.5);
			}
		fclose(fp);
	}
}

void initW_Random(double *wm, double *wp)
{
	int i, j, k;
	long ijk;

	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny; j++)
		{
			ijk = (long)(i * Ny + j);

			wm[ijk] = 0.1 * (drand48() - 0.5);
			wp[ijk] = 0.1 * (drand48() - 0.5);
		}
	}
}

void write_ph(int t, double *phA, double *phB, double *wm, double *wp)
{
	int i, j, k;
	long ijk;
	char phname[30];
	sprintf(phname, "pha%d.dat", t);
	FILE *fp = fopen(phname, "w");

	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny; j++)
		{
			ijk = (long)(i * Ny + j);
			fprintf(fp, "%lf %lf %lf %lf\n",
					phA[ijk], phB[ijk], wp[ijk] - wm[ijk], wp[ijk] + wm[ijk]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void freeE(double *wm, double *wp, double *phA, double *phB)
{
	int i, j, k, t, iter;
	long ijk;
	double freeEnergy, freeOld;
	double freeW, freeAB, freeS, freeDiff;
	char poname[30];
	double sm;
	// FILE *fp;

	sm = 1e-10;
	t = 0;

	freeEnergy = 0.0;

	do
	{
		sprintf(poname, "printout%d.dat", t);
		// printf("t = %f\n", t * dt);

		wmtowp(wm, wp, phA, phB);

		getConc(phA, phB, 1.0, wm, wp);

		freeW = 0.0;
		freeAB = 0.0;
		freeS = 0.0;

		for (ijk = 0; ijk < NxNy; ijk++)
		{
			freeW += ((2 * fA - 1) * wm[ijk] - wp[ijk]);
			freeAB += wm[ijk] * wm[ijk] / hAB;
		}

		freeW /= NxNy;
		freeAB /= NxNy;

		freeS = -log(q);

		freeOld = freeEnergy;
		freeEnergy = freeW + freeAB + freeS;

		freeDiff = fabs(freeEnergy - freeOld);
		printf(" %5d : %.8e, %.8e\n", t, freeEnergy, freeDiff);

		if (t % 1000 == 0)
			write_ph(t, phA, phB, wm, wp);
		update_wm(phA, phB, wm);

		t += 1;

	} while (freeDiff > sm);
}

void update_wm(double *phA, double *phB, double *wm)
{
	double *dfdwm;
	double dfdx, dfdy;
	double flow;
	int i, j, ijk;

	dfdwm = (double *)malloc(sizeof(double) * NxNy);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		dfdwm[ijk] = (2 * fA - 1) + 2 / hAB * wm[ijk] + (phB[ijk] - phA[ijk]);
	}
	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny; j++)
		{
			ijk = (long)(i * Ny + j);
			dfdx = (dfdwm[((i + 1) % Nx) * Ny + j] + dfdwm[((i - 1 + Nx) % Nx) * Ny + j] - 2 * dfdwm[i * Ny + j]) / dx / dx;
			dfdy = (dfdwm[i * Ny + ((j + 1) % Ny)] + dfdwm[i * Ny + ((j - 1 + Ny) % Ny)] - 2 * dfdwm[i * Ny + j]) / dy / dy;

			flow = dfdx + dfdy;
			wm[ijk] += rdt * flow;
		}
	}
	free(dfdwm);
}

double getConc(double *phlA, double *phlB, double phs0, double *wm, double *wp)
{
	double *wA, *wB;
	int i, j, k, iz;
	long ijk, ijkiz;
	double *qA, *qcA, *qB, *qcB;
	double ffl, *qInt, qtmp;

	wA = (double *)malloc(sizeof(double) * NxNy);
	wB = (double *)malloc(sizeof(double) * NxNy);
	qA = (double *)malloc(sizeof(double) * NxNy * (NsA + 1));
	qcA = (double *)malloc(sizeof(double) * NxNy * (NsA + 1));
	qB = (double *)malloc(sizeof(double) * NxNy * (NsB + 1));
	qcB = (double *)malloc(sizeof(double) * NxNy * (NsB + 1));
	qInt = (double *)malloc(sizeof(double) * NxNy);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		wA[ijk] = wp[ijk] - wm[ijk];
		wB[ijk] = wp[ijk] + wm[ijk];
	}

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = 1.0;
	}

	sovDifFft(qA, wA, qInt, fA, NsA, 1, epA);
	sovDifFft(qcB, wB, qInt, fB, NsB, -1, epB);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = qA[ijk * (NsA + 1) + NsA];
	}

	sovDifFft(qB, wB, qInt, fB, NsB, 1, epB);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = qcB[ijk * (NsB + 1)];
	}

	sovDifFft(qcA, wA, qInt, fA, NsA, -1, epA);

	q = 0.0;
	for (ijk = 0; ijk < NxNy; ijk++)
	{
		q += qcA[ijk * (NsA + 1)];
	}

	q /= NxNy;

	ffl = phs0 / q * ds0;

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
			kxyzdz[ijk] = exp(-dzc * kxyz[ijk] * epk);
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

void fftw(double *in, fftw_complex *out, int sign)
{
	long ijk;
	fftw_plan p_forward, p_backward;

	p_forward = fftw_plan_dft_r2c_2d(Nx, Ny, in, out, FFTW_ESTIMATE);
	p_backward = fftw_plan_dft_c2r_2d(Nx, Ny, out, in, FFTW_ESTIMATE);

	if (sign == 1)
		fftw_execute(p_forward);
	else
	{
		fftw_execute(p_backward);
		for (ijk = 0; ijk < NxNy; ijk++)
			in[ijk] /= NxNy;
	}

	fftw_destroy_plan(p_forward);
	fftw_destroy_plan(p_backward);
}

void wmtowp(double *wm, double *wp, double *phA, double *phB)
{
	int i, j, k, iter, maxIter;
	long ijk, ijkr;
	double inCompMax, fpsuC, qC;
	double *dwp;
	double sm;
	fftw_complex *wpk, *dwpk;

	dwp = (double *)malloc(sizeof(double) * NxNy);

	wpk = (fftw_complex *)malloc(sizeof(fftw_complex) * NxNy1);
	dwpk = (fftw_complex *)malloc(sizeof(fftw_complex) * NxNy1);

	iter = 0;
	maxIter = 200;
	sm = 5.0e-6;

	fftw(wp, wpk, 1);

	do
	{
		iter += 1;

		getConc(phA, phB, 1.0, wm, wp);

		inCompMax = 0.0;

		for (ijk = 0; ijk < NxNy; ijk++)
		{
			dwp[ijk] = phA[ijk] + phB[ijk] - 1.0;
			fpsuC = fabs(dwp[ijk]);
			if (fpsuC > inCompMax)
				inCompMax = fpsuC;
		}

		fftw(dwp, dwpk, 1);

		for (i = 0; i < Nx; i++)
			for (j = 0; j < Nyh1; j++)
			{
				ijk = i * Nyh1 + j;
				ijkr = i * Ny + j;

				dwpk[ijk][0] *= gfact[ijkr];
				dwpk[ijk][1] *= gfact[ijkr];
				dwpk[ijk][0] += wpk[ijk][0];
				dwpk[ijk][1] += wpk[ijk][1];

				wpk[ijk][0] = dwpk[ijk][0];
				wpk[ijk][1] = dwpk[ijk][1];

				if (ijk == 0)
				{
					dwpk[ijk][0] = 0.0;
					dwpk[ijk][1] = 0.0;
				}
			}

		fftw(wp, dwpk, -1);

		// printf(" %5d : %.8e\n", iter, inCompMax);

	} while (iter < maxIter && inCompMax > sm);

	free(dwp);
	free(wpk);
	free(dwpk);
}
