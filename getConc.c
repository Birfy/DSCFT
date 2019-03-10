#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include </export/home/cdy/fftw3/include/fftw3.h>
#include </data2/opt/include/gsl/gsl_blas.h>
#include </data2/opt/include/gsl/gsl_linalg.h>

#define Pi 3.141592653589

double getConc(double *ph, double phs0, double *wA, double *wB);
void sovDifFft(double *g, double *w, double *qInt, double z, int ns, int sign);
void write_ph(double *ph);

void init(double *wA, double *wB);

int ZDIMM, NsA1, NsB, NsA2;
double lx, ly, ds0;
double hAB, fA, fB, fA1, fA2;
int Nx, Ny, NxNy, Nyh1, NxNy1;
double *kxyz, dx, dy;
double q1, q2;

int main(int argc, char **argv)
{
	double *wA, *wB, *ph;
	int i, j, k, in, iseed = -3;
	double lylx;
	long ijk;
	double ksq;

	FILE *fp;

	fp = fopen("para_conc", "r");
	fscanf(fp, "%lf", &lylx);
	fscanf(fp, "%lf", &hAB);
	fscanf(fp, "%lf, %lf", &fA1, &fB);
	fscanf(fp, "%lf, %lf", &lx, &ly);
	fscanf(fp, "%d, %d", &Nx, &Ny);
	fscanf(fp, "%lf", &ds0);
	fclose(fp);

	NxNy = Nx * Ny;
	Nyh1 = Ny / 2 + 1;
	NxNy1 = Nx * Nyh1;
	double kx[Nx], ky[Ny];

	wA = (double *)malloc(sizeof(double) * NxNy);
	wB = (double *)malloc(sizeof(double) * NxNy);
	ph = (double *)malloc(sizeof(double) * NxNy);
	kxyz = (double *)malloc(sizeof(double) * NxNy);

	if (lylx != 0)
		ly = lx * sqrt(lylx);

	dx = lx / Nx;
	dy = ly / Ny;

	fA2 = 1.0 - fA1 - fB;
	fA = fA1 + fA2;

	printf("hAB = %.3lf\n", hAB);
	printf("fA1 = %.3lf, fB = %.3lf, fA2 = %.3lf\n", fA1, fB, fA2);
	printf("dx = %.3lf, dy = %.3lf\n", dx, dy);

	NsA1 = ((int)(fA1 / ds0 + 1.0e-6));
	NsB = ((int)(fB / ds0 + 1.0e-6));
	NsA2 = ((int)(fA2 / ds0 + 1.0e-6));

	fp = fopen("fet.dat", "w");
	fprintf(fp, "Nx = %d, Ny = %d\n", Nx, Ny);
	fprintf(fp, "hAB = %lf\n", hAB);
	fprintf(fp, "fA1 = %lf, fB = %lf, fA2 = %lf\n", fA1, fB, fA2);
	fprintf(fp, "NsA1 = %d, NsB = %d, NsA2 = %d\n", NsA1, NsB, NsA2);
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

	
	init(wA, wB);

    getConc(ph, 1.0, wA, wB);

    write_ph(ph);

	free(wA);
	free(wB);
	free(ph);
	free(kxyz);

	return 1;
}

void init(double *wA, double *wB)
{
	FILE *fp;
	long ijk;
	double e1, e2, e3, e4;

	
	fp = fopen("in.d", "r");
	for (ijk = 0; ijk < NxNy; ijk++)
	{
		fscanf(fp, "%lf %lf %lf %lf", &e1, &e2, &e3, &e4);
		wA[ijk] = e3;
		wB[ijk] = e4;
	}
	fclose(fp);
}

void write_ph(double *ph)
{
	int i, j, k;
	long ijk;
	FILE *fp = fopen("pha.dat", "w");

	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny; j++)
		{
			ijk = (long)(i * Ny + j);
			fprintf(fp, "%lf\n", ph[ijk]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

double getConc(double *ph, double phs0, double *wA, double *wB)
{
    int i, j, k, iz;
	long ijk, ijkiz;
	double *qA1, *qcA1, *qB, *qcB, *qA2, *qcA2;
	double ffl1, ffl2, *qInt, qtmp;

	qA1 = (double *)malloc(sizeof(double) * NxNy * (NsA1 + 1));
	qcA1 = (double *)malloc(sizeof(double) * NxNy * (NsA1 + 1));
	qB = (double *)malloc(sizeof(double) * NxNy * (NsB + 1));
	qcB = (double *)malloc(sizeof(double) * NxNy * (NsB + 1));
	qA2 = (double *)malloc(sizeof(double) * NxNy * (NsA2 + 1));
	qcA2 = (double *)malloc(sizeof(double) * NxNy * (NsA2 + 1));
	qInt = (double *)malloc(sizeof(double) * NxNy);


	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = 1.0;
	}

	sovDifFft(qA1, wA, qInt, fA1, NsA1, 1);
	sovDifFft(qcB, wB, qInt, fB, NsB, -1);
	sovDifFft(qA2, wA, qInt, fA2, NsA2, 1);
	sovDifFft(qcA2, wA, qInt, fA2, NsA2, -1);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = qA1[ijk * (NsA1 + 1) + NsA1];
	}

	sovDifFft(qB, wB, qInt, fB, NsB, 1);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = qcB[ijk * (NsB + 1)];
	}

	sovDifFft(qcA1, wA, qInt, fA1, NsA1, -1);

	q1 = 0.0;
	q2 = 0.0;
	for (ijk = 0; ijk < NxNy; ijk++)
	{
		q1 += qcA1[ijk * (NsA1 + 1)];
		q2 += qcA2[ijk * (NsA2 + 1)];
	}

	q1 /= NxNy;
	q2 /= NxNy;

	ffl1 = phs0 / q1 * ds0;
	ffl2 = phs0 / q2 * ds0;

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		ph[ijk] = 0.0;

		ZDIMM = NsA2 + 1;
		for (iz = NsA2; iz <= NsA2; iz++)
		{
			ijkiz = ijk * ZDIMM + iz;
			if (iz == 0 || iz == NsA2)
				ph[ijk] += (0.50 * qA1[ijkiz] * qcA1[ijkiz] * ffl1);
			else
				ph[ijk] += (qA1[ijkiz] * qcA1[ijkiz] * ffl1);
		}
	}
	free(qA1);
	free(qA2);
	free(qcA1);
	free(qcA2);
	free(qB);
	free(qcB);
	free(qInt);
    free(wA);
    free(wB);
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
	for (ijk = 0; ijk < NxNy; ijk++)
	{
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
