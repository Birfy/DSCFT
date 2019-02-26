#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include </export/home/cdy/fftw3/include/fftw3.h>
#include </data2/opt/include/gsl/gsl_blas.h>
#include </data2/opt/include/gsl/gsl_linalg.h>
#define MaxIT 50000

#define Pi 3.141592653589

double freeE(double *wm, double *wp, double *phA, double *phB, double *eta);
double getConc(double *phlA, double *phlB, double phs0, double *wm, double *wp);
void sovDifFft(double *g, double *w, double *qInt, double z, int ns, int sign, double epK);
void write_ph(int t, double *phA, double *phB, double *wm, double *wp);
double update_wm(double *phA, double *phB, double *wm);

void init(int in, double *wm, double *wp);
void initW_Random(double *wm, double *wp);

int ZDIMM, NsA1, NsB, NsA2;
double lx, ly, lz, ds0;
double hAB, fA1, fB, fA2;
double epA, epB;
int Nx, Ny, Nz, NxNyNz, Nzh1, NxNyNz1;
double *kxyz, dx, dy, dz;
double wopt, wcmp;
double fAinit, fBinit;
int aismatrix;
double timestep;
int timerepeat;

int main(int argc, char **argv)
{
    double *wm, *wp, *eta, *phA, *phB;
    int i, j, k, in, iseed = -3;
    double lylx, lzlx;
    long ijk;

    FILE *fp;
    time_t ts;
    iseed = time(&ts);
    srand48(iseed);

    fp = fopen("para", "r");
    fscanf(fp, "%d", &in);
	fscanf(fp, "%lf, %lf", &wopt, &wcmp);
	fscanf(fp, "%lf, %lf", &lylx, &lzlx);
	fscanf(fp, "%d", &aismatrix);
	fscanf(fp, "%lf", &hAB);
	fscanf(fp, "%lf, %lf", &fA1, &fB);
	fscanf(fp, "%lf, %lf, %lf", &lx, &ly, &lz);
	fscanf(fp, "%d, %d, %d", &Nx, &Ny, &Nz);
	fscanf(fp, "%lf", &ds0);
	fscanf(fp, "%lf, %lf", &epA, &epB);
    fscanf(fp, "%lf, %d", &timestep, &timerepeat);
	fclose(fp);

    NxNyNz = Nx * Ny * Nz;
    Nzh1 = Nz / 2 + 1;
    NxNyNz1 = Nx * Ny * Nzh1;
    double kx[Nx], ky[Ny], kz[Nz];

    wm = (double *)malloc(sizeof(double) * NxNyNz);
    wp = (double *)malloc(sizeof(double) * NxNyNz);
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
    
    fA2 = 1.0 - fB - fA1;

	printf("hAB = %.3lf\n", hAB);
	printf("fA1 = %.3lf, fB = %.3lf, fA2 = %.3lf\n", fA1, fB, fA2);
	printf("dx = %.3lf, dy = %.3lf, dz = %.3lf\n", dx, dy, dz);

    NsA1 = ((int)(fA1 / ds0 + 1.0e-6));
	NsB = ((int)(fB / ds0 + 1.0e-6));
	NsA2 = ((int)(fA2 / ds0 + 1.0e-6));

	fp = fopen("fet.dat", "w");
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
    
    if (aismatrix == 0)
    {
        fAinit = fA1 + fA2;
        fBinit = fB;
        init(in, wm, wp);
    }
    else if (aismatrix == 1)
    {
        fAinit = fB;
        fBinit = fA1 + fA2;
        init(in, wm, wp);
    }
    
    freeE(wm, wp, phA, phB, eta);

    free(wm);
    free(wp);
    free(phA);
    free(phB);
    free(eta);
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
				for (k = 0; k < Nz; k++)
				{
					fscanf(fp, "%lf %lf %lf %lf", &e1, &e2, &e3, &e4);
					ijk = (long)((i * Ny + j) * Nz + k);
					wm[ijk] = (e3 - e4) / 2;
					wp[ijk] = (e3 + e4) / 2;
				}
		fclose(fp);
    }
}

void initW_Random(double *wm, double *wp)
{
    int i, j, k;
    long ijk;

    for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				ijk = (long)((i * Ny + j) * Nz + k);

				wm[ijk] = hAB * (fA1 + fA2 - fB) + 0.10 * (drand48() - 0.5);
				wp[ijk] = hAB + 0.10 * (drand48() - 0.5);
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
			for (k = 0; k < Nz; k++)
			{
				ijk = (long)((i * Ny + j) * Nz + k);
				fprintf(fp, "%lf %lf %lf %lf\n",
						phA[ijk], phB[ijk], wm[ijk] + wp[ijk], wp[ijk] - wm[ijk]);
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

double freeE(double *wm, double *wp, double *phA, double *phB, double *eta)
{
    int i, j, k, t, iter, maxIter;
    long ijk;
    double freeEnergy, freeOld, qC;
    double freeH, freeS, freeDiff;
    double Sm1, Sm2, beta, psum, fpsum;
    double *wpDiff, inCompMax;
    double *wpnew;
    char poname[30];
    // FILE *fp;

    wpDiff = (double *)malloc(sizeof(double) * NxNyNz);
	wpnew = (double *)malloc(sizeof(double) * NxNyNz);
	
    Sm1 = 5.0e-6;
    Sm2 = 1.0e-8;
    maxIter = MaxIT;
    beta = 1.0;

    for (t = 0; t < timerepeat; t++)
    {
        sprintf(poname, "printout%d.dat", t);
        printf("t = %f\n", t * timestep);

        iter = 0;
		freeEnergy = 0.0;
        do
        {
            iter = iter + 1;
            qC = getConc(phA, phB, 1.0, wm, wp);

            freeH = 0.0;
            freeS = 0.0;

            inCompMax = 0.0;

            for (ijk = 0; ijk < NxNyNz; ijk++)
            {
                eta[ijk] = wp[ijk] - hAB / 2.0;
                psum = 1.0 - phA[ijk] - phB[ijk];
                fpsum = fabs(psum);
                if (fpsum > inCompMax)
                    inCompMax = fpsum;
                wpnew[ijk] = hAB * (phA[ijk] + phB[ijk]) / 2 + eta[ijk];
                wpDiff[ijk] = wpnew[ijk] - wp[ijk];
                wpDiff[ijk] -= wcmp * psum;

                freeH += 1.0 / hAB * wm[ijk] * wm[ijk] - wp[ijk];
            }

            freeH /= NxNyNz;
            
            freeS = -log(qC);

            freeOld = freeEnergy;
            freeEnergy = freeH + freeS;

            for (ijk = 0; ijk < NxNyNz; ijk++)
            {
                wp[ijk] += wopt * wpDiff[ijk];
            }

            // if (iter == 1 || iter % 20 == 0 || iter >= maxIter)
            // {
            //     if (iter == 1)
            //         fp = fopen(poname, "w");
            //     else
            //         fp = fopen(poname, "a");
            //     fprintf(fp, "%d\n", iter);
            //     fprintf(fp, "%10.8e, %10.8e, %10.8e, %e\n", freeEnergy, freeH, freeS, inCompMax);
            //     fclose(fp);
            // }
            printf(" %5d : %.8e, %.8e\n", iter, freeEnergy, inCompMax);
		    freeDiff = fabs(freeEnergy - freeOld);
        } while (iter < maxIter && (inCompMax > Sm1 || freeDiff > Sm2));
        
        // fp = fopen(poname, "a");
        // fprintf(fp, "%d\n", iter);
        // fprintf(fp, "%10.8e, %10.8e, %10.8e, %e\n", freeEnergy, freeH, freeS, inCompMax);
        // fclose(fp);
        

        // for (i = 0; i < Nx; i++)
        // {
        //     for (j = 0; j < Ny; j++)
        //     {
        //         ijk = (long)(i * Ny + j);
        //         wm[ijk] -= timestep * getFlow(i, j, wm, wp, phA, phB);
        //     }
        // }
        if (t % 50 == 0)
			write_ph(t, phA, phB, wm, wp);
        update_wm(phA, phB, wm);
		
        
    }
	free(wpDiff);
	free(wpnew);

	return freeDiff;
}

double update_wm(double *phA, double *phB, double *wm) 
{
    double *dfdwm;
    double dfdx, dfdy, dfdz;
    double flow;
    int i, j, k, ijk;

    dfdwm = (double *)malloc(sizeof(double) * NxNyNz);

    for (ijk = 0; ijk < NxNyNz; ijk++)
    {
        dfdwm[ijk] = 2 / hAB * wm[ijk] + (phA[ijk] - phB[ijk]);
    }
    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
        {
			for (k = 0; k < Nz; k++)
			{
				ijk = (long)((i * Ny + j) * Nz + k);
				dfdx = (dfdwm[(((i + 1) % Nx) * Ny + j) * Nz + k] + dfdwm[(((i - 1 + Nx) % Nx) * Ny + j) * Nz + k] - 2 * dfdwm[(i * Ny + j) * Nz + k]) / dx / dx;
				dfdy = (dfdwm[(i * Ny + ((j + 1) % Ny)) * Nz + k] + dfdwm[(i * Ny + ((j - 1 + Ny) % Ny)) * Nz + k] - 2 * dfdwm[(i * Ny + j) * Nz + k]) / dy / dy;
				dfdz = (dfdwm[(i * Ny + j) * Nz + ((k + 1) % Nz)] + dfdwm[(i * Ny + j) * Nz + ((k - 1 + Nz) % Nz)] - 2 * dfdwm[(i * Ny + j) * Nz + k]) / dz / dz;
				flow = dfdx + dfdy + dfdz;

				// if (flow > 100) 
				// {
				//     flow = 100;
				// }
				// if (flow < 100)
				// {
				//     flow = 100;
				// }

				wm[ijk] += timestep * flow;
			}
        }
    }
    free(dfdwm);
    return flow;
}

double getConc(double *phlA, double *phlB, double phs0, double *wm, double *wp)
{
    double *wA, *wB;
    int i, j, k, iz;
	long ijk, ijkiz;
	double *qA1, *qcA1, *qB, *qcB, *qA2, *qcA2;
	double ql, ffl, *qInt, qtmp;

    wA = (double *)malloc(sizeof(double) * NxNyNz);
    wB = (double *)malloc(sizeof(double) * NxNyNz);
	qA1 = (double *)malloc(sizeof(double) * NxNyNz * (NsA1 + 1));
	qcA1 = (double *)malloc(sizeof(double) * NxNyNz * (NsA1 + 1));
	qB = (double *)malloc(sizeof(double) * NxNyNz * (NsB + 1));
	qcB = (double *)malloc(sizeof(double) * NxNyNz * (NsB + 1));
	qA2 = (double *)malloc(sizeof(double) * NxNyNz * (NsA2 + 1));
	qcA2 = (double *)malloc(sizeof(double) * NxNyNz * (NsA2 + 1));
	qInt = (double *)malloc(sizeof(double) * NxNyNz);

    for (ijk = 0; ijk < NxNyNz; ijk++){
        wA[ijk] = wp[ijk] + wm[ijk];
        wB[ijk] = wp[ijk] - wm[ijk];
    }

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		qInt[ijk] = 1.0;
	}

	sovDifFft(qA1, wA, qInt, fA1, NsA1, 1, epA);
	sovDifFft(qcA2, wA, qInt, fA2, NsA2, -1, epA);

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		qInt[ijk] = qcA2[ijk * (NsA2 + 1)];
	}

	sovDifFft(qcB, wB, qInt, fB, NsB, -1, epB);

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		qInt[ijk] = qcB[ijk * (NsB + 1)];
	}

	sovDifFft(qcA1, wA, qInt, fA1, NsA1, -1, epA);

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		qInt[ijk] = qA1[ijk * (NsA1 + 1) + NsA1];
	}
	sovDifFft(qB, wB, qInt, fB, NsB, 1, epB);

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		qInt[ijk] = qB[ijk * (NsB + 1) + NsB];
	}

	sovDifFft(qA2, wA, qInt, fA2, NsA2, 1, epA);

	ql = 0.0;
	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		ql += qcA1[ijk * (NsA1 + 1)];
	}

	ql /= NxNyNz;

	ffl = phs0 / ql * ds0;

	for (ijk = 0; ijk < NxNyNz; ijk++)
	{
		phlA[ijk] = 0.0;
		phlB[ijk] = 0.0;

		ZDIMM = NsA1 + 1;
		for (iz = 0; iz <= NsA1; iz++)
		{
			ijkiz = ijk * ZDIMM + iz;
			if (iz == 0 || iz == NsA1)
				phlA[ijk] += (0.50 * qA1[ijkiz] * qcA1[ijkiz]);
			else
				phlA[ijk] += (qA1[ijkiz] * qcA1[ijkiz]);
		}

		ZDIMM = NsB + 1;
		for (iz = 0; iz <= NsB; iz++)
		{
			ijkiz = ijk * ZDIMM + iz;
			if (iz == 0 || iz == NsB)
				phlB[ijk] += (0.50 * qB[ijkiz] * qcB[ijkiz]);
			else
				phlB[ijk] += (qB[ijkiz] * qcB[ijkiz]);
		}

		ZDIMM = NsA2 + 1;
		for (iz = 0; iz <= NsA2; iz++)
		{
			ijkiz = ijk * ZDIMM + iz;
			if (iz == 0 || iz == NsA2)
				phlA[ijk] += (0.50 * qA2[ijkiz] * qcA2[ijkiz]);
			else
				phlA[ijk] += (qA2[ijkiz] * qcA2[ijkiz]);
		}

		phlA[ijk] *= ffl;
		phlB[ijk] *= ffl;
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

	return ql;
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
