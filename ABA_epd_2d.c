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
double lx, ly, ds0;
double hAB, fA1, fB, fA2;
double epA, epB;
int Nx, Ny, NxNy, Nyh1, NxNy1;
double *kxyz, dx, dy;
double wopt, wcmp;
double fAinit, fBinit;
int aismatrix;
double timestep;
int timerepeat;

int main(int argc, char **argv)
{
    double *wm, *wp, *eta, *phA, *phB;
    int i, j, k, in, iseed = -3;
    double lylx;
    long ijk;

    FILE *fp;
    time_t ts;
    iseed = time(&ts);
    srand48(iseed);

    fp = fopen("para", "r");
    fscanf(fp, "%d", &in);
	fscanf(fp, "%lf, %lf", &wopt, &wcmp);
	fscanf(fp, "%lf", &lylx);
	fscanf(fp, "%d", &aismatrix);
	fscanf(fp, "%lf", &hAB);
	fscanf(fp, "%lf, %lf", &fA1, &fB);
	fscanf(fp, "%lf, %lf", &lx, &ly);
	fscanf(fp, "%d, %d", &Nx, &Ny);
	fscanf(fp, "%lf", &ds0);
	fscanf(fp, "%lf, %lf", &epA, &epB);
    fscanf(fp, "%lf, %d", &timestep, &timerepeat);
	fclose(fp);

    NxNy = Nx * Ny;
    Nyh1 = Ny / 2 + 1;
    NxNy1 = Nx * Nyh1;
    double kx[Nx], ky[Ny];

    wm = (double *)malloc(sizeof(double) * NxNy);
    wp = (double *)malloc(sizeof(double) * NxNy);
	phA = (double *)malloc(sizeof(double) * NxNy);
	phB = (double *)malloc(sizeof(double) * NxNy);
	eta = (double *)malloc(sizeof(double) * NxNy);
	kxyz = (double *)malloc(sizeof(double) * NxNy);

    if (lylx != 0)
        ly = lx * sqrt(lylx);
    
    dx = lx / Nx;
    dy = ly / Ny;
    
    fA2 = 1.0 - fB - fA1;

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
				{
					fscanf(fp, "%lf %lf %lf %lf", &e1, &e2, &e3, &e4);
					ijk = (long)(i * Ny + j);
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
	{
		for (j = 0; j < Ny; j++)
			{
				ijk = (long)(i * Ny + j);

				wm[ijk] = hAB * (fA1 + fA2 - fB) + 0.10 * (drand48() - 0.5);
				wp[ijk] = hAB + 0.10 * (drand48() - 0.5);
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
						phA[ijk], phB[ijk], wm[ijk] + wp[ijk], wp[ijk] - wm[ijk]);
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

    wpDiff = (double *)malloc(sizeof(double) * NxNy);
	wpnew = (double *)malloc(sizeof(double) * NxNy);
	
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

            for (ijk = 0; ijk < NxNy; ijk++)
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

            freeH /= NxNy;
            
            freeS = -log(qC);

            freeOld = freeEnergy;
            freeEnergy = freeH + freeS;

            for (ijk = 0; ijk < NxNy; ijk++)
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
    double dfdx, dfdy;
    double flow;
    int i, j, ijk;

    dfdwm = (double *)malloc(sizeof(double) * NxNy);

    for (ijk = 0; ijk < NxNy; ijk++)
    {
        dfdwm[ijk] = 2 / hAB * wm[ijk] + (phA[ijk] - phB[ijk]);
    }
    for (i = 0; i < Nx; i++)
    {
        for (j = 0; j < Ny; j++)
        {
            ijk = (long)(i * Ny + j);
            dfdx = (dfdwm[((i + 1) % Nx) * Ny + j] + dfdwm[((i - 1 + Nx) % Nx) * Ny + j] - 2 * dfdwm[i * Ny + j]) / dx / dx;
            dfdy = (dfdwm[i * Ny + ((j + 1) % Ny)] + dfdwm[i * Ny + ((j - 1 + Ny) % Ny)] - 2 * dfdwm[i * Ny + j]) / dy / dy;

            flow = dfdx + dfdy;

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

    wA = (double *)malloc(sizeof(double) * NxNy);
    wB = (double *)malloc(sizeof(double) * NxNy);
	qA1 = (double *)malloc(sizeof(double) * NxNy * (NsA1 + 1));
	qcA1 = (double *)malloc(sizeof(double) * NxNy * (NsA1 + 1));
	qB = (double *)malloc(sizeof(double) * NxNy * (NsB + 1));
	qcB = (double *)malloc(sizeof(double) * NxNy * (NsB + 1));
	qA2 = (double *)malloc(sizeof(double) * NxNy * (NsA2 + 1));
	qcA2 = (double *)malloc(sizeof(double) * NxNy * (NsA2 + 1));
	qInt = (double *)malloc(sizeof(double) * NxNy);

    for (ijk = 0; ijk < NxNy; ijk++){
        wA[ijk] = wp[ijk] + wm[ijk];
        wB[ijk] = wp[ijk] - wm[ijk];
    }

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = 1.0;
	}

	sovDifFft(qA1, wA, qInt, fA1, NsA1, 1, epA);
	sovDifFft(qcA2, wA, qInt, fA2, NsA2, -1, epA);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = qcA2[ijk * (NsA2 + 1)];
	}

	sovDifFft(qcB, wB, qInt, fB, NsB, -1, epB);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = qcB[ijk * (NsB + 1)];
	}

	sovDifFft(qcA1, wA, qInt, fA1, NsA1, -1, epA);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = qA1[ijk * (NsA1 + 1) + NsA1];
	}
	sovDifFft(qB, wB, qInt, fB, NsB, 1, epB);

	for (ijk = 0; ijk < NxNy; ijk++)
	{
		qInt[ijk] = qB[ijk * (NsB + 1) + NsB];
	}

	sovDifFft(qA2, wA, qInt, fA2, NsA2, 1, epA);

	ql = 0.0;
	for (ijk = 0; ijk < NxNy; ijk++)
	{
		ql += qcA1[ijk * (NsA1 + 1)];
	}

	ql /= NxNy;

	ffl = phs0 / ql * ds0;

	for (ijk = 0; ijk < NxNy; ijk++)
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
