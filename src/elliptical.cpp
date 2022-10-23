#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 1e-7
#define PI 3.1415926535
#define N1 100       /* ÓÎ-‚Ó ˇ˜ÂÂÍ ÔÓ ÓÒË r*/
#define N2 500      /* ÓÎ-‚Ó ˇ˜ÂÂÍ ÔÓ ÓÒË phi*/
#define l1 1.0        /*ƒÎËÌ‡ ÒÚÓÓÌ˚ ÔˇÏÓÛ„ÓÎ¸ÌËÍ‡ ÔÓ ÓÒË x*/
#define l2 2.0 * PI

double function(double *, double *, int, int);
double function(double *r, double *phi, int k, int j)
{
	return -(1.0 / 12.0)*pow(r[k], 2.0)*sin(phi[j])*cos(phi[j])*(pow(r[k], 2.0) - 1);
}


/* ¬€◊»—À≈Õ»≈ ÕŒ–Ã€  ||x||_inf = max|x_i| */

double norm_inf(double *, int);
double norm_inf(double *a, int n)
{
	double maximum = 0;
	int i = 0;

	maximum = fabs(a[0]);
	for (i = 1; i < n; i++)
	{
		if (fabs(a[i]) > maximum)
		{
			maximum = fabs(a[i]);
		}
	}

	return maximum;
}

int main(void)
{
	int i = 0,
		j = 0,
		k = 1,
		param_iter = 0;

	double h1 = 0,
		h2 = 0,
		temp1 = 0,
		temp2 = 0,
		temp3 = 0,
		temp4 = 0;

	double* u, * u_next, * r, * phi, * norm, * u_r, *u_phi;

	u = (double*)calloc((N1 - 1)*(N2 + 1), sizeof(double));
	u_next = (double*)calloc((N1 - 1)*(N2 + 1), sizeof(double));
	norm = (double*)calloc((N1 - 1)*(N2 + 1), sizeof(double));
	r = (double*)calloc(N1 + 1, sizeof(double));
	phi = (double*)calloc(N2 + 1, sizeof(double));
	u_r = (double*)calloc(N1 + 1, sizeof(double));
	u_phi = (double*)calloc(N2 + 1, sizeof(double));


	if ((u == NULL) || (u_next == NULL) || (r == NULL) || (phi == NULL) || (norm == NULL) || (u_r == NULL) || (u_phi == NULL))
	{
		printf("\nMemory allocation error!\n");
		return -1;
	}

	/********* »Õ»÷»¿À»«¿÷»ﬂ Ã¿——»¬Œ¬ r Ë phi, ¿ “¿ ∆≈ «¿ƒ¿Õ»≈ Õ¿◊¿À‹ÕŒ√Œ «Õ¿◊≈Õ»ﬂ ¬≈ “Œ–¿ u *********/

	h1 = l1 / N1;
	h2 = l2 / N2;

	u_r[0] = u_r[N1] = 0;

	r[0] = 0.0;
	phi[0] = 0.0;

	k = 1;
	j = 0;

	for (i = 0; i < (N1 - 1)*(N2 + 1); i++)
	{
		u[i] = 0.0;
		
		if (i < N1)
		{
			r[i + 1] = r[i] + h1;
		}
		if (i < N2)
		{
			phi[i + 1] = phi[i] + h2;
		}
	}

	/*********  ŒÕ≈÷ »Õ»÷»¿À»«¿÷»» *********/


	/********* Õ¿◊¿ÀŒ »“≈–¿÷»ŒÕÕŒ√Œ œ–Œ÷≈——¿ *********/
	
	do
	{
		j = 0;
		k = 1;

		for (i = 0; i < (N1 - 1)*(N2 + 1); i++)
		{
			if (k > N1 - 1)
			{
				j++;
				k = 1;
			}

			if (j == 0 || j == N2)
			{
				if (j == 0)
				{
					if (k == 1 || k == N1 - 1)
					{
						if (k == 1)
						{
							temp1 = 2.0 * h2 * h2 * r[k] * r[k] * (u[i + 1]);
							temp2 = h1 * h2 * h2 * r[k] * (u[i + 1]);
							temp3 = 2.0 * h1 * h1 * (u[i + N1 - 1] + u[i + (N2 - 1)*(N1 - 1)]);
							temp4 = 2.0 * h1 * h1 * h2 * h2 * pow(r[k], 4.0)*sin(phi[j])*cos(phi[j]);
							u_next[i] = (temp1 + temp2 + temp3 + temp4) / (4.0*h2*h2*pow(r[k], 2.0) + 4.0*h1*h1);
						}
						if (k == N1 - 1)
						{
							temp1 = 2.0 * h2 * h2 * r[k] * r[k] * (u[i - 1]);
							temp2 = h1 * h2 * h2 * r[k] * (-u[i - 1]);
							temp3 = 2.0 * h1 * h1 * (u[i + N1 - 1] + u[i + (N2-1) * (N1 - 1)]);
							temp4 = 2.0 * h1 * h1 * h2 * h2 * pow(r[k], 4.0)*sin(phi[j])*cos(phi[j]);
							u_next[i] = (temp1 + temp2 + temp3 + temp4) / (4.0*h2*h2*pow(r[k], 2.0) + 4.0*h1*h1);
						}
					}
					else
					{
						temp1 = 2.0 * h2 * h2 * r[k] * r[k] * (u[i + 1] + u[i - 1]);
						temp2 = h1 * h2 * h2 * r[k] * (u[i + 1] - u[i - 1]);
						temp3 = 2.0 * h1 * h1 * (u[i + N1 - 1] + u[i + (N2-1) * (N1 - 1)]);
						temp4 = 2.0 * h1 * h1 * h2 * h2 * pow(r[k], 4.0)*sin(phi[j])*cos(phi[j]);
						u_next[i] = (temp1 + temp2 + temp3 + temp4) / (4.0*h2*h2*pow(r[k], 2.0) + 4.0*h1*h1);
						/*if (fabs(r[k] - 0.5) < 0.000001)
						{
							u_phi[j] = u_next[i];
						}*/
					}
				}
				if (j == N2)
				{
					u_next[i] = u_next[i - (N2) * (N1 - 1)];
					/*if (fabs(r[k] - 0.5) < 0.000001)
					{
						u_phi[j] = u_next[i];
					}*/
					/*if (k == 1 || k == N1 - 1)
					{
						if (k == 1)
						{
							temp1 = 2.0 * h2 * h2 * r[k] * r[k] * (u[i + 1]);
							temp2 = h1 * h2 * h2 * r[k] * (u[i + 1]);
							temp3 = 2.0 * h1 * h1 * (u[i - N1 + 1] + u[i - (N2-1) * (N1 - 1)]);
							temp4 = 2.0 * h1 * h1 * h2 * h2 * pow(r[k], 4.0)*sin(phi[j])*cos(phi[j]);
							u_next[i] = (temp1 + temp2 + temp3 + temp4) / (4.0*h2*h2*pow(r[k], 2.0) + 4.0*h1*h1);
						}
						if (k == N1 - 1)
						{
							temp1 = 2.0 * h2 * h2 * r[k] * r[k] * (u[i - 1]);
							temp2 = h1 * h2 * h2 * r[k] * (-u[i - 1]);
							temp3 = 2.0 * h1 * h1 * (u[i - N1 + 1] + u[i - (N2-1) * (N1 - 1)]);
							temp4 = 2.0 * h1 * h1 * h2 * h2 * pow(r[k], 4.0)*sin(phi[j])*cos(phi[j]);
							u_next[i] = (temp1 + temp2 + temp3 + temp4) / (4.0*h2*h2*pow(r[k], 2.0) + 4.0*h1*h1);
						}
					}
					else
					{
						temp1 = 2.0 * h2 * h2 * r[k] * r[k] * (u[i - 1] + u[i+1]);
						temp2 = h1 * h2 * h2 * r[k] * (u[i + 1] - u[i-1]);
						temp3 = 2 * h1 * h1 * (u[i - N1 + 1] + u[i - (N2-1) * (N1 - 1)]);
						temp4 = 2 * h1 * h1 * h2 * h2 * pow(r[k], 4.0)*sin(phi[j])*cos(phi[j]);
						u_next[i] = (temp1 + temp2 + temp3 + temp4) / (4.0*h2*h2*pow(r[k], 2.0) + 4.0*h1*h1);
					}*/
				}
			}
			else
			{
				if (k == 1 || k == N1 - 1)
				{
					if (k == 1)
					{
						temp1 = 2.0 * h2 * h2 * r[k] * r[k] * (0 + u[i + 1]);
						temp2 = h1 * h2 * h2 * r[k] * (0 + u[i + 1]);
						temp3 = 2 * h1 * h1 * (u[i - N1 + 1] + u[i + N1 - 1]);
						temp4 = 2 * h1 * h1 * h2 * h2 * pow(r[k], 4.0)*sin(phi[j])*cos(phi[j]);
						u_next[i] = (temp1 + temp2 + temp3 + temp4) / (4.0*h2*h2*pow(r[k], 2.0) + 4.0*h1*h1);
						/*if (fabs(phi[j] - (7.0 * PI / 4.0)) < 0.01)
						{
							u_r[k] = u_next[i];
						}*/
					}
					if (k == N1 - 1)
					{
						temp1 = 2.0 * h2 * h2 * r[k] * r[k] * (u[i - 1]);
						temp2 = h1 * h2 * h2 * r[k] * (-u[i - 1] );
						temp3 = 2 * h1 * h1 * (u[i - N1 + 1] + u[i + N1 -1]);
						temp4 = 2 * h1 * h1 * h2 * h2 * pow(r[k], 4.0)*sin(phi[j])*cos(phi[j]);
						u_next[i] = (temp1 + temp2 + temp3 + temp4) / (4.0*h2*h2*pow(r[k], 2.0) + 4.0*h1*h1);	
						/*if (fabs(phi[j] - (7.0 * PI / 4.0)) < 0.01)
						{
							u_r[k] = u_next[i];
							//printf("OK");
						}*/
					}
				}
				else
				{
					temp1 = 2.0 * h2 * h2 * r[k] * r[k] * (u[i + 1] + u[i - 1]);
					temp2 = h1 * h2*h2*r[k] * (u[i + 1] - u[i - 1]);
					temp3 = 2 * h1*h1*(u[i + N1 - 1] + u[i - N1 + 1]);
					temp4 = 2 * h1*h1*h2*h2*pow(r[k], 4.0)*sin(phi[j])*cos(phi[j]);
					u_next[i] = (temp1 + temp2 + temp3 + temp4) / (4.0*h2*h2*pow(r[k], 2.0) + 4.0*h1*h1);
					/*if (fabs(phi[j] - (7.0 * PI / 4.0)) < 0.01)
					{
						u_r[k] = u_next[i];
					}*/
					/*if (fabs(r[k] - 0.5) < 0.000001)
					{
						u_phi[j] = u_next[i];
					}*/
				}
			}

			norm[i] = u_next[i] - u[i];
			//u[i] = u_next[i];
			k++;
		}

		for (i = 0; i < (N1 - 1)*(N2 + 1); i++)
		{
			u[i] = u_next[i];
		}

		param_iter++;

		//printf("%lf ", norm_inf(norm, (N1 - 1)*(N2 + 1)));
	} while (fabs(norm_inf(norm, (N1 - 1) * (N2 + 1)))>EPS);

	/*********  ŒÕ≈÷ »“≈–¿÷»ŒÕÕŒ√Œ œ–Œ÷≈——¿ *********/
	//printf("%d", param_iter);
	FILE* g;
	//fopen_s(&g, "output_iter8000_cells_u(r)7pi4_200na1000.txt", "w");
	fopen_s(&g, "exact.txt", "w");
	k = 1;
	j = 0;

	/*for (i = 0; i < N1 + 1; i++)
	{
		fprintf(g, "%lf %lf\n",r[i], u_r[i]);
	}*/
	for (int i = 0; i < (N1 - 1) * (N2 + 1); i++)
	{
		norm[i] = function(r, phi, k, j);
		if (k > N1 - 1)
		{
			k = 1;
			j++;
			fprintf(g, " \n");
		}

		fprintf(g, "%lf ", norm[i]);
		k++;
	}
	fclose(g);


	system("pause");

	free(u);
	free(u_next);
	free(r);
	free(phi);
	free(norm);
	free(u_r);

	return 0;
}