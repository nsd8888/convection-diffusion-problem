#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#define grid 50


// This is the SIMPLE algorithm
/* temporal discritizaton
 Using 1st Order Explicit Scheme (Backward interpolation profile for the temporal)
 (phi(n+1) - phi(n))/delT + L(phi(n)) = 0;

 Note: As it is the explicit scheme check the courant number for the current mesh
 */

// Guess the Pressure 1st
// Initialize the Velocity Field in the domain
// Using the guess pressure and the initialized velocity calculate the guess (uncorrect) velocity
/*--> discrtization schemes for the velocity
  --> central difference scheme / second order upwind / upwing scheme / Quick scheme 
  --> as it is the structured mesh.......no need to worry about the 


  openfoam terms
gradient scheme

gauss linear........contins celllimited gauss linear 1;....1 indicates the fully bounded and o is not;
least square
gauss cubic

  CONVECTION SCHEMES/ DIVERGENCE SCHEMES
  --> linear interpolation........means central difference scheme.....this is the 2nd order accurate
  --> upwind interpolation scheme ........ 1st order accurate......too inaccurate
  --> linear upwing interpolation ........MEANS SECOND ORDER UPWIND 2nd order acurate
  --> LUST....blended 75% linear and 25% lineAR UPWIND......REQUIRES VELOCITY GRADIENT
  -->limited linear.....that limit towards upwind scheme in region of rapidly changig gradient; requires a coefficient, 1 is strongest limiting, tending towards linear as coefficient tends to zero
  --> V- schemes....

*/

// initialize pressure in the field






using namespace std;
int main()
{
    double u[grid][grid+1], un[grid][grid+1], ustr[grid][grid+1], u_E[grid][grid+1],u_P[grid][grid+1], u_N[grid][grid+1], u_S[grid][grid+1], a_e[grid][grid+1], uc[grid][grid];
    double a_E[grid][grid+1],a_P[grid][grid+1], a_N[grid][grid+1], a_S[grid][grid+1];

	double v[grid+1][grid], vn[grid+1][grid], vstr[grid+1][grid],v_E[grid+1][grid],v_W[grid+1][grid],v_N[grid+1][grid],v_S[grid+1][grid], b_e[grid+1][grid], vc[grid][grid];
	double b_E[grid+1][grid],b_W[grid+1][grid],b_N[grid+1][grid],b_S[grid+1][grid];

	double p[grid+1][grid+1], pn[grid+1][grid+1], pstr[grid+1][grid+1], p_e[grid+1][grid+1],pcor[grid+1][grid+1],c_p[grid+1][grid+1],pc[grid][grid];
	//double m[grid+1][grid+1];
	int i, j, step;
	double dx, dy, alpha, error, Re;
	double length= 1.0;
	step =1;
	dx = length/(grid-1);
	dy = length/(grid-1);
	alpha= 0.8;
	error = 1.0;
	Re = 100.0;


for (i = 1; i <= (grid - 1); i++)
{
	for (j = 1; j <= (grid - 1); j++)
	{
		p[i][j]= 1 + ((9810*((grid-1)-((j-1)+(1/2))*dy)* 0.00001));
	}
}

for (i = 1; i <= (grid - 1); i++)
{
	p[i][grid] = p[i][grid - 1];
}

for (j = 0; j <= (grid); j++)
{
	p[0][j] = p[1][j];
	p[grid][j] = p[grid-1][j];
}

//initialize the velocity field
for (i = 0; i <= (grid - 1); i++)
{
	for (j = 0; j <= (grid); j++)
	{
		u[i][j] = 0.0;
		u[i][grid] = 1.0;
		u[i][grid - 1] = 1.0;
	}
}

for (i = 0; i <= (grid); i++)
{
	for (j = 0; j <= (grid - 1); j++)
	{
		v[i][j] = 0.0;
	}
}



while (error<=0.000001)
{

// guess the presure

for (i = 0; i <= (grid); i++)
{
	for (j = 0; j <= (grid); j++)
	{
		pstr[i][j] = p[i][j];
	}
}

for(i = 0; i <= (grid - 1); i++)
{
	for (j = 0; j <= (grid); j++)
	{
		ustr[i][j] = u[i][j];
	}
}

for (i = 0; i <= (grid); i++)
{
	for (j = 0; j <= (grid - 1); j++)
	{
		vstr[i][j] = v[i][j];
	}
}



// calcu;ate uncorrect velocity
//using expicit method.....that one is backward interpolation for the time step

for (i = 1; i < (grid - 2); i++)
	{
		for (j = 1; j <= (grid - 1); j++)
		{
			u_E[i][j] = (ustr[i + 1][j] + ustr[i][j]) / 2;
			u_P[i][j] = (ustr[i][j] + ustr[i - 1][j]) / 2;
			u_N[i][j] = (vstr[i][j] + vstr[i + 1][j]) / 2;
			u_S[i][j] = (vstr[i][j - 1] + vstr[i + 1][j - 1]) / 2;

			a_E[i][j] = ((dy / dx * Re) - (u_E[i][j] * dy / 2));
			a_P[i][j] = ((dy / Re * dx) + (u_P[i][j] * dy / 2));
			a_N[i][j] = ((dx / dy * Re) - (u_N[i][j] * dx / 2));
			a_S[i][j] = ((dx / Re * dy) + (u_S[i][j] * dx / 2));

			a_e[i][j] = ((0.5 * u_E[i][j]) - (0.5 * u_P[i][j])) * dy + ((0.5 * u_N[i][j]) - (0.5 * u_S[i][j])) * dx + ((dy / dx)*(2 / Re)) + ((dx / dy)*(2 / Re));

			ustr[i][j] = ((a_E[i][j] / a_e[i][j]) * ustr[i + 1][j]) + ((a_P[i][j] / a_e[i][j]) * ustr[i - 1][j]) + ((a_N[i][j] / a_e[i][j]) * ustr[i][j - 1]) +
						 ((dy / a_e[i][j]) * (pstr[i + 1][j] - pstr[i][j]));
		}
	}

for (j = 1; j <= (grid - 1); j++)
{
	ustr[0][j] = 0;
	ustr[grid][j] = 0;
}

for (i = 0; i <= (grid - 1); i++)
{
	ustr[i][0] = -ustr[i][1];
	ustr[i][grid] = 2 - ustr[i][grid - 1];
}

// guess v velocity
// v velocity coreection
for (i = 1; i <= (grid - 1); i++)
{
	for (j = 1; j <= (grid - 2); j++)
	{
		v_E[i][j] = (ustr[i][j + 1] + ustr[i][j]) / 2;
		v_W[i][j] = (ustr[i - 1][j] + ustr[i - 1][j + 1]) / 2;
		v_N[i][j] = (vstr[i][j + 1] + vstr[i][j]) / 2;
		v_S[i][j] = (vstr[i][j] + vstr[i][j - 1]) / 2;

		b_W[i][j] = ((v_W[i][j] * (dy / 2)) + (dy / (dx * Re)));
		b_E[i][j] = ((dy / (dx * Re)) - (v_E[i][j] * (dy / 2)));
		b_N[i][j] = ((dx / (Re * dy)) - (v_N[i][j] * (dx / 2)));
	    b_S[i][j]= ((dx/(Re*dy))+ (v_S[i][j])*(dx/2));

	 b_e[i][j] = ((v_E[i][j] - v_W[i][j]) * (dy / 2)) + ((v_N[i][j] - v_S[i][j]) * (dx / 2)) + ((2 / Re) * (dy / dx)) + ((2 / Re) * (dx / dy));

	 vstr[i][j]= (b_W[i][j]/b_e[i][j])* vstr[i-1][j] + (b_E[i][j]/b_e[i][j])* vstr[i+1][j] + (b_N[i][j]/b_e[i][j])* vstr[i][j+1]
	               + (b_S[i][j]/b_e[i][j])* vstr[i][j-1] + (pstr[i][j+1]-pstr[i][j])*(dx/b_e[i][j]);
	}
}

for (j = 1; j <= (grid - 2); j++)
{
	vstr[0][j] = -vstr[1][j];
	vstr[grid][j] = -vstr[grid - 1][j];
}

for (i = 0; i <= (grid); i++)
{
	vstr[i][0] = 0.0;
	vstr[i][grid - 1] = 0.0;
}

//pressure correction
//first initialie them as the zero
for (i = 0; i <= (grid); i++)
{
	for (j = 0; j <= (grid); j++)
	{
		pcor[i][j] = 0;
	}
}

//calculate the pressure corection at interior points
for (i = 1; i <= (grid - 1); i++)
{
	for (j = 1; j <= (grid - 1); j++)
	{
		p_e[i][j] = (-(dy * dy / a_e[i - 1][j]) - (dy * dy / a_e[i][j]) - (dx * dx / b_e[i][j]) - (dx * dx / b_e[i][j - 1]));
	}
}

for (i = 1; i <= (grid - 1); i++)
{
	for (j = 1; j <= (grid - 1); j++)
	{

		c_p[i][j] = (dy * (ustr[i - 1][j] - ustr[i][j])) + (dx * (vstr[i][j - 1] - vstr[i][j]));
		pcor[i][j] = (-dy * dy / (a_e[i][j] * p_e[i][j])) * pcor[i + 1][j] + (-dy * dy / (a_e[i - 1][j] * p_e[i][j]) * pcor[i - 1][j]) + (-dx * dx / (b_e[i][j] * p_e[i][j])) * pcor[i][j + 1] + (-dx * dx / (b_e[i][j - 1] * p_e[i][j])) * pcor[i][j - 1] + c_p[i][j];
	}
}

//interior pressure boundaru condition

for (i = 1; i <= (grid - 1); i++)
{
	for (j = 1; j <= (grid - 1); j++)
	{
		pn[i][j] = pstr[i][j] + (alpha * pcor[i][j]);
	}
}

//boundary p
for (i = 1; i <= (grid - 1); i++)
{
	pn[i][0] = pn[i][1];
	pn[i][grid] = pn[i][grid - 1];
}

for (j = 0; j <= (grid); j++)
{
	pn[0][j] = pn[1][j];
	pn[grid][j] = pn[grid - 1][j];
}

//update u and v values
for (i = 1; i <= (grid - 2); i++)
{
	for (j = 1; j <= (grid - 1); j++)
	{
		un[i][j] = ustr[i][j] + (alpha * dy * (pcor[i + 1][j] - pcor[i][j]));
	}
}

for (i = 1; i <= (grid - 1); i++)
{
	for (j = 1; j <= (grid - 2); j++)
	{
		vn[i][j] = vstr[i][j] + (alpha * dx * (pcor[i][j + 1] - pcor[i][j]));
	}
}
//update the boundary conditions

for (j = 1; j <= (grid - 1); j++)
{
	un[0][j] = 0;
	un[grid][j] = 0;
}

for (i = 0; i <= (grid - 1); i++)
{
	un[i][0] = -un[i][1];
	un[i][grid] = 2 - un[i][grid - 1];
}

for (j = 1; j <= (grid - 2); j++)
{
	vn[0][j] = -vn[1][j];
	vn[grid][j] = -vn[grid - 1][j];
}

for (i = 0; i <= (grid); i++)
{
	vn[i][0] = 0.0;
	vn[i][grid - 1] = 0.0;
}

//error
for(i = 1; i <= (grid - 1); i++)
{
	for (j = 1; j <= (grid - 1); j++)
	{
		error = 0;
		error +=  fabs(c_p[i][j]);
	}
}

cout<<"the error at the step "<<step<<" is "<<error;

// Iterating u
for(i = 0; i <= (grid - 1); i++)
{
	for (j = 0; j <= (grid); j++)
	{
		u[i][j] = un[i][j];
	}
}

// Iterating v
for (i = 0; i <= (grid); i++)
{
	for (j = 0; j <= (grid - 1); j++)
	{
		v[i][j] = vn[i][j];
	}
}

// Iterating p
for (i = 0; i <= (grid); i++)
{
	for (j = 0; j <= (grid); j++)
	{
		p[i][j] = pn[i][j];
	}
}

step= step+1;

}

for (i=0; i<=(grid-1); i++)
	{
		for (j=0; j<=(grid-1); j++)
		{	
			uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
			vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
			pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
		}
	}

}






