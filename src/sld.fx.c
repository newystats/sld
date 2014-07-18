/* sld.fx.c - Part of Robert King's sld package for the R statistical
 * language, implementing the skew logistic distribution
 *
 * Copyright (C) Robert King 2014
 * robert.king@newcastle.edu.au
 * or robert.king.newcastle@gmail.com
 * http://tolstoy.newcastle.edu.au/~rking/publ/software.html
 *
 * This program is free software; you can redistribute it and/or modify it 
 * under the terms of the GNU General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) 
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA  02111-1307  USA
 *
 * These functions calculate F(x) for the quantile-based Skew Logistic distribution
 * It is adapted from gld.rs.fx.c, (c) Robert King 1993,2000,2001,2006 in the gld package.
 * Because of that, this is not very good C.  The original code was my first 
 * real program in C and I haven't tidied it up very much for release with the R package.
 * However, it has been used for some time and works fine.
 * 
 * This is code is loaded into R.
 * The two include lines are handled by R, but need to be here to keep the compile step happy
 */

#include <stdio.h>
#include <math.h>
#include <R.h> 

/* proto */

void sld_funcd( double , double , double *, double *, double *, double *, double *);

/* the function that finds the root */

void sld_distfunc( double *pa,double *pb,double *pd, 
double *pu1,double *pu2,double *pxacc, int *max_it,
double *ecks, double *u, int *lengthofdata)
{

/* pa, pb and pd:    pointers to the values of the parameters of the sld
 * pu1:         minimum value of u, should be zero
 * pu2:         maximum value of u, should be 1
 * pxacc:       desired accuracy of the calculation
 * max_it:      maximum iterations for N-R root finder
 * ecks:        the quantiles of the gld given
 * u:           array to put the calculated depths
 * pl:          length of the data
 */

double  u1, u2, xacc; 		
int i,j; /* indices */
double df,dx,dxold,f,fh,fl;
double x;
double temp,xh,xl,rts;
/* Initialising things */
i=0; j=0;
df=0.0;dx=0.0;dxold=0.0;f=0.0;fh=0.0;fl=0.0;
x=0.0;
temp=0.0;xh=0.0;xl=0.0;rts=0.0;

u1 = *pu1; u2 = *pu2; xacc = *pxacc;

for (i=0;i<*lengthofdata;i++)
{
    x = ecks[i];
	u[i] = 0.0;
	sld_funcd(u1,x,&fl,&df,pa,pb,pd);
	sld_funcd(u2,x,&fh,&df,pa,pb,pd);
	if (fl*fh >= 0.0) 
	{
		/* This is suggested in writing R extensions, but still gives the warning (not anymore) */
		error("Program aborted at parameter values %f, %f, %f\n The data value being investigated was index %d, value: %f\n", *pa, *pb, *pd, i, x);
	}
	if (fl < 0.0) {
		xl = u1;
		xh = u2;
		}
	else {
		xh = u1;
		xl = u2;
		}
	rts = 0.5*(u1+u2);
	dxold = fabs(u2-u1);
	dx = dxold;
	sld_funcd(rts,x,&f,&df,pa,pb,pd);
	for (j=1;j<=*max_it;j++) {
		if ((((rts - xh)*df - f)* ( (rts-xl)*df - f) >= 0.0 ) ||
( fabs(2.0*f) > fabs (dxold*df))) {
			dxold = dx;
			dx = .5* (xh - xl);
			rts = xl +dx;
			if (xl == rts ) { 
				u[i] = rts; 
				break; }
			}
		else {
			dxold = dx;
			dx = f/df;
			temp = rts;
			rts -= dx;
			if (temp == rts) { 
				u[i] = rts; 
				break; }
			}
		if (fabs(dx) < xacc) { 
			u[i] = rts; 
			break; }
		sld_funcd(rts,x,&f,&df,pa,pb,pd);
		if (f < 0.0)
			xl =rts;
		else 
			xh =rts;
		}
}
}

void sld_funcd( double u, double x, double *F, double *dFdu, double *pa, double *pb, double *pd)
{

/* *F is the sld Q(u)  */            
/* *dFdu is dQ(u)/du	*/
/* NOTE: *dFdu is 1/f(x) a.k.a. 1/f(Q(u)) - the opposite of what's required for the pdf */

if ( *pd == 0.0 ) {
  /* reflected exponential */
  *F = *pa + *pb * log(u) - x;
  }
if (*pd == 1.0) {
  /* exponential */
  *F = *pa + *pb *(1.0 - log(u)) -x;
  }
if ((*pd > 0.0) && (*pd < 1.0)) {
  *F = *pa + *pb * ( (1.0 - *pd)*log(u) - *pd * log(1.0 - u)) - x;
  }  
*dFdu = (*pb * (*pd * u + (1.0 - *pd)*(1.0 - u))) / (u * (1.0 - u));
}