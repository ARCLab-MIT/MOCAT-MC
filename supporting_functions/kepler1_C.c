#include <math.h> /*Include library */
#include "mex.h" /*Always include this */

void kepler1_C(double *manom,double *ecc,double pi,int n_sats,double *eanom){

int i, niter;
double ktol, pi2, xmai, eanomi, ecci, toli, si, ci, f, fp, fpp, fppp, delta, deltastar, deltak, sta, cta;
ktol = 1.0e-10;
pi2 = 2.*pi;
for(i=0; i<n_sats; i++){
    xmai = manom[i] - pi2*trunc(manom[i]/pi2);
    ecci = ecc[i];
    if (ecci == 0.0) {
        // circular orbit
//         tanom[i] = xmai;
        eanom[i] = xmai;
        continue;
    } else if (ecci<1.0) { /* more than one coefficient */
        // elliptic orbit
        eanomi = xmai + 0.85 * ecci*copysign(1,sin(xmai));
    } else {
        // hyperbolic orbit
        eanomi = log(2.0 * xmai / ecci + 1.8);
    }
    toli = 2.*ktol;
    niter = 0;
    while (toli>ktol){
        if (ecci<1.) {
            // elliptic orbit
            si = ecci * sin(eanomi);
            ci = ecci * cos(eanomi);
            
            f = eanomi - si - xmai;
            fp = 1. - ci;
            fpp = si;
            fppp = ci;

        } else {
            // hyperbolic orbit
            si = ecci * sinh(eanomi);
            ci = ecci * cosh(eanomi);

            f = si - eanomi - xmai;
            fp = ci - 1.;
            fpp = si;
            fppp = ci;
        }
        niter++;
        if (niter==20) {
            printf("Number of iterations exceeded limit of 20 in kepler1_C.");
            return;
        }
        
        toli = fabs(f);
        
        delta = -f / fp;
        deltastar = -f / (fp + 0.5 * delta * fpp);
        deltak = -f / (fp + 0.5 * deltastar * fpp + deltastar * deltastar * fppp / 6.);
        eanomi += deltak;
    }
    eanom[i] = eanomi;
//     if (ecci<1.) {
//         // elliptic orbit
//         sta = sqrt(1. - ecci * ecci) * sin(eanomi);
//         cta = cos(eanomi) - ecci;
//     }
//     else {
//         // hyperbolic orbit
//         sta = sqrt(ecci * ecci - 1.) * sinh(eanomi);
//         cta = ecci - cosh(eanomi);
//     }
//     tanom[i] = atan2(sta, cta);
}    
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
/* Macros for the ouput and input arguments */
int n_sats;
double *manom, *ecc, pi, *eanom;

/* Inputs:  kepler1_C(double *manom,double *ecc,double pi, int n_sats,double *eanom) */

/* create a pointer to the real data in the input vectors/matrices  */
    manom = mxGetPr(prhs[0]);
    ecc = mxGetPr(prhs[1]);

    /* get the value of the scalar inputs  */
    pi = mxGetScalar(prhs[2]);
    n_sats = mxGetScalar(prhs[3]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(n_sats,1,mxREAL);
//     plhs[1] = mxCreateDoubleMatrix(n_sats,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    eanom = mxGetPr(plhs[0]);
//     tanom = mxGetPr(plhs[1]);

    /* call the computational routine */
    kepler1_C(manom,ecc,pi,n_sats,eanom);
}