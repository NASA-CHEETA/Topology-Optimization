#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "CalculiX.h"

// Apply heaviside filter
void apply_heavisidefilter(double eta, double beta , double* x, double* xh, double ne )
{
    int i;

    for(i=0;i<ne; i++){
        xh[i] = (tanh(beta*eta) + tanh(beta*(x[i] - eta))) / (tanh(beta*eta) + tanh(beta*(1.0 - eta)));
    }

    return;
}

//sum x-xh
double sum_projection(double eta, double beta, double* x, double* xh, double ne){
    double sumx = 0, sumxh = 0;

    apply_heavisidefilter(eta, beta, x, xh, ne);

    for(int i = 0; i<ne;i++){
        sumx+= x[i];
        sumxh += xh[i];
    }

    return sumx - sumxh;
}

double eta_bisection(double a, double b, double tol, double* x, double* xh, double beta, double ne){
    
    double m;
    m= (a+b)/2.0;

    double fm = 0, fa = 0, fb = 0;


    fm = sum_projection(m, beta, x, xh, ne);
    fa = sum_projection(a, beta, x, xh, ne);
    fb = sum_projection(b, beta, x, xh, ne);

    if (fabs(fm) < tol){
        return m;
    } else if( fa*fm > 0){
        // make recursive call with a = m
        return eta_bisection(m, b, tol, x, xh, beta, ne);
    } else if( fb*fm > 0){
        // make recursive call with b = m
        return eta_bisection(a, m, tol, x, xh, beta, ne);
    } else{
        printf("None condition satisfied in eta_bisection \n");
        return -1;
    }

}


void calculate_eta(double* designFiltered, ITG ne, double* beta_H, double* eta_H){

    double a, b, beta, eta;
    a = 0;
    b = 1;
    beta = *beta_H;

    double tol = 0.0001;

    double* rhoH = NULL;
    rhoH = (double*)malloc(ne * sizeof(double));

    printf("Calculating volume preserving Eta threshold \n");
    eta = eta_bisection(a, b, tol, designFiltered, rhoH, beta, ne);
    *eta_H = eta;

    free(rhoH); 

}