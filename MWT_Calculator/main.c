/*
 *  main.c
 *  MWT_Calculator
 *
 *  Created by Claudio Pica on 06/09/2010.
 *  Modified by Eugenio Del Nobile on 01/11/2010.
 *  Copyright 2010 All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "read_input.h"

/* Mathematica stuff */
#define Power(x, y)     (cpow((double complex)(x), (double complex)(y)))
#define Sqrt(x)         (csqrt((double complex)(x)))

#define Pi 3.1415926535897932384626433832795


/* EXTERNAL PARAMETERS */

/* CKMBLOCK Parameters */
double cabi= 0.227736; /* Cabibbo angle */
/* SMInput Parameters */
double aS= 0.118; /* alpha strong */
double EE= 0.313429; /*  */
double GF= 0.0000116637; /* Fermi coupling constant */
double MZ= 91.1876; /* Z0 mass */
/* TCInput Parameters */
double gt= 2; /* g^tilde */
double MA= 500; /* mass of the axial */
double PS= 0.3; /* S parameter */
double rs= 0.; /* C-M parameter */
double MH= 200.; /* Higgs mass */
/* YUKAWA Parameters */
double MC= 1.3; /* Charm Yukawa mass */
double MB= 4.2; /* Bottom Yukawa mass */
double MT= 172; /* Top Yukawa mass */
double MTA= 1.777; /* Tau Yukawa mass */

/* MASSES */
double HM=INFINITY;
double ZM=INFINITY;
double TAM=INFINITY;
double CM=INFINITY;
double TM=INFINITY;
double BM=INFINITY;

/* WIDTHS */
double wZ=2.4952;
double wW=2.141;
double wT=1.50833649;

double wH=INFINITY;
double w1N=INFINITY;
double w2N=INFINITY;
double w1C=INFINITY;
double w2C=INFINITY;


/* INTERNAL PARAMETERS */

double gs=INFINITY; /* Strong coupling constant */
double v= INFINITY; /* SM VEV */
double r3= INFINITY; /* C-M parameter */
double r2= INFINITY; /* C-M parameter */
double f= INFINITY; /* Vector meson mass scale */
double MV= INFINITY; /* Strong vector meson mass */
double FV= INFINITY; /* Strong vector decay constant */
double FPi= INFINITY; /* Strong pion decay constant */
double FA= INFINITY; /* Strong axial decay constant */
double g2= INFINITY; /* electroweak SU2L gauge coupling */
double g1= INFINITY; /* electroweak U1Y gauge coupling */
double M1N= INFINITY; /* Neutral heavy vector meson R1 mass */
double M2N= INFINITY; /* Neutral heavy vector meson R2 mass */
double ThetaC= INFINITY; /* Charged vector meson mass angle */
double MW= INFINITY; /* W mass */
double M1C= INFINITY; /* Charged heavy vector meson R1 mass */
double M2C= INFINITY; /* Charged heavy vector meson R2 mass */
double VC11= INFINITY; /* VC11 */
double VC12= INFINITY; /* VC12 */
double VC13= INFINITY; /* VC13 */
double VC21= INFINITY; /* VC21 */
double VC22= INFINITY; /* VC22 */
double VC23= INFINITY; /* VC23 */
double VC31= INFINITY; /* VC31 */
double VC32= INFINITY; /* VC32 */
double VC33= INFINITY; /* VC33 */
double VN11= INFINITY; /* VN11 */
double VN12= INFINITY; /* VN12 */
double VN13= INFINITY; /* VN13 */
double VN14= INFINITY; /* VN14 */
double VN21= INFINITY; /* VN21 */
double VN22= INFINITY; /* VN22 */
double VN23= INFINITY; /* VN23 */
double VN24= INFINITY; /* VN24 */
double VN31= INFINITY; /* VN31 */
double VN32= INFINITY; /* VN32 */
double VN33= INFINITY; /* VN33 */
double VN34= INFINITY; /* VN34 */
double VN41= INFINITY; /* VN41 */
double VN42= INFINITY; /* VN42 */
double VN43= INFINITY; /* VN43 */
double VN44= INFINITY; /* VN44 */
double CN1= INFINITY; /* Charged vector meson normalization factor */
double CN2= INFINITY; /* Charged vector meson normalization factor */
double CN3= INFINITY; /* Charged vector meson normalization factor */
double C11= INFINITY; /* Charged vector meson mixing matrix element in VA base */
double C12= INFINITY; /* Charged vector meson mixing matrix element in VA base */
double C13= INFINITY; /* Charged vector meson mixing matrix element in VA base */
double C21= INFINITY; /* Charged vector meson mixing matrix element in VA base */
double C22= INFINITY; /* Charged vector meson mixing matrix element in VA base */
double C23= INFINITY; /* Charged vector meson mixing matrix element in VA base */
double C31= INFINITY; /* Charged vector meson mixing matrix element in VA base */
double C32= INFINITY; /* Charged vector meson mixing matrix element in VA base */
double C33= INFINITY; /* Charged vector meson mixing matrix element in VA base */
double NN2= INFINITY; /* Neutral vector meson normalization factor */
double NN3= INFINITY; /* Neutral vector meson normalization factor */
double NN4= INFINITY; /* Neutral vector meson normalization factor */
double N11= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N12= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N13= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N14= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N21= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N22= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N23= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N24= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N31= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N32= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N33= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N34= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N41= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N42= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N43= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double N44= INFINITY; /* Neutral vector meson mixing matrix element in VA base */
double ytau= INFINITY; /* Lepton Yukawa coupling ( ytau ) */
double yc= INFINITY; /* U-quark Yukawa coupling ( yc ) */
double yt= INFINITY; /* U-quark Yukawa coupling ( yt ) */
double yb= INFINITY; /* D-quark Yukawa coupling ( yb ) */
double CKM11= INFINITY; /* CKM-Matrix ( CKM11 ) */
double CKM12= INFINITY; /* CKM-Matrix ( CKM12 ) */
double CKM21= INFINITY; /* CKM-Matrix ( CKM21 ) */
double CKM22= INFINITY; /* CKM-Matrix ( CKM22 ) */


double UnitStep(double x) {
	return (x<0.)?0.:1.;
}

void calculate_widths() {
    
    double complex cwH, cw1N, cw2N, cw1C, cw2C;
    
#define _warn_if_imag(x) if(fabs(cimag((x)))>1.e-15) fprintf(stderr, "warning: Width with imaginary part! [%e]\n",cimag(x))

    cwH = (-3*(4*Power(BM,2) - Power(HM,2))*Sqrt(-4*Power(BM,2)*Power(HM,2) + Power(HM,4))*Power(yb,2)*UnitStep(-2*Sqrt(Power(BM,2)) + Sqrt(Power(HM,2))))/(16.*Power(Power(HM,2),1.5)*Pi);
    cwH += (3*(Power(HM,2) - 4*Power(TM,2))*Sqrt(Power(HM,4) - 4*Power(HM,2)*Power(TM,2))*Power(yt,2)*UnitStep(Sqrt(Power(HM,2)) - 2*Sqrt(Power(TM,2))))/(16.*Power(Power(HM,2),1.5)*Pi);
    cwH += (-3*(4*Power(CM,2) - Power(HM,2))*Sqrt(-4*Power(CM,2)*Power(HM,2) + Power(HM,4))*Power(yc,2)*UnitStep(-2*Sqrt(Power(CM,2)) + Sqrt(Power(HM,2))))/(16.*Power(Power(HM,2),1.5)*Pi);
    cwH += ((Power(HM,2) - 4*Power(TAM,2))*Sqrt(Power(HM,4) - 4*Power(HM,2)*Power(TAM,2))*Power(ytau,2)*UnitStep(Sqrt(Power(HM,2)) - 2*Sqrt(Power(TAM,2))))/(16.*Power(Power(HM,2),1.5)*Pi);
    cwH += (Sqrt(Power(HM,4) - 4*Power(HM,2)*Power(M2C,2))*(Power(HM,4) - 4*Power(HM,2)*Power(M2C,2) + 12*Power(M2C,4))*Power(Power(C13,2)*Power(g2,2)*(1 - r3 + rs) + Power(gt,2)*(Power(C33,2)*(-r2 + rs) + Power(C23,2)*(r2 + rs)) - Sqrt(2)*C13*g2*gt*(C33*(-r2 + rs) + C23*(r2 - r3 + rs)),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - 2*Sqrt(Power(M2C,2))))/(256.*Power(Power(HM,2),1.5)*Power(M2C,4)*Pi);
    cwH += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(HM,2)*Power(M2C,2) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4))*(Power(HM,4) + Power(M1C,4) + 10*Power(M1C,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(HM,2)*(Power(M1C,2) + Power(M2C,2)))*Power(gt*(-2*gt*(C32*C33*(-r2 + rs) + C22*C23*(r2 + rs)) + Sqrt(2)*C13*g2*(C32*(-r2 + rs) + C22*(r2 - r3 + rs))) + C12*g2*(2*C13*g2*(-1 + r3 - rs) + Sqrt(2)*gt*(C33*(-r2 + rs) + C23*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(M1C,2)) - Sqrt(Power(M2C,2))))/(1024.*Power(Power(HM,2),1.5)*Power(M1C,2)*Power(M2C,2)*Pi);
    cwH += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(HM,2)*Power(MW,2) - 2*Power(M2C,2)*Power(MW,2) + Power(MW,4))*(Power(HM,4) + Power(M2C,4) + 10*Power(M2C,2)*Power(MW,2) + Power(MW,4) - 2*Power(HM,2)*(Power(M2C,2) + Power(MW,2)))*Power(gt*(-2*gt*(C31*C33*(-r2 + rs) + C21*C23*(r2 + rs)) + Sqrt(2)*C13*g2*(C31*(-r2 + rs) + C21*(r2 - r3 + rs))) + C11*g2*(2*C13*g2*(-1 + r3 - rs) + Sqrt(2)*gt*(C33*(-r2 + rs) + C23*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(M2C,2)) - Sqrt(Power(MW,2))))/(1024.*Power(Power(HM,2),1.5)*Power(M2C,2)*Power(MW,2)*Pi);
    cwH += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(HM,2)*Power(M2C,2) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4))*(Power(HM,4) + Power(M1C,4) + 10*Power(M1C,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(HM,2)*(Power(M1C,2) + Power(M2C,2)))*Power(gt*(-2*gt*(C32*C33*(-r2 + rs) + C22*C23*(r2 + rs)) + Sqrt(2)*C13*g2*(C32*(-r2 + rs) + C22*(r2 - r3 + rs))) + C12*g2*(2*C13*g2*(-1 + r3 - rs) + Sqrt(2)*gt*(C33*(-r2 + rs) + C23*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(M1C,2)) - Sqrt(Power(M2C,2))))/(1024.*Power(Power(HM,2),1.5)*Power(M1C,2)*Power(M2C,2)*Pi);
    cwH += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(HM,2)*Power(MW,2) - 2*Power(M2C,2)*Power(MW,2) + Power(MW,4))*(Power(HM,4) + Power(M2C,4) + 10*Power(M2C,2)*Power(MW,2) + Power(MW,4) - 2*Power(HM,2)*(Power(M2C,2) + Power(MW,2)))*Power(gt*(-2*gt*(C31*C33*(-r2 + rs) + C21*C23*(r2 + rs)) + Sqrt(2)*C13*g2*(C31*(-r2 + rs) + C21*(r2 - r3 + rs))) + C11*g2*(2*C13*g2*(-1 + r3 - rs) + Sqrt(2)*gt*(C33*(-r2 + rs) + C23*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(M2C,2)) - Sqrt(Power(MW,2))))/(1024.*Power(Power(HM,2),1.5)*Power(M2C,2)*Power(MW,2)*Pi);
    cwH += (Sqrt(Power(HM,4) - 4*Power(HM,2)*Power(M1C,2))*(Power(HM,4) - 4*Power(HM,2)*Power(M1C,2) + 12*Power(M1C,4))*Power(Power(C12,2)*Power(g2,2)*(1 - r3 + rs) + Power(gt,2)*(Power(C32,2)*(-r2 + rs) + Power(C22,2)*(r2 + rs)) - Sqrt(2)*C12*g2*gt*(C32*(-r2 + rs) + C22*(r2 - r3 + rs)),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - 2*Sqrt(Power(M1C,2))))/(256.*Power(Power(HM,2),1.5)*Power(M1C,4)*Pi);
    cwH += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(HM,2)*Power(MW,2) - 2*Power(M1C,2)*Power(MW,2) + Power(MW,4))*(Power(HM,4) + Power(M1C,4) + 10*Power(M1C,2)*Power(MW,2) + Power(MW,4) - 2*Power(HM,2)*(Power(M1C,2) + Power(MW,2)))*Power(gt*(-2*gt*(C31*C32*(-r2 + rs) + C21*C22*(r2 + rs)) + Sqrt(2)*C12*g2*(C31*(-r2 + rs) + C21*(r2 - r3 + rs))) + C11*g2*(2*C12*g2*(-1 + r3 - rs) + Sqrt(2)*gt*(C32*(-r2 + rs) + C22*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(M1C,2)) - Sqrt(Power(MW,2))))/(1024.*Power(Power(HM,2),1.5)*Power(M1C,2)*Power(MW,2)*Pi);
    cwH += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(HM,2)*Power(MW,2) - 2*Power(M1C,2)*Power(MW,2) + Power(MW,4))*(Power(HM,4) + Power(M1C,4) + 10*Power(M1C,2)*Power(MW,2) + Power(MW,4) - 2*Power(HM,2)*(Power(M1C,2) + Power(MW,2)))*Power(gt*(-2*gt*(C31*C32*(-r2 + rs) + C21*C22*(r2 + rs)) + Sqrt(2)*C12*g2*(C31*(-r2 + rs) + C21*(r2 - r3 + rs))) + C11*g2*(2*C12*g2*(-1 + r3 - rs) + Sqrt(2)*gt*(C32*(-r2 + rs) + C22*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(M1C,2)) - Sqrt(Power(MW,2))))/(1024.*Power(Power(HM,2),1.5)*Power(M1C,2)*Power(MW,2)*Pi);
    cwH += (Sqrt(Power(HM,4) - 4*Power(HM,2)*Power(MW,2))*(Power(HM,4) - 4*Power(HM,2)*Power(MW,2) + 12*Power(MW,4))*Power(Power(C11,2)*Power(g2,2)*(1 - r3 + rs) + Power(gt,2)*(Power(C31,2)*(-r2 + rs) + Power(C21,2)*(r2 + rs)) - Sqrt(2)*C11*g2*gt*(C31*(-r2 + rs) + C21*(r2 - r3 + rs)),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - 2*Sqrt(Power(MW,2))))/(256.*Power(Power(HM,2),1.5)*Power(MW,4)*Pi);
    cwH += (Sqrt(Power(HM,4) - 4*Power(HM,2)*Power(M2N,2))*(Power(HM,4) - 4*Power(HM,2)*Power(M2N,2) + 12*Power(M2N,4))*Power(Power(g1,2)*Power(N14,2)*(1 - r3 + rs) + Power(g2,2)*Power(N24,2)*(1 - r3 + rs) + Power(gt,2)*(Power(N44,2)*(-r2 + rs) + Power(N34,2)*(r2 + rs)) - Sqrt(2)*g2*gt*N24*(N44*(-r2 + rs) + N34*(r2 - r3 + rs)) + g1*N14*(-2*g2*N24*(1 + r2 - r3) + Sqrt(2)*gt*(N44*(r2 - rs) + N34*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - 2*Sqrt(Power(M2N,2))))/(512.*Power(Power(HM,2),1.5)*Power(M2N,4)*Pi);
    cwH += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1N,2) + Power(M1N,4) - 2*Power(HM,2)*Power(M2N,2) - 2*Power(M1N,2)*Power(M2N,2) + Power(M2N,4))*(Power(HM,4) + Power(M1N,4) + 10*Power(M1N,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(HM,2)*(Power(M1N,2) + Power(M2N,2)))*Power(2*Power(g1,2)*N13*N14*(1 - r3 + rs) + 2*Power(g2,2)*N23*N24*(1 - r3 + rs) + 2*Power(gt,2)*(N43*N44*(-r2 + rs) + N33*N34*(r2 + rs)) - Sqrt(2)*g2*gt*(N24*(N43*(-r2 + rs) + N33*(r2 - r3 + rs)) + N23*(N44*(-r2 + rs) + N34*(r2 - r3 + rs))) + g1*(-2*g2*(N14*N23 + N13*N24)*(1 + r2 - r3) + Sqrt(2)*gt*(N14*(N43*(r2 - rs) + N33*(r2 - r3 + rs)) + N13*(N44*(r2 - rs) + N34*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(M1N,2)) - Sqrt(Power(M2N,2))))/(1024.*Power(Power(HM,2),1.5)*Power(M1N,2)*Power(M2N,2)*Pi);
    cwH += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(HM,2)*Power(ZM,2) - 2*Power(M2N,2)*Power(ZM,2) + Power(ZM,4))*(Power(HM,4) + Power(M2N,4) + 10*Power(M2N,2)*Power(ZM,2) + Power(ZM,4) - 2*Power(HM,2)*(Power(M2N,2) + Power(ZM,2)))*Power(2*Power(g1,2)*N12*N14*(1 - r3 + rs) + 2*Power(g2,2)*N22*N24*(1 - r3 + rs) + 2*Power(gt,2)*(N42*N44*(-r2 + rs) + N32*N34*(r2 + rs)) - Sqrt(2)*g2*gt*(N24*(N42*(-r2 + rs) + N32*(r2 - r3 + rs)) + N22*(N44*(-r2 + rs) + N34*(r2 - r3 + rs))) + g1*(-2*g2*(N14*N22 + N12*N24)*(1 + r2 - r3) + Sqrt(2)*gt*(N14*(N42*(r2 - rs) + N32*(r2 - r3 + rs)) + N12*(N44*(r2 - rs) + N34*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(M2N,2)) - Sqrt(Power(ZM,2))))/(1024.*Power(Power(HM,2),1.5)*Power(M2N,2)*Power(ZM,2)*Pi);
    cwH += (3*Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M2N,2) + Power(M2N,4))*Power(2*Power(g1,2)*N11*N14*(1 - r3 + rs) + 2*Power(g2,2)*N21*N24*(1 - r3 + rs) + 2*Power(gt,2)*(N41*N44*(-r2 + rs) + N31*N34*(r2 + rs)) - Sqrt(2)*g2*gt*(N24*(N41*(-r2 + rs) + N31*(r2 - r3 + rs)) + N21*(N44*(-r2 + rs) + N34*(r2 - r3 + rs))) + g1*(-2*g2*(N14*N21 + N11*N24)*(1 + r2 - r3) + Sqrt(2)*gt*(N14*(N41*(r2 - rs) + N31*(r2 - r3 + rs)) + N11*(N44*(r2 - rs) + N34*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(M2N,2))))/(256.*Power(Power(HM,2),1.5)*Pi);
    cwH += (Sqrt(Power(HM,4) - 4*Power(HM,2)*Power(M1N,2))*(Power(HM,4) - 4*Power(HM,2)*Power(M1N,2) + 12*Power(M1N,4))*Power(Power(g1,2)*Power(N13,2)*(1 - r3 + rs) + Power(g2,2)*Power(N23,2)*(1 - r3 + rs) + Power(gt,2)*(Power(N43,2)*(-r2 + rs) + Power(N33,2)*(r2 + rs)) - Sqrt(2)*g2*gt*N23*(N43*(-r2 + rs) + N33*(r2 - r3 + rs)) + g1*N13*(-2*g2*N23*(1 + r2 - r3) + Sqrt(2)*gt*(N43*(r2 - rs) + N33*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - 2*Sqrt(Power(M1N,2))))/(512.*Power(Power(HM,2),1.5)*Power(M1N,4)*Pi);
    cwH += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1N,2) + Power(M1N,4) - 2*Power(HM,2)*Power(ZM,2) - 2*Power(M1N,2)*Power(ZM,2) + Power(ZM,4))*(Power(HM,4) + Power(M1N,4) + 10*Power(M1N,2)*Power(ZM,2) + Power(ZM,4) - 2*Power(HM,2)*(Power(M1N,2) + Power(ZM,2)))*Power(2*Power(g1,2)*N12*N13*(1 - r3 + rs) + 2*Power(g2,2)*N22*N23*(1 - r3 + rs) + 2*Power(gt,2)*(N42*N43*(-r2 + rs) + N32*N33*(r2 + rs)) - Sqrt(2)*g2*gt*(N23*(N42*(-r2 + rs) + N32*(r2 - r3 + rs)) + N22*(N43*(-r2 + rs) + N33*(r2 - r3 + rs))) + g1*(-2*g2*(N13*N22 + N12*N23)*(1 + r2 - r3) + Sqrt(2)*gt*(N13*(N42*(r2 - rs) + N32*(r2 - r3 + rs)) + N12*(N43*(r2 - rs) + N33*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(M1N,2)) - Sqrt(Power(ZM,2))))/(1024.*Power(Power(HM,2),1.5)*Power(M1N,2)*Power(ZM,2)*Pi);
    cwH += (3*Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1N,2) + Power(M1N,4))*Power(2*Power(g1,2)*N11*N13*(1 - r3 + rs) + 2*Power(g2,2)*N21*N23*(1 - r3 + rs) + 2*Power(gt,2)*(N41*N43*(-r2 + rs) + N31*N33*(r2 + rs)) - Sqrt(2)*g2*gt*(N23*(N41*(-r2 + rs) + N31*(r2 - r3 + rs)) + N21*(N43*(-r2 + rs) + N33*(r2 - r3 + rs))) + g1*(-2*g2*(N13*N21 + N11*N23)*(1 + r2 - r3) + Sqrt(2)*gt*(N13*(N41*(r2 - rs) + N31*(r2 - r3 + rs)) + N11*(N43*(r2 - rs) + N33*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(M1N,2))))/(256.*Power(Power(HM,2),1.5)*Pi);
    cwH += (Sqrt(Power(HM,4) - 4*Power(HM,2)*Power(ZM,2))*(Power(HM,4) - 4*Power(HM,2)*Power(ZM,2) + 12*Power(ZM,4))*Power(Power(g1,2)*Power(N12,2)*(1 - r3 + rs) + Power(g2,2)*Power(N22,2)*(1 - r3 + rs) + Power(gt,2)*(Power(N42,2)*(-r2 + rs) + Power(N32,2)*(r2 + rs)) - Sqrt(2)*g2*gt*N22*(N42*(-r2 + rs) + N32*(r2 - r3 + rs)) + g1*N12*(-2*g2*N22*(1 + r2 - r3) + Sqrt(2)*gt*(N42*(r2 - rs) + N32*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - 2*Sqrt(Power(ZM,2))))/(512.*Power(Power(HM,2),1.5)*Power(ZM,4)*Pi);
    cwH += (3*Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(ZM,2) + Power(ZM,4))*Power(2*Power(g1,2)*N11*N12*(1 - r3 + rs) + 2*Power(g2,2)*N21*N22*(1 - r3 + rs) + 2*Power(gt,2)*(N41*N42*(-r2 + rs) + N31*N32*(r2 + rs)) - Sqrt(2)*g2*gt*(N22*(N41*(-r2 + rs) + N31*(r2 - r3 + rs)) + N21*(N42*(-r2 + rs) + N32*(r2 - r3 + rs))) + g1*(-2*g2*(N12*N21 + N11*N22)*(1 + r2 - r3) + Sqrt(2)*gt*(N12*(N41*(r2 - rs) + N31*(r2 - r3 + rs)) + N11*(N42*(r2 - rs) + N32*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2)) - Sqrt(Power(ZM,2))))/(256.*Power(Power(HM,2),1.5)*Pi);
    cwH += (Sqrt(Power(HM,4))*Power(Power(g1,2)*Power(N11,2)*(1 - r3 + rs) + Power(g2,2)*Power(N21,2)*(1 - r3 + rs) + Power(gt,2)*(Power(N41,2)*(-r2 + rs) + Power(N31,2)*(r2 + rs)) - Sqrt(2)*g2*gt*N21*(N41*(-r2 + rs) + N31*(r2 - r3 + rs)) + g1*N11*(-2*g2*N21*(1 + r2 - r3) + Sqrt(2)*gt*(N41*(r2 - rs) + N31*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(Sqrt(Power(HM,2))))/(32.*Power(Power(HM,2),1.5)*Pi);    
    wH=creal(cwH);
    _warn_if_imag(cwH);
    
    cw1N = -(Sqrt(-4*Power(BM,2)*Power(M1N,2) + Power(M1N,4))*(-(Power(M1N,2)*((5*Power(g1,2)*Power(N13,2))/144. - (g1*g2*N13*N23)/24. + (Power(g2,2)*Power(N23,2))/16.)) + Power(BM,2)*((5*Power(g1,2)*Power(N13,2))/144. - (g1*g2*N13*N23)/4. + (Power(g2,2)*Power(N23,2))/16. + (g1*N13*((g1*N13)/2. - (g2*N23)/4.))/6.))*UnitStep(-2*Sqrt(Power(BM,2)) + Sqrt(Power(M1N,2))))/(2.*Power(Power(M1N,2),1.5)*Pi);
    cw1N += (Sqrt(Power(M1N,4))*((5*Power(g1,2)*Power(N13,2))/144. - (g1*g2*N13*N23)/24. + (Power(g2,2)*Power(N23,2))/16.)*UnitStep(Sqrt(Power(M1N,2))))/(2.*Sqrt(Power(M1N,2))*Pi);
    cw1N += (Sqrt(Power(M1N,4))*((5*Power(g1,2)*Power(N13,2))/144. - (g1*g2*N13*N23)/24. + (Power(g2,2)*Power(N23,2))/16.)*UnitStep(Sqrt(Power(M1N,2))))/(2.*Sqrt(Power(M1N,2))*Pi);
    cw1N += (Sqrt(Power(M1N,4) - 4*Power(M1N,2)*Power(TM,2))*(Power(M1N,2)*((17*Power(g1,2)*Power(N13,2))/144. + (g1*g2*N13*N23)/24. + (Power(g2,2)*Power(N23,2))/16.) - Power(TM,2)*((17*Power(g1,2)*Power(N13,2))/144. - (g1*g2*N13*N23)/2. + (Power(g2,2)*Power(N23,2))/16. + (g1*N13*(-(g1*N13) + (g2*N23)/4.))/6.))*UnitStep(Sqrt(Power(M1N,2)) - 2*Sqrt(Power(TM,2))))/(2.*Power(Power(M1N,2),1.5)*Pi);
    cw1N += -(Sqrt(-4*Power(CM,2)*Power(M1N,2) + Power(M1N,4))*(-(Power(M1N,2)*((17*Power(g1,2)*Power(N13,2))/144. + (g1*g2*N13*N23)/24. + (Power(g2,2)*Power(N23,2))/16.)) + Power(CM,2)*((17*Power(g1,2)*Power(N13,2))/144. - (g1*g2*N13*N23)/2. + (Power(g2,2)*Power(N23,2))/16. + (g1*N13*(-(g1*N13) + (g2*N23)/4.))/6.))*UnitStep(-2*Sqrt(Power(CM,2)) + Sqrt(Power(M1N,2))))/(2.*Power(Power(M1N,2),1.5)*Pi);
    cw1N += (Sqrt(Power(M1N,4))*((17*Power(g1,2)*Power(N13,2))/144. + (g1*g2*N13*N23)/24. + (Power(g2,2)*Power(N23,2))/16.)*UnitStep(Sqrt(Power(M1N,2))))/(2.*Sqrt(Power(M1N,2))*Pi);
    cw1N += (Sqrt(Power(M1N,4) - 4*Power(M1N,2)*Power(TAM,2))*(Power(M1N,2)*((5*Power(g1,2)*Power(N13,2))/16. + (g1*g2*N13*N23)/8. + (Power(g2,2)*Power(N23,2))/16.) - Power(TAM,2)*((5*Power(g1,2)*Power(N13,2))/16. - (3*g1*g2*N13*N23)/4. + (Power(g2,2)*Power(N23,2))/16. - (g1*N13*((3*g1*N13)/2. - (g2*N23)/4.))/2.))*UnitStep(Sqrt(Power(M1N,2)) - 2*Sqrt(Power(TAM,2))))/(6.*Power(Power(M1N,2),1.5)*Pi);
    cw1N += (Sqrt(Power(M1N,4))*((5*Power(g1,2)*Power(N13,2))/16. + (g1*g2*N13*N23)/8. + (Power(g2,2)*Power(N23,2))/16.)*UnitStep(Sqrt(Power(M1N,2))))/(6.*Sqrt(Power(M1N,2))*Pi);
    cw1N += (Sqrt(Power(M1N,4))*((5*Power(g1,2)*Power(N13,2))/16. + (g1*g2*N13*N23)/8. + (Power(g2,2)*Power(N23,2))/16.)*UnitStep(Sqrt(Power(M1N,2))))/(6.*Sqrt(Power(M1N,2))*Pi);
    cw1N += (Sqrt(Power(M1N,4))*Power(-(g1*N13)/4. + (g2*N23)/4.,2)*UnitStep(Sqrt(Power(M1N,2))))/(6.*Sqrt(Power(M1N,2))*Pi);
    cw1N += (Sqrt(Power(M1N,4))*Power(-(g1*N13)/4. + (g2*N23)/4.,2)*UnitStep(Sqrt(Power(M1N,2))))/(6.*Sqrt(Power(M1N,2))*Pi);
    cw1N += (Sqrt(Power(M1N,4))*Power(-(g1*N13)/4. + (g2*N23)/4.,2)*UnitStep(Sqrt(Power(M1N,2))))/(6.*Sqrt(Power(M1N,2))*Pi);
    cw1N += (Sqrt(Power(M1N,4) - 4*Power(M1N,2)*Power(M2C,2))*(Power(M1N,6) + 16*Power(M1N,4)*Power(M2C,2) - 68*Power(M1N,2)*Power(M2C,4) - 48*Power(M2C,6))*Power(-(Power(C13,2)*g2*N23) - (gt*(2*C23*C33*N33 + Power(C23,2)*N43 + Power(C33,2)*N43))/Sqrt(2),2)*UnitStep(Sqrt(Power(M1N,2)) - 2*Sqrt(Power(M2C,2))))/(192.*Power(Power(M1N,2),1.5)*Power(M2C,4)*Pi);
    cw1N += -(Sqrt(Power(M1N,4) - 2*Power(M1N,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(M1N,2)*Power(M2C,2) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4))*(Power(M1N,8) + 8*Power(M1N,6)*(Power(M1C,2) + Power(M2C,2)) + Power(Power(M1C,2) - Power(M2C,2),2)*(Power(M1C,4) + 10*Power(M1C,2)*Power(M2C,2) + Power(M2C,4)) - 2*Power(M1N,4)*(9*Power(M1C,4) + 16*Power(M1C,2)*Power(M2C,2) + 9*Power(M2C,4)) + 8*Power(M1N,2)*(Power(M1C,6) - 4*Power(M1C,4)*Power(M2C,2) - 4*Power(M1C,2)*Power(M2C,4) + Power(M2C,6)))*(-(C12*C13*g2*N23) - (gt*(C23*C32*N33 + C22*C33*N33 + C22*C23*N43 + C32*C33*N43))/Sqrt(2))*(2*C12*C13*g2*N23 + Sqrt(2)*gt*(C23*C32*N33 + C22*C33*N33 + C22*C23*N43 + C32*C33*N43))*UnitStep(Sqrt(Power(M1N,2)) - Sqrt(Power(M1C,2)) - Sqrt(Power(M2C,2))))/(384.*Power(M1N,2)*Power(Power(M1N,2),1.5)*Power(M1C,2)*Power(M2C,2)*Pi);
    cw1N += -(Sqrt(Power(M1N,4) - 2*Power(M1N,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(M1N,2)*Power(MW,2) - 2*Power(M2C,2)*Power(MW,2) + Power(MW,4))*(Power(M1N,8) + 8*Power(M1N,6)*(Power(M2C,2) + Power(MW,2)) + Power(Power(M2C,2) - Power(MW,2),2)*(Power(M2C,4) + 10*Power(M2C,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M1N,4)*(9*Power(M2C,4) + 16*Power(M2C,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M1N,2)*(Power(M2C,6) - 4*Power(M2C,4)*Power(MW,2) - 4*Power(M2C,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C13*g2*N23) - (gt*(C23*C31*N33 + C21*C33*N33 + C21*C23*N43 + C31*C33*N43))/Sqrt(2))*(2*C11*C13*g2*N23 + Sqrt(2)*gt*(C23*C31*N33 + C21*C33*N33 + C21*C23*N43 + C31*C33*N43))*UnitStep(Sqrt(Power(M1N,2)) - Sqrt(Power(M2C,2)) - Sqrt(Power(MW,2))))/(384.*Power(M1N,2)*Power(Power(M1N,2),1.5)*Power(M2C,2)*Power(MW,2)*Pi);
    cw1N += -(Sqrt(Power(M1N,4) - 2*Power(M1N,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(M1N,2)*Power(M2C,2) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4))*(Power(M1N,8) + 8*Power(M1N,6)*(Power(M1C,2) + Power(M2C,2)) + Power(Power(M1C,2) - Power(M2C,2),2)*(Power(M1C,4) + 10*Power(M1C,2)*Power(M2C,2) + Power(M2C,4)) - 2*Power(M1N,4)*(9*Power(M1C,4) + 16*Power(M1C,2)*Power(M2C,2) + 9*Power(M2C,4)) + 8*Power(M1N,2)*(Power(M1C,6) - 4*Power(M1C,4)*Power(M2C,2) - 4*Power(M1C,2)*Power(M2C,4) + Power(M2C,6)))*(-(C12*C13*g2*N23) - (gt*(C23*C32*N33 + C22*C33*N33 + C22*C23*N43 + C32*C33*N43))/Sqrt(2))*(2*C12*C13*g2*N23 + Sqrt(2)*gt*(C23*C32*N33 + C22*C33*N33 + C22*C23*N43 + C32*C33*N43))*UnitStep(Sqrt(Power(M1N,2)) - Sqrt(Power(M1C,2)) - Sqrt(Power(M2C,2))))/(384.*Power(M1N,2)*Power(Power(M1N,2),1.5)*Power(M1C,2)*Power(M2C,2)*Pi);
    cw1N += -(Sqrt(Power(M1N,4) - 2*Power(M1N,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(M1N,2)*Power(MW,2) - 2*Power(M2C,2)*Power(MW,2) + Power(MW,4))*(Power(M1N,8) + 8*Power(M1N,6)*(Power(M2C,2) + Power(MW,2)) + Power(Power(M2C,2) - Power(MW,2),2)*(Power(M2C,4) + 10*Power(M2C,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M1N,4)*(9*Power(M2C,4) + 16*Power(M2C,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M1N,2)*(Power(M2C,6) - 4*Power(M2C,4)*Power(MW,2) - 4*Power(M2C,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C13*g2*N23) - (gt*(C23*C31*N33 + C21*C33*N33 + C21*C23*N43 + C31*C33*N43))/Sqrt(2))*(2*C11*C13*g2*N23 + Sqrt(2)*gt*(C23*C31*N33 + C21*C33*N33 + C21*C23*N43 + C31*C33*N43))*UnitStep(Sqrt(Power(M1N,2)) - Sqrt(Power(M2C,2)) - Sqrt(Power(MW,2))))/(384.*Power(M1N,2)*Power(Power(M1N,2),1.5)*Power(M2C,2)*Power(MW,2)*Pi);
    cw1N += (Sqrt(Power(M1N,4) - 4*Power(M1N,2)*Power(M1C,2))*(Power(M1N,6) + 16*Power(M1N,4)*Power(M1C,2) - 68*Power(M1N,2)*Power(M1C,4) - 48*Power(M1C,6))*Power(-(Power(C12,2)*g2*N23) - (gt*(2*C22*C32*N33 + Power(C22,2)*N43 + Power(C32,2)*N43))/Sqrt(2),2)*UnitStep(Sqrt(Power(M1N,2)) - 2*Sqrt(Power(M1C,2))))/(192.*Power(Power(M1N,2),1.5)*Power(M1C,4)*Pi);
    cw1N += -(Sqrt(Power(M1N,4) - 2*Power(M1N,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(M1N,2)*Power(MW,2) - 2*Power(M1C,2)*Power(MW,2) + Power(MW,4))*(Power(M1N,8) + 8*Power(M1N,6)*(Power(M1C,2) + Power(MW,2)) + Power(Power(M1C,2) - Power(MW,2),2)*(Power(M1C,4) + 10*Power(M1C,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M1N,4)*(9*Power(M1C,4) + 16*Power(M1C,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M1N,2)*(Power(M1C,6) - 4*Power(M1C,4)*Power(MW,2) - 4*Power(M1C,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C12*g2*N23) - (gt*(C22*C31*N33 + C21*C32*N33 + C21*C22*N43 + C31*C32*N43))/Sqrt(2))*(2*C11*C12*g2*N23 + Sqrt(2)*gt*(C22*C31*N33 + C21*C32*N33 + C21*C22*N43 + C31*C32*N43))*UnitStep(Sqrt(Power(M1N,2)) - Sqrt(Power(M1C,2)) - Sqrt(Power(MW,2))))/(384.*Power(M1N,2)*Power(Power(M1N,2),1.5)*Power(M1C,2)*Power(MW,2)*Pi);
    cw1N += -(Sqrt(Power(M1N,4) - 2*Power(M1N,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(M1N,2)*Power(MW,2) - 2*Power(M1C,2)*Power(MW,2) + Power(MW,4))*(Power(M1N,8) + 8*Power(M1N,6)*(Power(M1C,2) + Power(MW,2)) + Power(Power(M1C,2) - Power(MW,2),2)*(Power(M1C,4) + 10*Power(M1C,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M1N,4)*(9*Power(M1C,4) + 16*Power(M1C,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M1N,2)*(Power(M1C,6) - 4*Power(M1C,4)*Power(MW,2) - 4*Power(M1C,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C12*g2*N23) - (gt*(C22*C31*N33 + C21*C32*N33 + C21*C22*N43 + C31*C32*N43))/Sqrt(2))*(2*C11*C12*g2*N23 + Sqrt(2)*gt*(C22*C31*N33 + C21*C32*N33 + C21*C22*N43 + C31*C32*N43))*UnitStep(Sqrt(Power(M1N,2)) - Sqrt(Power(M1C,2)) - Sqrt(Power(MW,2))))/(384.*Power(M1N,2)*Power(Power(M1N,2),1.5)*Power(M1C,2)*Power(MW,2)*Pi);
    cw1N += (Sqrt(Power(M1N,4) - 4*Power(M1N,2)*Power(MW,2))*(Power(M1N,6) + 16*Power(M1N,4)*Power(MW,2) - 68*Power(M1N,2)*Power(MW,4) - 48*Power(MW,6))*Power(-(Power(C11,2)*g2*N23) - (gt*(2*C21*C31*N33 + Power(C21,2)*N43 + Power(C31,2)*N43))/Sqrt(2),2)*UnitStep(Sqrt(Power(M1N,2)) - 2*Sqrt(Power(MW,2))))/(192.*Power(Power(M1N,2),1.5)*Power(MW,4)*Pi);
    cw1N += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1N,2) + Power(M1N,4) - 2*Power(HM,2)*Power(M2N,2) - 2*Power(M1N,2)*Power(M2N,2) + Power(M2N,4))*(Power(HM,4) + Power(M1N,4) + 10*Power(M1N,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(HM,2)*(Power(M1N,2) + Power(M2N,2)))*Power(2*Power(g1,2)*N13*N14*(1 - r3 + rs) + 2*Power(g2,2)*N23*N24*(1 - r3 + rs) + 2*Power(gt,2)*(N43*N44*(-r2 + rs) + N33*N34*(r2 + rs)) - Sqrt(2)*g2*gt*(N24*(N43*(-r2 + rs) + N33*(r2 - r3 + rs)) + N23*(N44*(-r2 + rs) + N34*(r2 - r3 + rs))) + g1*(-2*g2*(N14*N23 + N13*N24)*(1 + r2 - r3) + Sqrt(2)*gt*(N14*(N43*(r2 - rs) + N33*(r2 - r3 + rs)) + N13*(N44*(r2 - rs) + N34*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(-Sqrt(Power(HM,2)) + Sqrt(Power(M1N,2)) - Sqrt(Power(M2N,2))))/(3072.*Power(M1N,2)*Power(Power(M1N,2),1.5)*Power(M2N,2)*Pi);
    cw1N += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1N,2) + Power(M1N,4) - 2*Power(HM,2)*Power(ZM,2) - 2*Power(M1N,2)*Power(ZM,2) + Power(ZM,4))*(Power(HM,4) + Power(M1N,4) + 10*Power(M1N,2)*Power(ZM,2) + Power(ZM,4) - 2*Power(HM,2)*(Power(M1N,2) + Power(ZM,2)))*Power(2*Power(g1,2)*N12*N13*(1 - r3 + rs) + 2*Power(g2,2)*N22*N23*(1 - r3 + rs) + 2*Power(gt,2)*(N42*N43*(-r2 + rs) + N32*N33*(r2 + rs)) - Sqrt(2)*g2*gt*(N23*(N42*(-r2 + rs) + N32*(r2 - r3 + rs)) + N22*(N43*(-r2 + rs) + N33*(r2 - r3 + rs))) + g1*(-2*g2*(N13*N22 + N12*N23)*(1 + r2 - r3) + Sqrt(2)*gt*(N13*(N42*(r2 - rs) + N32*(r2 - r3 + rs)) + N12*(N43*(r2 - rs) + N33*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(-Sqrt(Power(HM,2)) + Sqrt(Power(M1N,2)) - Sqrt(Power(ZM,2))))/(3072.*Power(M1N,2)*Power(Power(M1N,2),1.5)*Power(ZM,2)*Pi);
    cw1N += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1N,2) + Power(M1N,4))*Power(2*Power(g1,2)*N11*N13*(1 - r3 + rs) + 2*Power(g2,2)*N21*N23*(1 - r3 + rs) + 2*Power(gt,2)*(N41*N43*(-r2 + rs) + N31*N33*(r2 + rs)) - Sqrt(2)*g2*gt*(N23*(N41*(-r2 + rs) + N31*(r2 - r3 + rs)) + N21*(N43*(-r2 + rs) + N33*(r2 - r3 + rs))) + g1*(-2*g2*(N13*N21 + N11*N23)*(1 + r2 - r3) + Sqrt(2)*gt*(N13*(N41*(r2 - rs) + N31*(r2 - r3 + rs)) + N11*(N43*(r2 - rs) + N33*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(-Sqrt(Power(HM,2)) + Sqrt(Power(M1N,2))))/(256.*Power(Power(M1N,2),1.5)*Pi);
	w1N=creal(cw1N);
    _warn_if_imag(cw1N);
    
    cw2N = -(Sqrt(-4*Power(BM,2)*Power(M2N,2) + Power(M2N,4))*(-(Power(M2N,2)*((5*Power(g1,2)*Power(N14,2))/144. - (g1*g2*N14*N24)/24. + (Power(g2,2)*Power(N24,2))/16.)) + Power(BM,2)*((5*Power(g1,2)*Power(N14,2))/144. - (g1*g2*N14*N24)/4. + (Power(g2,2)*Power(N24,2))/16. + (g1*N14*((g1*N14)/2. - (g2*N24)/4.))/6.))*UnitStep(-2*Sqrt(Power(BM,2)) + Sqrt(Power(M2N,2))))/(2.*Power(Power(M2N,2),1.5)*Pi);
    cw2N += (Sqrt(Power(M2N,4))*((5*Power(g1,2)*Power(N14,2))/144. - (g1*g2*N14*N24)/24. + (Power(g2,2)*Power(N24,2))/16.)*UnitStep(Sqrt(Power(M2N,2))))/(2.*Sqrt(Power(M2N,2))*Pi);
    cw2N += (Sqrt(Power(M2N,4))*((5*Power(g1,2)*Power(N14,2))/144. - (g1*g2*N14*N24)/24. + (Power(g2,2)*Power(N24,2))/16.)*UnitStep(Sqrt(Power(M2N,2))))/(2.*Sqrt(Power(M2N,2))*Pi);
    cw2N += -(Sqrt(-4*Power(TM,2)*Power(M2N,2) + Power(M2N,4))*(-(Power(M2N,2)*((17*Power(g1,2)*Power(N14,2))/144. + (g1*g2*N14*N24)/24. + (Power(g2,2)*Power(N24,2))/16.)) + Power(TM,2)*((17*Power(g1,2)*Power(N14,2))/144. - (g1*g2*N14*N24)/2. + (Power(g2,2)*Power(N24,2))/16. + (g1*N14*(-(g1*N14) + (g2*N24)/4.))/6.))*UnitStep(-2*Sqrt(Power(TM,2)) + Sqrt(Power(M2N,2))))/(2.*Power(Power(M2N,2),1.5)*Pi);
    cw2N += -(Sqrt(-4*Power(CM,2)*Power(M2N,2) + Power(M2N,4))*(-(Power(M2N,2)*((17*Power(g1,2)*Power(N14,2))/144. + (g1*g2*N14*N24)/24. + (Power(g2,2)*Power(N24,2))/16.)) + Power(CM,2)*((17*Power(g1,2)*Power(N14,2))/144. - (g1*g2*N14*N24)/2. + (Power(g2,2)*Power(N24,2))/16. + (g1*N14*(-(g1*N14) + (g2*N24)/4.))/6.))*UnitStep(-2*Sqrt(Power(CM,2)) + Sqrt(Power(M2N,2))))/(2.*Power(Power(M2N,2),1.5)*Pi);
    cw2N += (Sqrt(Power(M2N,4))*((17*Power(g1,2)*Power(N14,2))/144. + (g1*g2*N14*N24)/24. + (Power(g2,2)*Power(N24,2))/16.)*UnitStep(Sqrt(Power(M2N,2))))/(2.*Sqrt(Power(M2N,2))*Pi);
    cw2N += -(Sqrt(-4*Power(TAM,2)*Power(M2N,2) + Power(M2N,4))*(-(Power(M2N,2)*((5*Power(g1,2)*Power(N14,2))/16. + (g1*g2*N14*N24)/8. + (Power(g2,2)*Power(N24,2))/16.)) + Power(TAM,2)*((5*Power(g1,2)*Power(N14,2))/16. - (3*g1*g2*N14*N24)/4. + (Power(g2,2)*Power(N24,2))/16. - (g1*N14*((3*g1*N14)/2. - (g2*N24)/4.))/2.))*UnitStep(-2*Sqrt(Power(TAM,2)) + Sqrt(Power(M2N,2))))/(6.*Power(Power(M2N,2),1.5)*Pi);
    cw2N += (Sqrt(Power(M2N,4))*((5*Power(g1,2)*Power(N14,2))/16. + (g1*g2*N14*N24)/8. + (Power(g2,2)*Power(N24,2))/16.)*UnitStep(Sqrt(Power(M2N,2))))/(6.*Sqrt(Power(M2N,2))*Pi);
    cw2N += (Sqrt(Power(M2N,4))*((5*Power(g1,2)*Power(N14,2))/16. + (g1*g2*N14*N24)/8. + (Power(g2,2)*Power(N24,2))/16.)*UnitStep(Sqrt(Power(M2N,2))))/(6.*Sqrt(Power(M2N,2))*Pi);
    cw2N += (Sqrt(Power(M2N,4))*Power(-(g1*N14)/4. + (g2*N24)/4.,2)*UnitStep(Sqrt(Power(M2N,2))))/(6.*Sqrt(Power(M2N,2))*Pi);
    cw2N += (Sqrt(Power(M2N,4))*Power(-(g1*N14)/4. + (g2*N24)/4.,2)*UnitStep(Sqrt(Power(M2N,2))))/(6.*Sqrt(Power(M2N,2))*Pi);
    cw2N += (Sqrt(Power(M2N,4))*Power(-(g1*N14)/4. + (g2*N24)/4.,2)*UnitStep(Sqrt(Power(M2N,2))))/(6.*Sqrt(Power(M2N,2))*Pi);
    cw2N += (Sqrt(Power(M2N,4) - 4*Power(M2N,2)*Power(M2C,2))*(Power(M2N,6) + 16*Power(M2N,4)*Power(M2C,2) - 68*Power(M2N,2)*Power(M2C,4) - 48*Power(M2C,6))*Power(-(Power(C13,2)*g2*N24) - (gt*(2*C23*C33*N34 + Power(C23,2)*N44 + Power(C33,2)*N44))/Sqrt(2),2)*UnitStep(Sqrt(Power(M2N,2)) - 2*Sqrt(Power(M2C,2))))/(192.*Power(Power(M2N,2),1.5)*Power(M2C,4)*Pi);
    cw2N += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(M1C,2)*Power(M2C,2) - 2*Power(M2N,2)*Power(M2C,2) + Power(M2C,4))*(Power(M1C,8) + 8*Power(M1C,6)*(Power(M2N,2) + Power(M2C,2)) + Power(Power(M2N,2) - Power(M2C,2),2)*(Power(M2N,4) + 10*Power(M2N,2)*Power(M2C,2) + Power(M2C,4)) - 2*Power(M1C,4)*(9*Power(M2N,4) + 16*Power(M2N,2)*Power(M2C,2) + 9*Power(M2C,4)) + 8*Power(M1C,2)*(Power(M2N,6) - 4*Power(M2N,4)*Power(M2C,2) - 4*Power(M2N,2)*Power(M2C,4) + Power(M2C,6)))*(-(C12*C13*g2*N24) - (gt*(C23*C32*N34 + C22*C33*N34 + C22*C23*N44 + C32*C33*N44))/Sqrt(2))*(2*C12*C13*g2*N24 + Sqrt(2)*gt*(C23*C32*N34 + C22*C33*N34 + C22*C23*N44 + C32*C33*N44))*UnitStep(-Sqrt(Power(M1C,2)) + Sqrt(Power(M2N,2)) - Sqrt(Power(M2C,2))))/(384.*Power(M1C,2)*Power(M2N,2)*Power(Power(M2N,2),1.5)*Power(M2C,2)*Pi);
    cw2N += -(Sqrt(Power(M2N,4) - 2*Power(M2N,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(M2N,2)*Power(MW,2) - 2*Power(M2C,2)*Power(MW,2) + Power(MW,4))*(Power(M2N,8) + 8*Power(M2N,6)*(Power(M2C,2) + Power(MW,2)) + Power(Power(M2C,2) - Power(MW,2),2)*(Power(M2C,4) + 10*Power(M2C,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M2N,4)*(9*Power(M2C,4) + 16*Power(M2C,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M2N,2)*(Power(M2C,6) - 4*Power(M2C,4)*Power(MW,2) - 4*Power(M2C,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C13*g2*N24) - (gt*(C23*C31*N34 + C21*C33*N34 + C21*C23*N44 + C31*C33*N44))/Sqrt(2))*(2*C11*C13*g2*N24 + Sqrt(2)*gt*(C23*C31*N34 + C21*C33*N34 + C21*C23*N44 + C31*C33*N44))*UnitStep(Sqrt(Power(M2N,2)) - Sqrt(Power(M2C,2)) - Sqrt(Power(MW,2))))/(384.*Power(M2N,2)*Power(Power(M2N,2),1.5)*Power(M2C,2)*Power(MW,2)*Pi);
    cw2N += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(M1C,2)*Power(M2C,2) - 2*Power(M2N,2)*Power(M2C,2) + Power(M2C,4))*(Power(M1C,8) + 8*Power(M1C,6)*(Power(M2N,2) + Power(M2C,2)) + Power(Power(M2N,2) - Power(M2C,2),2)*(Power(M2N,4) + 10*Power(M2N,2)*Power(M2C,2) + Power(M2C,4)) - 2*Power(M1C,4)*(9*Power(M2N,4) + 16*Power(M2N,2)*Power(M2C,2) + 9*Power(M2C,4)) + 8*Power(M1C,2)*(Power(M2N,6) - 4*Power(M2N,4)*Power(M2C,2) - 4*Power(M2N,2)*Power(M2C,4) + Power(M2C,6)))*(-(C12*C13*g2*N24) - (gt*(C23*C32*N34 + C22*C33*N34 + C22*C23*N44 + C32*C33*N44))/Sqrt(2))*(2*C12*C13*g2*N24 + Sqrt(2)*gt*(C23*C32*N34 + C22*C33*N34 + C22*C23*N44 + C32*C33*N44))*UnitStep(-Sqrt(Power(M1C,2)) + Sqrt(Power(M2N,2)) - Sqrt(Power(M2C,2))))/(384.*Power(M1C,2)*Power(M2N,2)*Power(Power(M2N,2),1.5)*Power(M2C,2)*Pi);
    cw2N += -(Sqrt(Power(M2N,4) - 2*Power(M2N,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(M2N,2)*Power(MW,2) - 2*Power(M2C,2)*Power(MW,2) + Power(MW,4))*(Power(M2N,8) + 8*Power(M2N,6)*(Power(M2C,2) + Power(MW,2)) + Power(Power(M2C,2) - Power(MW,2),2)*(Power(M2C,4) + 10*Power(M2C,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M2N,4)*(9*Power(M2C,4) + 16*Power(M2C,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M2N,2)*(Power(M2C,6) - 4*Power(M2C,4)*Power(MW,2) - 4*Power(M2C,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C13*g2*N24) - (gt*(C23*C31*N34 + C21*C33*N34 + C21*C23*N44 + C31*C33*N44))/Sqrt(2))*(2*C11*C13*g2*N24 + Sqrt(2)*gt*(C23*C31*N34 + C21*C33*N34 + C21*C23*N44 + C31*C33*N44))*UnitStep(Sqrt(Power(M2N,2)) - Sqrt(Power(M2C,2)) - Sqrt(Power(MW,2))))/(384.*Power(M2N,2)*Power(Power(M2N,2),1.5)*Power(M2C,2)*Power(MW,2)*Pi);
    cw2N += -(Sqrt(-4*Power(M1C,2)*Power(M2N,2) + Power(M2N,4))*(48*Power(M1C,6) + 68*Power(M1C,4)*Power(M2N,2) - 16*Power(M1C,2)*Power(M2N,4) - Power(M2N,6))*Power(-(Power(C12,2)*g2*N24) - (gt*(2*C22*C32*N34 + Power(C22,2)*N44 + Power(C32,2)*N44))/Sqrt(2),2)*UnitStep(-2*Sqrt(Power(M1C,2)) + Sqrt(Power(M2N,2))))/(192.*Power(M1C,4)*Power(Power(M2N,2),1.5)*Pi);
    cw2N += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(M1C,2)*Power(MW,2) - 2*Power(M2N,2)*Power(MW,2) + Power(MW,4))*(Power(M1C,8) + 8*Power(M1C,6)*(Power(M2N,2) + Power(MW,2)) + Power(Power(M2N,2) - Power(MW,2),2)*(Power(M2N,4) + 10*Power(M2N,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M1C,4)*(9*Power(M2N,4) + 16*Power(M2N,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M1C,2)*(Power(M2N,6) - 4*Power(M2N,4)*Power(MW,2) - 4*Power(M2N,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C12*g2*N24) - (gt*(C22*C31*N34 + C21*C32*N34 + C21*C22*N44 + C31*C32*N44))/Sqrt(2))*(2*C11*C12*g2*N24 + Sqrt(2)*gt*(C22*C31*N34 + C21*C32*N34 + C21*C22*N44 + C31*C32*N44))*UnitStep(-Sqrt(Power(M1C,2)) + Sqrt(Power(M2N,2)) - Sqrt(Power(MW,2))))/(384.*Power(M1C,2)*Power(M2N,2)*Power(Power(M2N,2),1.5)*Power(MW,2)*Pi);
    cw2N += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(M1C,2)*Power(MW,2) - 2*Power(M2N,2)*Power(MW,2) + Power(MW,4))*(Power(M1C,8) + 8*Power(M1C,6)*(Power(M2N,2) + Power(MW,2)) + Power(Power(M2N,2) - Power(MW,2),2)*(Power(M2N,4) + 10*Power(M2N,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M1C,4)*(9*Power(M2N,4) + 16*Power(M2N,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M1C,2)*(Power(M2N,6) - 4*Power(M2N,4)*Power(MW,2) - 4*Power(M2N,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C12*g2*N24) - (gt*(C22*C31*N34 + C21*C32*N34 + C21*C22*N44 + C31*C32*N44))/Sqrt(2))*(2*C11*C12*g2*N24 + Sqrt(2)*gt*(C22*C31*N34 + C21*C32*N34 + C21*C22*N44 + C31*C32*N44))*UnitStep(-Sqrt(Power(M1C,2)) + Sqrt(Power(M2N,2)) - Sqrt(Power(MW,2))))/(384.*Power(M1C,2)*Power(M2N,2)*Power(Power(M2N,2),1.5)*Power(MW,2)*Pi);
    cw2N += (Sqrt(Power(M2N,4) - 4*Power(M2N,2)*Power(MW,2))*(Power(M2N,6) + 16*Power(M2N,4)*Power(MW,2) - 68*Power(M2N,2)*Power(MW,4) - 48*Power(MW,6))*Power(-(Power(C11,2)*g2*N24) - (gt*(2*C21*C31*N34 + Power(C21,2)*N44 + Power(C31,2)*N44))/Sqrt(2),2)*UnitStep(Sqrt(Power(M2N,2)) - 2*Sqrt(Power(MW,2))))/(192.*Power(Power(M2N,2),1.5)*Power(MW,4)*Pi);
    cw2N += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1N,2) + Power(M1N,4) - 2*Power(HM,2)*Power(M2N,2) - 2*Power(M1N,2)*Power(M2N,2) + Power(M2N,4))*(Power(HM,4) + Power(M1N,4) + 10*Power(M1N,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(HM,2)*(Power(M1N,2) + Power(M2N,2)))*Power(2*Power(g1,2)*N13*N14*(1 - r3 + rs) + 2*Power(g2,2)*N23*N24*(1 - r3 + rs) + 2*Power(gt,2)*(N43*N44*(-r2 + rs) + N33*N34*(r2 + rs)) - Sqrt(2)*g2*gt*(N24*(N43*(-r2 + rs) + N33*(r2 - r3 + rs)) + N23*(N44*(-r2 + rs) + N34*(r2 - r3 + rs))) + g1*(-2*g2*(N14*N23 + N13*N24)*(1 + r2 - r3) + Sqrt(2)*gt*(N14*(N43*(r2 - rs) + N33*(r2 - r3 + rs)) + N13*(N44*(r2 - rs) + N34*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(-Sqrt(Power(HM,2)) - Sqrt(Power(M1N,2)) + Sqrt(Power(M2N,2))))/(3072.*Power(M1N,2)*Power(M2N,2)*Power(Power(M2N,2),1.5)*Pi);
    cw2N += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(HM,2)*Power(ZM,2) - 2*Power(M2N,2)*Power(ZM,2) + Power(ZM,4))*(Power(HM,4) + Power(M2N,4) + 10*Power(M2N,2)*Power(ZM,2) + Power(ZM,4) - 2*Power(HM,2)*(Power(M2N,2) + Power(ZM,2)))*Power(2*Power(g1,2)*N12*N14*(1 - r3 + rs) + 2*Power(g2,2)*N22*N24*(1 - r3 + rs) + 2*Power(gt,2)*(N42*N44*(-r2 + rs) + N32*N34*(r2 + rs)) - Sqrt(2)*g2*gt*(N24*(N42*(-r2 + rs) + N32*(r2 - r3 + rs)) + N22*(N44*(-r2 + rs) + N34*(r2 - r3 + rs))) + g1*(-2*g2*(N14*N22 + N12*N24)*(1 + r2 - r3) + Sqrt(2)*gt*(N14*(N42*(r2 - rs) + N32*(r2 - r3 + rs)) + N12*(N44*(r2 - rs) + N34*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(-Sqrt(Power(HM,2)) + Sqrt(Power(M2N,2)) - Sqrt(Power(ZM,2))))/(3072.*Power(M2N,2)*Power(Power(M2N,2),1.5)*Power(ZM,2)*Pi);
    cw2N += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M2N,2) + Power(M2N,4))*Power(2*Power(g1,2)*N11*N14*(1 - r3 + rs) + 2*Power(g2,2)*N21*N24*(1 - r3 + rs) + 2*Power(gt,2)*(N41*N44*(-r2 + rs) + N31*N34*(r2 + rs)) - Sqrt(2)*g2*gt*(N24*(N41*(-r2 + rs) + N31*(r2 - r3 + rs)) + N21*(N44*(-r2 + rs) + N34*(r2 - r3 + rs))) + g1*(-2*g2*(N14*N21 + N11*N24)*(1 + r2 - r3) + Sqrt(2)*gt*(N14*(N41*(r2 - rs) + N31*(r2 - r3 + rs)) + N11*(N44*(r2 - rs) + N34*(r2 - r3 + rs)))),2)*Power(v,2)*UnitStep(-Sqrt(Power(HM,2)) + Sqrt(Power(M2N,2))))/(256.*Power(Power(M2N,2),1.5)*Pi);
	w2N=creal(cw2N);    
    _warn_if_imag(cw2N);
    
    cw1C = -(Power(C12,2)*Power(g2,2)*Sqrt(Power(BM,4) - 2*Power(BM,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(BM,2)*Power(TM,2) - 2*Power(M1C,2)*Power(TM,2) + Power(TM,4))*(Power(BM,4) - 2*Power(M1C,4) + Power(M1C,2)*Power(TM,2) + Power(TM,4) + Power(BM,2)*(Power(M1C,2) - 2*Power(TM,2)))*UnitStep(-Sqrt(Power(BM,2)) + Sqrt(Power(M1C,2)) - Sqrt(Power(TM,2))))/(32.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Pi);
    cw1C += -(Power(C12,2)*Power(CKM22,2)*Power(g2,2)*(Power(CM,4) + Power(CM,2)*Power(M1C,2) - 2*Power(M1C,4))*Sqrt(Power(CM,4) - 2*Power(CM,2)*Power(M1C,2) + Power(M1C,4))*UnitStep(-Sqrt(Power(CM,2)) + Sqrt(Power(M1C,2))))/(32.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Pi);
    cw1C += (Power(C12,2)*CKM12*CKM21*Power(g2,2)*Sqrt(Power(M1C,4))*UnitStep(Sqrt(Power(M1C,2))))/(16.*Sqrt(Power(M1C,2))*Pi);
    cw1C += -(Power(C12,2)*CKM12*CKM21*Power(g2,2)*(Power(CM,4) + Power(CM,2)*Power(M1C,2) - 2*Power(M1C,4))*Sqrt(Power(CM,4) - 2*Power(CM,2)*Power(M1C,2) + Power(M1C,4))*UnitStep(-Sqrt(Power(CM,2)) + Sqrt(Power(M1C,2))))/(32.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Pi);
    cw1C += (Power(C12,2)*Power(CKM11,2)*Power(g2,2)*Sqrt(Power(M1C,4))*UnitStep(Sqrt(Power(M1C,2))))/(16.*Sqrt(Power(M1C,2))*Pi);
    cw1C += (Power(C12,2)*Power(g2,2)*(2*Power(M1C,4) - Power(M1C,2)*Power(TAM,2) - Power(TAM,4))*Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(TAM,2) + Power(TAM,4))*UnitStep(Sqrt(Power(M1C,2)) - Sqrt(Power(TAM,2))))/(96.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Pi);
    cw1C += (Power(C12,2)*Power(g2,2)*Sqrt(Power(M1C,4))*UnitStep(Sqrt(Power(M1C,2))))/(48.*Sqrt(Power(M1C,2))*Pi);
    cw1C += (Power(C12,2)*Power(g2,2)*Sqrt(Power(M1C,4))*UnitStep(Sqrt(Power(M1C,2))))/(48.*Sqrt(Power(M1C,2))*Pi);
    cw1C += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(M1C,2)*Power(M2C,2) - 2*Power(M2N,2)*Power(M2C,2) + Power(M2C,4))*(Power(M1C,8) + 8*Power(M1C,6)*(Power(M2N,2) + Power(M2C,2)) + Power(Power(M2N,2) - Power(M2C,2),2)*(Power(M2N,4) + 10*Power(M2N,2)*Power(M2C,2) + Power(M2C,4)) - 2*Power(M1C,4)*(9*Power(M2N,4) + 16*Power(M2N,2)*Power(M2C,2) + 9*Power(M2C,4)) + 8*Power(M1C,2)*(Power(M2N,6) - 4*Power(M2N,4)*Power(M2C,2) - 4*Power(M2N,2)*Power(M2C,4) + Power(M2C,6)))*(-(C12*C13*g2*N24) - (gt*(C23*C32*N34 + C22*C33*N34 + C22*C23*N44 + C32*C33*N44))/Sqrt(2))*(2*C12*C13*g2*N24 + Sqrt(2)*gt*(C23*C32*N34 + C22*C33*N34 + C22*C23*N44 + C32*C33*N44))*UnitStep(Sqrt(Power(M1C,2)) - Sqrt(Power(M2N,2)) - Sqrt(Power(M2C,2))))/(384.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Power(M2N,2)*Power(M2C,2)*Pi);
    cw1C += -(Sqrt(Power(M1N,4) - 2*Power(M1N,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(M1N,2)*Power(M2C,2) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4))*(Power(M1N,8) + 8*Power(M1N,6)*(Power(M1C,2) + Power(M2C,2)) + Power(Power(M1C,2) - Power(M2C,2),2)*(Power(M1C,4) + 10*Power(M1C,2)*Power(M2C,2) + Power(M2C,4)) - 2*Power(M1N,4)*(9*Power(M1C,4) + 16*Power(M1C,2)*Power(M2C,2) + 9*Power(M2C,4)) + 8*Power(M1N,2)*(Power(M1C,6) - 4*Power(M1C,4)*Power(M2C,2) - 4*Power(M1C,2)*Power(M2C,4) + Power(M2C,6)))*(-(C12*C13*g2*N23) - (gt*(C23*C32*N33 + C22*C33*N33 + C22*C23*N43 + C32*C33*N43))/Sqrt(2))*(2*C12*C13*g2*N23 + Sqrt(2)*gt*(C23*C32*N33 + C22*C33*N33 + C22*C23*N43 + C32*C33*N43))*UnitStep(-Sqrt(Power(M1N,2)) + Sqrt(Power(M1C,2)) - Sqrt(Power(M2C,2))))/(384.*Power(M1N,2)*Power(M1C,2)*Power(Power(M1C,2),1.5)*Power(M2C,2)*Pi);
    cw1C += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(M1C,2)*Power(ZM,2) - 2*Power(M2C,2)*Power(ZM,2) + Power(ZM,4))*(Power(M1C,8) + 8*Power(M1C,6)*(Power(M2C,2) + Power(ZM,2)) + Power(Power(M2C,2) - Power(ZM,2),2)*(Power(M2C,4) + 10*Power(M2C,2)*Power(ZM,2) + Power(ZM,4)) - 2*Power(M1C,4)*(9*Power(M2C,4) + 16*Power(M2C,2)*Power(ZM,2) + 9*Power(ZM,4)) + 8*Power(M1C,2)*(Power(M2C,6) - 4*Power(M2C,4)*Power(ZM,2) - 4*Power(M2C,2)*Power(ZM,4) + Power(ZM,6)))*(-(C12*C13*g2*N22) - (gt*(C23*C32*N32 + C22*C33*N32 + C22*C23*N42 + C32*C33*N42))/Sqrt(2))*(2*C12*C13*g2*N22 + Sqrt(2)*gt*(C23*C32*N32 + C22*C33*N32 + C22*C23*N42 + C32*C33*N42))*UnitStep(Sqrt(Power(M1C,2)) - Sqrt(Power(M2C,2)) - Sqrt(Power(ZM,2))))/(384.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Power(M2C,2)*Power(ZM,2)*Pi);
    cw1C += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4))*(5*Power(M1C,6) - 17*Power(M1C,4)*Power(M2C,2) - 17*Power(M1C,2)*Power(M2C,4) + 5*Power(M2C,6))*(-(C12*C13*g2*N21) - (gt*(C23*C32*N31 + C22*C33*N31 + C22*C23*N41 + C32*C33*N41))/Sqrt(2))*(2*C12*C13*g2*N21 + Sqrt(2)*gt*(C23*C32*N31 + C22*C33*N31 + C22*C23*N41 + C32*C33*N41))*UnitStep(Sqrt(Power(M1C,2)) - Sqrt(Power(M2C,2))))/(192.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Power(M2C,2)*Pi);
    cw1C += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(HM,2)*Power(M2C,2) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4))*(Power(HM,4) + Power(M1C,4) + 10*Power(M1C,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(HM,2)*(Power(M1C,2) + Power(M2C,2)))*Power(gt*(-2*gt*(C32*C33*(-r2 + rs) + C22*C23*(r2 + rs)) + Sqrt(2)*C13*g2*(C32*(-r2 + rs) + C22*(r2 - r3 + rs))) + C12*g2*(2*C13*g2*(-1 + r3 - rs) + Sqrt(2)*gt*(C33*(-r2 + rs) + C23*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(-Sqrt(Power(HM,2)) + Sqrt(Power(M1C,2)) - Sqrt(Power(M2C,2))))/(3072.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Power(M2C,2)*Pi);
    cw1C += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(M1C,2)*Power(MW,2) - 2*Power(M2N,2)*Power(MW,2) + Power(MW,4))*(Power(M1C,8) + 8*Power(M1C,6)*(Power(M2N,2) + Power(MW,2)) + Power(Power(M2N,2) - Power(MW,2),2)*(Power(M2N,4) + 10*Power(M2N,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M1C,4)*(9*Power(M2N,4) + 16*Power(M2N,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M1C,2)*(Power(M2N,6) - 4*Power(M2N,4)*Power(MW,2) - 4*Power(M2N,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C12*g2*N24) - (gt*(C22*C31*N34 + C21*C32*N34 + C21*C22*N44 + C31*C32*N44))/Sqrt(2))*(2*C11*C12*g2*N24 + Sqrt(2)*gt*(C22*C31*N34 + C21*C32*N34 + C21*C22*N44 + C31*C32*N44))*UnitStep(Sqrt(Power(M1C,2)) - Sqrt(Power(M2N,2)) - Sqrt(Power(MW,2))))/(384.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Power(M2N,2)*Power(MW,2)*Pi);
    cw1C += -(Sqrt(Power(M1N,4) - 2*Power(M1N,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(M1N,2)*Power(MW,2) - 2*Power(M1C,2)*Power(MW,2) + Power(MW,4))*(Power(M1N,8) + 8*Power(M1N,6)*(Power(M1C,2) + Power(MW,2)) + Power(Power(M1C,2) - Power(MW,2),2)*(Power(M1C,4) + 10*Power(M1C,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M1N,4)*(9*Power(M1C,4) + 16*Power(M1C,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M1N,2)*(Power(M1C,6) - 4*Power(M1C,4)*Power(MW,2) - 4*Power(M1C,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C12*g2*N23) - (gt*(C22*C31*N33 + C21*C32*N33 + C21*C22*N43 + C31*C32*N43))/Sqrt(2))*(2*C11*C12*g2*N23 + Sqrt(2)*gt*(C22*C31*N33 + C21*C32*N33 + C21*C22*N43 + C31*C32*N43))*UnitStep(-Sqrt(Power(M1N,2)) + Sqrt(Power(M1C,2)) - Sqrt(Power(MW,2))))/(384.*Power(M1N,2)*Power(M1C,2)*Power(Power(M1C,2),1.5)*Power(MW,2)*Pi);
    cw1C += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(MW,2) + Power(MW,4) - 2*Power(M1C,2)*Power(ZM,2) - 2*Power(MW,2)*Power(ZM,2) + Power(ZM,4))*(Power(M1C,8) + 8*Power(M1C,6)*(Power(MW,2) + Power(ZM,2)) + Power(Power(MW,2) - Power(ZM,2),2)*(Power(MW,4) + 10*Power(MW,2)*Power(ZM,2) + Power(ZM,4)) - 2*Power(M1C,4)*(9*Power(MW,4) + 16*Power(MW,2)*Power(ZM,2) + 9*Power(ZM,4)) + 8*Power(M1C,2)*(Power(MW,6) - 4*Power(MW,4)*Power(ZM,2) - 4*Power(MW,2)*Power(ZM,4) + Power(ZM,6)))*(-(C11*C12*g2*N22) - (gt*(C22*C31*N32 + C21*C32*N32 + C21*C22*N42 + C31*C32*N42))/Sqrt(2))*(2*C11*C12*g2*N22 + Sqrt(2)*gt*(C22*C31*N32 + C21*C32*N32 + C21*C22*N42 + C31*C32*N42))*UnitStep(Sqrt(Power(M1C,2)) - Sqrt(Power(MW,2)) - Sqrt(Power(ZM,2))))/(384.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Power(MW,2)*Power(ZM,2)*Pi);
    cw1C += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(MW,2) + Power(MW,4))*(5*Power(M1C,6) - 17*Power(M1C,4)*Power(MW,2) - 17*Power(M1C,2)*Power(MW,4) + 5*Power(MW,6))*(-(C11*C12*g2*N21) - (gt*(C22*C31*N31 + C21*C32*N31 + C21*C22*N41 + C31*C32*N41))/Sqrt(2))*(2*C11*C12*g2*N21 + Sqrt(2)*gt*(C22*C31*N31 + C21*C32*N31 + C21*C22*N41 + C31*C32*N41))*UnitStep(Sqrt(Power(M1C,2)) - Sqrt(Power(MW,2))))/(192.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Power(MW,2)*Pi);
    cw1C += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(HM,2)*Power(MW,2) - 2*Power(M1C,2)*Power(MW,2) + Power(MW,4))*(Power(HM,4) + Power(M1C,4) + 10*Power(M1C,2)*Power(MW,2) + Power(MW,4) - 2*Power(HM,2)*(Power(M1C,2) + Power(MW,2)))*Power(gt*(-2*gt*(C31*C32*(-r2 + rs) + C21*C22*(r2 + rs)) + Sqrt(2)*C12*g2*(C31*(-r2 + rs) + C21*(r2 - r3 + rs))) + C11*g2*(2*C12*g2*(-1 + r3 - rs) + Sqrt(2)*gt*(C32*(-r2 + rs) + C22*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(-Sqrt(Power(HM,2)) + Sqrt(Power(M1C,2)) - Sqrt(Power(MW,2))))/(3072.*Power(M1C,2)*Power(Power(M1C,2),1.5)*Power(MW,2)*Pi);
	w1C=creal(cw1C);    
    _warn_if_imag(cw1C);
    
    cw2C = -(Power(C13,2)*Power(g2,2)*Sqrt(Power(BM,4) - 2*Power(BM,2)*Power(TM,2) + Power(TM,4) - 2*Power(BM,2)*Power(M2C,2) - 2*Power(TM,2)*Power(M2C,2) + Power(M2C,4))*(Power(BM,4) + Power(TM,4) + Power(TM,2)*Power(M2C,2) - 2*Power(M2C,4) + Power(BM,2)*(-2*Power(TM,2) + Power(M2C,2)))*UnitStep(-Sqrt(Power(BM,2)) - Sqrt(Power(TM,2)) + Sqrt(Power(M2C,2))))/(32.*Power(M2C,2)*Power(Power(M2C,2),1.5)*Pi);
    cw2C += -(Power(C13,2)*Power(CKM22,2)*Power(g2,2)*(Power(CM,4) + Power(CM,2)*Power(M2C,2) - 2*Power(M2C,4))*Sqrt(Power(CM,4) - 2*Power(CM,2)*Power(M2C,2) + Power(M2C,4))*UnitStep(-Sqrt(Power(CM,2)) + Sqrt(Power(M2C,2))))/(32.*Power(M2C,2)*Power(Power(M2C,2),1.5)*Pi);
    cw2C += (Power(C13,2)*CKM12*CKM21*Power(g2,2)*Sqrt(Power(M2C,4))*UnitStep(Sqrt(Power(M2C,2))))/(16.*Sqrt(Power(M2C,2))*Pi);
    cw2C += -(Power(C13,2)*CKM12*CKM21*Power(g2,2)*(Power(CM,4) + Power(CM,2)*Power(M2C,2) - 2*Power(M2C,4))*Sqrt(Power(CM,4) - 2*Power(CM,2)*Power(M2C,2) + Power(M2C,4))*UnitStep(-Sqrt(Power(CM,2)) + Sqrt(Power(M2C,2))))/(32.*Power(M2C,2)*Power(Power(M2C,2),1.5)*Pi);
    cw2C += (Power(C13,2)*Power(CKM11,2)*Power(g2,2)*Sqrt(Power(M2C,4))*UnitStep(Sqrt(Power(M2C,2))))/(16.*Sqrt(Power(M2C,2))*Pi);
    cw2C += -(Power(C13,2)*Power(g2,2)*(Power(TAM,4) + Power(TAM,2)*Power(M2C,2) - 2*Power(M2C,4))*Sqrt(Power(TAM,4) - 2*Power(TAM,2)*Power(M2C,2) + Power(M2C,4))*UnitStep(-Sqrt(Power(TAM,2)) + Sqrt(Power(M2C,2))))/(96.*Power(M2C,2)*Power(Power(M2C,2),1.5)*Pi);
    cw2C += (Power(C13,2)*Power(g2,2)*Sqrt(Power(M2C,4))*UnitStep(Sqrt(Power(M2C,2))))/(48.*Sqrt(Power(M2C,2))*Pi);
    cw2C += (Power(C13,2)*Power(g2,2)*Sqrt(Power(M2C,4))*UnitStep(Sqrt(Power(M2C,2))))/(48.*Sqrt(Power(M2C,2))*Pi);
    cw2C += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(M2N,2) + Power(M2N,4) - 2*Power(M1C,2)*Power(M2C,2) - 2*Power(M2N,2)*Power(M2C,2) + Power(M2C,4))*(Power(M1C,8) + 8*Power(M1C,6)*(Power(M2N,2) + Power(M2C,2)) + Power(Power(M2N,2) - Power(M2C,2),2)*(Power(M2N,4) + 10*Power(M2N,2)*Power(M2C,2) + Power(M2C,4)) - 2*Power(M1C,4)*(9*Power(M2N,4) + 16*Power(M2N,2)*Power(M2C,2) + 9*Power(M2C,4)) + 8*Power(M1C,2)*(Power(M2N,6) - 4*Power(M2N,4)*Power(M2C,2) - 4*Power(M2N,2)*Power(M2C,4) + Power(M2C,6)))*(-(C12*C13*g2*N24) - (gt*(C23*C32*N34 + C22*C33*N34 + C22*C23*N44 + C32*C33*N44))/Sqrt(2))*(2*C12*C13*g2*N24 + Sqrt(2)*gt*(C23*C32*N34 + C22*C33*N34 + C22*C23*N44 + C32*C33*N44))*UnitStep(-Sqrt(Power(M1C,2)) - Sqrt(Power(M2N,2)) + Sqrt(Power(M2C,2))))/(384.*Power(M1C,2)*Power(M2N,2)*Power(M2C,2)*Power(Power(M2C,2),1.5)*Pi);
    cw2C += -(Sqrt(Power(M1N,4) - 2*Power(M1N,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(M1N,2)*Power(M2C,2) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4))*(Power(M1N,8) + 8*Power(M1N,6)*(Power(M1C,2) + Power(M2C,2)) + Power(Power(M1C,2) - Power(M2C,2),2)*(Power(M1C,4) + 10*Power(M1C,2)*Power(M2C,2) + Power(M2C,4)) - 2*Power(M1N,4)*(9*Power(M1C,4) + 16*Power(M1C,2)*Power(M2C,2) + 9*Power(M2C,4)) + 8*Power(M1N,2)*(Power(M1C,6) - 4*Power(M1C,4)*Power(M2C,2) - 4*Power(M1C,2)*Power(M2C,4) + Power(M2C,6)))*(-(C12*C13*g2*N23) - (gt*(C23*C32*N33 + C22*C33*N33 + C22*C23*N43 + C32*C33*N43))/Sqrt(2))*(2*C12*C13*g2*N23 + Sqrt(2)*gt*(C23*C32*N33 + C22*C33*N33 + C22*C23*N43 + C32*C33*N43))*UnitStep(-Sqrt(Power(M1N,2)) - Sqrt(Power(M1C,2)) + Sqrt(Power(M2C,2))))/(384.*Power(M1N,2)*Power(M1C,2)*Power(M2C,2)*Power(Power(M2C,2),1.5)*Pi);
    cw2C += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(M1C,2)*Power(ZM,2) - 2*Power(M2C,2)*Power(ZM,2) + Power(ZM,4))*(Power(M1C,8) + 8*Power(M1C,6)*(Power(M2C,2) + Power(ZM,2)) + Power(Power(M2C,2) - Power(ZM,2),2)*(Power(M2C,4) + 10*Power(M2C,2)*Power(ZM,2) + Power(ZM,4)) - 2*Power(M1C,4)*(9*Power(M2C,4) + 16*Power(M2C,2)*Power(ZM,2) + 9*Power(ZM,4)) + 8*Power(M1C,2)*(Power(M2C,6) - 4*Power(M2C,4)*Power(ZM,2) - 4*Power(M2C,2)*Power(ZM,4) + Power(ZM,6)))*(-(C12*C13*g2*N22) - (gt*(C23*C32*N32 + C22*C33*N32 + C22*C23*N42 + C32*C33*N42))/Sqrt(2))*(2*C12*C13*g2*N22 + Sqrt(2)*gt*(C23*C32*N32 + C22*C33*N32 + C22*C23*N42 + C32*C33*N42))*UnitStep(-Sqrt(Power(M1C,2)) + Sqrt(Power(M2C,2)) - Sqrt(Power(ZM,2))))/(384.*Power(M1C,2)*Power(M2C,2)*Power(Power(M2C,2),1.5)*Power(ZM,2)*Pi);
    cw2C += -(Sqrt(Power(M1C,4) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4))*(5*Power(M1C,6) - 17*Power(M1C,4)*Power(M2C,2) - 17*Power(M1C,2)*Power(M2C,4) + 5*Power(M2C,6))*(-(C12*C13*g2*N21) - (gt*(C23*C32*N31 + C22*C33*N31 + C22*C23*N41 + C32*C33*N41))/Sqrt(2))*(2*C12*C13*g2*N21 + Sqrt(2)*gt*(C23*C32*N31 + C22*C33*N31 + C22*C23*N41 + C32*C33*N41))*UnitStep(-Sqrt(Power(M1C,2)) + Sqrt(Power(M2C,2))))/(192.*Power(M1C,2)*Power(M2C,2)*Power(Power(M2C,2),1.5)*Pi);
    cw2C += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M1C,2) + Power(M1C,4) - 2*Power(HM,2)*Power(M2C,2) - 2*Power(M1C,2)*Power(M2C,2) + Power(M2C,4))*(Power(HM,4) + Power(M1C,4) + 10*Power(M1C,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(HM,2)*(Power(M1C,2) + Power(M2C,2)))*Power(gt*(-2*gt*(C32*C33*(-r2 + rs) + C22*C23*(r2 + rs)) + Sqrt(2)*C13*g2*(C32*(-r2 + rs) + C22*(r2 - r3 + rs))) + C12*g2*(2*C13*g2*(-1 + r3 - rs) + Sqrt(2)*gt*(C33*(-r2 + rs) + C23*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(-Sqrt(Power(HM,2)) - Sqrt(Power(M1C,2)) + Sqrt(Power(M2C,2))))/(3072.*Power(M1C,2)*Power(M2C,2)*Power(Power(M2C,2),1.5)*Pi);
    cw2C += -(Sqrt(Power(M2N,4) - 2*Power(M2N,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(M2N,2)*Power(MW,2) - 2*Power(M2C,2)*Power(MW,2) + Power(MW,4))*(Power(M2N,8) + 8*Power(M2N,6)*(Power(M2C,2) + Power(MW,2)) + Power(Power(M2C,2) - Power(MW,2),2)*(Power(M2C,4) + 10*Power(M2C,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M2N,4)*(9*Power(M2C,4) + 16*Power(M2C,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M2N,2)*(Power(M2C,6) - 4*Power(M2C,4)*Power(MW,2) - 4*Power(M2C,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C13*g2*N24) - (gt*(C23*C31*N34 + C21*C33*N34 + C21*C23*N44 + C31*C33*N44))/Sqrt(2))*(2*C11*C13*g2*N24 + Sqrt(2)*gt*(C23*C31*N34 + C21*C33*N34 + C21*C23*N44 + C31*C33*N44))*UnitStep(-Sqrt(Power(M2N,2)) + Sqrt(Power(M2C,2)) - Sqrt(Power(MW,2))))/(384.*Power(M2N,2)*Power(M2C,2)*Power(Power(M2C,2),1.5)*Power(MW,2)*Pi);
    cw2C += -(Sqrt(Power(M1N,4) - 2*Power(M1N,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(M1N,2)*Power(MW,2) - 2*Power(M2C,2)*Power(MW,2) + Power(MW,4))*(Power(M1N,8) + 8*Power(M1N,6)*(Power(M2C,2) + Power(MW,2)) + Power(Power(M2C,2) - Power(MW,2),2)*(Power(M2C,4) + 10*Power(M2C,2)*Power(MW,2) + Power(MW,4)) - 2*Power(M1N,4)*(9*Power(M2C,4) + 16*Power(M2C,2)*Power(MW,2) + 9*Power(MW,4)) + 8*Power(M1N,2)*(Power(M2C,6) - 4*Power(M2C,4)*Power(MW,2) - 4*Power(M2C,2)*Power(MW,4) + Power(MW,6)))*(-(C11*C13*g2*N23) - (gt*(C23*C31*N33 + C21*C33*N33 + C21*C23*N43 + C31*C33*N43))/Sqrt(2))*(2*C11*C13*g2*N23 + Sqrt(2)*gt*(C23*C31*N33 + C21*C33*N33 + C21*C23*N43 + C31*C33*N43))*UnitStep(-Sqrt(Power(M1N,2)) + Sqrt(Power(M2C,2)) - Sqrt(Power(MW,2))))/(384.*Power(M1N,2)*Power(M2C,2)*Power(Power(M2C,2),1.5)*Power(MW,2)*Pi);
    cw2C += -(Sqrt(Power(M2C,4) - 2*Power(M2C,2)*Power(MW,2) + Power(MW,4) - 2*Power(M2C,2)*Power(ZM,2) - 2*Power(MW,2)*Power(ZM,2) + Power(ZM,4))*(Power(M2C,8) + 8*Power(M2C,6)*(Power(MW,2) + Power(ZM,2)) + Power(Power(MW,2) - Power(ZM,2),2)*(Power(MW,4) + 10*Power(MW,2)*Power(ZM,2) + Power(ZM,4)) - 2*Power(M2C,4)*(9*Power(MW,4) + 16*Power(MW,2)*Power(ZM,2) + 9*Power(ZM,4)) + 8*Power(M2C,2)*(Power(MW,6) - 4*Power(MW,4)*Power(ZM,2) - 4*Power(MW,2)*Power(ZM,4) + Power(ZM,6)))*(-(C11*C13*g2*N22) - (gt*(C23*C31*N32 + C21*C33*N32 + C21*C23*N42 + C31*C33*N42))/Sqrt(2))*(2*C11*C13*g2*N22 + Sqrt(2)*gt*(C23*C31*N32 + C21*C33*N32 + C21*C23*N42 + C31*C33*N42))*UnitStep(Sqrt(Power(M2C,2)) - Sqrt(Power(MW,2)) - Sqrt(Power(ZM,2))))/(384.*Power(M2C,2)*Power(Power(M2C,2),1.5)*Power(MW,2)*Power(ZM,2)*Pi);
    cw2C += -(Sqrt(Power(M2C,4) - 2*Power(M2C,2)*Power(MW,2) + Power(MW,4))*(5*Power(M2C,6) - 17*Power(M2C,4)*Power(MW,2) - 17*Power(M2C,2)*Power(MW,4) + 5*Power(MW,6))*(-(C11*C13*g2*N21) - (gt*(C23*C31*N31 + C21*C33*N31 + C21*C23*N41 + C31*C33*N41))/Sqrt(2))*(2*C11*C13*g2*N21 + Sqrt(2)*gt*(C23*C31*N31 + C21*C33*N31 + C21*C23*N41 + C31*C33*N41))*UnitStep(Sqrt(Power(M2C,2)) - Sqrt(Power(MW,2))))/(192.*Power(M2C,2)*Power(Power(M2C,2),1.5)*Power(MW,2)*Pi);
    cw2C += (Sqrt(Power(HM,4) - 2*Power(HM,2)*Power(M2C,2) + Power(M2C,4) - 2*Power(HM,2)*Power(MW,2) - 2*Power(M2C,2)*Power(MW,2) + Power(MW,4))*(Power(HM,4) + Power(M2C,4) + 10*Power(M2C,2)*Power(MW,2) + Power(MW,4) - 2*Power(HM,2)*(Power(M2C,2) + Power(MW,2)))*Power(gt*(-2*gt*(C31*C33*(-r2 + rs) + C21*C23*(r2 + rs)) + Sqrt(2)*C13*g2*(C31*(-r2 + rs) + C21*(r2 - r3 + rs))) + C11*g2*(2*C13*g2*(-1 + r3 - rs) + Sqrt(2)*gt*(C33*(-r2 + rs) + C23*(r2 - r3 + rs))),2)*Power(v,2)*UnitStep(-Sqrt(Power(HM,2)) + Sqrt(Power(M2C,2)) - Sqrt(Power(MW,2))))/(3072.*Power(M2C,2)*Power(Power(M2C,2),1.5)*Power(MW,2)*Pi);
	w2C=creal(cw2C);    
    _warn_if_imag(cw2C);
    
#undef _warn_if_imag
}


void calculate() {
    
    /* Masses */
    HM = MH;
    ZM = MZ;
    CM = MC;
    BM = MB;
    TM = MT;
    TAM = MTA;

    /* "Strong coupling constant" */
    gs=sqrt(4*Pi*aS);

    /* "SM VEV" */
    v=(0.282094791773878*sqrt(pow(gt,2)*(8.88576587631673 - 1.*GF*pow(MA,2)*PS) + GF*pow(MA,2)*(50.265482457437 - 10.026513098524*sqrt(25.1327412287183 - 1.*pow(gt,2)*PS))))/sqrt(GF*pow(gt,2));
    
    /* "C-M parameter" */
    r3=(-4.*GF*pow(MA,2)*(-12.5663706143592 + 2.506628274631*sqrt(25.1327412287183 - 1.*pow(gt,2)*PS)))/(pow(gt,2)*(8.88576587631673 - 1.*GF*pow(MA,2)*PS) + GF*pow(MA,2)*(50.2654824574367 - 10.026513098524*sqrt(25.1327412287183 - 1.*pow(gt,2)*PS)));
    
    /* "C-M parameter" */
    r2=-1. + r3;
    
    /* "Vector meson mass scale" */
    f=0.282094791773878*sqrt(8.88576587631673/GF + pow(MA,2)*(50.2654824574367/pow(gt,2) - 1.*PS));
    
    /* "Strong vector meson mass" */
    MV=sqrt((0.353553390593274*pow(gt,2))/GF + pow(MA,2)*(1. - 0.0397887357729738*pow(gt,2)*PS));
    
    /* "Strong vector decay constant" */
    FV=(1.4142135623731*MV)/gt;
    
    /* "Strong pion decay constant" */
    FPi=0.840896415253715/sqrt(GF);
    
    /* "Strong axial decay constant" */
    FA=sqrt(-1.*pow(FPi,2) + pow(FV,2));
    
    /* "electroweak SU2L gauge coupling" */
    g2=1.4142135623731/sqrt(pow(EE,-2) - 2./pow(gt,2) + sqrt(pow(ZM,2)*(pow(MA,2) - 1.*pow(ZM,2))*(pow(MV,2) - 1.*pow(ZM,2))*(pow(FV,2) + (pow(EE,-2) - 2./pow(gt,2))*(pow(MV,2) - 1.*pow(ZM,2)))*(pow(FA,2)*pow(MA,2) + (pow(MA,2) - 1.*pow(ZM,2))*(-1.*pow(FV,2) + (pow(EE,-2) - 2./pow(gt,2))*pow(ZM,2))))/(pow(ZM,2)*(pow(MA,2) - 1.*pow(ZM,2))*(-1.*pow(MV,2) + pow(ZM,2))));
    
    /* "electroweak U1Y gauge coupling" */
    g1=1/sqrt(pow(EE,-2) - 1./pow(g2,2) - 2./pow(gt,2));
    
    /* "Neutral heavy vector meson R1 mass" */
    M1N=(0.5*sqrt((pow(FV,2)*(pow(g1,2) + pow(g2,2)) + 2.*(pow(MA,2) + pow(MV,2)))*pow(ZM,2) - 2.*pow(ZM,4) - 0.5*sqrt(-16.*pow(FPi,2)*pow(MA,2)*(pow(FV,2)*pow(g1,2)*pow(g2,2) + (pow(g1,2) + pow(g2,2))*pow(MV,2))*pow(ZM,2) + 4.*pow(ZM,4)*pow(pow(FV,2)*(pow(g1,2) + pow(g2,2)) + 2.*(pow(MA,2) + pow(MV,2) - 1.*pow(ZM,2)),2))))/ZM;
    
    /* "Neutral heavy vector meson R2 mass" */
    M2N=(0.5*sqrt((pow(FV,2)*(pow(g1,2) + pow(g2,2)) + 2.*(pow(MA,2) + pow(MV,2)))*pow(ZM,2) - 2.*pow(ZM,4) + 0.5*sqrt(-16.*pow(FPi,2)*pow(MA,2)*(pow(FV,2)*pow(g1,2)*pow(g2,2) + (pow(g1,2) + pow(g2,2))*pow(MV,2))*pow(ZM,2) + 4.*pow(ZM,4)*pow(pow(FV,2)*(pow(g1,2) + pow(g2,2)) + 2.*(pow(MA,2) + pow(MV,2) - 1.*pow(ZM,2)),2))))/ZM;
    
    /* "Charged vector meson mass angle" */
    ThetaC=acos((0.5*(2.*pow(FV,6)*pow(g2,6) + 3.*pow(FV,4)*pow(g2,4)*(pow(MA,2) + pow(MV,2)) + 2.*(pow(MA,2) - 2.*pow(MV,2))*(-9.*pow(FPi,2)*pow(g2,2)*pow(MA,2) + 8.*pow(MA,4) + 4.*pow(MA,2)*pow(MV,2) - 4.*pow(MV,4)) + 3.*pow(FV,2)*pow(g2,2)*(-3.*pow(FPi,2)*pow(g2,2)*pow(MA,2) + 2.*(pow(MA,4) - 4.*pow(MA,2)*pow(MV,2) + pow(MV,4)))))/pow(pow(FV,4)*pow(g2,4) - 3.*pow(FPi,2)*pow(g2,2)*pow(MA,2) + pow(FV,2)*pow(g2,2)*(pow(MA,2) + pow(MV,2)) + 4.*(pow(MA,4) - 1.*pow(MA,2)*pow(MV,2) + pow(MV,4)),1.5));
    
    /* "W mass" */
    MW=sqrt(0.333333333333333*(0.5*pow(FV,2)*pow(g2,2) + pow(MA,2) + pow(MV,2)) - 0.333333333333333*sqrt(pow(FV,4)*pow(g2,4) + 3.*pow(FA,2)*pow(g2,2)*pow(MA,2) + pow(FV,2)*pow(g2,2)*(-2.*pow(MA,2) + pow(MV,2)) + 4.*(pow(MA,4) - 1.*pow(MA,2)*pow(MV,2) + pow(MV,4)))*sin(0.523598775598299 + 0.333333333333333*ThetaC));
    
    /* "Charged heavy vector meson R1 mass" */
    M1C=sqrt(0.333333333333333*(0.5*pow(FV,2)*pow(g2,2) + pow(MA,2) + pow(MV,2)) - 0.333333333333333*sqrt(pow(FV,4)*pow(g2,4) + 3.*pow(FA,2)*pow(g2,2)*pow(MA,2) + pow(FV,2)*pow(g2,2)*(-2.*pow(MA,2) + pow(MV,2)) + 4.*(pow(MA,4) - 1.*pow(MA,2)*pow(MV,2) + pow(MV,4)))*sin(0.523598775598299 - 0.333333333333333*ThetaC));
    
    /* "Charged heavy vector meson R2 mass" */
    M2C=sqrt(0.333333333333333*(0.5*pow(FV,2)*pow(g2,2) + pow(MA,2) + pow(MV,2)) + 0.333333333333333*sqrt(pow(FV,4)*pow(g2,4) + 3.*pow(FA,2)*pow(g2,2)*pow(MA,2) + pow(FV,2)*pow(g2,2)*(-2.*pow(MA,2) + pow(MV,2)) + 4.*(pow(MA,4) - 1.*pow(MA,2)*pow(MV,2) + pow(MV,4)))*cos(0.333333333333333*ThetaC));
    
    /* "VC11" */
    VC11=0.5*pow(FV,2)*pow(g2,2);
    
    /* "VC12" */
    VC12=-0.5*FA*g2*MA;
    
    /* "VC13" */
    VC13=-0.5*FV*g2*MV;
    
    /* "VC21" */
    VC21=VC12;
    
    /* "VC22" */
    VC22=pow(MA,2);
    
    /* "VC23" */
    VC23=0;
    
    /* "VC31" */
    VC31=VC13;
    
    /* "VC32" */
    VC32=0;
    
    /* "VC33" */
    VC33=pow(MV,2);
    
    /* "VN11" */
    VN11=0.5*pow(FV,2)*pow(g1,2);
    
    /* "VN12" */
    VN12=0;
    
    /* "VN13" */
    VN13=0.5*FA*g1*MA;
    
    /* "VN14" */
    VN14=-0.5*FV*g1*MV;
    
    /* "VN21" */
    VN21=0;
    
    /* "VN22" */
    VN22=0.5*pow(FV,2)*pow(g2,2);
    
    /* "VN23" */
    VN23=-0.5*FA*g2*MA;
    
    /* "VN24" */
    VN24=-0.5*FV*g2*MV;
    
    /* "VN31" */
    VN31=VN13;
    
    /* "VN32" */
    VN32=VN23;
    
    /* "VN33" */
    VN33=pow(MA,2);
    
    /* "VN34" */
    VN34=0;
    
    /* "VN41" */
    VN41=VN14;
    
    /* "VN42" */
    VN42=VN24;
    
    /* "VN43" */
    VN43=0;
    
    /* "VN44" */
    VN44=pow(MV,2);
    
    /* "Charged vector meson normalization factor" */
    CN1=1/sqrt(pow(VC13,2)*pow(-1.*pow(MW,2) + VC22,2) + pow(VC12,2)*pow(pow(MW,2) - 1.*VC33,2) + pow(-1.*pow(MW,2) + VC22,2)*pow(pow(MW,2) - 1.*VC33,2));
    
    /* "Charged vector meson normalization factor" */
    CN2=1/sqrt(pow(VC13,2)*pow(-1.*pow(M1C,2) + VC22,2) + pow(VC12,2)*pow(pow(M1C,2) - 1.*VC33,2) + pow(-1.*pow(M1C,2) + VC22,2)*pow(pow(M1C,2) - 1.*VC33,2));
    
    /* "Charged vector meson normalization factor" */
    CN3=1/sqrt(pow(VC13,2)*pow(-1.*pow(M2C,2) + VC22,2) + pow(VC12,2)*pow(pow(M2C,2) - 1.*VC33,2) + pow(-1.*pow(M2C,2) + VC22,2)*pow(pow(M2C,2) - 1.*VC33,2));
    
    /* "Charged vector meson mixing matrix element in VA base" */
    C11=CN1*(-1.*pow(MW,2) + VC22)*(-1.*pow(MW,2) + VC33);
    
    /* "Charged vector meson mixing matrix element in VA base" */
    C12=CN2*(-1.*pow(M1C,2) + VC22)*(-1.*pow(M1C,2) + VC33);
    
    /* "Charged vector meson mixing matrix element in VA base" */
    C13=CN3*(-1.*pow(M2C,2) + VC22)*(-1.*pow(M2C,2) + VC33);
    
    /* "Charged vector meson mixing matrix element in VA base" */
    C21=CN1*VC12*(pow(MW,2) - 1.*VC33);
    
    /* "Charged vector meson mixing matrix element in VA base" */
    C22=CN2*VC12*(pow(M1C,2) - 1.*VC33);
    
    /* "Charged vector meson mixing matrix element in VA base" */
    C23=CN3*VC12*(pow(M2C,2) - 1.*VC33);
    
    /* "Charged vector meson mixing matrix element in VA base" */
    C31=CN1*VC13*(pow(MW,2) - 1.*VC22);
    
    /* "Charged vector meson mixing matrix element in VA base" */
    C32=CN2*VC13*(pow(M1C,2) - 1.*VC22);
    
    /* "Charged vector meson mixing matrix element in VA base" */
    C33=CN3*VC13*(pow(M2C,2) - 1.*VC22);
    
    /* "Neutral vector meson normalization factor" */
    NN2=1/sqrt(pow(g1 - 1.*g2,2)*pow(g1 + g2,2)*pow(gt,2)*pow(VN24,2)*pow(pow(ZM,2) - 1.*VN33,2) + pow(g1,2)*pow(g2,2)*pow(pow(ZM,2) - 1.*VN33,2)*pow(1.4142135623731*g2*VN24 + gt*(pow(ZM,2) - 1.*VN44),2) + pow(g2,2)*pow(pow(ZM,2) - 1.*VN33,2)*pow(1.4142135623731*pow(g1,2)*VN24 + g2*gt*(pow(ZM,2) - 1.*VN44),2) + pow(VN23,2)*pow(2.82842712474619*pow(g1,2)*g2*VN24 + (pow(g1,2) + pow(g2,2))*gt*(pow(ZM,2) - 1.*VN44),2));
    
    /* "Neutral vector meson normalization factor" */
    NN3=1/sqrt(pow(g1 - 1.*g2,2)*pow(g1 + g2,2)*pow(gt,2)*pow(VN24,2)*pow(pow(M1N,2) - 1.*VN33,2) + pow(g1,2)*pow(g2,2)*pow(pow(M1N,2) - 1.*VN33,2)*pow(1.4142135623731*g2*VN24 + gt*(pow(M1N,2) - 1.*VN44),2) + pow(g2,2)*pow(pow(M1N,2) - 1.*VN33,2)*pow(1.4142135623731*pow(g1,2)*VN24 + g2*gt*(pow(M1N,2) - 1.*VN44),2) + pow(VN23,2)*pow(2.82842712474619*pow(g1,2)*g2*VN24 + (pow(g1,2) + pow(g2,2))*gt*(pow(M1N,2) - 1.*VN44),2));
    
    /* "Neutral vector meson normalization factor" */
    NN4=1/sqrt(pow(g1 - 1.*g2,2)*pow(g1 + g2,2)*pow(gt,2)*pow(VN24,2)*pow(pow(M2N,2) - 1.*VN33,2) + pow(g1,2)*pow(g2,2)*pow(pow(M2N,2) - 1.*VN33,2)*pow(1.4142135623731*g2*VN24 + gt*(pow(M2N,2) - 1.*VN44),2) + pow(g2,2)*pow(pow(M2N,2) - 1.*VN33,2)*pow(1.4142135623731*pow(g1,2)*VN24 + g2*gt*(pow(M2N,2) - 1.*VN44),2) + pow(VN23,2)*pow(2.82842712474619*pow(g1,2)*g2*VN24 + (pow(g1,2) + pow(g2,2))*gt*(pow(M2N,2) - 1.*VN44),2));
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N11=EE/g1;
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N12=g1*g2*NN2*(-1.*pow(ZM,2) + VN33)*(1.4142135623731*g2*VN24 + gt*(pow(ZM,2) - 1.*VN44));
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N13=g1*g2*NN3*(-1.*pow(M1N,2) + VN33)*(1.4142135623731*g2*VN24 + gt*(pow(M1N,2) - 1.*VN44));
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N14=g1*g2*NN4*(-1.*pow(M2N,2) + VN33)*(1.4142135623731*g2*VN24 + gt*(pow(M2N,2) - 1.*VN44));
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N21=EE/g2;
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N22=-1.*g2*NN2*(-1.*pow(ZM,2) + VN33)*(1.4142135623731*pow(g1,2)*VN24 + g2*gt*(pow(ZM,2) - 1.*VN44));
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N23=-1.*g2*NN3*(-1.*pow(M1N,2) + VN33)*(1.4142135623731*pow(g1,2)*VN24 + g2*gt*(pow(M1N,2) - 1.*VN44));
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N24=-1.*g2*NN4*(-1.*pow(M2N,2) + VN33)*(1.4142135623731*pow(g1,2)*VN24 + g2*gt*(pow(M2N,2) - 1.*VN44));
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N31=0;
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N32=NN2*VN23*(2.82842712474619*pow(g1,2)*g2*VN24 + (pow(g1,2) + pow(g2,2))*gt*(pow(ZM,2) - 1.*VN44));
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N33=NN3*VN23*(2.82842712474619*pow(g1,2)*g2*VN24 + (pow(g1,2) + pow(g2,2))*gt*(pow(M1N,2) - 1.*VN44));
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N34=NN4*VN23*(2.82842712474619*pow(g1,2)*g2*VN24 + (pow(g1,2) + pow(g2,2))*gt*(pow(M2N,2) - 1.*VN44));
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N41=(1.4142135623731*EE)/gt;
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N42=(pow(g1,2) - 1.*pow(g2,2))*gt*NN2*VN24*(-1.*pow(ZM,2) + VN33);
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N43=(pow(g1,2) - 1.*pow(g2,2))*gt*NN3*VN24*(-1.*pow(M1N,2) + VN33);
    
    /* "Neutral vector meson mixing matrix element in VA base" */
    N44=(pow(g1,2) - 1.*pow(g2,2))*gt*NN4*VN24*(-1.*pow(M2N,2) + VN33);
    
    /* "Lepton Yukawa coupling ( ytau )" */
    ytau=(1.4142135623731*MTA)/v;
    
    /* "U-quark Yukawa coupling ( yc )" */
    yc=(1.4142135623731*MC)/v;
    
    /* "U-quark Yukawa coupling ( yt )" */
    yt=(1.4142135623731*MT)/v;
    
    /* "D-quark Yukawa coupling ( yb )" */
    yb=(1.4142135623731*MB)/v;
    
    /* "CKM-Matrix ( CKM11 )" */
    CKM11=cos(cabi);
    
    /* "CKM-Matrix ( CKM12 )" */
    CKM12=sin(cabi);
    
    /* "CKM-Matrix ( CKM21 )" */
    CKM21=-1.*sin(cabi);
    
    /* "CKM-Matrix ( CKM22 )" */
    CKM22=cos(cabi);
    
    calculate_widths();
        
}

void write_param_card() {

    printf("Block CKMBLOCK\n");
    printf("  1  %.13e # cabi\n",cabi);
    
    printf("Block SMInput\n");
    printf("  1  %.13e  # aS\n",aS);
    printf("  2  %.13e  # EE\n",EE);
    printf("  3  %.13e  # GF\n",GF);
    printf("  4  %.13e  # MZ\n",MZ);
    printf("Block TCInput\n");
    printf("  1  %.13e  # gt\n",gt);
    printf("  2  %.13e  # MA\n",MA);
    printf("  3  %.13e  # PS\n",PS);
    printf("  4  %.13e  # rs\n",rs);
    printf("  5  %.13e  # MH\n",MH);
    printf("Block YUKAWA\n");
    printf("  4  %.13e  # MC\n",MC);
    printf("  5  %.13e  # MB\n",MB);
    printf("  6  %.13e  # MT\n",MT);
    printf("  15  %.13e  # MTA\n",MTA);
    printf("Block MASS\n");
    printf("  4  %.13e  # CM\n",MC);
    printf("  5  %.13e  # BM\n",MB);
    printf("  6  %.13e  # TM\n",MT);
    printf("  15  %.13e  # TAM\n",MTA);
    printf("  23  %.13e  # ZM\n",MZ);
    printf("  24  %.13e  # MW\n",MW);
    printf("  25  %.13e  # HM\n",MH);
    printf("  50  %.13e  # M1N\n",M1N);
    printf("  51  %.13e  # M2N\n",M2N);
    printf("  52  %.13e  # M1C\n",M1C);
    printf("  53  %.13e  # M2C\n",M2C);
    printf("DECAY  6  %.13e  # wT\n",wT);
    printf("DECAY  23  %.13e  # wZ\n",wZ);
    printf("DECAY  24  %.13e  # wW\n",wW);
    
    printf("DECAY  25  %.13e  # wH\n",wH);
    printf("DECAY  50  %.13e  # w1N\n",w1N);
    printf("DECAY  51  %.13e  # w2N\n",w2N);
    printf("DECAY  52  %.13e  # w1C\n",w1C);
    printf("DECAY  53  %.13e  # w2C\n",w2C);
    
    
    printf("\n#===========================================================\n"
           "# QUANTUM NUMBERS OF NEW STATE(S) (NON SM PDG CODE) IF ANY\n"
           "#===========================================================\n\n"
           "Block QNUMBERS 50  # R1\n"
           "    1 0  # 3 times electric charge\n"
           "    2 3  # number of spin states (2S+1)\n"
           "    3 1  # colour rep (1: singlet, 3: triplet, 8: octet)\n"
           "    4 0  # Particle/Antiparticle distinction (0=own anti)\n"
           "Block QNUMBERS 51  # R2\n"
           "    1 0  # 3 times electric charge\n"
           "    2 3  # number of spin states (2S+1)\n"
           "    3 1  # colour rep (1: singlet, 3: triplet, 8: octet)\n"
           "    4 0  # Particle/Antiparticle distinction (0=own anti)\n"
           "Block QNUMBERS 52  # R1+\n"
           "    1 3  # 3 times electric charge\n"
           "    2 3  # number of spin states (2S+1)\n"
           "    3 1  # colour rep (1: singlet, 3: triplet, 8: octet)\n"
           "    4 1  # Particle/Antiparticle distinction (0=own anti)\n"
           "Block QNUMBERS 53  # R2+\n"
           "    1 3  # 3 times electric charge\n"
           "    2 3  # number of spin states (2S+1)\n"
           "    3 1  # colour rep (1: singlet, 3: triplet, 8: octet)\n"
           "    4 1  # Particle/Antiparticle distinction (0=own anti)\n"
           );
    
}
    

input_record_t input_var[] = {
    {"cabi", "cabi = %lf", DOUBLE_T, &cabi},  //      0                 
    {"aS", "aS = %lf", DOUBLE_T, &aS},        //      1       
    {"EE", "EE = %lf", DOUBLE_T, &EE},        //      2       
    {"GF", "GF = %lf", DOUBLE_T, &GF},        //      3       
    {"gt", "gt = %lf", DOUBLE_T, &gt},        //      4       
    {"MA", "MA = %lf", DOUBLE_T, &MA},        //      5       
    {"PS", "PS = %lf", DOUBLE_T, &PS},        //      6       
    {"rs", "rs = %lf", DOUBLE_T, &rs},        //      7       
    {"MC", "MC = %lf", DOUBLE_T, &MC},        //      8       
    {"MB", "MB = %lf", DOUBLE_T, &MB},        //      9       
    {"MT", "MT = %lf", DOUBLE_T, &MT},        //      10       
    {"MTA", "MTA = %lf", DOUBLE_T, &MTA},     //      11          
    {"MH", "MH = %lf", DOUBLE_T, &MH},        //      12       
    {"MZ", "MZ = %lf", DOUBLE_T, &MZ},        //      13       
    //{"TAM", "TAM = %lf", DOUBLE_T, &TAM},   //      14            
    //{"CM", "CM = %lf", DOUBLE_T, &CM},      //      15         
    //{"TM", "TM = %lf", DOUBLE_T, &TM},      //      16         
    //{"BM", "BM = %lf", DOUBLE_T, &BM},      //      17         
    {"wT", "wT = %lf", DOUBLE_T, &wT},        //      18       
    {"wZ", "wZ = %lf", DOUBLE_T, &wZ},        //      19       
    {"wW", "wW = %lf", DOUBLE_T, &wW},        //      20                   
    {NULL, NULL, 0, NULL}
};    

int main (int argc, const char * argv[])
{
    

	 *(double *) input_var[4].ptr  = (double) atof(argv[1]);  // "gt" : \tilde g         
	 *(double *) input_var[5].ptr  = (double) atof(argv[2]);  // "MA" : mass of the axial
	 *(double *) input_var[6].ptr  = (double) atof(argv[3]);  // "PS" : S parameter      
	 *(double *) input_var[7].ptr  = (double) atof(argv[4]);  // "rs" : C-M parameter    
	 *(double *) input_var[12].ptr = (double) atof(argv[5]);  // "MH" : Higgs mass       
	 
	 //void *ori_ptr = input_var[0].ptr;
	 //printf("[0].ptr: %x\n", input_var[0].ptr);
	 //printf("(*[0].ptr): %f\n", *(double*) input_var[0].ptr );


	 //input_var[0].ptr = &c;
	 //printf("(*[0].ptr): %f\n", *(double*) input_var[0].ptr );
	 //printf("[0].ptr: %x\n", input_var[0].ptr);


	 //printf("\n");
	 //printf("(*ori_ptr): %f\n", *(double*) ori_ptr );
	 //* (double*) ori_ptr = 999999999999999999.999;

	 //input_var[1].ptr = &c;
	 //input_var[2].ptr = &c;
	 //input_var[3].ptr = &c;
	 //input_var[4].ptr = &c;

    /* read input file */
    //read_input(input_var,argv[1]);
    
    calculate();
    
    /* write param_card.dat */
    write_param_card();
    
    return 0;
}




