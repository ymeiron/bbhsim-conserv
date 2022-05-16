#include "definitions.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#undef H

inline double sign(double x) { return (x>0) ? 1 : (-1); }

#define MAXIT 30
#define XACC 1e-8
double kepler_equation_E(double MA, double e)
// returns H that satisfies MA = E - e*sin(E).
//we use rtsec; must have e < 1
{
    int j;
    double fl,f,dx,swap,xl,rts;

    fl = MA;
    f  = MA - 2*M_PI;
    if (fabs(fl) < fabs(f)) {
        rts = 0;
        xl = 2*M_PI;
        swap = fl;
        fl = f;
        f = swap;
    } else {
        xl = 0;
        rts = 2*M_PI;
    }
    for (j = 1; j <= MAXIT; j++) {
        dx = (xl-rts)*f/(f-fl);
        xl = rts;
        fl = f;
        rts += dx;
        f = MA - rts + e*sin(rts);
//         fprintf(stderr,"%f\n",rts);
        if (fabs(dx) < XACC || f == 0.0) return rts;
    }
    
    // If you got here, then probably e is close to 1 and MA is close to either 0 or 2pi. Use an approximation.
    // TODO create another function "KeplerEquationE" that calculates the answer more carefully.
//     fprintf(stderr, "#####(1) MA=%e, e=%f\n",MA,e);
    swap = cbrt(MA);
    rts = (1.5874010519681996*(14 - 7*e + 2*e*e)*swap*swap - 2.8844991406148166*(5 - 7*e + 2*e*e))/(7.862224182626689*swap);
    f = MA - rts + e*sin(rts);
//     fprintf(stderr, "#####(1) rts=%f, f=%e\n",rts,f);
    if (f < 1e-6) return rts; // We must be more lenient here since we can't get arbitrarily close to the answer. 
//     fprintf(stderr, "#####(2) MA=%e, e=%f\n",MA,e);
    swap = cbrt(2*M_PI - MA);
    rts = 2*M_PI - (1.5874010519681996*(14 - 7*e + 2*e*e)*swap*swap - 2.8844991406148166*(5 - 7*e + 2*e*e))/(7.862224182626689*swap);
    f = MA - rts + e*sin(rts);
//     fprintf(stderr, "#####(2) rts=%f, f=%e\n",rts,f);
    if (f < 1e-6) return rts; // Ditto
    
    // that apparently didn't help.
    throw std::runtime_error("Maximum number of iterations exceeded in KeplerEquationE");
    return 0.0;
}

double kepler_equation_H(double MA, double e)
// returns H that satisfies MA = e*sinh(H) - H.
//we use rtsec; must have e > 1
{
    int j;
    double fl,f,dx,swap,xl,rts;
    double x1, x2;

    if (MA < 0) {
        x1 = cbrt(6*MA/e);
        x2 = std::min(0., -log(-2*MA/e));
    }
    else if (MA > 0) {
        x1 = std::max(0., log(2*MA/e));
        x2 = cbrt(6*MA/e);
    }
    else return 0;

    fl = MA - e*sinh(x1) + x1;
    f  = MA - e*sinh(x2) + x2;
    if (fabs(fl) < fabs(f)) {
        rts = x1;
        xl = x2;
        swap = fl;
        fl = f;
        f = swap;
    } else {
        xl = x1;
        rts = x2;
    }
    for (j = 1; j <= MAXIT; j++) {
        dx = (xl-rts)*f/(f-fl);
        xl = rts;
        fl = f;
        rts += dx;
        f = MA - e*sinh(rts) + rts;
        if (fabs(dx) < XACC || f == 0.0) return rts;
    }
    throw std::runtime_error("Maximum number of iterations exceeded in KeplerEquationH");
    return 0.0;
}
#undef MAXIT

int evolve_kepler_E(double CentralMass, double V[6], double DeltaT, double Rtarget, double *hdid, double Rin, double vin, double cosAlpha) {
    double E, L,                          // Energy and squared angular momentum
           a, e, phi0, n,                 // Parameters of the 2D hyperbola: semi-major axis, eccentricity, argument or pericentre & mean motion
           MA0, cosf0,      EA0, cosEA0,  // Initial mean anomany (MA), true anomaly (f) and hyperbolic eccentric anomaly (H)
           MA,  cosf, sinf, EA,  cosEA,   // See above
           Rout, vout,                    // Outgoing radious and velocity
           x, y, vx, vy,                  // Outgoing coordinates (cannonical)
           x1, y1, vx1, vy1,              // Actual outgoing coordinates
           tmp;                           // To store a square root used twice
    int Flag = 0;

    E = 0.5*vin*vin - CentralMass/Rin;
    L = V[0]*V[4] - V[1]*V[3];
    a = - CentralMass / (2*E);
    e = sqrt(1 + 2*E*L*L/(CentralMass*CentralMass));
    
    cosf0 = (a-a*e*e-Rin)/(e*Rin);
    if (cosf0 > 1)  cosf0 =  1;
    if (cosf0 < -1) cosf0 = -1;
    cosEA0 = (e + cosf0)/(1 + e*cosf0);
    
    if (sign(L) == sign(cosAlpha)) { // "Upper" side of the (cannonical) ellipse
        phi0 = atan2(V[1],V[0]) - acos(cosf0);
        EA0 =  acos(cosEA0);
    } else {                         // Lower side
        phi0 = atan2(V[1],V[0]) + acos(cosf0);
        EA0 =  2*M_PI-acos(cosEA0);
    }

    n  = sqrt(CentralMass/(a*a*a));
    fprintf(stderr, "a=%e,e=%e\n",a,e);
    MA0 = EA0 - e*sin(EA0);
    fprintf(stderr, "cosf0=%e,cosEA0=%e,phi0=%e,EA0=%e\n",cosf0,cosEA0,phi0,EA0);

    if (Rtarget > 0) {
        // check how long it will take to pull it out to Rtarget.... and find another name!!!
        cosEA = (a-Rtarget)/(a*e);
        if (fabs(cosEA) > 1) return ORBIT_INTEGRATOR_FAIL;
        EA = acos(cosEA);
        if (L < 0) EA = 2*M_PI - EA;
        MA = EA - e*sin(EA);
        if (((MA0<M_PI)&&(L<0)) || ((MA0>M_PI)&&(L>0))) *hdid = (2*M_PI - fabs(MA - MA0))/n;
        else *hdid = fabs(MA - MA0)/n;
        fprintf(stderr, "MA0=%e, MA=%e, time-it-willtake=%f, DeltaT=%f\n",MA0,MA,*hdid,DeltaT);
        if (*hdid <= DeltaT) {
            Rout = Rtarget;
            cosf = (e-cosEA) / (e*cosEA - 1);  // Cosine of target true anomaly
            sinf = sqrt(1 - cosf*cosf);
            if (L < 0) sinf = -sinf;
            vout = sqrt(2*E + 2*CentralMass/Rout);
            Flag = 1;
        }
        fprintf(stderr, "Rtarget=%e ==> Flag=%d\n",Rtarget, Flag);
    }

    if (Flag == 0) {
        MA = MA0 + n*DeltaT*sign(L);           // Target mean anomaly
        MA = MA - floor(MA/(2*M_PI))*2*M_PI;
        EA = kepler_equation_E(MA, e);           // Target hyperbolic eccentric anomaly
        cosEA = cos(EA);
        cosf = (e-cosEA) / (e*cosEA - 1);  // Cosine of target true anomaly
        sinf = sqrt(1 - cosf*cosf);
        if (MA > M_PI) sinf = -sinf;
        Rout = a*(1 - e*e) / (1 + e*cosf);
        vout = sqrt(2*E + 2*CentralMass/Rout);
        *hdid = DeltaT;
    }
    
    // Find coordinates on the cannonical ellipse
    x = Rout*cosf;
    y = Rout*sinf;
    tmp = sqrt(1+e*e+2*e*cosf); // Going to use this square root twice
    vx = -vout*sinf/tmp;
    //if (MA > M_PI) vx = -vx;
    vy = vout*(e + cosf)/tmp;
    if (L < 0) {
        vx = -vx;
        vy = -vy;
    }
    
    // Finally rotate hyperbola to original orientation
    x1 = x*cos(phi0) - y*sin(phi0);
    y1 = x*sin(phi0) + y*cos(phi0);
    vx1 = vx*cos(phi0) - vy*sin(phi0);
    vy1 = vx*sin(phi0) + vy*cos(phi0);
    
    V[0] = x1;
    V[1] = y1;
    V[3] = vx1;
    V[4] = vy1;
    return 0;
}

int evolve_kepler_H(double CentralMass, double V[6], double DeltaT, double Rtarget, double *hdid, double Rin, double vin, double cosAlpha) {
// This function advances a 2D hyperbolic orbit. This is a private function to be called only by AdvanceKepler,
// which determines whether the orbit is elliptical or hyperbolic, and transforms the initial conditions from 3D
// vectors to 2D vectors. The last three arguments are passed from the calling function to save double calculation.
    double E, L,                   // Energy and squared angular momentum
    a, e, phi0, n,                 // Parameters of the 2D hyperbola: semi-major axis, eccentricity, argument or pericentre & mean motion
    MA0, cosf0,      HA0, coshHA0, // Initial mean anomany (MA), true anomaly (f) and hyperbolic eccentric anomaly (HA)
    MA,  cosf, sinf, HA,  coshHA,  // See above
    Rout, vout,                    // Outgoing radious and velocity
    x, y, vx, vy,                  // Outgoing coordinates (cannonical)
    x1, y1, vx1, vy1,              // Actual outgoing coordinates
    tmp;                           // To store a square root used twice
    int Flag = 0;
    
    E = 0.5*vin*vin - CentralMass/Rin;
    L = V[0]*V[4] - V[1]*V[3];
    a = CentralMass / (2*E);
    e = sqrt(1 + 2*E*L*L/(CentralMass*CentralMass));
    
    cosf0 = (-a+a*e*e-Rin)/(e*Rin);
    coshHA0 = (e + cosf0)/(1 + e*cosf0);
    
    if (sign(L) == sign(cosAlpha)) { // "Upper" side of the (cannonical) hyperbola
        phi0 = atan2(V[1],V[0]) - acos(cosf0);
        HA0 =  acosh(coshHA0);
    } else {                         // Lower side
        phi0 = atan2(V[1],V[0]) + acos(cosf0);
        HA0 =  -acosh(coshHA0);
    }
    
    n  = sqrt(CentralMass/(a*a*a));
    MA0 = e*sinh(HA0) - HA0;
    
    if (Rtarget > 0) {
        // check how long it will take to pull it out to Rtarget.... and find another name!!!
        coshHA = (a+Rtarget)/(a*e);
        if (coshHA < 1) return ORBIT_INTEGRATOR_FAIL;
        HA = acosh(coshHA);
        if (L < 0) HA = - HA;
        MA = e*sinh(HA) - HA;
        *hdid = fabs(MA - MA0)/n;
//         fprintf(stderr, "MA0=%f, MA=%f, time-it-willtake=%f, DeltaT=%f\n",MA0,MA,*hdid,DeltaT);
        if (*hdid <= DeltaT) {
            Rout = Rtarget;
            cosf = (e-coshHA) / (e*coshHA - 1);  // Cosine of target true anomaly
            sinf = sqrt(1 - cosf*cosf);
            if (L < 0) sinf = -sinf;
            vout = sqrt(2*E + 2*CentralMass/Rout);
            Flag = 1;
        }
//         fprintf(stderr, "Rtarget=%f ==> Flag=%d\n",Rtarget, Flag);
    }
    
    if (Flag == 0) {
        MA = MA0 + n*DeltaT*sign(L);           // Target mean anomaly
        HA = kepler_equation_H(MA, e);           // Target hyperbolic eccentric anomaly
        coshHA = cosh(HA);
        cosf = (e-coshHA) / (e*coshHA - 1);  // Cosine of target true anomaly
        sinf = sqrt(1 - cosf*cosf);
        if (HA < 0) sinf = -sinf;
        Rout = a*(e*e - 1) / (1 + e*cosf);
        vout = sqrt(2*E + 2*CentralMass/Rout);
        *hdid = DeltaT;
    }
    
    // Find coordinates on the cannonical ellipse
    x = Rout*cosf;
    y = Rout*sinf;
    tmp = sqrt(1+e*e+2*e*cosf); // Going to use this square root twice
    vx = -vout*sinf/tmp;
    vy = vout*(e + cosf)/tmp;
    if (L < 0) {
        vx = -vx;
        vy = -vy;
    }
    
    // Finally rotate hyperbola to original orientation
    x1 = x*cos(phi0) - y*sin(phi0);
    y1 = x*sin(phi0) + y*cos(phi0);
    vx1 = vx*cos(phi0) - vy*sin(phi0);
    vy1 = vx*sin(phi0) + vy*cos(phi0);
    
    V[0] = x1;
    V[1] = y1;
    V[3] = vx1;
    V[4] = vy1;
    return 0;
}

int evolve_kepler(double m, double V[6], double dt, double r_target, double *hdid) {
    double Rin, vin,
    cosAlpha, sinAlpha,      // alpha is the angle between the radius and velocity vectors
    Z1, Z2, Z3, Z, cosEulerTheta, sinEulerTheta, phi, cosEulerPhi, sinEulerPhi,
           x1, y1, z1, vx1, vy1, vz1, // intermediate
    Vp[6]; // planar
    int Ret;
    Rin = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
    vin = sqrt(V[3]*V[3] + V[4]*V[4] + V[5]*V[5]);
    
//     fprintf(stderr,"Welcome to AdvanceKepler!\n");
//     fprintf(stderr,"Coordinates: %f,%f,%f,%f,%f,%f!\n",V[0],V[1],V[2],V[3],V[4],V[5]);

    cosAlpha = (V[0]*V[3] + V[1]*V[4] + V[2]*V[5])/Rin/vin;
    sinAlpha = sqrt(1 - cosAlpha*cosAlpha);
    
    // First, find the orbital plane and its Euler angles.
    Z1 = V[1]*V[5] - V[2]*V[4];
    Z2 = V[2]*V[3] - V[0]*V[5];
    Z3 = V[0]*V[4] - V[1]*V[3];
    Z = sqrt(Z1*Z1 + Z2*Z2 + Z3*Z3);
    cosEulerTheta = Z3/Z;
    sinEulerTheta = sqrt(1 - cosEulerTheta*cosEulerTheta);
    phi = atan2(Z1, -Z2);
    cosEulerPhi = cos(phi);
    sinEulerPhi = sin(phi);

    // Rotate the radius and velocity vectors to the xy plane.
    x1  =   cosEulerPhi*V[0] + sinEulerPhi*V[1];
    y1  = - sinEulerPhi*V[0] + cosEulerPhi*V[1];
    z1  =   V[2];
    vx1 =   cosEulerPhi*V[3] + sinEulerPhi*V[4];
    vy1 = - sinEulerPhi*V[3] + cosEulerPhi*V[4];
    vz1 =   V[5];

    Vp[0] =   x1;
    Vp[1] =   cosEulerTheta*y1  + sinEulerTheta*z1;
    Vp[2] =   0; // The matrix element is (- sinEulerTheta*y1  + cosEulerTheta*z1), but if we got the angles right it should be zero.
    Vp[3] =   vx1;
    Vp[4] =   cosEulerTheta*vy1 + sinEulerTheta*vz1;
    Vp[5] =   0; // Same comment

    if (vin*vin > 2*m/Rin) { // Unbound orbit / positive energy
        fprintf(stderr,"H; DeltaT=%f\n",dt);
        Ret = evolve_kepler_H(m, Vp, dt, r_target, hdid, Rin, vin, cosAlpha);
    } else {                           // Bound orbit / negative energy
        fprintf(stderr,"E; DeltaT=%f\n",dt);
        Ret = evolve_kepler_E(m, Vp, dt, r_target, hdid, Rin, vin, cosAlpha);
    }
    if (Ret != 0) return Ret; // If you failed, don't continue;
    
    // Now transform back to the original plane.
    x1  = Vp[0];
    y1  = cosEulerTheta*Vp[1];
    z1  = sinEulerTheta*Vp[1];
    vx1 = Vp[3];
    vy1 = cosEulerTheta*Vp[4];
    vz1 = sinEulerTheta*Vp[4];
    
    V[0] = cosEulerPhi*x1  - sinEulerPhi*y1;
    V[1] = sinEulerPhi*x1  + cosEulerPhi*y1;
    V[2] = z1;
    V[3] = cosEulerPhi*vx1 - sinEulerPhi*vy1;
    V[4] = sinEulerPhi*vx1 + cosEulerPhi*vy1;
    V[5] = vz1;
    
    return 0;
}