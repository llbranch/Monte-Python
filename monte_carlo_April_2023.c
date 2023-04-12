//////////////////////////////////
// AESOPE-Lite Monte Carlo
// Created by Liam Branch and Robert Johnson
// Copyright UCSC 2023
//////////////////////////////////

//////////////////////////////////
// IMPORTS
//////////////////////////////////

#include <math.h>
// #define pi 3.141592654
inline double radians(double radians) {
    return radians * (180.0 / M_PI);
}
typedef struct {
    float x, y, z;
} vector;
void normalize(vector* p) {
    float w = sqrt(p->x * p->x + p->y * p->y + p->z * p->z);
    p->x /= w; p->y /= w; p->z /= w;
}

//////////////////////////////////
// CONSTANTS
//////////////////////////////////
double c = 0.0299792; // Speed of Light in cm / ps
double q = 1.60217663e-19; // charge of electron columbs
double t_rise = 800; //ps 
double T3z=0; //cm is the bottom of T3
double T1z=33.782; //cm is the bottom of T1
double T4z=-28.07297; //cm is the bottom of T4
double T1_radius = 13; //cm
double T4_radius = 18; //cm 
double xPMT4 = 9.5*cos(110.)*2.54;
double yPMT4 = 9.5*sin(110.)*2.54;
double xPMT1 = 8.*cos(radians(-45.))*2.54;
double yPMT1 = 8.*np.sin(radians(-45.))*2.54;
double PMT1_radius = 4.6/2; //cm 
double PMT4_radius = 4.6/2; //cm 
double xPMT4=9.5*np.cos(110)*2.54;
double yPMT4=9.5*np.sin(110)*2.54;
int n_dynodes = 8;
// V = np.linspace(150,850,n_dynodes)
int V[] = {150,300,350,600,750,850};
int E_per_electron = 20
double QE = 0.23

//////////////////////////////////
// HELPER FUNCTIONS
//////////////////////////////////

// FIND SIGNIFICANT DIGIT POWER OF 10
def round_to_sig(x):
    return -int(np.floor(np.log10(np.abs(x))))

// NORMALIZE A VECTOR
def normalize(x):
    x /= np.linalg.norm(x)
    return x

// REDEF NORM FOR READABILITY
def mag(x):
    return np.linalg.norm(x)