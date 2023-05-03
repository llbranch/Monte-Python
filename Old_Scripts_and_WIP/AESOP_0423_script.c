//////////////////////////////////
// AESOP-Lite Monte Carlo
// AESOP_0423_script.c
// Created by Liam Branch and Robert Johnson
// Copyright UCSC 2023
//////////////////////////////////

//////////////////////////////////
// IMPORTS & DEFINITIONS
//////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// #include <AESOP_0423_script.h>
#define pi (double)(3.141592654)

//////////////////////////////////
// STRUCTS
//////////////////////////////////

typedef struct vector {
    double d[3];
}vector;


//////////////////////////////////
// HELPER FUNCTIONS
//////////////////////////////////


double radians(double radians) {
    return radians * (180.0 / pi);
}

// Normalize a double array of size 3
// Returns: double array
void normalize(double p[]) {
    double w = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    p[0] /= w; 
    p[1] /= w;
    p[2] /= w;
}
// Redefine norm for readability
// Returns: double scalar
double mag(vector p) {
    return sqrt(p.d[0]*p.d[0]+p.d[1]*p.d[1]+p.d[2]*p.d[2]);
}

// Find significant digit power of 10
double round_to_sig(double x) {
    return -1*((int) (floor(log(fabs(x)))));
}

// REF: https://stackoverflow.com/questions/20733590/dot-product-function-in-c-language
// Find dot prodcut of double arrays v and u of size 3
// Returns: double scalar
double dot_product(vector v, vector u) {
    int n = 3;
    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += v.d[i]*u.d[i];
    return result;
}
// Cross product across double arrays a and b of size 3 
// Returns: double array
vector cross_product(vector a, vector b) {
    vector product = {0.,0.,0.}; 
    product[0] = a.d[1] * b.d[2] - a.d[2] * b.d[1]; // error here where i left off
    product[1] = a.d[0] * b.d[2] - a.d[2] * b.d[0];
    product[2] = a.d[0] * b.d[1] - a.d[1] * b.d[0];
    return product;
}
// Multiply across double arrays a and b of size 3 
// Returns: double array
double* amult(double* a, double* b) {
    int n = 3;
    double result[n];
    for (int i = 0; i < n; i++) 
        result[n] = a[i]*b[i];
    return result;
}
// Add across double arrays a and b of size 3 
// Returns: double array
double* aadd(double* a, double* b) {
    int n = 3;
    double result[n];
    for (int i = 0; i < n; i++) 
        result[n] = a[i]+b[i];
    return result;
}
// Subtract across double arrays a and b of size 3 
// Returns: double array
double* asub(double* a, double* b) {
    int n = 3;
    double result[n];
    for (int i = 0; i < n; i++) 
        result[n] = a[i]-b[i];
    return result;
}
// Multiply double array a of size 3 and double scalar b
// Returns: double array 
double* smult(double* a, double b) {
    int n = 3;
    for (int i = 0; i < n; i++) 
        a[i] = a[i]*b;
    return a;
}
// Add double array a of size 3 and double scalar b
// Returns: double array
double* sadd(double* a, double b) {
    int n = 3;
    double result[n];
    for (int i = 0; i < n; i++)
        result[n] = a[i]+b;
    return result;
}
// Find minimum of array arr of integer size n 
// Returns: double scalar
double amin(double* arr, int n){
    double min=arr[0];
    for(int i=1;i<n;i++){
        if(arr[i]<min)
            min=arr[i];
    }
    return min;
}



// DISTANCE 2-DIM CIRCLE WITH LINE SEGMENT
// t = -D . ∆ ± √(D . ∆)^2 - |D|^2(|∆|^2 - R^2)
//     over |D|^2
// ARGUMENTS : # 3d directional vector, 3d point, center of scintillator, radius of scintillator, use corner circle boolean
double distance_circle(double* u, double* P, double* C, double radius, int quadrant) {
    double* D[3];
    if(dot_product(u,P)) { // does a normalized vector in 3d equate to normalized vector in 2d?
        double* D = smult(u,-1);
    } else {
        double* D = u;
    }
    double R = radius;
    double* bigDelta = asub(P,C);
    double magDsq = mag(D)*mag(D);
    double magDeltasq = mag(bigDelta)*mag(bigDelta);
    double DdotDelta = dot_product(D,bigDelta);
    if ((pow(DdotDelta,2) - magDsq * (magDeltasq - pow(R,2))) < 0){
        return 100.; // some large value that won't be chosen because photon has no intersection with circle
    }
    double sqrt_term = sqrt(pow(DdotDelta,2) - magDsq * (magDeltasq - pow(R,2)))/magDsq;
    double b_term = -DdotDelta/magDsq;
    double rootA = b_term - sqrt_term;
    double rootB = b_term + sqrt_term;
    if (quadrant != 0) { // if in corner don't use the other 3/4ths of the circle to find distance only 4th quadrant part
        if (fabs(rootA) > fabs(rootB)) return fabs(rootA);
        else return fabs(rootB);
    } else {
        if ((rootA < 0) && (dot_product(u,P) < 0)) return fbs(rootA);
        else return fabs(rootB);
    }
}

// ARGUMENTS : 3d directional vector, 3d point, z positions of planes bottom and top, plane dimension number
double distance_plane(double* u, double* P, double* plane, int dim) {                                     
    double d_plane = 100.;
    if (dim==2){
        if (u[dim] < 0) { // make sure direction matches location of plane 
            d_plane = plane[0];
        } else {
            d_plane = plane[1];
        }
    } else {
        d_plane = plane[0];
    }
    return abs((d_plane - P[dim])/u[dim]);
}

// SOLVE FOR DISTANCE LOGIC FUNCTION
double distance_solver(double* u, double* o, double* center, double* radius, double* plane_z, double* corner_center, double corner_radius, double* pmt_center, double pmt_radius, int *PMT_cond) {
    double dcircle = distance_circle(u,o,center,radius[0],0); // checks distance to circle boundary
    double dplane_z = distance_plane(u,o,plane_z,2); // checks distance to z boundary in general scint
    double dist = 100.; // initialize return variable
    if (dcircle > dplane_z) dist = dplane_z;
    else dist = dcircle;
    double* temp_o = amult(sadd(o,dist),u);
    int PMT_cond = 0;
    if ((temp_o[0] > 0) && (temp_o[1] < 0) && ((pow(temp_o[0],2)+pow(temp_o[1],2)) >= (pow(radius[0],2)-1))) {
        double check[4];
        check[0] = distance_plane(u,o,radius,0);                      // checks distance to x boundary
        check[1] = distance_plane(u,o,smult(radius,-1),1);                      // checks distance to y boundary
        check[2] = distance_plane(u,o,plane_z,2);                      // checks distance to z boundary inside light guide
        check[3] = distance_circle(u,o,corner_center, corner_radius, 1); // checks distance to corner boundary
        double light_guide_dist = amin(check,4);
        temp_o = aadd(o,smult(u,light_guide_dist));                                   // resuse this variable
        // if close to z = zero and within PMT circle
        if ((temp_o[2] < (plane_z[0]+0.01)) && (pow(temp_o[0]-pmt_center[0],2)+pow(temp_o[1]-pmt_center[1],2) <= pow(pmt_radius,2))) 
            PMT_cond = 1;
        return light_guide_dist;
    } else return dist;
}

//  PSEUDOCODE FOR EACH PHOTON INTERACTION WITH BOUNDARY
//         // if random number X_1 < mean ( Reflectance s_polarization + Reflectance s_polarization ):
//             // Reflect
//         // else if random number X_2 < absorbption into scintillator boundary probability:
//             // Absorbed and exit current particle simulation
//         // else if not absorbed:
//             // assume photon transmitted through boundary, 
//             // absorbed by white paint and reemmitted back 
//             // into scintillator with random direction given by random angles Phi_3, Theta_3
//             // with constraint of z coordinate entering
double* photon_interaction(double* o,double* u, double* n, int *notabsorbed, double n_1, double n_2) {
    double* u_r = aadd(u,smult(n,-2*dot_product(u, n))); // u_new = u - 2 (u . n)*n
    double* v[3];
    if(dot_product(u,n) < 0) {
        double* D = smult(u,-1);
    } else {
        double* D = u;
    }
    double theta = asin(mag(cross_product(v,n))/(mag(u)*mag(n)));
    double inside_sqrt = pow((n_1/n_2)*sin(theta),2);
    double sqrt_term = sqrt(1 - inside_sqrt);       // cos(theta)_transmission
    double Rs = pow(abs((n_1*cos(theta) - n_2*sqrt_term)/(n_1*cos(theta) + n_2*sqrt_term)),2);
    double Rp = pow(abs((n_1*sqrt_term - n_2*cos(theta))/(n_1*sqrt_term + n_2*cos(theta))),2);
    // Determine probability of reflectance
    if ((((double)rand())/RAND_MAX) < ((Rs+Rp)/2)) { 
        notabsorbed = 1;                            // if random chance is high enough reflect !
        return normalize(u_r);                      // return full internal reflection and not absorbed is True
    }                                               // else photon is transmitted to white paint
    if ((((double)rand())/RAND_MAX) < 0.80) {       // does it get absorbed? change probability when you get more data
        notabsorbed = 0;                            // not absorbed is False
        return normalize(u_r);                      
    } else {                                        // case of it didn't get absorbed!
        double theta_new = (((double)rand())/RAND_MAX)*2*pi; // new theta direction of photon
        double phi_new = (((double)rand())/RAND_MAX)*pi;     // new phi   direction of photon
        double new_u[3];
        u[0] = sin(phi_new)*cos(theta_new);
        u[1] = sin(phi_new)*sin(theta_new);
        u[2] = cos(phi_new);
        new_u = normalize(u);
        double change_factor = (((double)rand())/RAND_MAX)-0.5;      // chosen so new direction is at least...
        u_r = aadd(u_r, smult(new_u, change_factor));                // ...within bounds of scinatillator
        notabsorbed = 1;
        return normalize(u_r);
    }
}

// Calculate n vector for all planes and surfaces in apparatus
double* n_vec_calculate(double* o, double* scint_plane, double* light_guide_planes, double* corner_center, double corner_radius) {
    if (o[2] == scint_plane[0])                                    // bottom of scint
        return {0.,0.,+1.};
    if (o[2] == scint_plane[1])                                    // top of scint
        return {0.,0.,-1.};
    if (o[0] == light_guide_planes[0])                             // y plane of light guide 
        return {0.,+1.,0.};
    if (o[1] == light_guide_planes[1])                             // x plane of light guide
        return {-1.,0.,0.};
    if ((o[0] >= corner_center[0]) & (o[1] <= corner_center[1]))   // in corner
        return normalize(asub(o,corner_center));
    else {                                                         // else in main scintillator
        double* z = {0.,0.,0.};                        
        return normalize(asub(o,));
    }
}



int main() {
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
    double yPMT1 = 8.*sin(radians(-45.))*2.54;
    double PMT1_radius = 4.6/2; //cm 
    double PMT4_radius = 4.6/2; //cm 
    double xPMT4=9.5*cos(110)*2.54;
    double yPMT4=9.5*sin(110)*2.54;
    int n_dynodes = 8;
    int V[] = {150,300,350,600,750,850};
    int E_per_electron = 20;
    double QE = 1; // .23 is more realistic
    double artificial_gain = 10;
    double t_initial = 0; //ps
    double particle_init_angle_range = 40; //degrees
    double T1_width = 0.5; //cm
    double T4_width = 1; //cm
    double mean_free_path_scints = 0.00024; //cm
    double photons_produced_per_MeV = 10; //electrons
    double pr_of_scintillation = 0.8; 
    double max_simulated_reflections = 8;
    double pmt_electron_travel_time = 0; // approx 16 ns
    double artificial_gain = 10; // gain factor
}
