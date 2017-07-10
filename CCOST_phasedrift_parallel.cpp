/*Revised Source Code for Calculating the Phase Error in a Network of Crystal Oscillators. */
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <ctime>
#include <random>
#include <vector>
#include <algorithm>
#include <cstring>
#include "CCOST.h"

#define workers  5 /*How many threads to use, i.e. how many cores want to use*/

//Working Parameters
const double R1  = 30.9;
const double R2  = 1800;
const double L1  = 5.20e-4;
const double L2  = 2.60e-4;
const double C1  = 1.0e-10;
const double C2  = 2.5e-11;
const double a   = 939;
const double b   = 3e08;
const double lam = -0.90; //Coupling Parameter
const double h   = sqrt(L1*C2)*1.0e-2;
const double tol = 1.0e-11;
const double pi  = 4.0*atan(1.0);
const double B22 = sqrt(L1*C1);
const double B66 = sqrt(L2*C2);
int main (){
    
    int NIT = 9;   	/* Max osccilator = N_initial + NIT */
    int N_initial = 3;  /* The starting number of oscillators */
    int N_max = N_initial+NIT;
    int NUM=50;    /* How many samples to process, i.e. how many runs to average */
    double phased[NIT]; /* The function of average phase data*/
    double phasemin[NIT]; /*Maxes and mins for each pattern*/
    double phasemax[NIT];
	memset(&phased,0.0,sizeof(phased));
	memset(&phasemin,0.0,sizeof(phased));
	memset(&phasemax,0.0,sizeof(phased));
    double freq = 0;
    //saved integration time
    const double t1 = 0.0, t2 = 1000*B22;
    //number of saved integration points
    const long long points = ( (t2 - t1)/h ) + 1;
    //transient integration time
    const double t1_trans = 0.0, t2_trans = 10000*B22;//1e-04;
    const long long points_trans = ( (t2_trans - t1_trans)/h ) + 1;
    //build time vector
    std::vector<double> t(points, 0.0);
    //prepare time vector
#pragma omp parallel for num_threads(workers)
    for(int aa = 0; aa<points; aa++){
        t[aa] = t1 + aa*h;
    }
    
    //Open file to save time series
    std::ofstream myfile_tsN;
    myfile_tsN.open("ccost_phase_error.dat");
    
  //  int i=0;
    double par[9] = {R1,R2,L1,L2,C1,C2,a,b,lam};
    double intepar=h;
    
    //noise parameters
    double D =1.0e-07; //1.0e-04;//1.0e-3;//6.4e-9;//0.0;
    double tau_c = 1.0e-03;
    
    
    //EM-integration constants
    double htauc = h/tau_c, sqrtdht = sqrt(2*D*h)/tau_c;
    double mu = 0.0, sigma = 1.0; // for RNG in CCOST.h
    int seed = 37;//11//4
    int N = N_initial;
    //extra par arrays
    double noise_par[4] = {htauc,sqrtdht,mu,sigma};
    
    
    //Use Mersenne Twister to Generate Thread safe random numbers
    std::mt19937 generator( static_cast<unsigned int>( time( NULL ) ) );
    std::uniform_int_distribution<int> distribution( 1, 1000000 );
    int n=0;
    //Begin to run through the cases
    for(N=N_initial;N<N_max;N+=2){
        std::cout<<"N="<<N<<std::endl;
        /* The phase data that will be averaged for 22*/
        double phase1[NUM];
		memset(&phase1,0.0,sizeof(phase1));
#pragma omp parallel for num_threads(workers)
        for(int threadid=0; threadid<NUM;threadid++){
/*--------------------------- Initial Condition -------------------------------------------------------------*/
         double v_trans[4*N];
         int j;
            for(int ii=0; ii<N; ii++){
		j = (2*ii-1)%N;
		v_trans[4*ii]     = 5.0e-5*sin(2*pi*j/(N*B22));
		v_trans[4*ii + 1] = 5.0e-5*2*pi*j/(N*B22)*cos(2*pi*j/(N*B22));
		v_trans[4*ii + 2] = 5.0e-5*sin(2*pi*j/(N*B66));
		v_trans[4*ii + 3] = 5.0e-5*2*pi*j/(N*B66)*cos(2*pi*j/(N*B66));
	    }
/*--------------------------------Integrate CCOST------------------------------------------------------------*/
		std::vector<std::vector<double>> v;// set up solution object
	    v.resize(N, std::vector<double>(points, 0.0));//allocate space
            double  v_temp[4*N];
            double eta[2*N];
			memset(&eta,0.0, sizeof(eta));
			memset(&v_temp,0.0, sizeof(v_temp));
            //(v:= for transient run)
            ccost_trans(v_temp,v_trans,eta,par, intepar, noise_par, seed, points_trans, N);// Transient integration, saves memory cost
            //(v_i:= post-transient run)
            ccost_integration(v,v_temp,eta,par,intepar,noise_par,seed,points,N);//integration for saved time-series
                
            double phase = phasedrift(t,v,tol,points,N);
			phase1[threadid] = phase;
	    //deallocate memory from object. 
	    std::vector<std::vector<double>>().swap(v);
            //----------------------------------------------------------------------------------------
        }     
        //Calculate average over the oscillations 
        phased[n] = vave(phase1,NUM);
        phasemax[n] = dmax(phase1,NUM);
        phasemin[n] = dmin(phase1,NUM);
		//write summary statistics to file    
        myfile_tsN<<phased[n]<<'\t'; //mean
        myfile_tsN<<phasemin[n]<<'\t';//min
        myfile_tsN<<phasemax[n]<<'\t';//max
	myfile_tsN<<std::endl;
        n++;
 }
    //Close file
    myfile_tsN.close();
    //Program Done
    std::cout << "DONE!" << std::endl;
}
