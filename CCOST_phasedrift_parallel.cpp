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

#define workers  1 /*How many threads to use, i.e. how many cores want to use*/

//Working Parameters
const double R1  = 30.9;
const double R2  = 1000;
const double L1  = 520e-6;
const double L2  = 260e-6;
const double C1  = 1.0e-13;
const double C2  = 2.5e-14;
const double a   = 939;
const double b   = 3e08;
const double lam = 0.99; //Coupling Parameter
const double h   = 2.5e-12; // Euler-Maruyama is very sensitive for this problem, need time step to be small. 
const double tol = 1.0e-11;
int main (){
    
    int NIT = 6;   	/* How many oscillators to go through */
    int N_initial = 3;  /* The starting number of oscillators */
    int N_max = N_initial+NIT;
    int NUM=1;    /* How many samples to process, i.e. how many runs to average */
    double phased[NIT]={0}; /* The function of average phase data RW_1,2,3  and
                                Synchronized respectively*/
    double phasemin[NIT]={0}; /*Maxes and mins for each pattern*/
    double phasemax[NIT]={0};
    double freq = 0;
    //saved integration time
    const double t1 = 0.0, t2 = 1e-07;
    //number of saved integration points
    const long long points = ( (t2 - t1)/h ) + 1;
    //transient integration time
    const double t1_trans = 0.0, t2_trans = 1.5e-06;//1e-04;
    const long long points_trans = ( (t2_trans - t1_trans)/h ) + 1;
    double tol2 = 1e-11;
    //build time vector
    std::vector<double> t(points, 0.0);
    double pi = 4*atan(1.0);
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
        double phase1[NUM]={0};
#pragma omp parallel for num_threads(workers)
        for(int threadid=0; threadid<NUM;threadid++){
	 std::vector<std::vector<double>> v;// set up solution object
/*--------------------------- Initial Condition -------------------------------------------------------------*/
            double v_trans[4*N];
	    double B22 = 4.54545455e-8;
	    double B66 = 1.51515151e-8;
            for(int ii=0; ii<4; ii++){
            	v_trans[ii] =(distribution( generator )%100)*3.0e-04;
            }
	    for(int ii = 1; ii< N; ii++){
		v_trans[4*ii] = v_trans[0]*sin(2*pi*(ii/N)/B22);
		v_trans[4*ii + 1] = v_trans[1]*2*pi*ii/(N*B22)*cos(2*pi*(ii/N)/B22);
		v_trans[4*ii + 1] = v_trans[2]*sin(2*pi*(ii/N)/B66);
		v_trans[4*ii + 3] = v_trans[3]*2*pi*ii/(N*B66)*cos(2*pi*(ii/N)/B66);

	    }
/*--------------------------------Integrate CCOST------------------------------------------------------------*/
	    v.resize(N, std::vector<double>(points, 0.0));//allocate space
            double  v_temp[4*N]={0};
            double eta[2*N];
            //(v:= for transient run)
            ccost_trans(v_temp,v_trans,eta,par, intepar, noise_par, seed, points_trans, N);// Transient integration, saves memory cost
            //(v_i:= post-transient run)
            ccost_integration(v,v_temp,eta,par,intepar,noise_par,seed,points,N);//integration for saved time-series
                
            double phase = phasedrift(t,v,tol,points,N);
	    //deallocate memory from object. 
	    std::vector<std::vector<double>>().swap(v);
	    phase1[threadid] = phase;
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
