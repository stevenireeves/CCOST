/*CCOST_colored_noise.cpp: written by Steven Reeves 07/29/2015, C++ source code. Integrates the Coupled Crystal Oscillator System with the effects 
 of colored noise.*/
#include <iostream>
#include <stdio.h>
#include <ctime>
#include <random>
#include <vector>
#include <algorithm>
#include <cstring>
#include "CCOST.h"


//Working Parameters
const double R1  = 30.9;
const double R2  = 1000;
const double L1  = 520e-6;
const double L2  = 260e-6;
const double C1  = 1.0e-13;
const double C2  = 2.5e-14;
const double a   = 939;
const double b   = 3e08;
const double lam = 0.99;// 0.99;//0.5; //Coupling Parameter: can be varied for complete diagnostics
const double h   = 2.5e-12;
const double tol = 1.0e-11;
const double pi  = 4.0*atan(1.0);
int main (){
   	//number of oscillators
int N = 3;
int NUM=1;
    //integration time
const double t1 = 0.000, t2 = 1.0e-7; //5.0e-07;/*1.0e2;*/ //time scaled
    //number of saved integration points
const long long points = ( (t2 - t1)/h ) + 1;
	//transient integration time
const double t1_trans = 0.0, t2_trans = 3e-5; /*5.0e3*/ //time scaled
const long long points_trans = ( (t2_trans - t1_trans)/h ) + 1;   
//build time vector
 double t;
	//Open file to save time series
	std::ofstream myfile_tsN;
	myfile_tsN.open("CCOSTts.dat");
//prepare time vector
    t = t1;
	for(int aa = 1; aa<points; aa++){
		t += h;
		myfile_tsN << t << '\t';
		}
	myfile_tsN << std::endl;     

	int i=0;
   
    
    double par[9] = {R1,R2,L1,L2,C1,C2,a,b,lam};
    double intepar=h;
    //noise parameters
    double D = 1.0e-7; //1.0e-04;//1.0e-3;//6.4e-9;//0.0;
    double tau_c = 1.0e-01;
    
    
    //EM-Integration constants
    double htauc = h/tau_c, sqrtdht = sqrt(2*D*h)/tau_c;
    //Generate instance of object Normaldev_BM which generates normal deviates (Gaussian distribution) with mean mu and variance sigma
    //and seed     seed(Uses Box-Muller Transformation)
    double mu = 0.0, sigma = 1.0;
    srand(time(NULL));
    int seed = rand();
    //extra par arrays
    double noise_par[4] = {htauc,sqrtdht,mu,sigma};
    std::mt19937 generator(static_cast<unsigned int>( time( NULL ) ) );
    std::uniform_int_distribution<int> distribution( 1, 1000000 );
    double v_trans[4*N];
//----------------Integrate Crystal Oscillator-------------------------------------------------------
    //allocate arrays v (time series) dynamically to fix memory issues
   /*--------------------------- Initial Condition -------------------------------------------------------------*/
            double v_temp[4*N] = {0.0};
	    double B22 = 4.54545455e-8;
	    double B66 = 1.51515151e-8;
	    int j;
            for(int ii=0; ii<2; ii++){
            	v_trans[2*ii] 	  =(distribution( generator )%100)*3.0e-05;
            }
	    for(int ii = 1; ii< N; ii++){
		j = (2*ii - N)%N;
		v_trans[4*ii]     = v_trans[0]*sin(2*pi*(j/N)/B22);
		v_trans[4*ii + 1] = v_trans[0]*2*pi*j/(N*B22)*cos(2*pi*(j/N)/B22);
		v_trans[4*ii + 2] = v_trans[2]*sin(2*pi*(j/N)/B66);
		v_trans[4*ii + 3] = v_trans[2]*2*pi*j/(N*B66)*cos(2*pi*(j/N)/B66);

	    }
	double eta[2*N] = {0.0};
	std::vector<std::vector<double>> v(N, std::vector<double>(points, 0.0));// set up solution object
/*--------------------------------Integrate CCOST------------------------------------------------------------*/
//	    v.resize(N, std::vector<double>(points, 0.0));//allocate space
	std::cout<<"Burning Transients"<<std::endl;
 	//(v:= for transient run)
            ccost_trans(v_temp,v_trans,eta,par, intepar, noise_par, seed, points_trans, N);// Transient integration, saves memory cost
        //(v_i:= post-transient run)
            ccost_integration(v,v_temp,eta,par,intepar,noise_par,seed,points,N);//integration for saved time-series
/*---------------------------- Write time series to file ----------------------------------------------------*/
	for(int ii=0;ii<N;ii++){
        	for(int kk=0;kk<points;kk++){
			//std::cout<< v[ii][kk] << std::endl;
			myfile_tsN << v[ii][kk] << '\t';
		} 
		myfile_tsN << std::endl;
	}//*/
	//Close file
	myfile_tsN.close();
	//Program Done
	std::cout << "DONE!" << std::endl;
}
