/*CCOST_phasedrift_totalvariation.cpp: written by Steven Reeves 03/11/2016, C++ source code. Averages phasedrift over NIT runs, as N and Lambda Vary.*/
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <ctime>
#include <random>
#include "CCOST.h"
#define index(i,j,N)   i*N + j
using namespace std;

//function declaration
void ccost_integration(double v1[], double v2[], double v3[], double v4[], double v_trans[], double par[9],double intepar,double noise_par[4],int seed, long long points, int N);
void ccost_trans(double v[], double v_trans[], double par[9],double intepar,double noise_par[4],int seed, long long points,int N);
double phasedrift(double t[], double v[],double& freq, double pd[], double& T, double tol,long long points, int N);
void add(double v1[],double v3[], double u[], long long points, const int N);

//Working Parameters
const Doub R1  = 30.9;
const Doub R2  = 1000;
const Doub L1  = 520e-6;
const Doub L2  = 260e-6;
const Doub C1  = 1.0e-13;
const Doub C2  = 2.5e-14;
const Doub a   = 939;
const Doub b   = 3e08;
const Doub h   = 2.5e-12;
const Doub tol = 1.0e-11;
int main (){

int nlam = 10; /* how many points in the Lambda vector*/
int NIT = 5; //13;   	/* How many oscillators to go through */
int N_initial = 4;  /* The starting number of oscillators */
int N_max = N_initial+NIT;
int NUM=4;    /* How many threads to process, i.e. how many runs to average */
int num1;   
double phased[NIT][nlam+1]; /* The function of average phase data RW_1,2,3  and Synchronized respectively*/ 
double freq = 0;
	//Set up Lambda vector; 
	double lam1 = -0.99, lam2 = 0;
	double delta_lam = (lam2-lam1)/nlam;
	double lam[nlam+1]; 
	for(int ii = 0; ii<nlam+1 ; ii++) lam[ii]=lam1 + ii*delta_lam;
    //integration time
const Doub t1 = 0.000, t2 = 3e-05;
    //number of saved integration points
const long long points = ( (t2 - t1)/h ) + 1;
	//transient integration time
const double t1_trans = 0.0;
const double t2_trans = 3e-04;
long long points_trans = ( (t2_trans - t1_trans)/h ) + 1; 
double tol2 = 1e-11;
//build time vector
    double* t;
    t = new double[points];

//prepare time vector
    t [0] = t1;
	for(int aa = 1; aa<points; aa++){
		t[aa] = t[aa-1] + h;
		}
    double intepar=h;  
//noise parameters
    Doub D =1.0e-07;
    Doub tau_c = 1.0e-03;
    
    
//EM-Integration constants
    double htauc = h/tau_c, sqrtdht = sqrt(2*D*h)/tau_c;
//Generate instance of object Normaldev_BM which generates normal deviates (Gaussian distribution) with mean mu and variance sigma
//and seed     seed(Uses Box-Muller Transformation)
	double mu = 0.0, sigma = 1.0;
	int seed = 37;//11//4
    Normaldev_BM GaussRNG(mu,sigma,seed);
	int N = N_initial;
//extra par arrays
    Doub noise_par[4] = {htauc,sqrtdht,mu,sigma};     
//Open file to save time series
	ofstream myfile_tsN;
	myfile_tsN.open("CCOSTts.txt");
	
	int i=0;
for(int jj = 0; jj< nlam+1; jj++){
    double par[9] = {R1,R2,L1,L2,C1,C2,a,b,lam[jj]};


//Use Mersenne Twister to Generate Thread safe random numbers
	    std::mt19937 generator( static_cast<unsigned int>( time( NULL ) ) );
		std::uniform_int_distribution<int> distribution( 1, 1000000 );
    int n=0;
//Begin to run through the cases	  
//#pragma omp parallel for  
for(N=N_initial;N<N_max;N++){
//	if(N%2==0) N = N+1;
    cout<<"N="<<N<<endl;
    double phase1[NUM];// 
#pragma omp parallel for
 for(i=0;i<NUM;i++){
 //	cout<<" "<<i<<" "<<endl;
//randomize initial conditions
		double v_trans[4*N];
    for(Int ii=0; ii<4*N; ii++){	
        v_trans[ii] =(distribution( generator )%100)*1.0e-05;//1e-04;
       // cout<< "  v_trans =  " << v_trans[ii]<<endl;
    	} 
//--------------------------Integrate Crystal Oscillator-----------------------------------------
    //allocate arrays v (time series) dynamically to fix memory issues
    double* v1; double* v2; double* v3; double* v4;
    	v1 = new double[N*points]();  // v1 and v3 are the vectors needed. 
    	v2 = new double[N](); //v2 and v4 need not be saved, saves memory. 
    	v3 = new double[N*points]();
    	v4 = new double[N]();
	double  v[4*N];
		num1=0;
		double eta1[4];
//(v:= for transient run)	
	ccost_trans(v,v_trans,eta1,par, intepar, noise_par, seed, points_trans, N);// Transient integration, saves memory cost
//(v_i:= post-transient run)
	ccost_integration(v1,v2,v3,v4,v,eta1,par,intepar,noise_par,seed,points,N);// Integration for saved time-series

//Prepare for phasedrift calculation
	 double pd[N];
//	double pd = 0;
	 double* u1;
	 u1 = new double[N*points]; //v1 and v3 summed together
	 add(v1,v3,u1,points,N);
//Calculate Phase Drift	on i_k1+i_k3 for summed solution
	double T = 0;
	 double phase = phasedrift(t,u1,lam[jj],freq,pd,T,tol,points,N);
	 
//cout<<"Phase Error" << phase <<endl;
			phase1[i]=phase;
//unallocate the dynamic variables   
	delete[] v1; delete[] v2; delete[] v3; delete[] v4; delete[] u1;
//-----------------------------------------------------------------------------------------------------------------
}

//Calculate average over the oscillations 
 phased[n][jj] = vave(phase1,NUM);
 n=n+1;
}
}
//write phase error summary statistics to file
for(int j = 0; j<nlam +1; j++){
myfile_tsN<<lam[j]<<'\t';
for(int k=0;k<NIT;k++){
 myfile_tsN<<phased[k][j]<<'\t';
}
myfile_tsN<<endl;
}	
	delete[] t;	
//Close file
	myfile_tsN.close();
//Program Done
	cout << "DONE!" << endl;
//Manually cause execution window to remain on screen until a keystroke
	cin.get();	
}
