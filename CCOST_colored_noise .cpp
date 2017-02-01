/*CCOST_colored_noise.cpp: written by Steven Reeves 07/29/2015, C++ source code. Integrates the Coupled Crystal Oscillator System with the effects 
 of colored noise.*/
#include <iostream>
#include "CCOST.h"
#include <stdio.h>
#include <ctime>
#include <random>
#define index(i,j,N)   i*N + j
using namespace std;

//function declaration
void ccost_integration(double v1[], double v2[], double v3[], double v4[], double v_trans[],double eta1[],double noise[], double par[9],double intepar,double noise_par[4],int seed, long long points, int N);
void ccost_trans(double v[], double v_trans[], double eta1[], double par[9],double intepar,double noise_par[4],int seed, long long points,int N);
void phasedrift(double t[], double v[],double& freq,double sum[],double tol,long long points,int N);
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
const Doub lam = 0.99;// 0.99;//0.5; //Coupling Parameter: can be varied for complete diagnostics
const Doub h   = 2.5e-12;
const Doub tol = 1.0e-11;
int main (){
   	//number of oscillators
int N = 3;
int NUM=1;
double freq;
double phase1[NUM];
double phase2[NUM];
    //integration time
const Doub t1 = 0.000, t2 = 1.0e-07; //5.0e-07;/*1.0e2;*/ //time scaled
    //number of saved integration points
const long long points = ( (t2 - t1)/h ) + 1;
	//transient integration time
const double t1_trans = 0.0, t2_trans = 5e-4; /*5.0e3*/ //time scaled
const long long points_trans = ( (t2_trans - t1_trans)/h ) + 1;   
//build time vector
    double* t;
    t = new double[points];
//prepare time vector
    t [0] = t1;
	for(int aa = 1; aa<points; aa++){
		t[aa] = t[aa-1] + h;
		}
     
	//Open file to save time series
	ofstream myfile_tsN;
	myfile_tsN.open("CCOSTts.txt");

//write time to file
	for(int kk=0;kk<points;kk++){
	 //if(kk%1==0){
	 myfile_tsN << t[kk] << '\t';
     //		}
 		} 
 		myfile_tsN << endl;
 //*/
int i=0;
   
    
    double par[9] = {R1,R2,L1,L2,C1,C2,a,b,lam};
    double intepar=h;
    //noise parameters
    Doub D = 1.0e-07; //1.0e-04;//1.0e-3;//6.4e-9;//0.0;
    Doub tau_c = 1.0e-03;
    
    
    //EM-Integration constants
    double htauc = h/tau_c, sqrtdht = sqrt(2*D*h)/tau_c;
    //Generate instance of object Normaldev_BM which generates normal deviates (Gaussian distribution) with mean mu and variance sigma
    //and seed     seed(Uses Box-Muller Transformation)
    double mu = 0.0, sigma = 1.0;
    srand(time(NULL));
    
    int seed = rand()%100;//11//4
    Normaldev_BM GaussRNG(mu,sigma,seed);
    //extra par arrays
    Doub noise_par[4] = {htauc,sqrtdht,mu,sigma};
  //Use Mersenne Twister to Generate Thread safe random numbers
    std::mt19937 generator( static_cast<unsigned int>( time( NULL ) ) );
    std::uniform_int_distribution<int> distribution( 1, 1000000 );
	cout<<"N="<<N<<endl;
double v_trans[4*N];

// Doub time[2]      = {t1, t2};
 //randomize initial conditions

   for(int ii=0; ii<4*N; ii++){
 			v_trans[ii] =(distribution( generator )%100)*1.0e-04;
    //    cout<<v_trans[ii]<<endl;
		} 
//----------------Integrate Crystal Oscillator-------------------------------------------------------
    //allocate arrays v (time series) dynamically to fix memory issues
        double* v1; double* v2; double* v3; double* v4;
    v1 = new double[N*points];
    v2 = new double[N];
    v3 = new double[N*points];
    v4 = new double[N];
	double  v[4*N];
	double eta[N];
	cout<<"Burning Transients"<<endl;
  //(v:= for transient run)	
ccost_trans(v,v_trans,eta,par, intepar, noise_par, seed, points_trans, N);// Transient integration, saves memory cost
  //(v_i:= post-transient run)
ccost_integration(v1,v2,v3,v4,v,eta,par,intepar,noise_par,seed,points,N);// Integration for saved time-series
//for(int ii=0;ii<4*N;ii++) cout<<v[ii]<<endl;
double peak[N]={0};
	 double pd = 0;
	 double T = 0;
	 double* u1;
	 double* u2;
	 u1 = new double[N*points];
	 u2 = new double[points];
	 add(v1,v3,u1,points,N);
	//vsum(u1, u2,points,N);
	for(int aa = 0;aa<points;aa++ ) u2[aa] = u1[index(0,aa,points)];
	 //cout<<"Calculating Phase Drift" << endl;	 
	 phasedrift_kernel(t,u2,freq,pd,T,tol,points); //Phase Drift Calculation
	 cout<<"Phase Drift = "<< pd <<endl;
	 peaks(t,u1,peak,points,N);
	 for(int kk=0;kk<N;kk++){
	 	cout<<"Peak ="<<peak[kk]<<endl;
	 }
	 double delT1 = dabs(peak[0]-peak[1]);
	 double delT2 = dabs(peak[0]-peak[2]);
	 //double T = 4.5454e-08;
	 cout<< "T/N= " << T/N <<" delT1 = "<<delT1<<" delT2 = "<<delT2<< endl;
	 cout<<"T/N mod delT1 = "<< fmod(T/N,delT1)<<" delT2 mod T/N = "<<fmod(delT2,T/N)<<endl;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Write important data to file

//Write time series to file
for(int ii=0;ii<N;ii++){
        for(int kk=0;kk<points;kk++){
		//if(kk%1==0){
		myfile_tsN << u1[index(ii,kk,points)] << '\t';
//		}
} 
}//*/
myfile_tsN << endl;
// */
    //unallocate the dynamic variables   
	delete[] v1; delete[] v2; delete[] v3; delete[] v4; delete[] u1; delete[] u2;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*	// Write phase drift to file
		for(int kk=0;kk<NUM;kk++){
		myfile_tsN<<phase1[kk]<<'\t';
		}
		myfile_tsN<<endl; 
		
		for(int kk=0;kk<NUM;kk++){
		myfile_tsN<<phase2[kk]<<'\t';
		}
		myfile_tsN<<endl; 
		
		 //*/

	delete[] t;	
	//Close file
	myfile_tsN.close();
	//Program Done
	cout << "DONE!" << endl;
	//Manually cause execution window to remain on screen until a keystroke
	cin.get();	
}
