/*CCOST_phasedrift_parallel.cpp: written by Steven Reeves 11/04/2015, C++ source
 code. Averages phasedrift over nit runs, for N oscillators.*/
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <ctime>
#include <random>
#include "CCOST.h"
#define index(i,j,N)   i*N + j
using namespace std;

//function declaration
void ccost_integration(double v1[], double v2[], double v3[], double v4[], double
                       v_trans[],double noise[], double par[9],double intepar,double noise_par[4],int seed, long
                       long points, int N);
void ccost_trans(double v[], double v_trans[], double par[9],double intepar,double
                 noise_par[4],int seed, long long points,int N);
double phasedrift(double t[], double v[],double lam,double& freq, double pd[], double& T, double
                  tol,long long points, int N);
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
const Doub lam = 0.99; //Coupling Parameter
const Doub h   = 2.5e-12; // Euler-Maruyama is very sensitive for this problem, need time step to be small. 
const Doub tol = 1.0e-11;
int main (){
    
    int NIT = 6;   	/* How many oscillators to go through */
    int N_initial = 3;  /* The starting number of oscillators */
    int N_max = N_initial+NIT;
    int workers = 5; /*How many threads to use, i.e. how many cores want to use*/
    int NUM=40;    /* How many samples to process, i.e. how many runs to average */
    int num1;
    double phased[4][NIT]={0}; /* The function of average phase data RW_1,2,3  and
                                Synchronized respectively*/
    double phasemin[4][NIT]={0}; /*Maxes and mins for each pattern*/
    double phasemax[4][NIT]={0};
    double freq = 0;
    //saved integration time
    const Doub t1 = 0.000, t2 = 3e-05;
    //number of saved integration points
    const long long points = ( (t2 - t1)/h ) + 1;
    //transient integration time
    const double t1_trans = 0.0, t2_trans = 2.5e-04;//1e-04;
    const long long points_trans = ( (t2_trans - t1_trans)/h ) + 1;
    double tol2 = 1e-11;
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
    
  //  int i=0;
    double par[9] = {R1,R2,L1,L2,C1,C2,a,b,lam};
    double intepar=h;
    
    //noise parameters
    Doub D =1.0e-07; //1.0e-04;//1.0e-3;//6.4e-9;//0.0;
    Doub tau_c = 1.0e-03;
    
    
    //EM-Integration constants
    double htauc = h/tau_c, sqrtdht = sqrt(2*D*h)/tau_c;
    double mu = 0.0, sigma = 1.0; // for RNG in CCOST.h
    int seed = 37;//11//4
    int N = N_initial;
    //extra par arrays
    Doub noise_par[4] = {htauc,sqrtdht,mu,sigma};
    
    
    //Use Mersenne Twister to Generate Thread safe random numbers
    std::mt19937 generator( static_cast<unsigned int>( time( NULL ) ) );
    std::uniform_int_distribution<int> distribution( 1, 1000000 );
    int n=0;
    //Begin to run through the cases
    //#pragma omp parallel for
    for(N=N_initial;N<N_max;N++){
        cout<<"N="<<N<<endl;
        /* The phase data that will be averaged for 22 and 66 MHz respectively*/
        double phase1[NUM]={0};// RW_1
        double phase2[NUM]={0};// RW_2
        double phase3[NUM]={0};// RW_3
        double phase4[NUM]={0};// SYNCHRONIZED
        int it = 0; //Iteration number
while(it<NUM/workers){
	cout<<"iteration = "<<it<<endl;
#pragma omp parallel for
        for(int threadid=0; threadid<workers;threadid++){
    	int	i = workers*it+threadid;
          	//cout<<" "<<i<<" "<<endl;
            //randomize initial conditions
            double v_trans[4*N];
            for(Int ii=0; ii<4*N; ii++){
            	if(lam == 0) v_trans[ii] =(distribution( generator )%100)*3.0e-04;
            	if(lam != 0) v_trans[ii] =(distribution( generator )%100)*1.0e-05;
            }
            //--------------------------Integrate Crystal Oscillator-----------------------------------------
            //allocate arrays v (time series) dynamically
            double* v1; double* v2; double* v3; double* v4;
            v1 = new double[N*points]();
            v2 = new double[N*points]();
            v3 = new double[N*points]();
            v4 = new double[N*points]();
            double  v[4*N]={0};
            double* noise;
            noise = new double[points]();
            num1=0;
            double eta[N];
            //(v:= for transient run)
            ccost_trans(v,v_trans,eta,par, intepar, noise_par, seed, points_trans, N);// Transient integration, saves memory cost
            //(v_i:= post-transient run)
            ccost_integration(v1,v2,v3,v4,v,eta,noise,par,intepar,noise_par,seed,points,N);//Integration for saved time-series
                
                //Prepare for phasedrift calculation
                double pd[N]={0};
            //	double pd = 0;
            double* u1;
            double* u3; //first oscillator for sorting purposes
            u1 = new double[N*points];
            u3 = new double[points];
            add(v1,v3,u1,points,N);
            // vsum(u1, u2,points,N);
            for(int hh = 0; hh<points; hh++) u3[hh] = u1[index(0,hh,points)];
            //Calculate Phase Drift	on i_k1+i_k3 for summed solution
            double T = 0;
            double phase = phasedrift(t,u1,lam,freq,pd,T,tol,points,N);
            
            //Sort Pattern
            double T1 = periodfinder(t,u3,tol,points);
            //	 cout<<"Phase Error" << phase <<endl;
            int type = 0;
            // Determine Pattern Type;
            if(lam != 0){
            type = sort(t,u1,points,lam,T1,N); //Uncoupled oscillators will not yield a pattern, therefore only calculate when lambda is nonzero
        	}
           // cout<<" Pattern " << type << endl;
            if(type == 1){
                phase1[i]=phase;
                //	cout<<phase1[i]<<endl;
            }
            else if(type == 2){
                phase2[i]=phase;
                		cout<<phase2[i]<<endl;
            }
            else if(type == 0){
                phase3[i]=phase;
               // cout<<"phase = "<< phase <<endl;
            }
            else if (type == -1){
                phase4[i]=phase;
            }
            //unallocate the dynamic variables
            delete[] v1; delete[] v2; delete[] v3; delete[] v4; delete[] u1; delete[] u3; delete[]
            noise;
            //----------------------------------------------------------------------------------------
        }
        it = it+1; 
}
        
        //Insure NUM nonzero samples in each Array where necessary
            int nz;
         // if(lam == 0) nz = numzero(phase3,NUM);
          if(lam > 0){
            if(N%2 ==1) nz = numzero(phase1,NUM);
            if(N%2 ==0) nz = numzero(phase2,NUM);
        }//*/
        
            //cout<<"nz = " << nz << endl;
            while(nz>0){
            	cout<<"nz = "<<nz<<endl;
                int indz[nz];
                indexz(indz, phase1,NUM);
#pragma omp parallel for
                for(int aa = 0; aa<workers; aa++){
                    		int j = indz[aa%nz];
                     //cout<<j<<endl; // */
                    double v_trans[4*N];
                    for(Int ii=0; ii<4*N; ii++){
                        v_trans[ii] =(distribution( generator )%100)*1.0e-05;//1e-04;
                        // cout<< "  v_trans =  " << v_trans[ii]<<endl;
                    }
                    //--------------------------Integrate Crystal Oscillator-----------------------------------------
                    //allocate arrays v (time series) dynamically to fix memory issues
                    double* v1; double* v2; double* v3; double* v4;
                    v1 = new double[N*points]();
                    v2 = new double[N*points]();
                    v3 = new double[N*points]();
                    v4 = new double[N*points]();
                    double eta[N];
                    double* noise;
                    noise = new double[points];
                    double  v[4*N]={0};
                    num1=0;
                    //(v:= for transient run)
                    ccost_trans(v,v_trans,eta,par, intepar, noise_par, seed, points_trans, N);// Transient integration, saves memory cost
                    //(v_i:= post-transient run)
                    ccost_integration(v1,v2,v3,v4,v,eta,noise,par,intepar,noise_par,seed,points,N);// Integration for saved time-series
                        
                        //Prepare for phasedrift calculation
                        double pd[N]={0};
                    double* u1;
                    double* u2; //first oscillator for sorting purposes
                    u1 = new double[N*points];
                    u2 = new double[points];
                    add(v1,v3,u1,points,N);
                    // vsum(u1, u2,points,N);
                    for(int hh = 0; hh<points; hh++) u2[hh] = u1[index(0,hh,points)];
                    //Calculate Phase Drift	on i_k1+i_k3 for summed solution
                    double T = 0;
                    double phase = phasedrift(t,u1,lam,freq,pd,T,tol,points,N);
                    //Sort Pattern
                    double T1 = periodfinder(t,u2,tol,points);
                    int type = 0;
                    type = sort(t,u1,points,lam,T1,N);
                    if(N%2 == 1){
                        //if(phase1[aa]==0){
                        if(phase1[j]==0){
                            if(type == 1){
                                phase1[j]=phase;
                                if(phase>1e-10) phase1[j]=0;
                                //	cout<<phase1[aa]<<endl;
                            }	
                        }
                    }
                    if(N%2 == 0){
                        //if(phase2[aa]==0){
                        if(phase2[j]==0){
                            if(type == 2){
                                phase2[j]=phase;
                                //		cout<<phase2[aa]<<endl;
                            }		
                        }
                    }
                    //unallocate the dynamic variables   
                    delete[] v1; delete[] v2; delete[] v3; delete[] v4; delete[] u1; delete[] u2; delete[] 
                    noise;
                }
                //if(lam == 0) nz = numzero(phase3,NUM);
                
                if(lam >0){
                if(N%2 ==1) nz = numzero(phase1,NUM);
                if(N%2 ==0) nz = numzero(phase2,NUM);
            }
                //cout<<"nz = "<<nz<<endl;
            }
//	}
        
        
        
        
        //Calculate average over the oscillations 
        phased[0][n] = vave(phase1,NUM);
        phased[1][n] = vave(phase2,NUM);
        phased[2][n] = vave(phase3,NUM);
        phased[3][n] = vave(phase4,NUM);
        phasemax[0][n] = vmax(phase1,NUM);
        phasemin[0][n] = nzvmin(phase1,NUM);
        phasemax[1][n] = vmax(phase2,NUM);
        phasemin[1][n] = nzvmin(phase2,NUM);
        phasemax[2][n] = vmax(phase3,NUM);
        phasemin[2][n] = nzvmin(phase3,NUM);
        phasemax[3][n] = vmax(phase4,NUM);
        phasemin[3][n] = nzvmin(phase4,NUM);
	
		//write summary statistics to file    
    for(int i = 0;i<4;i++){
            myfile_tsN<<phased[i][n]<<'\t'; //mean
    	}
    	 myfile_tsN<<endl;
    for(int i = 0;i<4;i++){
            myfile_tsN<<phasemin[i][n]<<'\t';//min
    	}
    	 myfile_tsN<<endl;
    for(int i = 0;i<4;i++){
            myfile_tsN<<phasemax[i][n]<<'\t';//max
    	}	
    	 myfile_tsN<<endl;
        n=n+1;
    }
   
    delete[] t;	
    //Close file
    myfile_tsN.close();
    //Program Done
    cout << "DONE!" << endl;
    //Manually cause execution window to remain on screen until a keystroke
    cin.get();	
}
