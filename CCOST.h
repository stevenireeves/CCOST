#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <math.h>
#include "nr3.h"
#include "ran.h"
#include "deviates.h"
#include "gamma.h"

#define index(i,j,N)   i*N + j


void ccost_integration(double v1[], double v2[], double v3[], double v4[], double v_trans[], double eta1[], double par[9],double intepar,double noise_par[4],int seed, long long points, int N){
    
    //store initial conditions (constant history) in the solution array
    for(int jj=0; jj<N; jj++){ v1[index(jj,0,points)] = v_trans[4*jj]; v2[jj] = v_trans[4*jj+1];
        v3[index(jj,0,points)] = v_trans[4*jj+2]; v4[jj] = v_trans[4*jj+3];
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //noise parameters
    //EM-Integration constants
    //const Doub sdh = sqrt(2*D*h);
    
    double htauc = noise_par[0], sqrtdht = noise_par[1];
    //Generate instance of object Normaldev_BM which generates normal deviates (Gaussian distribution) with mean mu and variance sigma
    //and seed seed(Uses Box-Muller Transformation)
    double mu = noise_par[2], sigma = noise_par[3];
    //int seed = 2;//37
    std::random_device rd;
    std::mt19937 generator(rd());
	std::normal_distribution<double> dist1(mu,sigma);
    
    //parameters
    Doub R1  = par[0];
    Doub R2  = par[1];
    Doub L1  = par[2];
    Doub L2  = par[3];
    Doub C1  = par[4];
    Doub C2  = par[5];
    Doub a   = par[6];
    Doub b   = par[7];
    Doub lam = par[8];
    
    //Integration Parameter
    double h=intepar;
    double randn[N];
    double v2new[N];
    double v4new[N];
    //Euler-Maruyama method
    int j=0;
    //uncoupled case
    if(lam == 0){
    for (int ii=0; ii<points-1; ii++){
        //Generate normal deviates and store in array randn
        for(int rr = 0; rr < N; rr++){
		randn[rr] = dist1(generator);
		dist1.reset(); //forces eta to be uncorrelated.
    }//*/
        for(int i = 0; i < N; i++){
            //Oscillator Array
            v1[index(i,ii+1,points)] = v1[index(i,ii,points)] + h*v2[i];
           
            v2new[i] = v2[i] + h*(1/L1*((a-3*b*pow(v1[index(i,ii,points)]+v3[index(i,ii,points)],2))*(v2[i]+v4[i])-(R1*v2[i]+v1[index(i,ii,points)]/C1))) + eta1[i];
            
            v3[index(i,ii+1,points)] = v3[index(i,ii,points)] + h*v4[i];
            
            v4new[i] = v4[i] + h*(1/L2*((a-3*b*pow(v1[index(i,ii,points)]+v3[index(i,ii,points)],2))*(v2[i]+v4[i])-(R2*v4[i]+v3[index(i,ii,points)]/C2))) + eta1[i];
            
            //integrate colored noise array "eta"
           eta1[i] = eta1[i] - htauc*eta1[i] + sqrtdht*randn[i];
           
          //  noise[index(i,ii,points)] = eta1[i]; 
        }
        for(int i=0;i<N;i++){
		v2[i] = v2new[i];
    	v4[i] = v4new[i];
	}
    }
}
    //coupled case
    if(lam!=0){
    for (int ii=0; ii<points-1; ii++){
        //Generate normal deviates and store in array randn
       for(int rr = 0; rr < N; rr++){
		randn[rr] = dist1(generator);
		dist1.reset(); //forces eta to be uncorrelated.
    }//*/
        for(int i = 0; i < N; i++){
            j = (i+1)%N;
            //Oscillator Array
            v1[index(i,ii+1,points)] = v1[index(i,ii,points)] + h*v2[i];
            
            v2new[i] = v2[i] + h*(1/L1*((a-3*b*pow(v1[index(i,ii,points)]+v3[index(i,ii,points)]-lam*(v1[index(j,ii,points)]+v3[index(j,ii,points)]),2))*(v2[i]+v4[i]-lam*(v2[j]+v4[j]))-(R1*v2[i]+v1[index(i,ii,points)]/C1))) + eta1[i];
            
            v3[index(i,ii+1,points)] = v3[index(i,ii,points)] + h*v4[i];
            
            v4new[i] = v4[i] + h*(1/L2*((a-3*b*pow(v1[index(i,ii,points)]+v3[index(i,ii,points)]-lam*(v1[index(j,ii,points)]+v3[index(j,ii,points)]),2))*(v2[i]+v4[i]-lam*(v2[j]+v4[j]))-(R2*v4[i]+v3[index(i,ii,points)]/C2))) + eta1[i];
            
            //integrate colored noise array "eta"
        	int k = (i+1);
            eta1[i] = eta1[i] - htauc*eta1[i] + sqrtdht*randn[i];
          //  eta2[i] = eta2[i] - htauc*eta2[i] + sqrtdht*randn[2*i+1];
          //  noise[index(i,ii,points)] = eta1[i];
        }
        for(int i=0;i<N;i++){
		v2[i] = v2new[i];
    	v4[i] = v4new[i];
	}
    }
}
}

//distance modulus
double dmod(double a, double b){
	double temp;
	if(fmod(a,b)<fmod(b,a)) temp = fmod(a,b);
	else temp = fmod (b,a);
	return temp;
}

//calculates absolute value of num1
double dabs(double num1){
    
    if(num1 >= 0.0){
        return num1;
    }
    else{ return -num1;}
    
}

//finds index of zeros
void indexz(int ind[], double v[], int N){
	int k=0; 
	for(int ii = 0; ii<N; ii++){
		if(v[ii]==0){
		 ind[k]=ii;
		 k++;
	}
	}
}
//calculates the number of zeros in vector
int numzero(double v[], int N){
	int temp=0; 
	for(int ii=0; ii<N;ii++){
		if(v[ii]==0) temp++;
	}
	return temp;
}

int intmax(int v[], int N){
	 int temp=v[0];
	 for(int ii=0;ii<N;ii++){
		if(v[ii]>temp){
			temp=v[ii];
		 }
	 }
	 return temp;
}

double vmax(double v[], int N){
	 double temp=v[0];
	 for(int ii=0;ii<N;ii++){
		if(v[ii]>temp){
			temp=v[ii];
		 }
	 }
	 return temp;
}

double vmin(double v[], int N){
	 double temp=v[0];
	 for(int ii=0;ii<N;ii++){
		if(v[ii]<temp){
			temp=v[ii];
		 }
	 }
	 return temp;
}

double nzvmin(double v[], int N){
	double temp;
	for (int ii = 0; ii<N;ii++){
		if(v[ii]>0){
			temp = v[ii];
			break;
		}
	}
	
	for(int ii=0;ii<N;ii++){
		if(v[ii]<=temp && v[ii]>0){
			temp=v[ii];
		 }
	 }
	 return temp;
}

void ccost_trans(double v[], double v_trans[], double eta1[], double par[9],double intepar,double noise_par[4],int seed, long long points, const int N){
    
    double v1[4*N];
    //store initial conditions (constant history) in the solution array
    for(int jj=0; jj<4*N; jj++){ v1[jj] = v_trans[jj];
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //noise parameters
    //EM-Integration constants
    //const Doub sdh = sqrt(2*D*h);
    
    double htauc = noise_par[0], sqrtdht = noise_par[1];
    //Generate instance of object Normaldev_BM which generates normal deviates (Gaussian distribution) with mean mu and variance sigma
    //and seed seed(Uses Box-Muller Transformation)
    double mu = noise_par[2], sigma = noise_par[3];
   std::random_device rd;
    std::mt19937 generator(rd());
   std::normal_distribution<double> dist1(0,1);
    
    //parameters
    Doub R1  = par[0];
    Doub R2  = par[1];
    Doub L1  = par[2];
    Doub L2  = par[3];
    Doub C1  = par[4];
    Doub C2  = par[5];
    Doub a   = par[6];
    Doub b   = par[7];
    Doub lam = par[8];
    
    //Integration Parameter
    double h=intepar;
    double randn[N];
    
    //Euler-Maruyama method
    //uncoupled case
if(lam == 0){
        for(long long int ii=0;ii<points;ii++){
        //Generate normal deviates and store in array randn
       for(int rr = 0; rr < N; rr++){
		randn[rr] = dist1(generator);
		dist1.reset(); //forces eta to be uncorrelated.
    }//*/
        
        for(int i = 0; i < N; i++){
        
        	
            //Oscillator Array
            v[4*i]   = v1[4*i]   + h*v1[4*i+1];
            
            v[4*i+1] = v1[4*i+1] + h*(1/L1*((a-3*b*pow(v1[4*i]+v1[4*i+2],2))*(v1[4*i+1]+v1[4*i+3])-(R1*v1[4*i+1]+v1[4*i]/C1))) + eta1[i];
            
            v[4*i+2] = v1[4*i+2] + h*v1[4*i+3];
            
            v[4*i+3] = v1[4*i+3] + h*(1/L2*((a-3*b*pow(v1[4*i]+v1[4*i+2],2))*(v1[4*i+1]+v1[4*i+3])-(R2*v1[4*i+3]+v1[4*i+2]/C2))) + eta1[i];
            
            //integrate colored noise array "eta"
        	eta1[i] = eta1[i] - htauc*eta1[i] + sqrtdht*randn[i];
           // eta1[i] = eta1[i] + sqrtdht*randn[i];
           // eta2[i] = eta2[i] - htauc*eta2[i] + sqrtdht*randn[2*i+1];
            
        }
        for(int i=0;i<4*N;i++) v1[i]=v[i];
}
}

if(lam != 0){

    for(long long int ii=0;ii<points;ii++){
        //Generate normal deviates and store in array randn
        for(int i = 0; i < N; i++){
        		int j = (i+1)%N;
        	randn[i] = dist1(generator);
	 	   dist1.reset(); //forces eta to be uncorrelated.*/
            //Oscillator Array
            v[4*i]   = v1[4*i]   + h*v1[4*i+1];
            
            v[4*i+1] = v1[4*i+1] + h*(1/L1*((a-3*b*pow(v1[4*i]+v1[4*i+2]-lam*(v1[4*j]+v1[4*j+2]),2))*(v1[4*i+1]+v1[4*i+3]-lam*(v1[4*j+1]+v1[4*j+3]))-(R1*v1[4*i+1]+v1[4*i]/C1)))+ eta1[i];
            
            v[4*i+2] = v1[4*i+2] + h*v1[4*i+3];
            
            v[4*i+3] = v1[4*i+3] + h*(1/L2*((a-3*b*pow(v1[4*i]+v1[4*i+2]-lam*(v1[4*j]+v1[4*j+2]),2))*(v1[4*i+1]+v1[4*i+3]-lam*(v1[4*j+1]+v1[4*j+3]))-(R2*v1[4*i+3]+v1[4*i+2]/C2)))+ eta1[i];
            
            //integrate colored noise array "eta"
           eta1[i] = eta1[i] - htauc*eta1[i] + sqrtdht*randn[i];
           // eta2[i] = eta2[i] - htauc*eta2[i] + sqrtdht*randn[2*i+1];
        }
        for(int i=0;i<4*N;i++) v1[i]=v[i];
    }
}
}

void peaks(double t[],double v[], double peak[],long long points, int N){
	int ii = 1;
	double temp=0;
	while(v[index(0,ii-1,points)]>v[index(0,ii,points)] || v[index(0,ii,points)]<v[index(0,ii+1,points)]){
		temp = t[ii];
		ii++;
		if(ii == points)break;
	}
	peak[0]=temp;
	for(int i=1;i<N;i++){	
	ii =1; 
	temp = 0;
while(t[ii]< peak[0]-1.0e-09 ||(v[index(i,ii-1,points)]>v[index(i,ii,points)] || v[index(i,ii,points)]<v[index(i,ii+1,points)])){
		temp = t[ii];
		ii++;
		if(ii == points)break;
	}
	peak[i]=temp;
}
}

int sort(double t[],double v[], long long points,double lam, double T,int N){
int type;
double peak[N];
double delT1;
double delT2;

peaks(t,v,peak, points, N);
//for(int ii=0;ii<N;ii++) cout<<"peak"<<ii<<"="<< peak[ii]<< endl;

double del1 =peak[0]-peak[1]; 
double del2 =peak[0]-peak[2]; 
delT1 = dabs(del1);
delT2 = dabs(del2);
//cout<<"delT1 = "<< delT1 <<"  delT2= "<< delT2 << "  T/N= "<< T/N <<endl;
if((dmod(T/N,delT1)<9e-09 || dmod(delT2,T/N)<9e-09) && N%2!=0 && lam>0){ //RW_1
    type = 1;}
else if(delT1<1e-09 && delT2<1e-09 || lam<0){//SYNCHRONIZED
    type = -1;}
else if(dabs(delT2)<1e-09 || (N%2==0 && lam>0)){//RW_2 && RW_3
    type = 2;}
else type = 0;   //Unidentified, generally from numerical error.  

return type;
}

double vave(double v[],int N){
	double sum=0;
	double ave=0;
	for(int ii=0;ii<N;ii++){
		if(std::isnan(v[ii])==0) sum = sum+ v[ii];
	}
	ave = sum/N;
	return ave;
}

double nzvave(double v[], int N){
	double ave;
	int k = 0;
	for(int ii = 0; ii<N ; ii++){
		if(v[ii]>0 && std::isnan(v[ii])==0) k++;
	}
	if(k == 0){
//		cout<<"vector is all zero"<<endl;
		ave = 0; 
		return ave;
	}
	double nzave[k]; 
	int kk = 0; 
	for(int ii = 0; ii < N; ii++){
		if(v[ii]>0){
		 nzave[kk] = v[ii];
		 kk++;
		}
	}
	
	ave = vave(nzave,k);
	return ave;
}

double argmin(double t[],double v[],int N){
	double temp=v[0];
	double temp1 = t[0];
	 for(int ii=0;ii<N;ii++){
		if(v[ii]<temp){
			temp = v[ii];
			temp1 = t[ii];
		 }
	 }
	 return temp1;	
	}


void add(double v1[], double v3[], double u[], long long points, const int N){
    for (int ii=0;ii<N;ii++){
        for(int jj=0;jj<points;jj++){
            u[index(ii,jj,points)]=v1[index(ii,jj,points)]+v3[index(ii,jj,points)];
        }
    }
}

void vsum(double v[], double u[], long long points, const int N){
	for (int ii=0;ii<N;ii++){
        for(int jj=0;jj<points;jj++){
            u[jj]=u[jj]+v[index(ii,jj,points)];
        }
    }
}

int count(double t[],double v[],int N, long long points){
	int num1=0;
    int n[N];
    //Locating times such that x(t)=0
    for(int jj=0; jj<N; jj++){
    	for(int ii=1; ii<points-1;ii++){
            if((v[index(jj,ii-1,points)]<0) && (0<v[index(jj,ii+1,points)])){
               n[jj]=n[jj]+1;
             //  cout<<n[jj]<<endl;
            }
        }
    } 
    //num1:=number of - to + crossings
    num1=intmax(n,N);
    return num1;
}

void averagesignal(double v[],double u[], long long points, int N){
	vsum(v,u,points,N);
	for(int ii = 0; ii<points ;ii++) u[ii] = u[ii]/N;
}

double vstd(double v[], double mu, int N){
	if(N==1){
		double sigma = 0;
		return sigma;
	}
	
	double sigma = 0; 
	double delta[N];
	for(int ii=0; ii<N;ii++){
		delta[ii]=pow(v[ii]-mu,2);
		sigma = sigma +delta[ii];
	}
	sigma = sqrt(sigma/(N-1));
	return sigma;
}

void phasedrift1(double t[], double v[],double& freq, double& av, double tol,long long points, const int N){
    //cout<<N<<endl;
    double* p;
    p= new double[points]; 
        for(int jj=0;jj<points-1;jj++){
            p[jj]=0;
        }
    double slope; 
    double a0;
    double a1;
    double a2;
    double A, B, C;
    double temp1, temp2;
    int num1=0;
    //Locating times such that x(t)=0
		for(int ii=1; ii<points-1;ii++){
            if((v[ii]<0) && (0<v[ii+1])){
                //Quadratic Interpolation
                a0 = v[ii-1]/((t[ii-1]-t[ii])*(t[ii-1]-t[ii+1]));
                a1 = v[ii]/((t[ii]-t[ii-1])*(t[ii]-t[ii+1]));
                a2 = v[ii+1]/((t[ii+1]-t[ii])*(t[ii+1]-t[ii-1]));
 				 A = a0 + a1 + a2;
				 B = -((t[ii]+t[ii+1])*a0 + (t[ii-1]+t[ii+1])*a1 + (t[ii-1]+t[ii])*a2);
				 C = a0*t[ii]*t[ii+1] + a1*t[ii-1]*t[ii+1] + a2*t[ii-1]*t[ii];              
                temp1 = (-B - sqrt(pow(B,2)-4*A*C))/(2*A);
                temp2 = (-B + sqrt(pow(B,2)-4*A*C))/(2*A);
                if(temp1 < t[ii+1] && temp1 > t[ii-1]){
                	p[ii] = temp1;
				}
				else if(temp2 < t[ii+1] && temp2 > t[ii-1]){
					p[ii] = temp2;
				}//*/
				else{
						slope = (v[ii+1]-v[ii-1])/(t[ii+1]-t[ii-1]);
                		p[ii]=t[ii]-v[ii]/slope;//linear interpolation  */
				}
                num1=num1+1;
               // cout<<num1<<endl;
            }        
        }
   
    
    //num1:=number of - to + crossings
    //Break if num1=0
//	cout<<num1<<endl;
    if(num1==0) {
        cout<<"Phase Drift Failed"<<endl;
        return;}
    
    
    //remove the zero vaules in the p vector
    double* tp;
    tp = new double[num1];
    
     int kk=0;
         for(int ii=1;ii<points-1;ii++){
                if(p[ii]!=0){                    
                   tp[kk]=p[ii];
                    kk++;
                }
            }

    
    //Build periods
    double stan_period1 = 4.5454e-8;
    //Standard Period of a 22 MHz oscillation 
    double stan_period2 = 1.5152e-8;
    //Standard Period of a 66 MHz oscillation
    double* period;
    double temp;
    int num;
    period = new double[num1-2];
    
        for(int ii=0;ii<num1-2;ii++){
            temp = tp[ii+2]-tp[ii];
            if(tp[ii+2]<1e-12) temp = 0; 
          //  cout<<tp[jj][ii+2]<< " " <<tp[jj][ii] << " " << temp <<endl;
			period[ii]= temp;
            if(period[ii]>0){
            	num=num+1;
			}
        }

    int num2=num;
    //cout<<num2<<endl;
	
    int k;
    double* phasetemp;
    //Eliminating Zeros from period vector
    	k=0;
    	phasetemp=new double[num2];
        for(int ii=0;ii<num1-2;ii++){
                if(dabs(period[ii])>tol) {
				phasetemp[k]=period[ii];
				 k++;
			}
        }

    //Subtract periods to find change in period
    double* phase;
	phase= new double[num2-2];
        for(int ii=0;ii<num2-3;ii++){
            phase[ii]=dabs(phasetemp[ii+2]-phasetemp[ii]);
            if(dabs(phasetemp[ii+2]<1e-8)) phase[ii] = 0; 
            if(std::isnan(phase[ii])==1) phase[ii] = 0;
         //   cout<<phasetemp[ii+2]<< " " <<phasetemp[ii] << " " << phase[ii] <<endl;
        }

//*/ 

//cout<<"Phase temp 00  "<<phasetemp[0][num2]<<endl;
double stan_period; 
if(dabs(phasetemp[0]-stan_period1)<4e-09){
	stan_period = stan_period1;
	freq = 22; 
//	cout<< " 22 Mhz " << endl;
}
if(dabs(phasetemp[0]-stan_period2)<4e-09){
	stan_period = stan_period2;
	freq = 66; 
//	cout<< " 66 MHz " << endl;
}
//      cout<<"length phase ="<<num2-2<<endl;
    //Average over all measurements to find phase drift for each oscillator
 
        av = vave(phase,num2-2);
 		delete[] phasetemp;
    	delete[] phase; 
	//*/
    //delete dynamic allocations from memory
    delete[] p;
    delete[] tp;
    delete[] period;
   
}

double periodfinder(double t[], double v[], double tol,long long points){
    double* p;
    p= new double[points]; 
        for(int jj=0;jj<points-1;jj++){
            p[jj]=0;
        }
    
    double a0;
    double a1;
    double a2;
    double A, B, C;
    double temp1, temp2;
    int num1=0;
    //Locating times such that x(t)=0
		for(int ii=1; ii<points-1;ii++){
            if((v[ii]<0) && (0<v[ii+1])){
            	/*slope = (v[index(jj,ii+1,points)]-v[index(jj,ii-1,points)])/(t[ii+1]-t[ii-1]);
                p[jj][ii]=t[ii]-v[index(jj,ii,points)]/slope;//linear interpolation  */
                
                //Quadratic Interpolation
                a0 = v[ii-1]/((t[ii-1]-t[ii])*(t[ii-1]-t[ii+1]));
                a1 = v[ii]/((t[ii]-t[ii-1])*(t[ii]-t[ii+1]));
                a2 = v[ii+1]/((t[ii+1]-t[ii])*(t[ii+1]-t[ii-1]));
 				 A = a0 + a1 + a2;
				 B = -((t[ii]+t[ii+1])*a0 + (t[ii-1]+t[ii+1])*a1 + (t[ii-1]+t[ii])*a2);
				 C = a0*t[ii]*t[ii+1] + a1*t[ii-1]*t[ii+1] + a2*t[ii-1]*t[ii];              
                temp1 = (-B - sqrt(pow(B,2)-4*A*C))/(2*A);
                temp2 = (-B + sqrt(pow(B,2)-4*A*C))/(2*A);
                if(temp1 < t[ii+1] && temp1 > t[ii+1]){
                	p[ii] = temp1;
				}
				else{
					p[ii] = temp2;
				}
                num1=num1+1;
               // cout<<num1<<endl;
            }        
        }
   
    
    //num1:=number of - to + crossings
    
    //remove the zero vaules in the p vector
    double* tp;
    tp = new double[num1];
    
     int kk=0;
         for(int ii=1;ii<points-1;ii++){
                if(p[ii]!=0){                    
                   tp[kk]=p[ii];
                    kk++;
                }
            }

    
    //Build periods
    double* period;
    double temp;
    int num;
    period = new double[num1];
    
        for(int ii=0;ii<num1-2;ii++){
            temp = tp[ii+1]-tp[ii];
           // cout<<tp[ii+1]<< " " <<tp[ii] << " " << temp <<endl;
            if(tp[ii+2]<1e-12) temp = 0; 
			period[ii]= temp;
            if(period[ii]>0){
            	num=num+1;
			}
        }

    int num2=num;
    //cout<<num2<<endl;
	
    int k;
    double* phasetemp;
    //Eliminating Zeros from period vector
    	k=0;
    	phasetemp=new double[num2];
        for(int ii=0;ii<num1-2;ii++){
                if(dabs(period[ii])>tol) {
				phasetemp[k]=period[ii];
				 k++;
			}
        }

double T = vave(phasetemp,num2);
	delete[] phasetemp;
	//*/
    //delete dynamic allocations from memory
    delete[] p;
    delete[] tp;
    delete[] period;
return T;   
}

void phasedrift_kernel(double t[], double v[],double& freq, double& phase_error, double& T, double tol,long long points){
    double* p;
    p= new double[points]; 
        for(int jj=0;jj<points-1;jj++){
            p[jj]=0;
        }
    double slope; 
    double a0;
    double a1;
    double a2;
    double A, B, C;
    double temp1, temp2;
    int num1=0;
    //Locating times such that x(t)=0
		for(int ii=1; ii<points-1;ii++){
            if((v[ii]<0) && (0<v[ii+1])){
                
                //Quadratic Interpolation
               a0 = v[ii-1]/((t[ii-1]-t[ii])*(t[ii-1]-t[ii+1]));
                a1 = v[ii]/((t[ii]-t[ii-1])*(t[ii]-t[ii+1]));
                a2 = v[ii+1]/((t[ii+1]-t[ii])*(t[ii+1]-t[ii-1]));
 				 A = a0 + a1 + a2;
				 B = -((t[ii]+t[ii+1])*a0 + (t[ii-1]+t[ii+1])*a1 + (t[ii-1]+t[ii])*a2);
				 C = a0*t[ii]*t[ii+1] + a1*t[ii-1]*t[ii+1] + a2*t[ii-1]*t[ii];              
                temp1 = (-B - sqrt(pow(B,2)-4*A*C))/(2*A);
                temp2 = (-B + sqrt(pow(B,2)-4*A*C))/(2*A);
                if(temp1 < t[ii+1] && temp1 > t[ii-1]){
                	p[ii] = temp1;
				}
				else if(temp2 < t[ii+1] && temp2 > t[ii-1]){
					p[ii] = temp2;
				}
				else{
						slope = (v[ii+1]-v[ii-1])/(t[ii+1]-t[ii-1]);
                		p[ii]=t[ii]-v[ii]/slope;//linear interpolation  */
				} //*/
                num1=num1+1;
               // cout<<num1<<endl;
            }        
        }
   
    
    //num1:=number of - to + crossings
    //Break if num1=0
//	cout<<num1<<endl;
    if(num1==0) {
        cout<<"Phase Drift Failed"<<endl;
        return;}
    
    
    //remove the zero vaules in the p vector
    double* tp;
    tp = new double[num1];
    
     int kk=0;
         for(int ii=1;ii<points-1;ii++){
                if(p[ii]!=0){                    
                   tp[kk]=p[ii];
                    kk++;
                }
            }

    
    //Build periods
    double stan_period1 = 4.7e-8;
    //Standard Period of a 22 MHz oscillation 
    double stan_period2 = 1.5152e-8;
    //Standard Period of a 66 MHz oscillation
    double* period;
    double temp;
    int num;
    period = new double[num1];
    
        for(int ii=0;ii<num1-2;ii++){
            temp = tp[ii+1]-tp[ii];
           // cout<<tp[ii+1]<< " " <<tp[ii] << " " << temp <<endl;
            if(tp[ii+2]<1e-12) temp = 0; 
			period[ii]= temp;
            if(period[ii]>0){
            	num=num+1;
			}
        }

    int num2=num;
    //cout<<num2<<endl;
	
    int k;
    double* phasetemp;
    //Eliminating Zeros from period vector
    	k=0;
    	phasetemp=new double[num2];
        for(int ii=0;ii<num1-2;ii++){
                if(dabs(period[ii])>tol) {
				phasetemp[k]=period[ii];
				 k++;
			}
        }

 T = vave(phasetemp,num2);
 phase_error = vstd(phasetemp,T,num2);
//cout<<"Average period = "<<period_aver<<endl;
double dif1 = dabs(stan_period1-T);
double dif2 = dabs(stan_period2-T);
//cout<<"dif1 = "<< dif1<< " dif2 = "<<dif2 <<endl;
if(dif1< 1.0e-08){
	freq = 22; 
//	cout<< " 22 Mhz " << endl;
}
if(dif2< 1.0e-08){
	freq = 66; 
//	cout<< " 66 MHz " << endl;
}

	delete[] phasetemp;
	//*/
    //delete dynamic allocations from memory
    delete[] p;
    delete[] tp;
    delete[] period;
   
}

double phasedrift(double t[], double v[],double lam,double& freq, double pd[], double& T, double tol,long long points, int N){
	double phase;
	double* u;
	u = new double[points];
	if(lam != 0 ){
	//	vsum(v,u,points,N);
	for(int i =0; i<N; i++){
		phasedrift_kernel(t,u,freq,pd[i],T,tol,points);} //*/
	//	cout<<"phase error"<<ii<<" = "<<pd[ii]<<endl;
	phase = vave(pd,N);
	}//*/
	if(lam == 0){
			averagesignal(v,u,points,N);
			phasedrift_kernel(t,u,freq,phase,T,tol,points);
	}
	delete[] u;
	if(phase>1e-10) phase=0;
	return phase;
}
