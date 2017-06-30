#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <vector>
#include <cmath>


void ccost_integration(std::vector<std::vector<double> > &v_out, double v[], double eta1[], const double par[9],
		       const double intepar, const double noise_par[4], const int seed, const long long points, const int N){
    
   double v1[4*N];
    //store initial conditions (constant history) in the solution array
    for(int jj=0; jj<4*N; jj++){ 
	v1[jj] = v[jj];
    }
    
    /*Generate instance of object Normaldev_BM which generates normal deviates (Gaussian distribution) with mean mu and variance sigma
 		and seed(Uses Box-Muller Transformation)*/
    double htauc = noise_par[0], sqrtdht = noise_par[1];
    double mu = noise_par[2], sigma = noise_par[3];
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> dist1(mu,sigma);
    
    //parameters
    double R1  = par[0];
    double R2  = par[1];
    double L1  = par[2];
    double L2  = par[3];
    double C1  = par[4];
    double C2  = par[5];
    double a   = par[6];
    double b   = par[7];
    double lam = par[8];
    
    //Integration Parameter
    double h=intepar;
    double randn[2*N];
    //Euler-Maruyama method
    int j=0;
    //uncoupled case
    if(lam == 0){
    for (int ii=0; ii<points-1; ii++){
        //Generate normal deviates and store in array randn
        for(int rr = 0; rr < 2*N; rr++){
		randn[rr] = dist1(generator);
		dist1.reset(); //forces eta to be uncorrelated.
    }//*/
         for(int i = 0; i < N; i++){
        
        	
            //Oscillator Array
       v[4*i]   = v1[4*i]   + h*v1[4*i+1];
            
       v[4*i+1] = v1[4*i+1] + h*(1/L1*((a-3*b*pow(v1[4*i]+v1[4*i+2],2))
		  *(v1[4*i+1]+v1[4*i+3])-(R1*v1[4*i+1]+v1[4*i]/C1))) + eta1[2*i];
            
       v[4*i+2] = v1[4*i+2] + h*v1[4*i+3];
            
       v[4*i+3] = v1[4*i+3] + h*(1/L2*((a-3*b*pow(v1[4*i]+v1[4*i+2],2))
		  *(v1[4*i+1]+v1[4*i+3])-(R2*v1[4*i+3]+v1[4*i+2]/C2))) + eta1[2*i+1];
            
	//soln update
	v_out[i][ii] = v[4*i] + v[4*i + 2];
            //update colored noise array "eta"
        	eta1[2*i] += -htauc*eta1[2*i] + sqrtdht*randn[2*i];
        	eta1[2*i+1] += -htauc*eta1[2*i+1] + sqrtdht*randn[2*i+1];
        }
        for(int i=0;i<4*N;i++) v1[i]=v[i];
            
    }
}
    //coupled case
    if(lam!=0){
    for (int ii=0; ii<points-1; ii++){
        //Generate normal deviates and store in array randn
       for(int rr = 0; rr < 2*N; rr++){
		randn[rr] = dist1(generator);
		dist1.reset(); //forces eta to be uncorrelated.
    		}//*/
        for(int i = 0; i < N; i++){
            j = (i+1)%N;
                   //Oscillator Array
          v[4*i]   = v1[4*i]   + h*v1[4*i+1];
            
          v[4*i+1] = v1[4*i+1] + h*(1/L1*((a-3*b*pow(v1[4*i]+v1[4*i+2]-lam*(v1[4*j]+v1[4*j+2]),2))
		     *(v1[4*i+1]+v1[4*i+3]-lam*(v1[4*j+1]+v1[4*j+3]))-(R1*v1[4*i+1]+v1[4*i]/C1)))+ eta1[2*i];
            
          v[4*i+2] = v1[4*i+2] + h*v1[4*i+3];
            
          v[4*i+3] = v1[4*i+3] + h*(1/L2*((a-3*b*pow(v1[4*i]+v1[4*i+2]-lam*(v1[4*j]+v1[4*j+2]),2))
		     *(v1[4*i+1]+v1[4*i+3]-lam*(v1[4*j+1]+v1[4*j+3]))-(R2*v1[4*i+3]+v1[4*i+2]/C2)))+ eta1[2*i+1];
         
  	//soln update
	v_out[i][ii] = v[4*i] + v[4*i + 2]; 
            //update colored noise array "eta"
        	eta1[2*i] += -htauc*eta1[2*i] + sqrtdht*randn[2*i];
        	eta1[2*i+1] += -htauc*eta1[2*i+1] + sqrtdht*randn[2*i+1];
        }
        for(int i=0;i<4*N;i++){ 
		if(std::isnan(v1[i])) {
			std::cout<<"Failed to Converge"<<std::endl;
			return;
		}
		v1[i]=v[i];
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

double vmax(std::vector<std::vector<double> > &vec){
	int n = vec[0].size();
	double k = vec[0][0];
	for(int i = 0; i<n; i++){
		double temp = *std::min_element(std::begin(vec[i]),std::end(vec[i]));
		if(temp > k) k = temp; 
	}
	return k; 
}

double dmin(const double vec[], const int N){
	double k = vec[0];
	for (int i = 0; i<N; i++)
	{	
		if( k > vec[i]) k = vec[i];
	}
	return k;
}

double dmax(const double vec[], const int N){
	double k = vec[0];
	for (int i = 0; i<N; i++)
	{	
		if( k < vec[i]) k = vec[i];
	}
	return k;
}

double vmin(std::vector<std::vector<double> > &vec){
	int n = vec[0].size();
	double k = vec[0][0];
	for(int i = 0; i<n; i++){
		double temp = *std::min_element(std::begin(vec[i]),std::end(vec[i]));
		if(temp < k) k = temp; 
	}
	return k; 
}
void ccost_trans(double v[], double v_trans[], double eta1[], double par[9],double intepar,double noise_par[4],int seed, long long points, const int N){
    
    double v1[4*N];
    //store initial conditions (constant history) in the solution array
    for(int jj=0; jj<4*N; jj++) v1[jj] = v_trans[jj];
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //noise parameters
    //EM-Integration constants
    double htauc = noise_par[0], sqrtdht = noise_par[1];
    //Generate instance of object Normaldev_BM which generates normal deviates (Gaussian distribution) with mean mu and variance sigma
    //and seed seed(Uses Box-Muller Transformation)
    double mu = noise_par[2], sigma = noise_par[3];
   std::random_device rd;
    std::mt19937 generator(rd());
   std::normal_distribution<double> dist1(0,1);
    
    //parameters
    double R1  = par[0];
    double R2  = par[1];
    double L1  = par[2];
    double L2  = par[3];
    double C1  = par[4];
    double C2  = par[5];
    double a   = par[6];
    double b   = par[7];
    double lam = par[8];
    
    //Integration Parameter
    double h=intepar;
    double randn[2*N];
    
    //Euler-Maruyama method
    //uncoupled case
if(lam == 0){
        for(long long int ii=0;ii<points;ii++){
        //Generate normal deviates and store in array randn
       for(int rr = 0; rr < 2*N; rr++){
		randn[rr] = dist1(generator);
		dist1.reset(); //forces eta to be uncorrelated.
    }//*/
        
        for(int i = 0; i < N; i++){
        
        	
            //Oscillator Array
       v[4*i]   = v1[4*i]   + h*v1[4*i+1];
            
       v[4*i+1] = v1[4*i+1] + h*(1/L1*((a-3*b*pow(v1[4*i]+v1[4*i+2],2))
		  *(v1[4*i+1]+v1[4*i+3])-(R1*v1[4*i+1]+v1[4*i]/C1))) + eta1[2*i];
            
       v[4*i+2] = v1[4*i+2] + h*v1[4*i+3];
            
       v[4*i+3] = v1[4*i+3] + h*(1/L2*((a-3*b*pow(v1[4*i]+v1[4*i+2],2))
		  *(v1[4*i+1]+v1[4*i+3])-(R2*v1[4*i+3]+v1[4*i+2]/C2))) + eta1[2*i+1];
            
            //update colored noise array "eta"
        	eta1[2*i] += -htauc*eta1[2*i] + sqrtdht*randn[2*i];
        	eta1[2*i+1] += -htauc*eta1[2*i+1] + sqrtdht*randn[2*i+1];
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
            
          v[4*i+1] = v1[4*i+1] + h*(1/L1*((a-3*b*pow(v1[4*i]+v1[4*i+2]-lam*(v1[4*j]+v1[4*j+2]),2))
		     *(v1[4*i+1]+v1[4*i+3]-lam*(v1[4*j+1]+v1[4*j+3]))-(R1*v1[4*i+1]+v1[4*i]/C1))) + eta1[2*i];
            
          v[4*i+2] = v1[4*i+2] + h*v1[4*i+3];
            
          v[4*i+3] = v1[4*i+3] + h*(1/L2*((a-3*b*pow(v1[4*i]+v1[4*i+2]-lam*(v1[4*j]+v1[4*j+2]),2))
		     *(v1[4*i+1]+v1[4*i+3]-lam*(v1[4*j+1]+v1[4*j+3]))-(R2*v1[4*i+3]+v1[4*i+2]/C2))) + eta1[2*i+1];
            
            //update colored noise array "eta"
        	eta1[2*i] += -htauc*eta1[2*i] + sqrtdht*randn[2*i];
        	eta1[2*i+1] += -htauc*eta1[2*i+1] + sqrtdht*randn[2*i+1];
        }
        for(int i=0;i<4*N;i++) v1[i]=v[i];
    }
}
}


double vave(double v[],int N){
	double sum=0.0;
	double ave=0.0;
	for(int ii=0;ii<N;ii++){
		if(std::isnan(v[ii])==0) sum += v[ii];
	}
	ave = sum/N;
	return ave;
}


void vsum(std::vector<std::vector<double> > &v, std::vector<double> &u, const long long points, const int N){
	for (int ii=0;ii<N;ii++){
        for(int jj=0;jj<points;jj++){
            u[jj]+=v[ii][jj];
        }
    }
}


void averagesignal(std::vector<std::vector<double>> &v,std::vector<double> &u,const long long points,const int N){
	vsum(v,u,points,N);
	for(int ii = 0; ii<points ;ii++) u[ii] /= N;
}

double vstd(std::vector<double> &v, double mu){
	int N = v.size();
	double sigma = 0;
	if(N>1){
	double delta;
	sigma = 0; 
	for(int ii=0; ii<N;ii++){
		delta=abs(v[ii]-mu);
		sigma+=delta;
	}
	sigma = sqrt(sigma/(N-1));
	}
	return sigma;
}


void phasedrift_kernel(std::vector<double> &t, std::vector<std::vector<double>> &v,
		       double& phase_error, double tol,const long long points){
    std::vector<double> p(points, 0.0);
    int N = v.size();
    double slope; 
    double a0;
    double a1;
    double a2;
    double A, B, C;
    double temp1, temp2;
    int num1=0;
    double u1 = 0.0, u2 = 0.0, u0 = 0.0;
    //Locating times such that x(t)=0
	for(int ii=1; ii<points-1;ii++){
		for(int k = 0; k < N; k++){
			u1 += v[k][ii];
			u2 += v[k][ii+1];
			u0 += v[k][ii-1];
		}
    		if((u1<0) && (0<u2)){      
                //Quadratic Interpolation
	                a0 = u0/((t[ii-1]-t[ii])*(t[ii-1]-t[ii+1]));
	                a1 = u1/((t[ii]-t[ii-1])*(t[ii]-t[ii+1]));
       			a2 = u2/((t[ii+1]-t[ii])*(t[ii+1]-t[ii-1]));
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
				slope = (u2-u0)/(t[ii+1]-t[ii-1]);
                		p[ii]=t[ii]-u1/slope;//linear interpolation  */
			} //*/
                num1++;
            }        
        }

    if(num1==0) {
        std::cout<<"Solution Failed to Converge, outputting solution"<<std::endl;
	std::ofstream myfile_err;
    	myfile_err.open("CCOST_ERROR.dat");
	for(int j = 0; j < N; j++){
		for(int i = 0; i < points; i++) myfile_err<<v[j][i]<<'\t';
			myfile_err<<std::endl;
	}
        return;
	}
    
    
    //remove the zero vaules in the p vector
    std::vector<double> tp(num1,0.0);
    
     int kk=0;
         for(int ii=1;ii<points-1;ii++){
                if(p[ii]!=0){                    
                   tp[kk]=p[ii];
                    kk++;
                }
            }
    //remove p
    std::vector<double>().swap(p);
    
    std::vector<double> period(num1,0.0);
    double temp;
    int num;
        for(int ii=0;ii<num1-2;ii++){
            temp = tp[ii+1]-tp[ii];
            if(tp[ii+2]<1e-12) temp = 0; 
		period[ii]= temp;
            if(period[ii]>0){
            	num++;
		}
        }
    //remove tp
    std::vector<double>().swap(tp);	
    int k;
    double T = 0.0;
    std::vector<double> phasetemp(num,0.0);
    //Eliminating Zeros from period vector
    	k=0;
        for(int ii=0;ii<num1;ii++){
                if(std::abs(period[ii])>tol) {
			phasetemp[k]=period[ii];
			T += period[ii];
			std::cout<<"period = "<<period[ii]<<std::endl;
			k++;
		}
        }
    //remove period
    std::vector<double>().swap(period);
 T /= (num1-2);
 phase_error = vstd(phasetemp,T);
}

double phasedrift(std::vector<double> &t, std::vector<std::vector<double>> &v, const double tol, const long long points,const int N){
	double phase;
	phasedrift_kernel(t,v,phase,tol,points);
	return phase;
	std::cout<<"phase" << phase << std::endl;
}
