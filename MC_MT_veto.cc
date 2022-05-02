#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <chrono>
#include <bits/stdc++.h>

#include "DGS.h"
#include "dShowerStdlib.h"
#include "mstwpdf.h"
#include "PDF.h"

using namespace dShower;

double integrand_PDGS(shared_ptr<DPDF> dpdf,double x1,double x2, double xbar1,double xbar2, double yr, double q,double flv1, double flv2, double flv3,double flv4){
	//function to pull integrand for pure DGS dPD, i.e. sqaure of no polarisation
	dpdf->setYParam(yr);
	
	double delta_1=dpdf->intrinsic(x1,xbar1,q,flv1,flv2,0,0);
	double delta_2=dpdf->intrinsic(x2,xbar2,q,flv3,flv4,0,0);
		
	return delta_1*delta_2;
}

double terms_mixpol(shared_ptr<DPDF> dpdf,double x1, double x2, double yr, double q, int flv1, int flv2,int n){
	//function to pull one term of integrand for mixed polarisation DGS, i.e. non polarised DGS dPD +/- aligned or anti-aligned dPD
	dpdf->setYParam(yr);
	
	double delta;
	
	if(n==0){
		delta=2.*dpdf->intrinsic(x1,x2,q,flv1,flv2,1,1);
	}else{
		delta=2.*dpdf->intrinsic(x1,x2,q,flv1,flv2,1,-1);
	}
	
	/*
	delta=dpdf->intrinsic(x1,x2,q,flv1,flv2,0,0)+1.*pow(-1,n)*(dpdf->intrinsic(x1,x2,q,flv1,flv2,1,1)-dpdf->intrinsic(x1,x2,q,flv1,flv2,-1,1));
	*/
	return delta;
}

double integrand_mixpol(shared_ptr<DPDF> dpdf,double x1, double x2, double x3, double x4, double yr, double q, int flv1, int flv2,int flv3,int flv4,int n){
	// function to 'square' terms to produce integrand
	
	double delta_1=terms_mixpol(dpdf,x1,x3,yr,q,flv1,flv2,n);
	
	double delta_2=terms_mixpol(dpdf,x2,x4,yr,q,flv3,flv4,n);
	
	return delta_1*delta_2;
}

double terms_MSTW(shared_ptr<SPDF> spdf,double x1,double x2,double yr, double q,int flv1, int flv2,double sigma_eff){
	//function to pull one term of integrand for the MSTW product ansatz, i.e. dPDF=sPDF*sPDF*f(y,effective cross-section)
	
	double delta_1=spdf->value(x1,q,flv1,false);
	double delta_2=spdf->value(x2,q,flv2,false);
	double gfn=exp(-1.*(yr-4.)*(yr-4.)/(2*sigma_eff*sigma_eff))*pow(sqrt(2*M_PI)*sigma_eff,-1);
	return delta_1*delta_2*gfn;
}

double integrand_MSTW(shared_ptr<SPDF> spdf,double x1, double x2, double x3, double x4, double yr, double q, int flv1, int flv2,int flv3,int flv4,double sigma_eff){
	//squares MSTW terms to produce itnegrand analogue
	double delta=terms_MSTW(spdf,x1,x3,yr,q,flv1,flv2,sigma_eff)*terms_MSTW(spdf,x2,x4,yr,q,flv3,flv4,sigma_eff);
	return delta;
}

vector<double> lim_finder_positive(double iters,double com_energy,double mass,double rapprod){
	/*
	Method to determine the high and low limits on Y1/Y2 that satisfy
	x1+x3 = e^Y1+e^Y2<q/Mw for a given Y1*Y2>0. x1+x3 has more stringent
	limits than x2+x4, which is trivially<q/Mw for all Y1/Y2.
	*/
	
	//set initial values for the low and high limits, the lowest 
	double frac_l=1.;
	double frac_h=1.;
	
	vector<double> lims; //vector to store limits
	lims.push_back(frac_l); //fill vector so if method fails there will still exist a value
	lims.push_back(frac_h);
	
	//Find low limit
	for(int i=0;i<iters;i++){
		
		//starting at 1, test whether frac_l breaks the condition
		double x1=(mass/com_energy)*exp(sqrt(rapprod*frac_l));
		//double x2=(mass/com_energy)*exp(-sqrt(rapprod*frac_l));
		double x3=(mass/com_energy)*exp(sqrt(rapprod/frac_l));
		//double x4=(mass/com_energy)*exp(-sqrt(rapprod/frac_l));
			
		
		if (x1+x3>1){
			//save first value that breaks condition, better to veto
			//than to miss some phase-space
			lims[0]=frac_l;
			//cout<<"low "<<i<<endl;
			break;
		}if(i==iters-1){
			lims[0]=frac_l;
		}
		//decrease frac_l towards zero until condition broken
		frac_l+=-1/iters;
	}
	
	//Find high limit
	for(int i=0;i<iters;i++){

		//starting at 1, test whether frac_h breaks the condition
		double x1=(mass/com_energy)*exp(sqrt(rapprod*frac_h));
		//double x2=(mass/com_energy)*exp(-sqrt(rapprod*frac_h));
		double x3=(mass/com_energy)*exp(sqrt(rapprod/frac_h));
		//double x4=(mass/com_energy)*exp(-sqrt(rapprod/frac_h));
		
		if (x1+x3>1){
			lims[1]=frac_h;//save first value that breaks condition
			//cout<<"high "<<i<<endl;
			break;
		}
		//increase frac_h until condition broken
		frac_h+=pow(log(com_energy/mass),2)/(rapprod*iters);
	}
	return lims;		
}

vector<double> lim_finder_negative(double iters,double com_energy,double mass,double rapprod){
	/*
	Method to determine the high and low limits on Y1/Y2 that satisfy
	x1+x4 = e^Y1+e^-Y2<q/Mw and x2+x3=e^-Y1+e^Y2<q/Mw for a given 
	Y1*Y2<0. here we must consider both x1+x4 and x2+x3 as they give
	the higher and lower limits respectively.
	*/
	
	//set initial values for the low and high limits, the lowest 
	double frac_l=-1.;
	double frac_h=-1.;
	
	vector<double> lims; //vector to store limits
	lims.push_back(frac_l); //fill vector so if method fails there will still exist a value
	lims.push_back(frac_h);
	
	//Find low limit
	for(int i=0;i<iters;i++){

		//starting at -1, test whether frac_l breaks the condition
		double x1 = (mass/com_energy)*exp(sqrt(rapprod*frac_l));
		//double x2 = (mass/com_energy)*exp(-sqrt(rapprod*frac_l));
		//double x3 = (mass/com_energy)*exp(sqrt(rapprod/frac_l));
		double x4 = (mass/com_energy)*exp(-sqrt(rapprod/frac_l));
		
		if (x1+x4>1){
			lims[0]=frac_l;//save first value that breaks condition
			//cout<<"low "<<i<<endl;
			break;
		}if(i==iters-1){
			lims[0]=frac_l;
		}
		//decrease frac_l until condition broken
		frac_l-=pow(log(com_energy/(2*mass)+sqrt(pow(com_energy/(2*mass),2)-1)),2)/(-rapprod*iters);
	}
	
	
	//Find high limit
	for(int i=0;i<iters;i++){
		
		//starting at -1, test whether frac_l breaks the condition
		//double x1 = (mass/com_energy)*exp(sqrt(rapprod*frac_h));
		double x2 = (mass/com_energy)*exp(-sqrt(rapprod*frac_h));
		double x3 = (mass/com_energy)*exp(sqrt(rapprod/frac_h));
		//double x4 = (mass/com_energy)*exp(-sqrt(rapprod/frac_h));
		
		if (x2+x3>1){
			//save first value that breaks condition, better to veto
			//than to miss some phase-space
			lims[1]=frac_h;
			//cout<<"high "<<i<<endl;
			break;
		}
		//increase frac_h towards zero until condition broken
		frac_h+=1/iters;	
	}
	
	return lims;		
}

double charge_dict(int flv){
	//charge dictionary to map quark labels from my internal notation used
	//in python script to determine relevant quark combos for W+W+/W-W- production
	//to notation in pDF script, which I should really just have used from beginning
	int flv_internal=flv+3;	
	vector<double> charges={1./3.,-2./3.,1./3.,0.,-1./3.,2./3.,-1./3.};
	return charges[flv_internal];
}

void integrator_1flv(shared_ptr<DPDF> dpdf,shared_ptr<SPDF> spdf,int iters, int mix,vector<double> maxes,vector<double> mins,double yrmax,double yrmin,int f1,int f2,int f3,int f4,double q,double com_energy,double Mw,double rap_prod_max,double rap_prod_min,vector<double> rapprods, double rapmin,int bins,double sigma_eff){
	
	//function that actually computes the MC integration for a given quark combo, rapidity product, minimum rapidity etc.
	
	//initialise momentum fractions outside of loop
	double x1;
	double x2;
	double xbar1;
	double xbar2;
	
	//string to log whether this process describes W+W+ or W-W- production.
	std::string Wpos;
	double chsum=charge_dict(f1)+charge_dict(f2)+charge_dict(f3)+charge_dict(f4);
	if(chsum>0){
		Wpos = "W+";
	}else{
		Wpos = "W-";
	}
	
	//create files in relevant directory for each dPDF choice, *will need to be changed for a different user/system*
	ofstream file_sums;
	bool asym=false;
	bool fixedy=true;
	if(asym){
		if(mix==0){
			std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Asym_plots/PDGS/Pure DGS_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
			file_sums.open(sumname);
		}if(mix==1){
			std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Asym_plots/Mixpol/Mixpol_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
			file_sums.open(sumname);
		}if(mix==2){
			std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Asym_plots/Pospol/Pospol_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
			file_sums.open(sumname);
		}if(mix==3){
			std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Asym_plots/MSTW/MSTW_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
			file_sums.open(sumname);
		}
	}else{
		if(fixedy){
			//create files in relevant directory for each dPDF choice, *will need to be changed for a different user/system*
			if(mix==0){
				std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Distributions/Fixed_Y/PDGS_noval_nomom/Pure DGS_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
				file_sums.open(sumname);
			}if(mix==1){
				std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Distributions/Fixed_Y/Mixpol_noval_nomom/Mixpol_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
				file_sums.open(sumname);
			}if(mix==2){
				std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Distributions/Fixed_Y/Pospol_noval_nomom/Pospol_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
				file_sums.open(sumname);
			}if(mix==3){
				std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Distributions/Fixed_Y/MSTW_noval_nomom/MSTW_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
				file_sums.open(sumname);
			}
		}else{
			//create files in relevant directory for each dPDF choice, *will need to be changed for a different user/system*
			if(mix==0){
				std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Distributions/PDGS_noval_nomom/Pure DGS_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
				file_sums.open(sumname);
			}if(mix==1){
				std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Distributions/Mixpol_noval_nomom/Mixpol_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
				file_sums.open(sumname);
			}if(mix==2){
				std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Distributions/Pospol_noval_nomom/Pospol_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
				file_sums.open(sumname);
			}if(mix==3){
				std::string sumname = "/home/lthsmith/Desktop/dShowerOL/src/Integrals/Distributions/MSTW_noval_nomom/MSTW_sums_"+Wpos+"_"+std::to_string(f1)+std::to_string(f2)+std::to_string(f3)+std::to_string(f4)+'_'+std::to_string(com_energy)+"_"+std::to_string(iters)+"_"+std::to_string(rap_prod_max)+" "+std::to_string(rap_prod_min)+" "+std::to_string(rapmin)+".txt";
				file_sums.open(sumname);
			}
		}
	}
	
	//write variable descriptions to file for readability
	file_sums<<"rapprod,sum,err\n";
	
	//random variables for y integration
	std::random_device re;
	std::default_random_engine eng_y(re());
	std::uniform_real_distribution<double> unif_y(yrmin,yrmax);
	
	//randomly seed the number generators for each iteration
	eng_y.seed(std::chrono::system_clock::now().time_since_epoch().count());
	//eng_y.seed(117); //fixed seeds for debugging
	
	double y_fixed=4.0;
	
	//vector to store each monte-carlo value of the integrand per iteration
	vector<double> partial_integral;
	partial_integral.resize(iters);
	
	//intialise variables for time logging/run duration prediction
	auto start= std::chrono::system_clock::now();
	auto end = std::chrono::system_clock::now();
	double diff_doub;
	
	double rapprod;
	
	for(int k=0;k<bins;k++){
		
		//use pre-calculated values of rapidity products and
		//rapidity ratio limits to iterate over product space and 
		//randomly generate rapidity ratio
		rapprod=rapprods[k];
		
		double ratmax=maxes[k];
		double ratmin=mins[k];
		
		//variable to track how many generated ratios we veto for any 
		//criteria break
		int excess=0.;
	
		
		//random variables for ratio Y1/Y2
		std::random_device rat;
		std::default_random_engine eng_rat(rat());
		std::uniform_real_distribution<double> unif_rat(ratmin,ratmax);
		
		//randomly seed the number generator for each iteration
		eng_rat.seed(std::chrono::system_clock::now().time_since_epoch().count());
		//eng_rat.seed(117); //fixed seed for debugging
		
		
		if(rapprod>0.){
			for (int i=0;i<iters;){

			double yr = unif_y(eng_y);
			yr=y_fixed;//quick 'n' dirty fix to fix y without accidentally breaking anything
			
			//randomly generate y1 & y2
			
			double raprat = unif_rat(eng_rat);
			
			double Y1=sqrt(rapprod*raprat);
			double Y2=sqrt(rapprod/raprat);
			
			//assign value to the momentum fraction variables
			x1=(Mw/com_energy)*exp(Y1);
			xbar1=(Mw/com_energy)*exp(-Y1);
			x2=(Mw/com_energy)*exp(Y2);
			xbar2=(Mw/com_energy)*exp(-Y2);
			
			//ensure x1+x3<1 & x2+x4<1 (same condition for each DPDF)
			if(x1+x2<1. && xbar1+xbar2<1. && Y1>=rapmin && Y2>=rapmin){
				//randomly sample from the dPDF distribution according to the generated parameters, and log each value
				//into the 'partial integral' holding vector
				if(mix==0){
					partial_integral[i]=2*yr*(1/(2.*raprat))*(integrand_PDGS(dpdf,x1,xbar1,x2,xbar2,yr,q,f1,f2,f3,f4)+integrand_PDGS(dpdf,xbar1,x1,xbar2,x2,yr,q,f1,f2,f3,f4));
					partial_integral[i]+=2*yr*(1/(2.*raprat))*(integrand_PDGS(dpdf,x1,xbar1,x2,xbar2,yr,q,f1,f4,f3,f2)+integrand_PDGS(dpdf,xbar1,x1,xbar2,x2,yr,q,f1,f4,f3,f2));
				}else if(mix==1 or mix ==2){
					partial_integral[i]=2*yr*(1/(2.*raprat))*(integrand_mixpol(dpdf,x1,xbar1,x2,xbar2,yr,q,f1,f2,f3,f4,0)+integrand_mixpol(dpdf,xbar1,x1,xbar2,x2,yr,q,f1,f2,f3,f4,0));
					partial_integral[i]+=2*yr*(1/(2.*raprat))*(integrand_mixpol(dpdf,x1,xbar1,x2,xbar2,yr,q,f1,f4,f3,f2,1)+integrand_mixpol(dpdf,xbar1,x1,xbar2,x2,yr,q,f1,f4,f3,f2,1));
				}else if(mix==3){
					partial_integral[i]=2*yr*(1/(2.*raprat))*(integrand_MSTW(spdf,x1,xbar1,x2,xbar2,yr,q,f1,f2,f3,f4,sigma_eff)+integrand_MSTW(spdf,xbar1,x1,xbar2,x2,yr,q,f1,f2,f3,f4,sigma_eff));
					partial_integral[i]+=2*yr*(1/(2.*raprat))*(integrand_MSTW(spdf,x1,xbar1,x2,xbar2,yr,q,f1,f4,f3,f2,sigma_eff)+integrand_MSTW(spdf,xbar1,x1,xbar2,x2,yr,q,f1,f4,f3,f2,sigma_eff));
				}				
				
				if((i+1)%iters==0){
					//each time we complete an integral over a specific rapprod,
					//print an estimation of remaining time for this quark combination
					end = std::chrono::system_clock::now();
					std::chrono::duration<double> diff=end-start;
					diff_doub=diff.count()*(1.*bins/(1.*k+1.)-1.);	
					std::cout<<"Estimated time remaining for Integral "<<f1<<f2<<f3<<f4<<" at rapmin: "<<rapmin<<" is: "<<diff_doub<<"s.\n";
					//I realise that this doesn't take into account the extra time (which can be very significant) from the vetos 
					//as we enforce a minimum rapidity value, but it makes me feel better to have a vague eta so I'm keeping it.
				}

				i++;
				
			}else{
				excess++;
				//if we reach a large number of exceptions (happens when close to the minimum rapidity, by a large degree)
				//print out that there has been a large number of vetos. This again has minimal utility, and will likely be removed later 
				if((excess+1)%100000==0){
					std::cout<<"This is exception: "<<excess+1<<" attempting Iter "<<i+1<<" for "<<f1<<f2<<f3<<f4<<" at rapmin: "<<rapmin<<"."<<endl;
				}
			}	
			}
		}else{
			for (int i=0;i<iters;){
				
				//this section is identical to the above, but with the veto conditions x1+xbar2<1 & xbar1+x2<1
				//due to the swapping of rapidity parameters in the integral
				
				double yr = unif_y(eng_y);
				yr=y_fixed;//quick 'n' dirty fix to fix y without accidentally breaking anything
				
				double raprat = unif_rat(eng_rat);
				
				double Y1=sqrt(rapprod*raprat);
				double Y2=sqrt(rapprod/raprat);
				
				
				x1=(Mw/com_energy)*exp(Y1);
				xbar1=(Mw/com_energy)*exp(-Y1);
				x2=(Mw/com_energy)*exp(Y2);
				xbar2=(Mw/com_energy)*exp(-Y2);
				
				//ensure x1+x4<1 & x2+x3<1 (same condition for each DPDF)
				if(x1+xbar2<1. && xbar1+x2<1. && Y1>=rapmin && Y2>=rapmin){
					if(mix==0){
						partial_integral[i]=-2*yr*(1/(2.*raprat))*(integrand_PDGS(dpdf,x1,xbar1,xbar2,x2,yr,q,f1,f2,f3,f4)+integrand_PDGS(dpdf,xbar1,x1,x2,xbar2,yr,q,f1,f2,f3,f4));
						partial_integral[i]-=2*yr*(1/(2.*raprat))*(integrand_PDGS(dpdf,x1,xbar1,xbar2,x2,yr,q,f1,f4,f3,f2)+integrand_PDGS(dpdf,xbar1,x1,x2,xbar2,yr,q,f1,f4,f3,f2));
					}else if(mix==1 or mix==2){
						partial_integral[i]=-2*yr*(1/(2.*raprat))*(integrand_mixpol(dpdf,x1,xbar1,xbar2,x2,yr,q,f1,f2,f3,f4,0)+integrand_mixpol(dpdf,xbar1,x1,x2,xbar2,yr,q,f1,f2,f3,f4,0));
						partial_integral[i]-=2*yr*(1/(2.*raprat))*(integrand_mixpol(dpdf,x1,xbar1,xbar2,x2,yr,q,f1,f4,f3,f2,1)+integrand_mixpol(dpdf,xbar1,x1,x2,xbar2,yr,q,f1,f4,f3,f2,1));
					}else if(mix==3){
						partial_integral[i]=-2*yr*(1/(2.*raprat))*(integrand_MSTW(spdf,x1,xbar1,xbar2,x2,yr,q,f1,f2,f3,f4,sigma_eff)+integrand_MSTW(spdf,xbar1,x1,x2,xbar2,yr,q,f1,f2,f3,f4,sigma_eff));
						partial_integral[i]-=2*yr*(1/(2.*raprat))*(integrand_MSTW(spdf,x1,xbar1,xbar2,x2,yr,q,f1,f4,f3,f2,sigma_eff)+integrand_MSTW(spdf,xbar1,x1,x2,xbar2,yr,q,f1,f4,f3,f2,sigma_eff));
					}	

					if((i+1)%iters==0){
						
						end = std::chrono::system_clock::now();
					
						std::chrono::duration<double> diff=end-start;
						
						diff_doub=diff.count()*(1.*bins/(1.*k+1.)-1.);
						std::cout<<"Estimated time remaining for Integral "<<f1<<f2<<f3<<f4<<" at rapmin: "<<rapmin<<" is: "<<diff_doub<<"s.\n";
					}

					i++;
					
				}else{
					excess++;
					if((excess+1)%100000==0){
						std::cout<<"This is exception: "<<excess+1<<" attempting Iter "<<i+1<<" for "<<f1<<f2<<f3<<f4<<" at rapmin: "<<rapmin<<"."<<endl;
					}
				}
				}
				

		}
		
		//gives a total number of exceptions for each rapidity product after completion of the integral
		std::cout<<"There were "<<excess<<" exceptions for "<<f1<<f2<<f3<<f4<<" at rapprod = "<<rapprod<<" at rapmin: "<<rapmin<<".\n";
		
		//this section calculates the values of the integrals from the MC-generated
		//'partial integrals' vector, and writes them to the relevant file
		double sumelem;
		double sumsqelem;
		double sqsumelem;
		
		sumelem=0.;
		sumsqelem=0.;
		sqsumelem=0.;
		
		//think error calculation might be off, check against notes
		for(int j=0;j<iters;j++){
		
			double elem=partial_integral[j];
		
			sumelem += elem;
			sumsqelem += elem*elem;
		}
		sumsqelem/=iters;
		
		double sqelem=sumelem/iters;
		
		sqsumelem=sqelem*sqelem;
		if(fixedy){
			sumelem*=(ratmax-ratmin)/(iters+excess);
		}else{
			sumelem*=(ratmax-ratmin)*(yrmax-yrmin)/(iters+excess);
		}
		double errelem;
		errelem=sqrt((sumsqelem-sqsumelem)/iters);
		
		file_sums<<rapprod<<","<<sumelem<<","<<errelem<<"\n";

	}
	
	file_sums.close();
}

int main(int argc, char *argv[]){	

	if (argc!=8){return -1;}
	//set mass of W boson
	double Mw=80.385;
	double sigma_eff;
	int sigtype=3;
	//set value of the effective cross-section for auxiliary f(y,sigma_eff)
	//in MSTW dPDF analogue
	if(sigtype==0){
		//valence option
		sigma_eff=1.62*sqrt(7.5*0.87)*20.;
	}else if(sigtype==1){
		//noval option
		sigma_eff=sqrt(380*17);
	}else if(sigtype==2){
		//nomom option
		sigma_eff=sqrt(330*17);
	}else if(sigtype==3){
		//noval and nomom option
		sigma_eff=sqrt(5100);
	}
	
	//select which dPDF to use, 0=PDGS, 1=Mixpol, 2=Pospol, 3=MSTW
	int mix=std::atoi(argv[1]);
	
	//take input number of iterations over y
	int iter=std::atoi(argv[2]);
	
	//take desired flavours 
	int f1=std::atoi(argv[3]);
	int f2=std::atoi(argv[4]);
	int f3=std::atoi(argv[5]);
	int f4=std::atoi(argv[6]);
	
	//take desired minimum rapidity
	double rapmin=std::atof(argv[7]);
	
	//Take hard scale equal to mass of W boson
	double q=Mw;
	
	//Set C.o.M energy to 14 TeV, 
	int com_energy=14000.;
	
	//assertations for kinematic & dShower limitations
	assert(com_energy > 2*Mw);
	assert(q<=172);
	
	//init. parameters for PDF initialisation
	ConstParam cstpar;
	Counter info;
	
	//set variety of important parameters
	cstpar.setRootS(com_energy);
	cstpar.setNFlv(3);
	cstpar.setMHatMax(172);
	cstpar.setUnequalScale(false);
	
	//minimum value of the impact parameter, calculated from constraints in DGS.cc
	double yrmin=1/sqrt(4*172*172/(M_B0*M_B0)-4);
	double yrmax=cstpar.getYCut();
	
	//declare varibles for s/dPDFs
	shared_ptr<DPDF> dpdf;
	shared_ptr<SPDF> spdf;
	//initialise s/dPDFs and set their parameters
	if(mix==0 or mix==1 or mix==2){
		dpdf.reset(new DGSLongPolDPDF(&cstpar,&info));
	}else if(mix==3){
		spdf.reset(new SMSTW2008(&cstpar,&info));
	}else{
		std::cout<<"Please select one of the given integral prescriptions"<<"\n.";
		return -1;
	}
	//bool to determine whether to use the full range or the curtailed
	//dist plot range
	bool fullrange=true;
	
	double rap_prod_max;
	double rap_prod_min;
	
	if(fullrange){
		//absolute maximum of the rapidity products, from noting that the minimum of
		//x1+x2 in +ve sector occurs at y1/y2=1 and that must always be <1 for mom. cons
		rap_prod_max=pow(log(com_energy/(2*Mw)),2);
		//absolute maximum of the rapidity products, from noting that x1+xbar2 intercepts x2+xbar1
		//i nthe -ve sector at y1/y2=-1 and that this value must be <1
		// n.b. this has modulus larger than the maximum, and in the interest of symmetry
		//I leave the full range as an option
		rap_prod_min=-rap_prod_max;//-pow(log(com_energy/(2*Mw)+sqrt(pow(com_energy/(2*Mw),2)-1)),2);
	}else{
		rap_prod_max=6.0;
		rap_prod_min=-rap_prod_max;
	}
	
	//move the rapprod limits slightly off the extreme to ensure convergence
	//rap_prod_max-=1e-6;
	//rap_prod_min+=1e-6;
	
	//set number of bins (number of different evenly-space rapprods to consider)
	int bins = 48;
	
	//set up vectors to hold values for these space rapprods and the max & min of
	//the rapidity ratio corresponding to them
	vector<double> rapprods;
	rapprods.resize(bins);
	
	vector<double> maxes;
	vector<double> mins;
	maxes.resize(bins);
	mins.resize(bins);
	
	//number of points to sample to find limits
	int precision=1e5;
	
	//variable to help exclusion of rapprod values that lie within ymin^2>u>-ymin^2
	bool rapshift=false;
	
	//method to generate vector of rapidity products
	double rapprod = rap_prod_max-(rap_prod_max-rap_prod_min-2*pow(rapmin,2))/(2.*bins); //starting value s.t. each bin is equally spaced over the cut range
	for (int k=0;k<bins;k++){
		rapprods[k]=rapprod;
		rapprod-=(rap_prod_max-rap_prod_min-2*pow(rapmin,2))/(bins);
		
		//jump over excluded values, maintaining spacing
		if(rapprod<=pow(rapmin,2) && rapshift==false){
			rapprod-=2*pow(rapmin,2);
			rapshift=true;
		}
	}
	
	for(int k=0;k<bins;k++){
		//loop to use limit finder functions to find max/min values of rapidity product
		if(rapprods[k]>0.){
			vector<double> lims=lim_finder_positive(precision,com_energy,Mw,rapprods[k]);
			maxes[k]=lims[1];
			mins[k]=lims[0];
		}else{
			
			vector<double> lims=lim_finder_negative(precision,com_energy,Mw,rapprods[k]);
			maxes[k]=lims[1];
			mins[k]=lims[0];
		}
	}
	
	//both informative and a warning
	std::cout<<"_______INTEGRATING__"<<f1<<f2<<f3<<f4<<"__HAVE__PATIENCE_______"<<endl;
	
	//call integrator function
	integrator_1flv(dpdf,spdf,iter,mix,maxes,mins,yrmax,yrmin,f1,f2,f3,f4,q,com_energy,Mw,rap_prod_max,rap_prod_min,rapprods,rapmin,bins,sigma_eff);
    
	return 0;
}
