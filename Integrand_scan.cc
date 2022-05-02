#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <chrono>
#include <bits/stdc++.h>

/*
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
	
	
	//delta=dpdf->intrinsic(x1,x2,q,flv1,flv2,0,0)+1.*pow(-1,n)*(dpdf->intrinsic(x1,x2,q,flv1,flv2,1,1)-dpdf->intrinsic(x1,x2,q,flv1,flv2,-1,1));
	
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
*/

int main(int argc, char *argv[]){	
	
	if (argc!=6){return -1;}
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

	//take desired flavours 
	int f1=std::atoi(argv[2]);
	int f2=std::atoi(argv[3]);
	int f3=std::atoi(argv[4]);
	int f4=std::atoi(argv[5]);
	
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
	
	int points=100;

	double x1_lo=0.01;
	double x1_hi=0.99;
	
	double x2_lo=0.01;
	double x2_hi=0.99;
	
	//loop to scan equally-spaced points in x1,x2-space obeying the 
	//x1+x2<1. condition
	for(int x1_frac=0;x1_frac<points+1;x1_frac++){
		
		double x1val=x1_lo+x1_frac*(x1_hi-x1_lo)/points;
		
		double x2_frac_lim=(1-x1val-x2_lo)*points/(x2_hi-x2_lo);
		
		
		for(int x2_frac=0;x2_frac<x2_frac_lim;x2_frac++){
			
			double x2val=x2_lo+x2_frac*(x2_hi-x2_lo)/points;
			
			
		}
	}
	
	return 0;
}
