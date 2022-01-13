#ifndef PARAMSIM_H
#define PARAMSIM_H

#include <iomanip>

//Class to read in the simulation parameters from a file
class PARAMS
{
  public:
  int L1_;
  int L2_; 
  double J_; // J for projection operator
  double Q_; // Q: product of two J's
  double Beta_;  //the inverse temperature
  int suN_; // N of SU(N) symmetry of model
  int EQL_; //the number of equilibriation MC steps
  int MCS_; //the number of production MC steps
  int nBin_; // # of Bins (for statistics)
  long SEED_; //Seed

  double hmin_; // min magnetic field biasing spin=0
  double hmax_; // min magnetic field biasing spin=0
  int Nh_; // number of hfield values
  double epsfac_; // constant shift factor
  vector<double> hvec;
  double h_;//the actual value used by this thread
  double eps_;//the actual value used by this thread


  double JbN_; // J/suN
  double QbNN_;// Q/(suN*suN)
  double hb_; //hb for a bond = h/4 (on the square lattice)

  PARAMS(int my_rank){
    //initializes commonly used parameters from a file
    ifstream pfin;
    pfin.open("param.dat");    
    pfin >> L1_;
    pfin >> L2_;
    pfin >> J_;
    pfin >> Q_;
    pfin >> hmin_;
    pfin >> hmax_;
    pfin >> Nh_;
    pfin >> epsfac_;
    pfin >> Beta_;
    pfin >> EQL_;
    pfin >> MCS_;
    pfin >> nBin_;
    pfin >> SEED_;
    pfin.close();


    suN_=2;
    JbN_=J_/suN_;
    QbNN_=Q_/(suN_*suN_);

    hvec.resize(Nh_,0.0);
    double tmph=hmin_;

    double fac;    
    if(Nh_>1)
    	fac = pow(hmax_/hmin_, 1.0/(1.0*Nh_ - 1.0));
    else
	fac =1.0;

    for(int n=0; n<Nh_; n++){
        hvec[n]=tmph;
        tmph=tmph*fac;
    }

    h_=hvec[my_rank%Nh_];
    eps_=h_*epsfac_;
    hb_=h_/4.0;



    // print out parameter file to make sure its read correctly

    cout <<"SSE simulations of SU(N) square lattice JQ simulations in a field" <<endl;  
    cout <<"(Lx, Ly) = ("<<L1_<<", "<<L2_<<")"<<endl;
    cout <<"J = "<<J_<<endl;
    cout <<"Q = "<<Q_<<endl;
    cout <<"J/Q = "<<J_/Q_<<endl;
    cout <<"h = "<<h_<<endl;
    cout <<"epsilon (shift) = "<<eps_<<endl;
    cout <<"Beta = "<<Beta_<<", SU(N), N = "<<suN_<<endl;
    cout <<"EQL steps = "<<EQL_<<endl;
    cout <<"number of steps in a bin = "<<MCS_<<", number of bins = "<<nBin_<<endl;
    cout <<"SEED = "<<SEED_<<endl;

    
  }//constructor

}; //PARAMS

#endif
