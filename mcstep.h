#ifndef MCSTEP_H 
#define MCSTEP_H 

//#include "RKK_determ.h"
#include "lattice.h"
#include "hamil.h"

// Class declarations which define the
// Monte Carlo Steps
//
class MCSTEP
{
  private:
    vector<int> OPLinkList;
    vector<int> LinkVXType;
    vector<int> First;
    vector<int> Last;
    vector<int> bop;
    vector<int> bnd1;
    vector<int> bnd2;
    vector<int> Bndry_Op; //0=middle, 1=start, 2=end (for t-winding #)


    long LLsize;
    double aa;  //factor that M>n
    int isEQ; //flag to determine if equilibriating
    long countVx;  //the length of each loop in #vertex
    double ll_uneq2;
    int numOperType;
    void ADJNL(const long tmpvtx);  //function to adjuct the number of loops
    double NLV; //number of vertices to cover each MCS loop update (see ADJNL)
    bool SW;

  public:
	// number of counted loops in SW algorithm OR estimated number of loops in the Wolff algorithm
	double numLoops;
	int m; // number of free spins
    // number of loops to do in each MCS for Wolff algorithm
    int numWOLFF;
    double eqlnumWOLFF; //for trying to estimate the best numWOLFF to use during equilibration
    //some thermodynamic estimators
    double nnT;  //for total energy

    double mag; // mag
    double magsq; // squared mag
    double N00; //number of diagonal J2 bonds acting on spins 0 0
    double N0; //Number of diagonal J2 bonds action on just one 0 spin
    double mN00; //same number multiplied by magnetization + for exponent
    double mN0;

  
    //things for measuring the x and y VBS order parameter
    vector<vector<vector<double> > > Nxyeo;  //Nxyeo[x or y][even or odd][00, 11, 01/10] = # of this operator
    
    MCSTEP(const LATTICE& latt, const PARAMS& p);//constructor
    
    void DIAGUPDATE(MTRand& ran, const LATTICE& latt, SSEDATAS& sdatas, MATRIXELEMS& mel, const PARAMS& p);
    void LINKOPERATOR(const LATTICE&,SSEDATAS&);// , const OPLOOP& );
    void OPERATORLOOP(MTRand& ran, const LATTICE& latt,SSEDATAS& sdatas, const MATRIXELEMS& mel, const PARAMS& p); //, const OPLOOP& loop);
    void UPDATEALPHA(MTRand& ran, const LATTICE& latt, SSEDATAS& sdatas, const PARAMS& p);
    void MEASURE(MTRand& ran, const PARAMS& p,const LATTICE& latt, const SSEDATAS& sdatas);
    int randomchoice(MTRand& ran, const vector<double> probabilities);//takes a vector of probabilities, gives me the choice (index).  Needed in below function.
    int vertexscatter(int time, int li, int flav, const SSEDATAS& sdatas, const MATRIXELEMS& mel, MTRand& ran);//gives me the outleg of a vertex

    //some measurement functions
    void MEASURE_clear(const PARAMS& p, const LATTICE& latt);

};
#endif
