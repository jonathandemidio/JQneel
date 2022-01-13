#ifndef HAMIL_H 
#define HAMIL_H 

#include "RKK_determ.h" 
#include "simparam.h"
#include "lattice.h"
#include "MersenneTwister.h"

// Class declarations which define the
// Matrix elements and vertices according to 
// the Hamiltonian
//

class MATRIXELEMS
{
 public:

    int SWRE[8];
    int SWCO[8];
    int BOUNCE[8];
    int COST[8];

    int SWREbnd[4];
    int SWCObnd[4];
    int COSTbnd[4];

    MATRIXELEMS(const PARAMS& param, const LATTICE& latt); //constructor;
    double plx(int s1, int s2, int s3, int s4); //
    double bnd(int s1, int s2);
    
    
    //for linear programming
    vector<vector<vector<vector<double> > > >  ham;// ham[s1][s2][s3][s4]=mel
    vector<vector<vector<vector<int> > > > meltype;// meltype[s1][s2][s3][s4]=mel number
    vector<vector<int> > typespin; // typespin[melnum][0-3] = spin;
    int Nmel; //number of matrix elements
    vector<vector<vector<double> > > vtxProbs; //vtxProbs[flavor][mel]=[Pbounce, Pswre, Pcost, Pswco]


    void constructprobs(vector<int> svec,int flav); //creates vtxprobs for a given Ham
    vector<vector<vector<vector<double> > > > initializeham(const PARAMS& p);//create ham with site coordination numbers
    vector<vector<double> > simplexsolve(vector<double> Wvec); //solves the directed loop equations for each table minimizing bounces. Returns the probability table
    int pivot(vector<vector<double> > &matrix, vector<int> &avec, vector<int> &yvec, int row, int col);//for pivots

    
 private:
    vector<vector<int> > visitedvtx;//tells me which vertices have been visited
    //visitedvtx[flav][mel#] = 0 or 1

}; 

class SSEDATAS
{
  public:
    long nn;  //number of non-zero operators = nnp + nnb
    long Ma;  //length of operator string (redundant)
    vector<int> oprtr; //Operator element, 1=diag bond, 2=diag plaq, 0=null, -1 = OD bond, -2=ODplaq
    vector<int> loc; // bond or plaquette number as the case maybe
    vector<vector <int> > bndvtx; // bond vertex type on bond 1:  0 -> suN
    vector<int> Spin;  //Sz basis state (base 0)!!!!
    vector <int> nnpos; //stores the imaginary time positions of the non-null operators (to easily find them for loop updates)

    SSEDATAS(MTRand& ran,const LATTICE& latt,const PARAMS& param);

    void INCREASEM(const int newM);
    void print(const LATTICE& latt);


}; //SSEDATAS 





#endif
