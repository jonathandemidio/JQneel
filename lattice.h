#ifndef LATTICE_H
#define LATTICE_H

#include "RKK_determ.h" 
#include "simparam.h"
// Class declarations which define the 
// square lattice using plaquettes
//
class LATTICE
{

 public:
  
  // NEW VARIABLES
  int Nsite; // number of sites
  int Nbond; // total bonds
  int Nplaq; // total plaqs

  int Lx;
  int Ly;

   vector<vector<int> >Bnd; 
   vector<vector<int> >Plq; 
   

   LATTICE(const PARAMS& param);  //takes all parameters as input
   void print();   //debugger


};

#endif
