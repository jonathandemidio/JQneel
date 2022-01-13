#ifndef FILEOPS_H
#define FILEOPS_H

#include "RKK_determ.h"
#include "hamil.h"
#include "mcstep.h"
#include <iomanip>

// Class declarations which define file
// names and operators for reading/writing
//
class FILEOPS
{
  public:
   char fname[8];  //configuration file
   char B1_name[11]; //houses x-type bond correlations
   char B2_name[11]; // "" y-type  ""

   char B1d_name[11]; //houses diagonal x-type bond correlations
   char B2d_name[11]; // "" y-type  ""

   char dname[8];  //data file

   FILEOPS(const int); 
   void WRITECONFIG(const PARAMS p,const SSEDATAS& sdatas, const MCSTEP& mcs);
   void READCONFIG(PARAMS& p,SSEDATAS& sdatas, MATRIXELEMS& mel ,MCSTEP& mcs, const LATTICE& latt);
   void TDprint(const PARAMS,const MCSTEP&, const LATTICE&);
   void CORRprint(const PARAMS p, const MCSTEP& mcs, const LATTICE& latt);

};

#endif
