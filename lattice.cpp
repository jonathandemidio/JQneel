#include "lattice.h"

LATTICE::LATTICE(const PARAMS& param){ //constructor

    //---------------------- Rectangular Lattice -------------------------      
    Nsite = param.L1_*param.L2_;
    Lx=param.L1_;
    Ly=param.L2_;
		
    Nbond=2*Nsite;// total number of bonds
    Bnd.resize(Nbond,vector<int>(2)); // each bond has exactly two sites

    int bond=0;
    //regular bonds
    for (int j=0; j< Ly; j++) 
        for (int i=0; i< Lx; i++) 
        {       
            Bnd.at(bond).at(0)=i+Lx*j;
            Bnd.at(bond).at(1)=(i+1)%Lx + Lx*j;// move +x
            bond++;
        }
		

	for (int j=0; j< Ly; j++) 
	    for (int i=0; i< Lx; i++) 
		{
			Bnd.at(bond).at(0)=i+Lx*j;
			Bnd.at(bond).at(1)=  i + Lx*((j+1)%Ly);// move +y
			bond++;
		}
      
    if (bond != Nbond)
      {cout <<"lattice: bug in code ... exiting"<<endl;exit(1);}

      //----------------------- PLAQS -----------------------------

     
    Nplaq = 2*Nsite;// total number of plaqs
    Plq.resize(Nplaq,vector<int>(4));  //redefine the bond/spin array size
    
    vector<vector<int> > blank;
    blank.resize(4,vector<int>(3));
    
    
    int plaq=0;
    // Q-interaction
    for (int j=0; j< Ly; j++) 
        for (int i=0; i< Lx; i++) 
	    { 
	        //xx
	        Plq.at(plaq).at(0)=i+Lx*j;
	        Plq.at(plaq).at(1)=(i+1)%Lx + Lx*j;// move +x
	        Plq.at(plaq).at(2)= i+ Lx*((j+1)%Ly);// move +y
	        Plq.at(plaq).at(3)=(i+1)%Lx + Lx*((j+1)%Ly); // move +x+y 
	        plaq++;
	 
        }

    for (int j=0; j< Ly; j++)
        for (int i=0; i< Lx; i++)
        {
	        //yy
	        Plq.at(plaq).at(0)=i+Lx*j;
	        Plq.at(plaq).at(2)=(i+1)%Lx + Lx*j;// move +x
	        Plq.at(plaq).at(1)= i+ Lx*((j+1)%Ly);// move +y
	        Plq.at(plaq).at(3)=(i+1)%Lx + Lx*((j+1)%Ly); // move +x+y 
	        plaq++;
	    }

   if (plaq != Nplaq)
      {cout <<"lattice: bug in code ... exiting"<<endl;exit(1);}
    
    
  cout <<"End of lattice.cpp constructor"<<endl;
  cout <<"Nsite = "<<Nsite<<endl;
  cout <<"Nbond = "<<Nbond<<"; Nplaq = "<<Nplaq<<endl;
  
}//constructor

void LATTICE::print(){ 

  cout << "bonds"<<endl;
  for (int bond=0;bond< Nbond;bond++) {       
    cout<<bond<<" : ";
    for(int i=0;i<Bnd.at(bond).size();i++)
      cout<<Bnd.at(bond).at(i)<<" ";
    
  }//bonds


  cout << "plaquettes"<<endl;
  for (int plaq=0; plaq < Nplaq; plaq++) {       
    cout<<plaq<<" : ";
    for (int i =0;i<Plq.at(plaq).size();i++)
      cout<<Plq.at(plaq).at(i)<<" ";

  }//plaq


}//print
