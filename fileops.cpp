
#include "RKK_determ.h"
#include "fileops.h"

FILEOPS::FILEOPS(const int N)
{ 
  int Ttemp=N;

  //need to establish conf file name, as well as data file name.
  sprintf(dname,"%d.data",Ttemp);
  sprintf(fname,"%d.conf",Ttemp);


}//constructor


void FILEOPS::WRITECONFIG(const PARAMS p,const SSEDATAS& sdatas, const MCSTEP& mcs)
{ 
  int i,j;
  ofstream cfout;

  cfout.open(fname);

  cfout<<mcs.numWOLFF<<"\n";
  cfout <<sdatas.Ma<<"\n"; 
  
  for (i=0 ; i<sdatas.Spin.size(); i++) 
    cfout<<sdatas.Spin[i]<<" ";
  cfout<<"\n";
  
  cfout<<"-77\n";  /*end of array tag*/

  for (i=0; i < sdatas.oprtr.size(); i++)  
    cfout<<sdatas.oprtr[i]<<" ";
  cfout<<"\n";

  cfout<<"-88\n";  /*end of array tag*/

  for (i=0; i < sdatas.loc.size(); i++)  
    cfout<<sdatas.loc[i]<<" ";
  cfout<<"\n";

  cfout<<"-99\n";  /*end of array tag*/

  for (i=0; i < sdatas.bndvtx.size(); i++)  
    {
      for(j=0;j<8;j++)
	cfout<<sdatas.bndvtx[i][j]<<" ";
      
      cfout<<"\n";
    }
  
  cfout<<"-909\n";  /*end of array tag*/
 
 
  cfout.close();
  
} /*end WRITECONFIG*/



void FILEOPS::READCONFIG(PARAMS& p,SSEDATAS& sdatas, MATRIXELEMS& mel ,MCSTEP& mcs, const LATTICE& latt)
{
  int opi, endOFf;
  int Ncount, Ncount2;

  Ncount=Ncount2=0;

  ifstream cfin;
  
  cfin.open(fname);


  cfin >> mcs.numWOLFF;
  if (mcs.numWOLFF< 1)
    cout<<"CHECK CONFIG FILE: Nl"<<"\n";

  cfin >> sdatas.Ma;
  if (sdatas.Ma< 1)
    cout<<"CHECK CONFIG FILE: Nl"<<"\n";

  ////////////////////////////////
  sdatas.Spin.clear();
  endOFf = 1;
  while(endOFf ==1){
    cfin >> opi;
    if (opi == -77)    /*EOF tag*/
      endOFf = 0;
    else sdatas.Spin.push_back(opi);
  }//end while
  if (sdatas.Spin.size() != latt.Nsite) cout<<"CONFIG Spin FUCKUP \n";
  ////////////////////////////////

  ////////////////////////////////
  sdatas.oprtr.clear();  //empty operator string
  endOFf = 1;
  while(endOFf ==1){
    cfin >> opi;
    if (opi == -88)    /*EOF tag*/
      endOFf = 0;
    else {
      sdatas.oprtr.push_back(opi);
      if (opi != 0) Ncount ++;
    } 
  }  /*end while*/
  if (sdatas.Ma != sdatas.oprtr.size())  cout<<"CONFIG FILE ERROR oprtr \n"; 
  ////////////////////////////////

  ////////////////////////////////
  sdatas.loc.clear();  //empty operator string
  endOFf = 1;
  while(endOFf ==1){
    cfin >> opi;
    if (opi == -99)    /*EOF tag*/
      endOFf = 0;
    else {
      sdatas.loc.push_back(opi);
      if (opi != -1) Ncount2++;
    } 
  }  /*end while*/
  ////////////////////////////////

  vector <int> blank;
  blank.resize(8);
  ////////////////////////////////
  sdatas.bndvtx.clear();  //empty operator string
  endOFf = 1;

  for (int i=0;i<sdatas.Ma;i++)
    {
      sdatas.bndvtx.push_back(blank);      
      for(int j=0;j<8;j++)
	cfin >>sdatas.bndvtx[i][j];
    }

  cfin >> opi;
  if (opi != -909)    /*EOF tag*/
    {cout<<"CONFIG FILE ERROR oprtr \n"; exit(1);}
  ////////////////////////////////


  if (Ncount != Ncount2)  {cout<<Ncount<<" "<<Ncount2<<"CONFIG FILE ERROR N \n"; exit(1);}

  sdatas.nn = Ncount;
  sdatas.nnpos.resize(sdatas.Ma,0);//this will be created in linkoperator, no worries

  cfin.close();

}//READCONFIG

//--------------------------------------------------------------------------------
void FILEOPS::TDprint(const PARAMS p, const MCSTEP& mcs, const LATTICE& latt)
//Print thermodynamic estimators to file.
{
  
  ofstream dout;
  dout.open(dname,ios::app);

  double energy;
  energy=-mcs.nnT/(1.0*p.MCS_*p.Beta_) + p.eps_*latt.Nbond;//shift back
  energy = energy/(1.0*latt.Nsite);


  double mag;
  mag=mcs.mag/(1.0*p.MCS_*latt.Nsite);

  double magsq;
  magsq=mcs.magsq/(1.0*p.MCS_*latt.Nsite*latt.Nsite);

  vector<double> xbonds;
  xbonds.resize(2,0.0);
  for(int b=0; b<2; b++){
    xbonds[b] = mcs.Nxyeo[0][b][0]/(1.0*p.MCS_*latt.Nsite*p.Beta_*(p.J_/2 + p.h_/2 + p.eps_))
              + mcs.Nxyeo[0][b][1]/(1.0*p.MCS_*latt.Nsite*p.Beta_*(p.J_/2 + p.eps_))
              + mcs.Nxyeo[0][b][2]/(1.0*p.MCS_*latt.Nsite*p.Beta_*(p.J_/2));
  }

  vector<double> ybonds;
  ybonds.resize(2,0.0);
  for(int b=0; b<2; b++){
    ybonds[b] = mcs.Nxyeo[1][b][0]/(1.0*p.MCS_*latt.Nsite*p.Beta_*(p.J_/2 + p.h_/2 + p.eps_))
              + mcs.Nxyeo[1][b][1]/(1.0*p.MCS_*latt.Nsite*p.Beta_*(p.J_/2 + p.eps_))
              + mcs.Nxyeo[1][b][2]/(1.0*p.MCS_*latt.Nsite*p.Beta_*(p.J_/2));
  }

  dout<<setprecision(8)<<energy<<" ";                 //energy
  dout<<setprecision(8)<<mag<<" ";                    //mag
  dout<<setprecision(8)<<magsq<<" ";                  //magsq
  dout<<setprecision(8)<<mcs.N00/(1.0*p.MCS_)<<" ";   //needed for exponent
  dout<<setprecision(8)<<mcs.N0/(1.0*p.MCS_)<<" ";    //needed for exponent
  dout<<setprecision(8)<<mcs.mN00/(1.0*p.MCS_*latt.Nsite)<<" ";   //needed for exponent
  dout<<setprecision(8)<<mcs.mN0/(1.0*p.MCS_*latt.Nsite)<<" ";    //needed for exponent
  dout<<setprecision(8)<<xbonds[0]-xbonds[1]<<" "; //Ovbs in x-direction
  dout<<setprecision(8)<<ybonds[0]-ybonds[1]<<" "; //Ovbs in y-direction


  // done
  dout<<endl;
  dout.flush();
  dout.close();

}//TDprint
