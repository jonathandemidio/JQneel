//functions for classes found in hamil.h

#include "hamil.h"

SSEDATAS::SSEDATAS(MTRand& ran,const LATTICE& latt,const PARAMS& param)
{//constructor
 

  nn = 0;  /*Operator string empty*/
  Ma = 20; /*Starting upper limit for OpString array*/
  nnpos.resize(Ma,0); //time locations of all non-nulls (only ever look at first nn elements)

  for(int i=0 ; i<latt.Nsite; i++)  //random starting spin configuration
    Spin.push_back(ran.randInt(param.suN_-1)); // randomly assign spins

  vector <int> blank;
  blank.resize(8); // eight integers that store vrtx type

  for (int i=0 ; i<Ma; i++)  /*initialize operator string to 0*/
    {
      // OLD NUMBERING      oprtr.push_back(-1);
      oprtr.push_back(0);// blanks
      loc.push_back(-1);
      bndvtx.push_back(blank);
    }

}//constructor

void SSEDATAS::INCREASEM(const int newM){
//******* Increases the expansion M
  int numFill;

  numFill = newM - Ma;  /*# of fill in operators needed*/
  
  vector <int> blank;
  blank.resize(8,-1); // eight integers that store vrtx type
  
  for(int i=0; i<numFill; i++)
    {
      // OLD NUMBERING oprtr.push_back(-1);  //pushes to end
      oprtr.push_back(0);  //pushes blank at end
      loc.push_back(-1);  //pushes to end
      bndvtx.push_back(blank);
    }
  
  Ma = newM;
  nnpos.resize(Ma,0);

  if (Ma != oprtr.size() ) {cout<<"ERROR_INCREASEM"<<endl;exit(1);}
  if (Ma != loc.size() ) {cout<<"ERROR_INCREASEM"<<endl; exit(1);}
  if (Ma != bndvtx.size() ) {cout<<"ERROR_INCREASEM"<<endl; exit(1);}
  
}//INCREASEM


void SSEDATAS::print(const LATTICE& latt)
{

    for(int i=0;i<Ma;i++)
    {
        cout<<oprtr[i]<<" "<<loc[i]<<" ";
        for (int l=0;l<8;l++)
	        cout<<bndvtx.at(i).at(l)<<" ";
        cout <<endl;
    }
    cout <<endl;
    for(int i=0;i<latt.Nsite;i++)
        cout<<Spin[i];
 
    cout << "|"<<nn<<" "<<Ma<<endl<<endl;

  
}



MATRIXELEMS::MATRIXELEMS(const PARAMS& p, const LATTICE& latt)
{ //constructor
  

    BOUNCE[0]=0;  BOUNCE[1]=1;
    BOUNCE[2]=2;  BOUNCE[3]=3;
    BOUNCE[4]=4;  BOUNCE[5]=5;
    BOUNCE[6]=6;  BOUNCE[7]=7;
    
    SWRE[0]=1;  SWRE[1]=0;
    SWRE[2]=3;  SWRE[3]=2;
    SWRE[4]=5;  SWRE[5]=4;
    SWRE[6]=7;  SWRE[7]=6;
    
    SWCO[0]=5;  SWCO[1]=4;
    SWCO[2]=7;  SWCO[3]=6;
    SWCO[4]=1;  SWCO[5]=0;
    SWCO[6]=3;  SWCO[7]=2;
    
    COST[0]=4;  COST[1]=5;
    COST[2]=6;  COST[3]=7;
    COST[4]=0;  COST[5]=1;
    COST[6]=2;  COST[7]=3;
    


    SWREbnd[0]=1;  SWREbnd[1]=0;
    SWREbnd[2]=3;  SWREbnd[3]=2;

    SWCObnd[0]=3;  SWCObnd[1]=2;
    SWCObnd[2]=1;  SWCObnd[3]=0;

    COSTbnd[0]=2;  COSTbnd[1]=3;
    COSTbnd[2]=0;  COSTbnd[3]=1;
    
   
    ///////////////   get all possible vertex update probabilites
    //get ham
    ham.resize(p.suN_,vector<vector<vector<double> > >(p.suN_,vector<vector<double> >(p.suN_,vector<double>(p.suN_,0.0))));
    ham=initializeham(p);


    //get nonzero matrix elements
    meltype.resize(p.suN_,vector<vector<vector<int> > >(p.suN_,vector<vector<int> >(p.suN_,vector<int>(p.suN_,-1))));
    Nmel=0;//count and enumerate matrix elements
    vector<int> tmpspins;
    tmpspins.resize(4,-1);
    for(int i=0; i<p.suN_; i++)
        for(int j=0; j<p.suN_; j++)
            for(int k=0; k<p.suN_; k++)
                for(int l=0; l<p.suN_; l++)//loop over all combos
                    if(ham[i][j][k][l]!=0.0){
                        meltype[i][j][k][l]=Nmel;
                        Nmel++;
                        tmpspins[0]=i;tmpspins[1]=j;tmpspins[2]=k;tmpspins[3]=l;
                        typespin.push_back(tmpspins);
                    }


    //the central object
    vtxProbs.resize(p.suN_,vector<vector<double> >(Nmel,vector<double>(4,0.0)));
    //vtxProbs[newflav][mel#] = [Pbounce,Pswre,Pcost,Pswco]
    
    //used for building: so I know what I already looked at
    visitedvtx.resize(p.suN_,vector<int>(Nmel,0));
    //visitedvtx[newflav][mel#] = 0 or 1

    //loop over mels and inleg flavors to generate all update probabilities
    for(int mel=0; mel<Nmel; mel++){
        for(int flav=0; flav<p.suN_; flav++){//loop over all inleg colors
            if(visitedvtx[flav][mel]==1 || flav==typespin[mel][0])
                continue;//matrix element already considered, or doesn't change mel
            else
                constructprobs(typespin[mel],flav); //get vtxProbs[flav][mel] & update visitedvtx[flav][mel]
        }
    }

    
    /*cout<<"let's be sure that update probabilities look reasonable..."<<endl;
    for(int mel=0; mel<Nmel; mel++){
        cout<<"___"<<endl;
        cout<<typespin[mel][2]<<" "<<typespin[mel][3]<<endl;
        cout<<typespin[mel][0]<<" "<<typespin[mel][1]<<endl;
        cout<<"___"<<endl;
        for(int flav=0; flav<p.suN_; flav++)
            if(flav!=typespin[mel][0])
                cout<<flav<<" -> "<<vtxProbs[flav][mel][0]<<","<<vtxProbs[flav][mel][1]<<","<<vtxProbs[flav][mel][2]<<","<<vtxProbs[flav][mel][3]<<endl;
        cout<<endl;
    }
    exit(1); //kill for time being
    */

}


void MATRIXELEMS::constructprobs(vector<int> svec, int newflav)
{

    int outflav;  //new color to change outleg
    vector<int> mvec; //vector holding the matrix element (s1,s2,s3,s4), svec=original, mvec=transformed
    mvec.resize(4);
    vector<double> Wvec; //weight vector for matrix elements defining a table
    Wvec.resize(4,-100.0); //always four matrix elements (four outlegs)
    vector<int> melset; //the matrix element numbers contained in the current table
    melset.resize(4,-1);
    vector<int> flavset; //the flavors that go along with melset (new color to change inleg)
    flavset.resize(4,-1);
    vector<vector<double> > Pij;//the probability table that simplex solve will spit out.    




    for(int ol=0; ol<4; ol++){//loop over possible outlegs
        for(int s=0; s<4; s++)
            mvec[s]=svec[s];
        mvec[0]=newflav; //transform inleg

        //determine how to change outleg
        if(ol==0 || ol==3)//bounce, or SWCO
            outflav=svec[0]; //dagger 
        else if(ol==1 || ol==2)//SWRE or COST
            outflav=newflav; //regular
        mvec[ol]=outflav; //transform outlet with regular or dagger operator

        //The weight of this new mel.
        Wvec[ol]=ham[mvec[0]][mvec[1]][mvec[2]][mvec[3]]; //--> All i need to solve DLE's ("A" matrix always same).

        //mark as visited and remember this mel for probabilities
        if(Wvec[ol]!=0.0){
            melset[ol]=meltype[mvec[ol]][mvec[SWREbnd[ol]]][mvec[COSTbnd[ol]]][mvec[SWCObnd[ol]]];
            if(ol==0)
                flavset[ol]=newflav; // the original process 
            else
                flavset[ol]=svec[ol];// the connected process (incoming color can be different)
            visitedvtx[flavset[ol]][melset[ol]]=1;//mark as visited
        }
    }//end loop over outlegs



    Pij=simplexsolve(Wvec); //use Wvec to solve for vtxmove probabilities, then enter them in table

    //now copy over Pij in the the vtx probabilities
    //the what was the outleg when the mels were generated is now the inleg
    for(int inleg=0; inleg<4; inleg++)
        if(melset[inleg]!=-1)
            for(int move=0; move<4; move++) //move = 0,1,2,3 = bounce, SWRE, COST, SWCO
                vtxProbs[flavset[inleg]][melset[inleg]][move]=Pij[inleg][move];


}



vector<vector<double> > MATRIXELEMS::simplexsolve(vector<double> Wvec)
{
    //******* BETTER SIMPLEX ALGORITHM ***********
    //OK I'm really in good shape now.  I just need to apply the simplex algorithm for the following tableau
    //
    //     | a12 a13 a14 a23 a24 a34 |   W
    // --------------------------------------------------------
    // a11 |  1   1   1   0   0   0  | Wvec[0]
    // a22 |  1   0   0   1   1   0  | Wvec[1]
    // a33 |  0   1   0   1   0   1  | Wvec[2]
    // a44 |  0   0   1   0   1   1  | Wvec[3]
    // ----------------------------------------------------
    // Z   | -2  -2  -2  -2  -2  -2  | - sum_i Wvec[i]


    vector<int> Avec; // keeps track of variables in top row (initially all non-bounces)
                      //(a11=1,a12=2,a13=3,a14=4,a22=5,a23=6,a24=7,a33=8,a34=9,a44=10)
    Avec.resize(6,0);
    Avec[0]=2; Avec[1]=3; Avec[2]=4; Avec[3]=6; Avec[4]=7; Avec[5]=9;

    vector<int> Yvec; // keeps track of variables in the left (initially the bounces)
    Yvec.resize(4,0);
    Yvec[0]=1; Yvec[1]=5; Yvec[2]=8; Yvec[3]=10;
    
    vector <double> tmpvec;
    tmpvec.resize(7,0.0); //to create the 7 columns (including W) of simplex tableau
    vector<vector<double> > STab; //simplex tableau, 5 rows (including Z-cost function)
    for(int k=0;k<5;k++)
        STab.push_back(tmpvec);

    //canonical form for the simplex tableau
    STab[0][0]=1.0; STab[0][1]=1.0; STab[0][2]=1.0; STab[0][3]=0.0; STab[0][4]=0.0; STab[0][5]=0.0; STab[0][6]=Wvec[0];
    STab[1][0]=1.0; STab[1][1]=0.0; STab[1][2]=0.0; STab[1][3]=1.0; STab[1][4]=1.0; STab[1][5]=0.0; STab[1][6]=Wvec[1];
    STab[2][0]=0.0; STab[2][1]=1.0; STab[2][2]=0.0; STab[2][3]=1.0; STab[2][4]=0.0; STab[2][5]=1.0; STab[2][6]=Wvec[2];
    STab[3][0]=0.0; STab[3][1]=0.0; STab[3][2]=1.0; STab[3][3]=0.0; STab[3][4]=1.0; STab[3][5]=1.0; STab[3][6]=Wvec[3];
    STab[4][0]=-2.; STab[4][1]=-2.; STab[4][2]=-2.; STab[4][3]=-2.; STab[4][4]=-2.; STab[4][5]=-2.; STab[4][6]=-Wvec[0]-Wvec[1]-Wvec[2]-Wvec[3];


    //I NEED TO THINK ABOUT WHAT IF WVEC[I] = 0, I think it's OK, I zero them out at the end


    //STEP 1: Pivot until all Z coefficients (the ones not associated with a slack variable) in the bottom row are positive.
    double tmpmin;
    int prow, pcol;
    int negcoeff=1;
    while(negcoeff==1){
        for(int l=0; l<Avec.size(); l++)
            if(STab[4][l]<0.0){
                tmpmin=10000.0;
                prow=-1; pcol=l;
                for(int k=0; k<Yvec.size(); k++)
                    if(STab[k][l]>0.0 && STab[k][6]/STab[k][l]<tmpmin){
                        tmpmin=STab[k][6]/STab[k][l];
                        prow=k;
                    }
                pivot(STab,Avec,Yvec,prow,pcol);
            }
        negcoeff=0;
        //check if there are still negative coefficients
        for(int l=0; l<Avec.size(); l++)
            if(STab[4][l]<0.0)
                negcoeff=1;
    }



    //load in the solution from the simplex algorithm
    vector<double> STsol;
    STsol.resize(10,0.0); //the values of the 10 aij's (including bounces)
    for(int k=0; k<Yvec.size(); k++)
        if(Yvec[k]>0)
            STsol[Yvec[k]-1]=STab[k][6];

    //finally, just need to create the probability tables.
    vector<vector<double> > Pij;
    tmpvec.resize(4,0.0);
    for(int k=0; k<4; k++)
        Pij.push_back(tmpvec);

    //   BOUNCE               SWRE               COST                SWCO
    Pij[0][0]=STsol[0]; Pij[0][1]=STsol[1]; Pij[0][2]=STsol[2]; Pij[0][3]=STsol[3]; //inleg = 0
    Pij[1][0]=STsol[4]; Pij[1][1]=STsol[1]; Pij[1][2]=STsol[6]; Pij[1][3]=STsol[5]; //inleg = 1
    Pij[2][0]=STsol[7]; Pij[2][1]=STsol[8]; Pij[2][2]=STsol[2]; Pij[2][3]=STsol[5]; //inleg = 2
    Pij[3][0]=STsol[9]; Pij[3][1]=STsol[8]; Pij[3][2]=STsol[6]; Pij[3][3]=STsol[3]; //inleg = 3

    //normalize to get probabilities
    for(int k=0; k<4; k++)
        for(int l=0; l<4; l++)
            if(Wvec[k]!=0.0)
                Pij[k][l]=Pij[k][l]/Wvec[k];
            else
                Pij[k][l]=0.0;



    return Pij;

}//end simplexsolve


//pivot function
int MATRIXELEMS::pivot(vector<vector<double> > &matrix, vector<int> &avec, vector<int> &yvec, int row, int col){

    int tmpval=yvec[row];
    yvec[row]=avec[col];
    avec[col]=tmpval;

    for(int i=0; i<matrix.size(); i++)
        for(int j=0; j<matrix[i].size(); j++)
            if(i!=row && j!=col)
                matrix[i][j]=matrix[i][j]-matrix[i][col]*matrix[row][j]/matrix[row][col];

    for(int i=0; i<matrix.size(); i++)
        if(i!=row && matrix[i][col]!=0.0)
            matrix[i][col]=-matrix[i][col]/matrix[row][col];

    for(int j=0; j<matrix[0].size(); j++)
        if(j!=col)
            matrix[row][j]=matrix[row][j]/matrix[row][col];

    matrix[row][col]=1.0/matrix[row][col];

    //matrix is now pivoted, just return zero

    return 0;

}


vector<vector<vector<vector<double> > > > MATRIXELEMS::initializeham(const PARAMS& p){

    vector<vector<vector<vector<double> > > >tmpham;
    tmpham.resize(p.suN_,vector<vector<vector<double> > >(p.suN_,vector<vector<double> >(p.suN_,vector<double>(p.suN_,0.0))));
    //ham[s1][s2][s3][s3] = value
    //   s3 s4
    //   =====   => matrix element
    //   s1 s2

    
    //take care of diagonal field
    for(int s1=0; s1<p.suN_; s1++)
        for(int s2=0; s2<p.suN_; s2++){
            tmpham[s1][s2][s1][s2] += p.eps_;
            if(s1==0)
                tmpham[s1][s2][s1][s2] += p.hb_;
            if(s2==0)
                tmpham[s1][s2][s1][s2] += p.hb_;
        }

    //take care of J projector
    for(int s1=0; s1<p.suN_; s1++)
        for(int s2=0; s2<p.suN_; s2++)
            tmpham[s1][s1][s2][s2] += p.JbN_;

    return tmpham;
}





double MATRIXELEMS::plx(int s1, int s2, int s3, int s4)
{

  if((s1==s2) && (s3==s4))
    return 1.0;
  else
    return 0.0;
}

double MATRIXELEMS::bnd(int s1, int s2)
{

 //not used
 if(s1==s2)
    return 1.0;
 else
   return 0.0;
}
