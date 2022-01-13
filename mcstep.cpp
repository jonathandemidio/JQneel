//functions for classes found in mcstep.h
#include "mcstep.h"

MCSTEP::MCSTEP(const LATTICE& latt, const PARAMS& p){

  isEQ = 1;
  numWOLFF = 5;
  NLV = 2.0;
  aa = 1.2;

  eqlnumWOLFF=0.0;

  Nxyeo.resize(2,vector<vector<double> >(2,vector<double>(3,0.0)));

}//constructor

void MCSTEP::DIAGUPDATE(MTRand& ran,const LATTICE& latt, SSEDATAS& sdatas, MATRIXELEMS& mel,const PARAMS& p)
{
  
    int time;
    long newM;
    int ranplx,ranbnd; /*Random plaqs for Heat-Bath attempt - called b in A.S. ref.*/
    double wplx,wbnd,wnull,wden;
    double dice;

    /**********************************************/
    /** GO THRU TIME SLICES ONCE ALLOWING CHANGES BETWEEN: NULL JD AND QDD **/
    /**********************************************/
    for (time=0; time < sdatas.Ma; time++) 
    {
        if (sdatas.oprtr[time] >=0) // this means its "fully" diagonal including null operator which is 0
	    {    
	        if(sdatas.oprtr[time]==1)
	            ranbnd = sdatas.loc[time];
	        else
	            ranbnd = ran.randInt(latt.Nbond-1);	 
            //wbnd=p.Beta_*latt.Nbond*p.JbN_*mel.bnd(sdatas.Spin[latt.Bnd[ranbnd][0]],sdatas.Spin[latt.Bnd[ranbnd][1]]);
            wbnd=p.Beta_*latt.Nbond*mel.ham[ sdatas.Spin[latt.Bnd[ranbnd][0]] ][ sdatas.Spin[latt.Bnd[ranbnd][1]] ][ sdatas.Spin[latt.Bnd[ranbnd][0]] ][ sdatas.Spin[latt.Bnd[ranbnd][1]] ];
        
	        if(sdatas.oprtr[time]==2)
	            ranplx=sdatas.loc[time];
	        else
	            ranplx = ran.randInt(latt.Nplaq-1);

            wplx=p.Beta_*latt.Nplaq*p.QbNN_*mel.plx(sdatas.Spin[latt.Plq[ranplx][0]],sdatas.Spin[latt.Plq[ranplx][1]],sdatas.Spin[latt.Plq[ranplx][2]],sdatas.Spin[latt.Plq[ranplx][3]]);
	  

	        if(sdatas.oprtr[time]==0)
	            wnull=(double)(sdatas.Ma-sdatas.nn);
	        else
	            wnull=(double)(sdatas.Ma-sdatas.nn+1);
	  
	        wden=wbnd+wplx+wnull;
	        dice = ran.randDblExc(); 
	  
	        if(dice > ((wbnd+wplx)/wden)){// INSERT NULL

	            if(sdatas.oprtr[time] > 0) // non-null J INTERACTION
	                sdatas.nn--; //reduce nn by 1
                
	            sdatas.oprtr[time]=0; // insert blank
	            sdatas.loc[time]=-1;
	        }
	        else if(dice > (wplx/wden)){// INSERT BND

	            if(sdatas.oprtr[time] == 0)
	                sdatas.nn++; //update # of non-zero operators in string
	            

	            sdatas.oprtr[time] = 1; //JD=1
	            sdatas.loc[time]=ranbnd;
	            sdatas.bndvtx[time][0]=sdatas.Spin[latt.Bnd[ranbnd][0]];
	            sdatas.bndvtx[time][1]=sdatas.Spin[latt.Bnd[ranbnd][1]]; 
	            sdatas.bndvtx[time][4]=sdatas.Spin[latt.Bnd[ranbnd][0]];// its diagonal
	            sdatas.bndvtx[time][5]=sdatas.Spin[latt.Bnd[ranbnd][1]]; // its diagonal
	            // dont need to update bndvtx[time][2 & 3 & 6 & 7] because its a bond
	        }
	        else{ //INSERT PLX
	            if(sdatas.oprtr[time] == 0)
	                sdatas.nn++; //update # of non-zero operators in string
                
	            sdatas.oprtr[time] = 2; // Plaquette diagonal =2
	            sdatas.loc[time]=ranplx;
	            sdatas.bndvtx[time][0]=sdatas.Spin[latt.Plq[ranplx][0]];
	            sdatas.bndvtx[time][1]=sdatas.Spin[latt.Plq[ranplx][1]]; 
	            sdatas.bndvtx[time][2]=sdatas.Spin[latt.Plq[ranplx][2]];
	            sdatas.bndvtx[time][3]=sdatas.Spin[latt.Plq[ranplx][3]]; 
	    
	            sdatas.bndvtx[time][4]=sdatas.Spin[latt.Plq[ranplx][0]];// its diagonal
	            sdatas.bndvtx[time][5]=sdatas.Spin[latt.Plq[ranplx][1]]; // its diagonal
	            sdatas.bndvtx[time][6]=sdatas.Spin[latt.Plq[ranplx][2]];// its diagonal
	            sdatas.bndvtx[time][7]=sdatas.Spin[latt.Plq[ranplx][3]]; // its diagonal
	        }
	    }// END DIAG OPERATORS
        ///////// THEN  IF OFF_DIAGONAL PROPOGATE STATE
        else if(sdatas.oprtr[time]==-1) { // JO
	        sdatas.Spin[latt.Bnd[sdatas.loc[time]][0]] = sdatas.bndvtx[time][4];
	        sdatas.Spin[latt.Bnd[sdatas.loc[time]][1]] = sdatas.bndvtx[time][5];
            //HERE I WILL NEED TO CATCH BROKEN LOOPS
        }
        else if (sdatas.oprtr[time] == -2){ // Plaq Off-diagonal
	        sdatas.Spin[latt.Plq[sdatas.loc[time]][0]] =  sdatas.bndvtx[time][4];
	        sdatas.Spin[latt.Plq[sdatas.loc[time]][1]] =  sdatas.bndvtx[time][5]; 
	        sdatas.Spin[latt.Plq[sdatas.loc[time]][2]] =  sdatas.bndvtx[time][6];
	        sdatas.Spin[latt.Plq[sdatas.loc[time]][3]] =  sdatas.bndvtx[time][7];
            //HERE I WILL NEED TO CATCH BROKEN LOOPS
        }
        else
	        {cout <<"TIME TRAVEL :: SHOULD NOT END UP HERE!!!"<<endl;exit(1);}
      
      
    }//END TIME LOOP 

   
    newM = (int)( ((double)sdatas.nn) * aa);
  
    if (isEQ==1)
    {
        if (newM  > sdatas.Ma)      /*Adjust the length of the OpString array only during EQ*/     
	    sdatas.INCREASEM(newM);
    }
    else if ( ((double)sdatas.nn)>(0.95*(double)sdatas.Ma)) // if op string gets too long during simiulation (after EQ), quit!
        {cout <<"INCREASE EQUILIBRIATION TIME ... EXITING ... SHOULDNT CHANGE M after Equilibiration"<<endl;exit(1);}  
  
}//DIAGUPDATE

void MCSTEP::LINKOPERATOR(const LATTICE& latt, SSEDATAS& sdatas)
{

  //initialize OPLinkList to -1: note because of bonds, this array is NOT fully packed!
  // some legs have -1 entries: this mean they dont exist in the "real" linked list
  // when starting a loop, should check that entered vertex is != -1
  if (isEQ == 1){
    OPLinkList.clear();
    OPLinkList.resize(8*sdatas.Ma,-1);// has to be initialized at each call
  }
  else{
    OPLinkList.assign(8*sdatas.Ma,-1);// seems to be faster to initialize this way
  }

  //has to be initialized at each call to -1
  First.assign(latt.Nsite,-1);  
  Last.assign(latt.Nsite,-1);   

  LLsize=0;  //linked list size (sdatas.oprtr[time] excluding -1 elements)
  int countnn=0;
  for(int time=0;time<sdatas.Ma;time++)
    if(sdatas.oprtr[time]!=0) {

     sdatas.nnpos[countnn]=time;
     countnn++;

     /*FIRST RECORD IN THE LINKED LIST*/
      if(abs(sdatas.oprtr[time])==1){//bnd operator
        for(int site=0;site<2;site++){
          if(First[latt.Bnd[sdatas.loc[time]][site]]==-1){
            First[latt.Bnd[sdatas.loc[time]][site]]=8*time+site;
            Last[latt.Bnd[sdatas.loc[time]][site]]=8*time+site+4;/*+4 to move up in time for a bnd*/
          }
          else{/*enter the link into the linked list*/
            OPLinkList[Last[latt.Bnd[sdatas.loc[time]][site]]]=8*time+site;		  
            OPLinkList[8*time+site]=Last[latt.Bnd[sdatas.loc[time]][site]];		  
            Last[latt.Bnd[sdatas.loc[time]][site]]=8*time+site+4;
            LLsize+=2;
          }
        }//end for site 0..2  
      }//end bnd operator
      else{  // its a plaquette operator
       //cout<<"SHOULD NOT GO HERE FOR Q=0!!!"<<endl;
        for(int site=0;site<4;site++){
          if(First[latt.Plq[sdatas.loc[time]][site]]==-1){
            First[latt.Plq[sdatas.loc[time]][site]]=8*time+site;
            Last[latt.Plq[sdatas.loc[time]][site]]=8*time+site+4;/*+4 to move up in time for a plaq*/
          }
          else{/*enter into the linked list*/
            OPLinkList[Last[latt.Plq[sdatas.loc[time]][site]]]=8*time+site;
            OPLinkList[8*time+site]=Last[latt.Plq[sdatas.loc[time]][site]];		  
            Last[latt.Plq[sdatas.loc[time]][site]]=8*time+site+4;
            LLsize+=2;
          }
        }//end for site 0..4
      }//end plaquette operator
  
    }// end if its a nullop, nothing to do

  
  /*join Last and First*/
  for(int spin=0;spin<latt.Nsite;spin++){
    if (First[spin] != -1){ 
      LLsize+=2;
      OPLinkList.at(Last[spin])=First[spin];
      OPLinkList.at(First[spin])=Last[spin];
    }
  }

  
} //LINKOPERATOR

void MCSTEP::OPERATORLOOP(MTRand& ran, const LATTICE& latt, SSEDATAS& sdatas, const MATRIXELEMS& mel, const PARAMS& param) 
//, const OPLOOP& loop);
{  

    int suNnew_flav;
    int suNold_flav; // loops updatate    new_flav <--> old_flav
    int tmp_flav;
    int Nlegmax;
    int startleg;
    int time,current,lo,li;
    int loop;
    bool exitOL;
    int ranpos;
    

    Nlegmax=8*sdatas.Ma; // max number of legs 	
    countVx = 0;

    for(int countW=0; countW<numWOLFF; countW++){  //number of Wolff clusters	
        if(sdatas.nn){// avoids crashing
		    //start at a random leg such that ==> OPLinkList[startleg] != -1
            //this start should be improved, it wastes time...
            //startleg=-1;
		    //while(startleg==-1){
		    //    startleg=ran.randInt(Nlegmax-1);
		    //    if(OPLinkList[startleg]==-1) 
			//        startleg=-1;
		    //}
	
            //I don't think this is a huge time savings here...
            ranpos = ran.randInt(sdatas.nn-1);
            if(abs(sdatas.oprtr[sdatas.nnpos[ranpos]])==1)
                startleg = 8*sdatas.nnpos[ranpos] + 4*ran.randInt(1)+ran.randInt(1);
            else
                startleg = 8*sdatas.nnpos[ranpos] + ran.randInt(7);


		    current=startleg;// start at startleg
		    exitOL = false;
		    loop = 0;
            
            //determine new and old colors
            suNold_flav = sdatas.bndvtx[current/8][current%8];
		    suNnew_flav=ran.randInt(param.suN_-1); // new suN index for startleg
            while(suNnew_flav == suNold_flav)
                suNnew_flav=ran.randInt(param.suN_-1);


		    while(loop<Nlegmax && exitOL==false)
		    {
			    loop++;
			    time=current/8; // index of interaction from 0 to Ma-1     
			    li=current%8;
		

			    if(abs(sdatas.oprtr[time])==1){// then this is a bond
                    lo=vertexscatter(time,li,suNnew_flav,sdatas,mel,ran);//scatter through vertex
                    //lo=mel.SWRE[li];
                }
			    else if(abs(sdatas.oprtr[time])==2)// then this is a plaquette interaction
			        lo=mel.SWRE[li];// out-leg //xx vs yy only matters when copying into Spin values
			    else
			        {cout<<"OPERATORLOOP part :: PROGRAM IS SCREWING UP:: Linklist problem"<<endl;exit(1);}


                //update the vertices to be suNnew_flav
                sdatas.bndvtx[time][li]=suNnew_flav;

                if(lo==mel.BOUNCE[li] || lo==mel.SWCO[li]){//need to switch colors
                    tmp_flav=suNnew_flav;
                    suNnew_flav = suNold_flav;
                    suNold_flav=tmp_flav;
                }

			    sdatas.bndvtx[time][lo]=suNnew_flav;
                countVx += 2;
			    current=8*time+lo;
                if(current==startleg) {exitOL = true; break;} // NEED THIS!!

			    current=OPLinkList[current];// connect thru the LinkList
			    if(current==startleg) {exitOL = true; break;} // NEED THIS!!
		    }
	  
		    if(loop==Nlegmax) 
		        cout<<"LOOP SIZES ARE GETTING TOO BIG!!!"<<endl;
		  
	    }// if nn!=0
    } //loop over wolff loops  

    //Equilibriation of loop length
    if (isEQ == 1){
        ADJNL(countVx);
        eqlnumWOLFF+=1.0*numWOLFF;//get running average during equilibration
    }
	

}//OPERATORLOOP

int MCSTEP::randomchoice(MTRand& ran, const vector<double> probabilities)
{
    double dice;
    dice=ran.randDblExc();
    for(int ind=0; ind<probabilities.size(); ind++){
        if(dice<=probabilities[ind]){
            return ind;
        }
        dice-=probabilities[ind];
    }
    cout<<"problem with randomchoice, exiting..."<<endl;exit(1);
    return -1;

}//randomchoice


int MCSTEP::vertexscatter(int time, int li, int flav, const SSEDATAS& sdatas, const MATRIXELEMS& mel, MTRand& ran)
{
    int type;
    int choseindex;
    //this needs to be changed for general lattices!

    type=mel.meltype[ sdatas.bndvtx[time][li] ][ sdatas.bndvtx[time][mel.SWRE[li]] ][ sdatas.bndvtx[time][mel.COST[li]] ][ sdatas.bndvtx[time][mel.SWCO[li]] ];

    choseindex=randomchoice(ran, mel.vtxProbs[flav][type]);

    if(choseindex==0)
        return li;
    else if(choseindex==1)
        return mel.SWRE[li];
    else if(choseindex==2)
        return mel.COST[li];
    else if(choseindex==3)
        return mel.SWCO[li];
    else
        {cout<<"COULDN'T FIND OUTLEG!! EXITING..."<<endl; exit(1); return 0;}


}//vertexscatter


void MCSTEP::UPDATEALPHA(MTRand& ran, const LATTICE& latt, SSEDATAS& sdatas, const PARAMS& param)
{

    for(int time=0;time<sdatas.Ma;time++)
        if(sdatas.oprtr[time]!=0){//if its not null, update the operator
            if(abs(sdatas.oprtr[time])==1) // J interaction
    	    {
	            if(sdatas.bndvtx[time][0]==sdatas.bndvtx[time][4])
	                sdatas.oprtr[time]=1; // bond diagonal
	            else
	                sdatas.oprtr[time]=-1; // bond off-diagonal
	        }
            else // Q interaction
	        {
                if((sdatas.bndvtx[time][0]==sdatas.bndvtx[time][4]) && (sdatas.bndvtx[time][3]==sdatas.bndvtx[time][7]))
                    sdatas.oprtr[time]=2; // plaquette diagonal
                else
                    sdatas.oprtr[time]=-2; // plaquette off-diagonal
	        }
        }//finished updating oprtr list



    for(int site=0; site<latt.Nsite; site++){
        if(First[site]!=-1)
            sdatas.Spin[site]=sdatas.bndvtx[First[site]/8][First[site]%8];
        else
            sdatas.Spin[site]=ran.randInt(param.suN_-1); //random between 0 and suN-1
    }// finish updating Spin vector

}//UPDATEALPHA


void MCSTEP::ADJNL(const long tempvtx)
//adjusts the number of simulation loops
{
   if (tempvtx < LLsize)
     numWOLFF++;  //incease the number of loops
                                                                                       
   else if(numWOLFF> 1) /*decrease # loops*/
     numWOLFF--;
                                                                                       
}//ADJNL


void MCSTEP::MEASURE_clear(const PARAMS& p, const LATTICE& latt){

  isEQ = 0;  //equilibriation is over
  nnT =0.0;
  mag=0.0;
  magsq=0.0;
  N00=0.0;
  N0=0.0;
  mN00=0.0;
  mN0=0.0;


  //all for x and y VBS order parameters

    for(int x=0; x<2; x++)
        for(int e=0; e<2; e++)
            for(int m=0; m<3; m++)
                Nxyeo[x][e][m]=0.0;  //Nxyeo[x or y][even or odd][00, 11, 01/10] = # of this operator

}

//---------------------------------------------------------------------------

void MCSTEP::MEASURE(MTRand& ran, const PARAMS& p, const LATTICE& latt, const SSEDATAS& sdatas)
{

    nnT += 1.0*sdatas.nn;

    double mag0 = 0.0; //magnetization at time slice 0
    double dmag=0.0;    //change in magnetization
    double sumdmag=0.0; //for measuring the magnetization
    double sumdmagsq=0.0; //for measuring the squared magnetization

    double tmpN00=0.0;
    double tmpN0=0.0;



    for(int site=0; site<latt.Nsite; site++)
        mag0 -= (1.0*sdatas.Spin[site]-0.5);

    for(int time=0; time<sdatas.Ma; time++){
        //measure number of operator for exponent
        if(sdatas.oprtr[time]==1) { // J-diagonal
            if(sdatas.bndvtx[time][0]==0 && sdatas.bndvtx[time][1]==0){
                tmpN00+=1.0;

                //measure mels for Ovbs
                if(sdatas.loc[time]/latt.Nsite ==0)
                    Nxyeo[0][sdatas.loc[time]%2][0]+=1.0;
                else
                    Nxyeo[1][((sdatas.loc[time]%latt.Nsite)/latt.Lx)%2][0]+=1.0;
            }
            else{
                if(sdatas.bndvtx[time][0]==0)
                    tmpN0+=1.0;
                if(sdatas.bndvtx[time][1]==0)
                    tmpN0+=1.0;
            }

            //measure mels for Ovbs
            if(sdatas.bndvtx[time][0]==1 && sdatas.bndvtx[time][1]==1){
                if(sdatas.loc[time]/latt.Nsite ==0)
                    Nxyeo[0][sdatas.loc[time]%2][1]+=1.0;
                else
                    Nxyeo[1][((sdatas.loc[time]%latt.Nsite)/latt.Lx)%2][1]+=1.0;
            }

        }



        // measure magnetization thoughout entire imaginary time history
        if(sdatas.oprtr[time]==-1) { // J-OD
            for(int l=0; l<2; l++)
                dmag-=1.0*(sdatas.bndvtx[time][l+4]-sdatas.bndvtx[time][l]);
        
            //Measure mels for Ovbs
            if(sdatas.loc[time]/latt.Nsite==0)
                Nxyeo[0][sdatas.loc[time]%2][2]+=1.0;
            else
                Nxyeo[1][((sdatas.loc[time]%latt.Nsite)/latt.Lx)%2][2]+=1.0;
        
        }
        else if(sdatas.oprtr[time]==-2){
            for(int l=0; l<4; l++)
                dmag-=1.0*(sdatas.bndvtx[time][l+4]-sdatas.bndvtx[time][l]);
        }
        sumdmag += (mag0+dmag)/(1.0*sdatas.Ma);
        sumdmagsq += ((mag0+dmag)*(mag0+dmag))/(1.0*sdatas.Ma);
    }//end loop over time

    mag+=sumdmag;
    magsq+=sumdmagsq;

    N00+=tmpN00;
    N0+=tmpN0;

    mN00+=sumdmag*tmpN00;
    mN0+=sumdmag*tmpN0;


}// end MEASUREMENT!
