#include "RKK_determ.h"
#include "fileops.h"
#include "lattice.h"
#include "simparam.h"
#include "hamil.h"
#include "mcstep.h"
#include "MersenneTwister.h"
#include <time.h>
//#include <mpi.h>



int main(int argc, char** argv){


     /*** Global Variables***/
    //int my_rank=0;  // MPI process rank
    //int numP=1;     // MPI # of processes
    int my_rank;
    int numP;

    //MPI initialize
    //MPI_Init(&argc, &argv);
    //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //MPI_Comm_size(MPI_COMM_WORLD, &numP);
	numP=1;
	my_rank=0; //(local) serial version


    PARAMS param(my_rank);
    clock_t start = clock();
    FILEOPS outpt(my_rank);
    MTRand ran(param.SEED_+my_rank);
    LATTICE sqlatt(param);  //lattice object
    //sqlatt.print(param);
    MATRIXELEMS mel(param,sqlatt);
    SSEDATAS ssekdatas(ran,sqlatt,param);
    MCSTEP mcstep(sqlatt,param);

    int decEQL;

    if(param.EQL_<=0){
        param.EQL_=-param.EQL_;
        ifstream cin;
        cin.open(outpt.fname);
        if(!cin.fail())
            outpt.READCONFIG(param,ssekdatas,mel,mcstep,sqlatt);
        else
            {cout<<"can't find file to read from, exiting..."<<endl;exit(1);}
        cin.close();
    }


    if(param.EQL_>0){
        decEQL=param.EQL_/10;
            for(int i=0;i<param.EQL_;i++){ //equilibration
                mcstep.DIAGUPDATE(ran,sqlatt,ssekdatas,mel,param);
                    mcstep.LINKOPERATOR(sqlatt,ssekdatas);
                    mcstep.OPERATORLOOP(ran,sqlatt,ssekdatas, mel,param);
                    mcstep.UPDATEALPHA(ran, sqlatt, ssekdatas, param);
                    if((i+1)%decEQL==0)
                        cout << my_rank << ": " << (i+1)*10/decEQL << "% done equilibrating..." << endl;
            }
            mcstep.numWOLFF = int(mcstep.eqlnumWOLFF/(1.0*param.EQL_));
    }
    outpt.WRITECONFIG(param,ssekdatas,mcstep) ;


    cout << "my_rank (numP)  seed = "<<my_rank<<" ("<<numP<<") "<<endl;


    // loop over bins
    for (int j=0; j<param.nBin_; j++) //measurements
    {
        mcstep.MEASURE_clear(param,sqlatt);
        for(int i=0;i<param.MCS_;i++)
            {
                mcstep.DIAGUPDATE(ran,sqlatt,ssekdatas,mel,param);
                mcstep.LINKOPERATOR(sqlatt,ssekdatas);
                mcstep.OPERATORLOOP(ran,sqlatt,ssekdatas, mel,param);
                mcstep.UPDATEALPHA(ran,sqlatt, ssekdatas, param);
                mcstep.MEASURE(ran,param,sqlatt,ssekdatas);
            }
        outpt.WRITECONFIG(param,ssekdatas,mcstep);   //write configuration to .conf file
        outpt.TDprint(param, mcstep, sqlatt);   //write data to .data file
    }

    cout << "clock " << double(clock()-start)/CLOCKS_PER_SEC << endl;

    //MPI_Finalize();

    ///////////////////

    return 0;

}
