#calculates energy, mag, and Neel exponent
#  python bip_obs.py JQvbs_example
import os
import time
import datetime
import sys
import numpy as np

#if (len(sys.argv) < 2):
#	sys.exit("please input: (1) dir name  ex: J1J2_example")

if(len(sys.argv) == 2):
	drop=sys.argv[2]
else:
	drop='0'


# used for my job structure, don't use this here ...
#dir_name = str(sys.argv[1])
#os.system('rm '+dir_name+'/obs 2>/dev/null')
#os.system('find `pwd`/'+dir_name+' -name *L*B*hmin* >tmp')
#f1array=open('tmp', 'r').read().split()
#N1=len(f1array)
N1 = 1 # assume only one simulation in current directory
f1array = [os.getcwd()]
dir_name = os.getcwd()

#read jobs in sorted deltamin order
hvals=np.zeros(N1)
for i in range(N1):
	print(f1array[i]+'/param.dat')
	hmin=float(open(f1array[i]+'/param.dat','r').read().split()[4])
	hvals[i]=hmin
ind=np.argsort(hvals)


for j in range(N1):
	i=ind[j]
	start_time = time.time()
	# get the parameters from param.dat
	Lx=int(open(f1array[i]+'/param.dat','r').read().split()[0])
	Ly=int(open(f1array[i]+'/param.dat','r').read().split()[1])
	J=float(open(f1array[i]+'/param.dat','r').read().split()[2])
	Q=float(open(f1array[i]+'/param.dat','r').read().split()[3])
	hmin=float(open(f1array[i]+'/param.dat','r').read().split()[4])
	hmax=float(open(f1array[i]+'/param.dat','r').read().split()[5])
	Nh=int(open(f1array[i]+'/param.dat','r').read().split()[6])
	epsfac=float(open(f1array[i]+'/param.dat','r').read().split()[7])
	beta=float(open(f1array[i]+'/param.dat','r').read().split()[8])
	eql=int(open(f1array[i]+'/param.dat','r').read().split()[9])
	mcs=int(open(f1array[i]+'/param.dat','r').read().split()[10])

	
	Nsite=1.0*Lx*Ly

	#get values of h
	tmph=hmin
	hvec=np.zeros(Nh)
	for n in range(Nh):
		hvec[n]=tmph
		if Nh != 1:
			tmph=tmph*np.power(hmax/hmin, 1.0/(1.0*Nh - 1.0))



	#now copy over the data from $DATA (scratch folder)
	#os.system('cp /gpfs/data/fs71452/jondemidio/JQstagfield/'+f1array[i].split('data/')[1]+'/*.data '+f1array[i]+' 2>/dev/null;')
    
	#Nproc=48
	Nproc=1

	for c in range(Nh):
		h=hvec[c]
		eps=h*epsfac

		#gather all data
		for proc in range(Nproc//Nh):
			if proc == 0:
				#os.system('cat '+f1array[i]+'/'+str(c)+'.data > tmp_data 2>/dev/null;')
				os.system("awk 'FNR>"+drop+"' "+f1array[i]+"/"+str(c)+".data > tmp_data 2>/dev/null;")
			else:
				#os.system('cat '+f1array[i]+'/'+str(c + proc*Nh)+'.data >> tmp_data 2>/dev/null;')
				os.system("awk 'FNR>"+drop+"' "+f1array[i]+"/"+str(c + proc*Nh)+".data >> tmp_data 2>/dev/null;")



		if np.size(np.fromfile('tmp_data'))==0:
			print('skipping ')
			continue #skip this value of i if data set is empty


		data=np.loadtxt('tmp_data',ndmin=2)	
		Nbin =len(data[:,0])
		print(Nbin, " bins")
		#print((1.0*Nbin)/(16.0*15), " bins per submit")

		Nfk = 100 #number of fake data sets
		
		fkenergy=np.zeros(Nfk)
		fkmag=np.zeros(Nfk)
		fkmagsq=np.zeros(Nfk)
		fkdmag=np.zeros(Nfk)
		fkexp=np.zeros(Nfk)	
	
		for f in range(Nfk):
		
			mybins = np.random.randint(Nbin,size=Nbin)
		   
			fkenergy[f]=np.mean(data[mybins,0])
			fkmag[f]=np.mean(data[mybins,1])
			fkmagsq[f]=np.mean(data[mybins,2])
			N00 = np.mean(data[mybins,3])
			N0 = np.mean(data[mybins,4])
			mN00 = np.mean(data[mybins,5])
			mN0 = np.mean(data[mybins,6])
			
			m = fkmag[f]
			On = (2*N00)/(h/2+eps+J/2) + (N0)/(h/4+eps)
			mOn= (2*mN00)/(h/2+eps+J/2) + (mN0)/(h/4+eps)
			
			fkdmag[f] = (mOn - m*On)/4
			fkexp[f] = h*(fkdmag[f]/m)
 
	
		#the statistical errors 
		errenergy = np.sqrt(np.var(fkenergy))
		errmag = np.sqrt(np.var(fkmag))
		errmagsq = np.sqrt(np.var(fkmagsq))
		errdmag = np.sqrt(np.var(fkdmag))
		errexp = np.sqrt(np.var(fkexp))
		
		
		#the mean values
		energy = np.mean(data[:,0])
		mag = np.mean(data[:,1])
		magsq = np.mean(data[:,2])
		N00 = np.mean(data[:,3])
		N0 = np.mean(data[:,4])
		mN00 = np.mean(data[:,5])
		mN0 = np.mean(data[:,6])
		
		m = mag
		On = (2*N00)/(h/2+eps+J/2) + (N0)/(h/4+eps)
		mOn= (2*mN00)/(h/2+eps+J/2) + (mN0)/(h/4+eps)
		dmag = (mOn - m*On)/4 #the 4 here comes from using hb not h
		exp = h*(dmag/m)

	
		
		
		outfile = open(dir_name+'/obs','a')
		outfile.write('%d %.6f %.6f %.10f %.6f ' % (Lx,J,Q,h,beta))
    		outfile.write('%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f ' % (energy,errenergy,mag,errmag,magsq,errmagsq,dmag,errdmag,exp,errexp))
    		outfile.write('\n')
    		outfile.close()

	
	
	print("--- %s seconds ---" % (time.time() - start_time))
	print("Finished "+f1array[i])


#end
