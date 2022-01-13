OBJS = sseJQ_determ.o fileops.o lattice.o hamil.o mcstep.o
CC = g++
#CC = mpiicc
CFLAGS = -O2 #-xHost
LIBS = -lm

a.out: $(OBJS)
	$(CC) $(OBJS) -o a.out $(LIBS)

sseJQ_determ_mpi.o: sseJQ_determ.cpp RKK_determ.h lattice.h fileops.h hamil.h MersenneTwister.h mcstep.h simparam.h
	$(CC) $(CFLAGS) -c sseJQ_determ.cpp

lattice.o: lattice.cpp lattice.h RKK_determ.h simparam.h
	$(CC) $(CFLAGS) -c lattice.cpp

hamil.o: hamil.cpp hamil.h lattice.h MersenneTwister.h RKK_determ.h simparam.h
	$(CC) $(CFLAGS) -c hamil.cpp

mcstep.o: mcstep.cpp mcstep.h hamil.h lattice.h MersenneTwister.h RKK_determ.h simparam.h
	$(CC) $(CFLAGS) -c mcstep.cpp

fileops.o: fileops.cpp fileops.h RKK_determ.h hamil.h mcstep.h
	$(CC) $(CFLAGS) -c fileops.cpp

clean :
	rm *.o

