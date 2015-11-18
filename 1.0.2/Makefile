objects = qm.o
FC = gfortran
switch = -O3
EXEC = qm.x
#LIB = /usr/lib64/atlas/liblapack.so.3.0
#LIB = -llapack  
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)
fft.x: fft.o 
	$(FC) -o fft.x $(LIB) $(switch) fft.o 

%.o: %.f
	$(FC) -c  $<
%.o: %.f90
	$(FC) -c  $<

#	$(FC) -c $<

clean:
	rm $(objects)
veryclean:
	rm *.o
