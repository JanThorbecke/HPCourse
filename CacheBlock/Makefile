# module swap PrgEnv-pgi PrgEnv-cray

F90 = gfortran
FFLAGS = -O3

LD = $(F90)
LDFLAGS = 

TIMELIBDIR = ./C_time_lib
TIMELIB = lib_timer.a

PRGS = stencil1 stencil2 tensor

all : $(PRGS)

$(TIMELIBDIR)/$(TIMELIB) : 
	make -C $(TIMELIBDIR)

% : %.o $(TIMELIBDIR)/$(TIMELIB)
	$(LD) $(LDFLAGS) -o $@ $^

%.o : %.f90
	$(F90) $(FFLAGS) -c -o $@ $<

clean :
	rm -f *.lst *.o $(PRGS)

distclean : clean
	make -C $(TIMELIBDIR) clean

