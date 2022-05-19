FFLAGS=
FC=ifort

all: clean compile run plot

plot:
	gnuplot plot.gnp

run:
	chmod +x prefetch
	./prefetch

compile:
	$(FC) $(FFLAGS) -o prefetch prefetch.f

clean:
	rm -f prefetch

slurm:
	rm -f job*
	sbatch run.slurm
