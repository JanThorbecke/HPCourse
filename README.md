# HPCourse
Exercises for the course 'Introduction to High Performance Computing'

Exercises:
+ Hardware:
Cycles
Load Store
Memory Latency

+ Operating System and tools:
IO
Floating point precision
Elastic Finite Difference, with some mistakes as exercise If you are looking for the correct latest version of the code click here.

+ Single core Optimization:
Effect of loop order on performance
Example of loop vectorization
Classic Matrix Multiplication
Classic Matrix Multiplication Loop Unrolling
Classic Matrix Multiplication cache blocking
1D Convolution
2D Convolution
Correlation
Three cache blocking exercises

+ Parallelization:
Compute min/max/sum parallel with OpenMP and try placement
OpenMP loop and cache trashing
OpenMP loop scheduling and load balancing
Compute PI parallel with OpenMP and MPI
Examples of MPI master-slave implementations
OpenMP exercises
MPI exercise to set up MPI-tasks for 2D domain decomposition

Note: for detailed information about the CPU in your system look at /proc/cpuinfo on your linux system or type: 'dmesg | grep CPU'. On windows there is a utility cpu-z which you can download and use to get information about your CPU.

In each directory you can find a file called README with detailed instructions how to compile and run the exercise.

