      program main 
************************************************************
* Serial program to solve a Helmholtz equation 
* (d2/dx2)u + (d2/dy2)u - alpha u = f  using Jacobi method
*  
* Input :  n - grid dimension in x direction 
*          m - grid dimension in y direction
*          tol   - error tolerance for iterative solver
*          relax - Successice over relaxation parameter
*          mits  - Maximum iterations for iterative solver
*          mtemp - For message passing 
* On output 
*       : u(n,m) - Dependent variable (solutions)
*       : f(n,m) - Right hand side function 
*************************************************************
      implicit none 

      integer n,m,mits,mtemp 
      double precision tol,relax,alpha 

      common /idat/ n,m,mits,mtemp
      common /fdat/tol,alpha,relax
     
* 
* Read info 
* 
      read(5,*) n,m          ! X & Y grid dimensions
      read(5,*) alpha        ! Helmholtz constant 
      read(5,*) relax        ! Relaxation factor for solver
      read(5,*) tol          ! Iterative solver tolerance 
      read(5,*) mits         ! Maximum number of iterations for it solver
      read(5,*) mtemp        ! dummy integer variable (not used here) 

*
* Calls a driver routine 
*  This allows dynamic allocation of arrays. 
*  Currently only works for F90 compiler 
* 
      call driver () 

      stop
      end 
