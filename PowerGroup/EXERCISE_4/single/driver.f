      subroutine driver ( ) 
*************************************************************
*
* This is where the arrays are allocated and initialzed. 
*
* Working varaibles/arrays 
*     dx  - grid spacing in x direction 
*     dy  - grid spacing in y direction 
*     
*************************************************************
      implicit none 

      integer n,m,mits,mtemp 
      double precision tol,relax,alpha 

      common /idat/ n,m,mits,mtemp
      common /fdat/tol,alpha,relax

* Allocate memory here 

      double precision u(n,m),f(n,m),dx,dy


* Initialize data

      call initialize (n,m,alpha,dx,dy,u,f)

* Solve Helmholtz equation

      call jacobi (n,m,dx,dy,alpha,relax,u,f,tol,mits)

* Check error between exact solution

      call  error_check (n,m,alpha,dx,dy,u,f)

      return 
      end 
