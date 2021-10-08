      subroutine jacobi (n,m,dx,dy,alpha,omega,u,f,tol,maxit)
******************************************************************
* Solves poisson equation on rectangular grid assuming : 
* (1) Uniform discretization in each direction, and 
* (2) Dirichlect boundary conditions 
* 
* Jacobi method is used in this routine 
*
* Input : n,m    number of grid points in x & y direction respectivly
*         dx,dy  Grid spacing in the x and y direction respectivly 
*         alpha  Helmholtz constant 
*         omega  Relaxation factor 
*         f(n,m) "Right hand side" of PDE 
*         u(n,m) Solution 
*         tol    Convergence tolerance fr l2 norm of residual
*         maxit  Maximum number of iterations for Jacobi Method 
*
* Output : u(n,m) - Solution 
*****************************************************************
      implicit none 
      integer n,m,maxit
      double precision dx,dy,f(n,m),u(n,m),alpha, tol,omega
*
* Local variables 
* 
      integer i,j,k
      double precision error,resid,rsum,ax,ay,b
      double precision uold(n,m)
*
* Initialize coefficients 
      ax = 1.0/(dx*dx) ! X-direction coef 
      ay = 1.0/(dy*dy) ! Y-direction vorf
      b  = -2.0/(dx*dx)-2.0/(dy*dy) - alpha ! Central coeff  


      k = 1 
      error = 10.0 * tol 

      do while (k.le.maxit .and. error.gt. tol) 
         error = 0.0 

* Copy new solution into old

         do j=1,m
            do i=1,n
               uold(i,j) = u(i,j) 
            enddo
         enddo

* Compute stencil, residual, & update 

         do j = 2,m-1
            do i = 2,n-1 

*     Evaluate residual 

               resid = (ax*(uold(i-1,j) + uold(i+1,j)) 
     &                + ay*(uold(i,j-1) + uold(i,j+1))
     &                 + b * uold(i,j) - f(i,j))/b 

* Update solution 

               u(i,j) = uold(i,j) - omega * resid

* Error accumulation 

               error = error + resid*resid 
            end do
         enddo

* Error check 

         error = sqrt(error)/dble(n*m) 
         k = k + 1
*
      enddo                     ! End iteration loop 

* Print information 

      print *, 'Total Its : ', k, '      l2 Error : ',error
      if (k.eq.maxit) print *, 'Exceeded maximum iterations'

      return 
      end 
