      subroutine jacobi (n,m,dx,dy,alpha,omega,u,f,tol,maxit)
******************************************************************
* Subroutine HelmholtzJ
* Solves poisson equation on rectangular grid assuming : 
* (1) Uniform discretization in each direction, and 
* (2) Dirichlect boundary conditions 
* 
* Jacobi method is used in this routine 
*
* Input : n,m 
*         dx,dy
*         alpha
*         omega
*         f(n,m)
*         u(n,m) 
*         tol
*         maxit 
*
* Output : u(n,m) - Solution 
*****************************************************************
      implicit none 
      integer n,m,maxit
      double precision dx,dy,f(n,m),u(n,m),alpha, tol,omega
*
* Local variables 
* 
      integer i,j,k,k_local 
      double precision error,resid,rsum,ax,ay,b
      double precision error_local, uold(n,m)

      real ta,tb,tc,td,te,ta1,ta2,tb1,tb2,tc1,tc2,td1,td2
      real te1,te2
      real second
      external second
*
* Initialize coefficients 
      ax = 1.0/(dx*dx) ! X-direction coef 
      ay = 1.0/(dy*dy) ! Y-direction coef
      b  = -2.0/(dx*dx)-2.0/(dy*dy) - alpha ! Central coeff  

      error = 10.0 * tol 
      k = 1
c$omp parallel 
c$omp& shared(k,maxit,omega,error,tol,n,m,ax,ay,b,alpha,uold,u,f,
c$omp&        ta,tb,tc,td,te)
c$omp& private(i,j,k_local,resid,error_local)

      error_local = error 
      k_local  = k

      do while (k.le.maxit .and. error.gt. tol) 

         error_local = 0.0 

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

* Accumulate residual error

               error_local = error_local + resid*resid 

            end do
         enddo

* Error check 

*  Add error from all processors 

         error = error + error_local

         k_local = k_local + 1

* Update shared data (use one proc only) 

         error = sqrt(error)/dble(n*m)
         k = k_local 

         error_local = error 
*
      enddo                     ! End iteration loop 
*
c$omp end parallel

      print *, 'Total Its : ', k-1, '      l2 Error : ',error
      if (k.eq.maxit) print *, 'Exceeded maximum iterations'

      return 
      end 
