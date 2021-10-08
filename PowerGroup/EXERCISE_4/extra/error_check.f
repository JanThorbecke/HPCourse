      subroutine error_check (n,m,alpha,dx,dy,u,f) 
************************************************************
* Checks error between numerical and exact solution 
*
************************************************************ 
      implicit none 
     
      integer n,m
      double precision u(n,m),f(n,m),dx,dy,alpha 
      
      integer i,j
      double precision xx,yy,temp,error 

      dx = 2.0 / (n-1)
      dy = 2.0 / (m-1)
      error = 0.0 
      do j = 1,m
         do i = 1,n
            xx = -1.0d0 + dx * dble(i-1)
            yy = -1.0d0 + dy * dble(j-1)
            temp  = u(i,j) - (1.0-xx*xx)*(1.0-yy*yy)
            error = error + temp*temp 
         enddo
      enddo
      
      error = sqrt(error)/dble(n*m)

      print *, 'Solution Error : ',error

      return 
      end 
