      program main
      double precision :: a(1024), b(1024), c(1024), f
      integer :: n, i
      double precision :: t0, t1, wallclock_time

      n = 1024;
      f = 1.0;
      do i=1,1024
          a(i) = i-1.0
          b(i) = i-1.0
          c(i) = 1.0
      end do
      t0=wallclock_time()
      do i=1,1024
          call addF(a, b, c, f, n)
      end do
      t1=wallclock_time()

      write(*,*),'F-function', t1-t0,'seconds'
      write(*,*),'c(1)=', c(1),'c(2)=', c(2)

      do i=1,1024
          c(i) = 1.0
      end do
      call addF(a, c(1), c(2), f, n-1)
      write(*,*),'c(1)=', c(1),'c(2)=', c(2)
      end

      subroutine addF(a, b, c, f, n)
      double precision a(n), b(n), c(n), f
      integer n
      integer k

      do k=1,n
          c(k) = a(k) + f*b(k)
      end do

      end subroutine



