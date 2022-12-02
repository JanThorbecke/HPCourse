! this shows the effect prefetching.
! Adapted by the Lukas van de Wiel from the example by Intel on:
! https://www.intel.com/content/www/us/en/develop/documentation/fortran-compiler-oneapi-dev-guide-and-reference/top/language-reference/a-to-z-reference/m-to-n/mm-prefetch.html


      program preFetch
      
      implicit none
      
      integer, parameter :: n            = 10000
      integer, parameter :: maxs         = 50
      integer, parameter :: dist         = 10
      integer, parameter :: maxScanRange = 200

      real, allocatable :: a(:,:), b(:,:)
      real :: start, finish, runtime
      integer :: i, j, s, sp1
      integer :: from, to
      integer :: prefetchParameter
      
      allocate(a(n,n))
      allocate(b(n,n)) 

      a = 0.0
      b = 1.0

!----- run without prefetch

      from = maxs + dist
      to   = n - maxs - dist

      call cpu_time(start)
      
      do j = from, to
         do i = from, to
           a(i, j) = b(i-dist, j) + b(i+dist, j)
         enddo
      enddo
            
      call cpu_time(finish)
      runtime = finish-start

      open(unit=22, file="noPrefetch.dat")
      do s = 1, maxScanRange
          write(22,*) s, runtime
      enddo
      close(22)

      write(*,*) "Time without prefetch = ", finish-start," seconds."

!----- run with prefetch 0

      open(unit=40, file="prefetchLVL0.dat")

      do s = 1, maxScanRange
          sp1 = s+1

          a = 0.0
          call cpu_time(start)

          do j = from, to
              do i = from, to
                  a(i, j) = b(i-dist, j) + b(i+dist, j)
                  call mm_prefetch (a(i+s, j), 0)
                  call mm_prefetch (b(i+sp1, j), 0)
              enddo
          enddo

          call cpu_time(finish)

          write(40,*) s, finish-start

      enddo
      close(40)


!----- run with prefetch 1

      open(unit=41, file="prefetchLVL1.dat")

      do s = 1, maxScanRange
          sp1 = s+1

          a = 0.0
          call cpu_time(start)

          do j = from, to
              do i = from, to
                  a(i, j) = b(i-dist, j) + b(i+dist, j)
                  call mm_prefetch (a(i+s, j), 1)
                  call mm_prefetch (b(i+sp1, j), 1)
              enddo
          enddo

          call cpu_time(finish)

          write(41,*) s, finish-start

      enddo
      close(41)


!----- run with prefetch 2

      open(unit=42, file="prefetchLVL2.dat")

      do s = 1, maxScanRange
          sp1 = s+1

          a = 0.0
          call cpu_time(start)

          do j = from, to
              do i = from, to
                  a(i, j) = b(i-dist, j) + b(i+dist, j)
                  call mm_prefetch (a(i+s, j), 2)
                  call mm_prefetch (b(i+sp1, j), 2)
              enddo
          enddo

          call cpu_time(finish)

          write(42,*) s, finish-start

      enddo
      close(42)


!----- run with prefetch 3

      open(unit=43, file="prefetchLVL3.dat")

      do s = 1, maxScanRange
          sp1 = s+1

          a = 0.0
          call cpu_time(start)

          do j = from, to
              do i = from, to
                  a(i, j) = b(i-dist, j) + b(i+dist, j)
                  call mm_prefetch (a(i+s, j), 3)
                  call mm_prefetch (b(i+sp1, j), 3)
              enddo
          enddo

          call cpu_time(finish)

          write(43,*) s, finish-start

      enddo
      close(43)



      end program
