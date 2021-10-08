      program test_omp

      implicit none
      integer i,j, size, chunk_size
      parameter (size=240000)
      parameter (chunk_size=size/24)

      double precision x_sum, x_max, x_min, x
      double precision rsize, x_psum(24), x_pmax(24), x_pmin(24)
      double precision t1, t2, omp_get_wtime

      t1 = omp_get_wtime()
!$OMP PARALLEL DO
!$OMP&   SHARED(x_psum, x_pmax, x_pmin)
!$OMP&   PRIVATE(i, j, x)
      do 200 j=1, 24
        call random_number(x)
        x_psum(j) = x
        x_pmax(j) = x
        x_pmin(j) = x
        do 100 i=2, chunk_size
          call random_number(x)
          x_psum(j) = x_psum(j) + x
          if (x .lt. x_pmin(j)) x_pmin(j) = x
          if (x .gt. x_pmax(j)) x_pmax(j) = x
  100   continue
  200 continue

      x_sum = x_psum(1)
      x_max = x_pmax(1)
      x_min = x_pmin(1)
      do 300 j=2, 24
        x_sum = x_sum + x_psum(j)
        if (x_pmin(j) .lt. x_min) x_min = x_pmin(j)
        if (x_pmax(j) .gt. x_max) x_max = x_pmax(j)
  300 continue
      t2 = omp_get_wtime()

      rsize=dble(size)
      write(*,'(6(1pe11.3))') rsize,x_sum,x_min,x_max,x_sum/rsize,t2-t1
      end
