program stencil2

  implicit none
  integer loop
  integer i, j, k, l
  integer jbegin, jend, kts, kte, ibegin, iend
  double precision startt, endt, wallclock_time

  double precision err
  double precision, parameter :: tol = 1d-12

  integer, parameter :: NLOOP = 1000

  real*4, allocatable :: c(:,:,:), a(:,:,:,:), b(:,:,:), d(:,:,:)

  ibegin = 1
  iend = 90

  kts = 1
  kte = 36

  jbegin = 1
  jend = 42

  allocate(c(ibegin:iend, kts:kte, jbegin:jend), &
       d(ibegin:iend, kts:kte, jbegin:jend), &
       a(19, ibegin-1:iend+1, kts-1:kte+1, jbegin-1:jend+1), &
       b(ibegin-1:iend+1, kts-1:kte+1, jbegin-1:jend+1))

  ! Initialize data
  do j = jbegin-1, jend+1
     do k = kts-1, kte+1
        do i = ibegin-1, iend+1
           do l = 1, 19
              a(l,i,k,j) = 3.0 - mod(l + i + j + k, 6)
           end do
           b(i,k,j) = sin(real(i,4))*cos(real(k,4))*real(j,4)
        end do
     end do
  end do

  ! Calculate correct answer
  do j = jbegin, jend
     do k = kts, kte
        do i = ibegin, iend
           d(i,k,j) = a(1 ,i,k,j)*b(i  ,k  ,j  )      &
                     +a(2 ,i,k,j)*b(i-1,k  ,j  )      &
                     +a(3 ,i,k,j)*b(i+1,k  ,j  )      &
                     +a(4 ,i,k,j)*b(i  ,k  ,j-1)      &
                     +a(5 ,i,k,j)*b(i  ,k  ,j+1)      &
                     +a(6 ,i,k,j)*b(i+1,k  ,j+1)      &
                     +a(7 ,i,k,j)*b(i+1,k  ,j-1)      &
                     +a(8 ,i,k,j)*b(i-1,k  ,j-1)      &
                     +a(9 ,i,k,j)*b(i-1,k  ,j+1)      &
                     +a(10,i,k,j)*b(i  ,k-1,j  )      &
                     +a(11,i,k,j)*b(i-1,k-1,j  )      &
                     +a(12,i,k,j)*b(i+1,k-1,j  )      &
                     +a(13,i,k,j)*b(i  ,k-1,j-1)      &
                     +a(14,i,k,j)*b(i  ,k-1,j+1)      &
                     +a(15,i,k,j)*b(i  ,k+1,j  )      &
                     +a(16,i,k,j)*b(i-1,k+1,j  )      &
                     +a(17,i,k,j)*b(i+1,k+1,j  )      &
                     +a(18,i,k,j)*b(i  ,k+1,j-1)      &
                     +a(19,i,k,j)*b(i  ,k+1,j+1)
        end do
     end do
  end do

  ! Timing loop
  do loop = 0, NLOOP
     if (loop == 1) startt = wallclock_time()
! BEGIN CACHE BLOCKING
     do j = jbegin, jend
        do k = kts, kte
           do i = ibegin, iend
              c(i,k,j) = a(1 ,i,k,j)*b(i  ,k  ,j  )      &
                        +a(2 ,i,k,j)*b(i-1,k  ,j  )      &
                        +a(3 ,i,k,j)*b(i+1,k  ,j  )      &
                        +a(4 ,i,k,j)*b(i  ,k  ,j-1)      &
                        +a(5 ,i,k,j)*b(i  ,k  ,j+1)      &
                        +a(6 ,i,k,j)*b(i+1,k  ,j+1)      &
                        +a(7 ,i,k,j)*b(i+1,k  ,j-1)      &
                        +a(8 ,i,k,j)*b(i-1,k  ,j-1)      &
                        +a(9 ,i,k,j)*b(i-1,k  ,j+1)      &
                        +a(10,i,k,j)*b(i  ,k-1,j  )      &
                        +a(11,i,k,j)*b(i-1,k-1,j  )      &
                        +a(12,i,k,j)*b(i+1,k-1,j  )      &
                        +a(13,i,k,j)*b(i  ,k-1,j-1)      &
                        +a(14,i,k,j)*b(i  ,k-1,j+1)      &
                        +a(15,i,k,j)*b(i  ,k+1,j  )      &
                        +a(16,i,k,j)*b(i-1,k+1,j  )      &
                        +a(17,i,k,j)*b(i+1,k+1,j  )      &
                        +a(18,i,k,j)*b(i  ,k+1,j-1)      &
                        +a(19,i,k,j)*b(i  ,k+1,j+1)
           end do
        end do
     end do
! END CACHE BLOCKING
  end do
  endt = wallclock_time()

  print '(a,i5,a,f10.6,a)', 'Elapsed time for ', NLOOP, ' loops: ', &
	dble(endt - startt), ' seconds'

  d = d - c
  err = 0d0
  do j = jbegin, jend
     do k = kts, kte
        do i = ibegin, iend
           err = err + d(i,k,j)**2
        end do
     end do
  end do
  err = sqrt(err)

  if (err < tol) then
     print *,'Check passed.'
  else
     print '(a,f20.16)','Check failed, error = ',err
  endif

  deallocate(c, a, b)

end program stencil2
