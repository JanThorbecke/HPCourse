program tensor

  implicit none

  double precision, allocatable :: A(:,:,:,:), B(:,:,:,:), C(:,:,:,:), D(:,:,:,:)
  integer i, j, k, l, m, n
  integer loop

  double precision err
  double precision, parameter :: tol = 1d-12
  double precision startt, endt, wallclock_time

  integer, parameter :: NLOOP = 2

  integer, parameter :: N1 = 10, M1 = 500
  integer, parameter :: N2 =  8, M2 = 454

  allocate(A(N1,N2,M1,M2), B(M1,M2,N1,N2), C(N1,N2,N1,N2), D(N1,N2,N1,N2))

  ! Initialize A and B
  do i = 1, N1
     do j = 1, N2
        do k = 1, M1
           do l = 1, M2
              A(i,j,k,l) = 1d0/(dble(i + j + k + l))
              B(k,l,i,j) = sin(dble(i + j + k + l))
           end do
        end do
     end do
  end do

  ! Calculate correct answer
  D = 0d0
  do n = 1, N2
     do m = 1, N1
        do j = 1, N2
           do i = 1, N1
              do l = 1, M2
                 do k = 1, M1
                    D(i,j,m,n) = D(i,j,m,n) + A(i,j,k,l)*B(k,l,m,n)
                 end do
              end do
           end do
        end do
     end do
  end do

  ! Timing loop
  do loop = 0, NLOOP
     if (loop == 1) startt = wallclock_time()
     C = 0d0
! BEGIN CACHE BLOCKING
     do n = 1, N2
        do m = 1, N1
           do j = 1, N2
              do i = 1, N1
                 do l = 1, M2
                    do k = 1, M1
                       C(i,j,m,n) = C(i,j,m,n) + A(i,j,k,l)*B(k,l,m,n)
                    end do
                 end do
              end do
           end do
        end do
     end do
! END CACHE BLOCKING
  end do
  endt = wallclock_time()

  print '(a,i5,a,f10.6,a)', 'Elapsed time for ', NLOOP, ' loops: ', &
	dble(endt - startt), ' seconds'

  D = D - C
  err = sqrt(dot_product( &
       reshape(D, (/ N1*N2*N1*N2 /)), &
       reshape(D, (/ N1*N2*N1*N2 /))))

  if (err < tol) then
     print *,'Check passed.'
  else
     print '(a,f20.16)','Check failed, error = ',err
  endif

  deallocate(A, B, C, D)

end program tensor
