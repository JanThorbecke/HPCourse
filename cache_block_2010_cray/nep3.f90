program tensor

  implicit none

  double precision, allocatable :: A(:,:,:,:), B(:,:,:,:), C(:,:,:,:), D(:,:,:,:)
  integer i, j, k, l, m, n
  integer loop

  double precision err
  double precision, parameter :: tol = 1d-12

  integer, parameter :: NLOOP = 2

  integer, parameter :: N1 = 10, M1 = 500
  integer, parameter :: N2 =  8, M2 = 454

  integer(kind=8) startt, endt
  double precision eclksec

  integer(kind=8) irtc, irtc_rate
  external irtc, irtc_rate

  eclksec = dble(irtc_rate())

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
              do l = 1, M2
                 do k = 1, M1
        do j = 1, N2
           do i = 1, N1
                    D(i,j,m,n) = D(i,j,m,n) + A(i,j,k,l)*B(k,l,m,n)
                 end do
              end do
           end do
        end do
     end do
  end do

  ! Timing loop
  do loop = 0, NLOOP
     if (loop == 1) startt = irtc()
     C = 0d0
! BEGIN CACHE BLOCKING
! dir$ blockable(n,m,l,k,j,i)
!dir$ blockingsize(555)
   do l = 1, M2
!dir$ blockingsize(444)
      do k = 1, M1
!dir$ noblocking
     do n = 1, N2
!dir$ noblocking
        do m = 1, N1
! dir$ unroll(8)
           do j = 1, N2
! dir$ unroll(10)
              do i = 1, N1
                       C(i,j,m,n) = C(i,j,m,n) + A(i,j,k,l)*B(k,l,m,n)
                    end do
                 end do
              end do
           end do
        end do
     end do
! END CACHE BLOCKING
  end do
  endt = irtc()

  print '(a,i5,a,f10.6,a)', 'Elapsed time for ', NLOOP, ' loops: ', &
    dble(endt - startt) / eclksec, ' seconds'

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
