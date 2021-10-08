program stencil1

  implicit none

  double precision, allocatable :: u(:,:,:), d(:,:,:), e(:,:,:)
  double precision err
  integer i, N
  double precision, parameter :: tol = 1d-16
  integer, parameter :: NLOOP = 10
  double precision startt, endt, wallclock_time

  N = 294

  allocate(u(N,N,N), d(N,N,N), e(N,N,N))

  call init(N, N, N, u, e)

  call diff(N, N, N, u, d)

  ! Timing loop
  startt = wallclock_time()
  do i = 1, NLOOP
     call diff(N, N, N, u, d)
  end do
  endt = wallclock_time()

  print '(a,i5,a,f10.6,a)', 'Elapsed time for ', NLOOP, ' loops: ', &
	dble(endt - startt), ' seconds'

  u = d - e
  err = sqrt(dot_product( &
       reshape(u(6:N-5, 6:N-5, 6:N-5), (/ (N-10)**3 /)), &
       reshape(u(6:N-5, 6:N-5, 6:N-5), (/ (N-10)**3 /)) )) / (N**3)

  if (err < tol) then
     print *,'Check passed.'
  else
     print '(a,f20.16)','Check failed, error = ',err
  endif

  deallocate(u, d, e)

contains

  subroutine init(nx, ny, nz, u, e)

    integer, intent(in) :: nx, ny, nz
    double precision, intent(out) :: u(nx, ny, nz)
    double precision, intent(out) :: e(nx, ny, nz)

    integer i, j, k
    double precision x, y, z

    if (nx < 1 .or. ny < 1 .or. nz < 1) then
       write(6, *) 'Improper array sizes passed to init().'
       stop
    end if

    do k = 1, nz
       z = (dble(k) - 1d0)/(dble(nz) - 1d0)
       do j = 1, ny
          y = (dble(j) - 1d0)/(dble(ny) - 1d0)
          do i = 1, nx
             x = (dble(i) - 1d0)/(dble(nx) - 1d0)
             u(i, j, k) = sin(x) + y**2 + cos(z)
             e(i, j, k) = cos(x) + 2*y  - sin(z)
          end do
       end do
    end do

  end subroutine init


  subroutine diff(nx, ny, nz, u, d)
    ! Only computes interior differences

    integer, intent(in) :: nx, ny, nz
    double precision, intent(in)  :: u(nx, ny, nz)
    double precision, intent(out) :: d(nx, ny, nz)

    double precision dxi, dyi, dzi
    integer i, j, k

    if (nx < 13 .or. ny < 13 .or. nz < 13) then
       write(6, *) 'Improper array sizes passed to diff().'
       stop
    end if

    d = 0d0

    ! Sixth-order centered differences in interior

    dxi = (dble(nx) - 1d0)/3840d0
    dyi = (dble(ny) - 1d0)/3840d0
    dzi = (dble(nz) - 1d0)/3840d0

! BEGIN CACHE BLOCKING
    do k = 6, nz-5
       do j = 6, ny-5
          do i = 6, nx-5
             d(i, j, k) = d(i, j, k) + &
                  (-   9d0*u(i-5, j, k) +  125d0*u(i-3, j, k) &
                   -2250d0*u(i-1, j, k) + 2250d0*u(i+1, j, k) &
                   - 125d0*u(i+3, j, k) +    9d0*u(i+5, j, k)) * dxi + &
                  (-   9d0*u(i, j-5, k) +  125d0*u(i, j-3, k) &
                   -2250d0*u(i, j-1, k) + 2250d0*u(i, j+1, k) &
                   - 125d0*u(i, j+3, k) +    9d0*u(i, j+5, k)) * dyi + &
                  (-   9d0*u(i, j, k-5) +  125d0*u(i, j, k-3) &
                   -2250d0*u(i, j, k-1) + 2250d0*u(i, j, k+1) &
                   - 125d0*u(i, j, k+3) +    9d0*u(i, j, k+5)) * dzi
          end do
       end do
    end do
! END CACHE BLOCKING

  end subroutine diff


end program stencil1
