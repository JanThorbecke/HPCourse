      program poisson

      real t0,t1,wclock


C     parameters:
C
C
      integer nnx,nny,nnz
C
      parameter (nnx=200, nny=200, nnz=32)
C
C
C     nnx is the maximum number of computational cells in x-direction;
C     nny is the maximum number of computational cells in y-direction;
C     nnz is the maximum number of computational cells in z-direction.
C     It is recommended to put the maximum problem sizes equal
C     to the actual sizes. That is: nnx=nx, nny=ny and nnz=nz.
C
C
C
C
      integer lognnz
C
      parameter (lognnz=5)
C
C     lognnz is the 2-log of the maximum number of computational cells
C     in z-direction.
C
C

C     global variables:
C
C

      common /fftdat/ pi,ib,c,s
      double precision pi
      integer ib(0:nnz-1)
      double precision c(0:nnz-1,1:lognnz),
     &                 s(0:nnz-1,1:lognnz)


C     The arrays in this common block are initiated by the subroutine
C     inifft, and are used for the fast fourier transforms, see fft and
C     ifft.

C
C
      double precision r,q,a,b
C
      common/workspace/r(0:nnx+1,0:nny+1,0:nnz+1),
     &                 q(0:nnx+1,0:nny+1,0:nnz+1),
     &                 a(0:nnx+1,0:nny+1,0:nnz+1),
     &                 b(0:nnx+1,0:nny+1,0:nnz+1)
C
C
C     The arrays that are defined here are used as workspace.
C
C
C
      double precision d,e,lx,ly
C
      common/laplac/ d(0:nnx+1,0:nny+1,0:nnz/2),
     &               e(0:nnx+1,0:nny+1,0:nnz/2),
     &              lx(0:nnx+1,0:nny+1,0:nnz/2),
     &              ly(0:nnx+1,0:nny+1,0:nnz/2)

C
C
C
C     The arrays lx, ly, e and d are initiated by "poissn" such
C     that lx(i,j) contains the western coefficient of the central
C     discretisation of the "normalised" Laplace operator
C
C                           2     2
C     -dx(i)*dy(j)*dz(k)*( d   + d   )
C                           xx    yy
C
C     arround cell (i,j,k).
C     ly(i,j,k) contains the southern coefficient;
C     d(i,j,k) contains the central coefficient.
C     The eastern and northern coefficients are given by
C     lx(i+1,j,k), ly(i,j+1,k)  respectively,
C     i.e., the "normalised" Laplace operator is symmetric.
C     Apart from the diagonal the Laplace operator given above
C     is equal to a  Laplace operator in a spectral space.
C     This spectral space is related to the 3D physical space
C     by a fourier transform in the z-direction.
C     This transform adds a value to the diagonal of the Laplacian.
C     This value is stored in e. See "poissn".
C     The arrays lx,ly, e and d are overwritten by the preconditioner.
C     Virtual elements of d, e, lx, and ly, that are elements with
C     i=0, or i=nx+1, or j=0 or j=ny+1 are introduced
C     to simplify the procedure for solving the Poisson equation
C     for the pressure.
C     The arrays e, lx, and ly are diagonally ordened (see "precon").
C
C
C
      double precision dx,dy,dz
C
      common/mesh/dx(0:nnx+1),dy(0:nny+1),dz(0:nnz+1)
C
C
C     dx(i) is the width of the computational cell (i,j,k)
C     dy(j) is the heigth of the computational cell (i,j,k)
C     dz(k) is the depth of the computational cell (i,j,k)

C     dx(0),dx(nx+1),dy(0),dy(ny+1),dz(0) and dz(nz+1) are virtual
C     meshsizes (see "grid" for details)
C
C
C
      double precision tol,alfa
C
C
C
C     tol is the tolerance used in the iterative Poisson solver
C     alfa controls the modification of the preconditioner
C
C
C
C
      integer nx,ny,nz
C
C
C
C     nx is the number of cells in x-direction
C     ny is the number of cells in y-direction
C     nz is the number of cells in z-direction
C
C
C
C
      double precision p
C
      common/press/p(0:nnx+1,0:nny+1,0:nnz+1)
C
C
C
C     The array p contains the pressure at time-level n+1,
C     The pressure is defined at the centres of the computational
C     cells, i.e., p(i,j,k) represents the pressure at x=dx(1)+..
C     +dx(i-1)+0.5*dx(i), y=dy(1)+..+dy(j-1)+0.5*dy(j), z=dz(1)+..
C     +dz(k-1)+0.5*dz(k), where i=1,..,nx, j=1,..,ny and k=1,2,..,nz.
C     The entries that fall outside this range are introduced for
C     convenience (see "iccg"). These entries do not have a physical
C     meaning.
C


C     local variables
      integer i,j,k





C     define the grid:
      call grid(nx,ny,nz)


C     discretise the Poisson equation:
      call poissn(nx,ny,nz)


C     define the solution:
      call solu(nx,ny,nz)


C     compute the right-handside of the poisson equation:
      call rhs(nx,ny,nz)


      alfa = 1.0D0

C     compute the preconditioner:
      call precon(nx,ny,nz,alfa)

C     initialise the fast fourier transform:
      call inifft(nz)


C     copy the right-handside of the poisson equation:

      do 10 k=0,nz-1
      do 11 j=1,ny
      do 12 i=1,nx
      r(i,j,k) = a(i,j,k)
 12   continue
 11   continue
 10   continue



C     compute the fast fourier transform:
      call fft(nx,ny,nz)


C     define the initial solution

      do 20 k=0,nz-1
      do 21 j=1,ny
      do 22 i=1,nx
      r(i,j,k) = 0.0D0
 22   continue
 21   continue
 20   continue




C     perform the iccg iterations

      tol = 1.0D-06

      t0 = wclock()
      call iccgp(nx,ny,nz,tol)
      t1 = wclock() - t0
      print *, ' Iccgp done, wall time is: ', t1
      print *, ' '
      print *, ' '


C     copy the result:


      do 30 k=0,nz-1
      do 31 j=1,ny
      do 32 i=1,nx
      q(i,j,k) = r(i,j,k)
 32   continue
 31   continue
 30   continue



C     compute the inverse fast fourier transform:
      call ifft(nx,ny,nz)


C     verify the result:
      call verify(nx,ny,nz)



      stop
      end
