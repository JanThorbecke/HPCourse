      subroutine rhs(nx,ny,nz)


      integer nx,ny,nz


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

C     global variables:
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


      integer i,j,k


C     compute the right-handside of the poisson equation:


      do 20 k=1,nz
      do 21 j=1,ny
      do 22 i=1,nx
      a(i,j,k) = lx(i,j,4)*p(i-1,j,k) + lx(i+1,j,4)*p(i+1,j,k)
     &         + ly(i,j,4)*p(i,j-1,k) + ly(i,j+1,4)*p(i,j+1,k)
     &         +  d(i,j,4)*p(i,j,k)
     &         + (dx(i)*dy(j)/dz(4))*
     &                 (-p(i,j,k-1) + 2*p(i,j,k) - p(i,j,k+1))
 22   continue
 21   continue
 20   continue


      do 31 j=1,ny
      do 32 i=1,nx
      a(i,j,0) = a(i,j,nz)
 32   continue
 31   continue


      return
      end
