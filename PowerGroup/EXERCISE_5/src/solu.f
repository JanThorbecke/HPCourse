      subroutine solu(nx,ny,nz)


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
      integer i,j,k,m
      double precision  x(0:nnx+1), y(0:nny+1)
      double precision fr(0:nnz+1),fi(0:nnz+1), pi


C     computation grid points:

      x(0) = 0.0D0
      do 10 i=1,nx+1
      x(i) = x(i-1) + dx(i)
 10   continue


      y(0) = 0.0D0
      do 11 j=1,ny+1
      y(j) = y(j-1) + dy(j)
 11   continue



C     For k=0, the solution becomes:


      do 20 j=0,ny+1
      do 21 i=0,nx+1
      p(i,j,0) = x(i)*x(i)*x(i)*y(j)*y(j)*y(j)
 21   continue
 20   continue




C     The k-dependent part of the solution reads


      pi = 4.0D0*atan(1.0D0)


      do 30 k=0,nz+1
      fr(k) = 0.0D0
      fi(k) = 0.0D0
      do 31 m=0,nz/2
      fr(k) = fr(k) + cos(2.0D0*m*pi*k*dz(k))
      fi(k) = fi(k) + sin(2.0D0*m*pi*k*dz(k))
 31   continue
 30   continue



C     For solution becomes:


      do 100 k=1,nz+1
      do 101 j=0,ny+1
      do 102 i=0,nx+1
      p(i,j,k) = p(i,j,0)*(fr(k)+fi(k))/nz
 102  continue
 101  continue
 100  continue



      do 120 j=0,ny+1
      do 121 i=0,nx+1
      p(i,j, 0) = 0.0D0
      p(i,j,nz) = 0.0D0
 121  continue
 120  continue


      return
      end
