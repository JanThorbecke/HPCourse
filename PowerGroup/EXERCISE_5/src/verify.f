      subroutine verify(nx,ny,nz)


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



C     local variables
      integer i,j,k
      double precision diff,diff0


C     compute the solution:
      call solu(nx,ny,nz)


C     compute the difference with r:

      diff=0.0D0


      diff0 = p(1,1,0) - r(1,1,0)

      do 80 k=0,nz-1
      do 81 j=1,ny
      do 82 i=1,nx
      if ( abs( p(i,j,k)-r(i,j,k)-diff0) .gt. diff) then
      diff = abs(p(i,j,k)-r(i,j,k)-diff0)
      endif
 82   continue
 81   continue
 80   continue

cp      write(*,*) diff



      return
      end
