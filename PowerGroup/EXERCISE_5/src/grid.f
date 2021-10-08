      subroutine grid(nx,ny,nz)


      integer nx,ny,nz


C     nx is the number of cells in x-direction
C     ny is the number of cells in y-direction
C     nz is the number of cells in z-direction


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


C     local variables:
      integer i,j,k
      double precision d,a,b,nx2,ny2


C     output:
C     This subroutine initiates the arrays dx, dy and dz.
C     dx(i) is overwritten for i=0,..,nx+1
C     dy(j) is overwritten for j=0,..,ny+1
C     dz(k) is overwritten for k=0,..,nz+1


C     description:
C     This subroutine defines the computational cells of a
C     stretched, staggered grid.
C     The physical coordinates of the centre of cell (i,j,k) are
C
C
C     x=dx(1)+dx(2)+dx(3)+...+dx(i-1)+0.5dx(i)
C
C     y=dy(1)+dy(2)+dy(3)+...+dy(j-1)+0.5dy(i)
C
C     z=dz(1)+dz(2)+dz(3)+...+dz(k-1)+0.5dz(k)
C
C
C     where i=1,2,..,nx, j=1,2,..,ny and k=1,..,nz.
C     Or, stated otherwise,
C     dx(i) is the width of the computational cell (i,j,k)
C     dy(j) is the heigth of the computational cell (i,j,k)
C     dz(k) is the depth of the computational cell (i,j,k)
C     For convience, the following virtual meshsizes are introduced:
C     dx(0), dx(nx+1), dy(0), dy(ny+1), dz(0) and dz(nz+1)






      nx = nnx
      ny = nny
      nz = nnz



C     The grid is stretched in the streamwise direction:


      nx2 = nx/2


      b=3.0E0
      a=0.5E0/SINH(0.5E0*b)
      d=1.0E0/(1.0E0*nx)

      do 11 i=1,nx2
      dx(i) = a*SINH(b*i*d)-a*SINH(b*(i-1)*d)
 11   continue


      do 12 i=nx2+1,nx
      dx(i) = dx(nx-i+1)
 12   continue


C     The grid is stretched in the wall-normal direction:


      ny2 = ny/2

      b=4.4678E0
      a=0.5E0/SINH(0.5E0*b)
      d=1.0E0/(1.0E0*ny)

      do 21 j=1,ny2
      dy(j) = a*SINH(b*j*d)-a*SINH(b*(j-1)*d)
 21   continue


      do 22 j=ny2+1,ny
      dy(j) = dy(ny-j+1)
 22   continue



C     The grid is uniform in the spanwise direction



      d=4.0E0/(1.0E0*nz)

      do 30 k=1,nz
      dz(k) = d
 30   continue







C     and, the virtual elements of dx, dy and dz are

      dx(0   ) = dx(1 )
      dx(nx+1) = dx(nx)

      dy(0   ) = dy(1)
      dy(ny+1) = dy(ny)

      dz(0   ) = dz(1)
      dz(nz+1) = dz(nz)


      return
      end
