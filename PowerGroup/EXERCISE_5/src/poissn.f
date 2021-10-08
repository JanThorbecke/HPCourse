      subroutine poissn(nx,ny,nz)


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

C     local variables:
      integer i,j,k,nz2
      double precision pi

C     output:
C     This routine overwrittes the arrays lx, ly, d and e.
C     lx(i,j,k) is overwritten for i=1,..,nx+1,  j=1,..,ny
C     and k=0,..,nz/2.
C     ly(i,j,k) is overwritten for i=1,..,nx, j=1,..,ny+1
C     and k=0,..,nz/2.
C     d(i,j,k) is overwritten for i=1,..,nx, j=1,..,ny
C     and k=0,..,nz/2.
C     e(i,j,k) is overwritten for i=1,..,nx, j=1,..,ny
C     and for k=0,..,nz/2.


C     description:
C     lx, ly, d and e are filled with the coefficients of the
C     discretisation of the fourier transformed Laplace operator
C     in the Poisson equation for the pressure p.
C     After a fourier-transform in the third direction the system of
C     equations reads
C
C     lx(i,j,k)*p(i-1,j,k) +
C
C     ly(i,j,k)*p(i,j-1,k) +
C
C     (d(i,j,k) + e(i,j,k))*p(i,j,k) +
C
C     lx(i+1,j,k)*p(i+1,j,k) +
C
C     ly(i,j+1,k)*p(i,j+1,k) =
C
C     r(i,j,k)
C
C     where i=1,2,..,nx, j=1,2,..,ny and k=1,...,nz/2.
C     Before the fourier transform, this equation reads
C
C                  2     2     2
C     -dx*dy*dz*( d   + d   + d   )(p)  = r
C                  xx    yy    zz
C
C
C     The symmetric, positive definite "normalised" Laplace
C     operator is discretised as follows
C
C
C     lx(i,j,k)*p(i-1,j,k)+lx(i+1,j,k)*p(i+1,j,k)+
C
C     ly(i,j,k)*p(i,j-1,k)+ly(i,j+1,k)*p(i,j+1,k)+
C
C     lz(i,j,k)*p(i,j,k-1)+lz(i,j,k+1)*p(i,j,k+1)+
C
C     diag(i,j,k)*p(i,j,k),
C
C
C     where,
C
C     lx(i,j,k)=-dx(i)*dy(j)*dz(k)/(dx(i)*0.5*(dx(i-1)+dx(i)))=
C
C            =-(2.0*dy(j)*dz(k))/(dx(i)+dx(i-1))
C
C     for 2,3,..,nx,  j=1,2,..,ny and k=1,..,nz,
C
C     and the entries that fall outside the range i=1,2,..,nx
C     j=1,2,..,ny and k=1,2..,nz  are set equal to zero:
C
C     lx(1,j,k)=lx(nx+1,j,k)=0
C
C     Likewise, the coefficients of the discretisation of the second
C     derivative with respect to y are given by
C
C     ly(i,j,k)=-(2.0*dx(i)*dz(k))/(dy(j)+dy(j-1))
C
C     for j=2,3,..,ny,  i=1,2,..,nx and k=1,..,nz,
C
C     and the entries that fall outside the range i=1,2,..,nx
C     j=1,2,..,ny and k=1,..,nz are set equal to zero:
C
C     ly(i,1,k)=ly(i,ny+1,k)=0



      pi =4.0D0*atan(1.0D0)

      nz2=nz/2


      do 10 k=0,nz2


      do 21 j=1,ny
      lx(1,j,k)=0.0D0
      do 22 i=2,nx
      lx(i,j,k)=-2.0D0*dy(j)*dz(k)/(dx(i)+dx(i-1))
 22   continue
      lx(nx+1,j,k)=0.0D0
 21   continue



      do 31 i=1,nx
      ly(i,1,k)=0.0D0
      do 32 j=2,ny
      ly(i,j,k)=-(2.0D0*dx(i)*dz(k))/(dy(j)+dy(j-1))
 32   continue
      ly(i,ny+1,k)=0.0D0
 31   continue




C     After the fourier transform, the diagonal coefficients become
C     equal to the minus the sum of the off-diagonal coefficients
C     lx and ly (stored in d) and the eigenvalues e(i,j,k):


      do 41 j=1,ny
      do 42 i=1,nx
      d(i,j,k) = - ( lx(i,j,k) + lx(i+1,j  ,k)
     &             + ly(i,j,k) + ly(i  ,j+1,k) )

      e(i,j,k) = dx(i)*dy(j)*
     &           (2.0D0-2.0D0*cos(2.0D0*pi*k/(1.0D0*nz)))/dz(k)
  42  continue
  41  continue





  10  continue


      return
      end
