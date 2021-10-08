      subroutine precon(nx,ny,nz,alfa)
C
C
      integer nx,ny,nz
C
C     nx is the number of cells in x-direction
C     ny is the number of cells in y-direction
C     nz is the number of cells in z-direction
C
C
      double precision alfa
C
C     alfa controls the modification of the preconditioner
C     see: Ivar Gustafson, A class of first order factorization
C     methods, BIT 18 (1978), pp. 142-156.
C
C
C
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
C     local variables:
      integer i,j,k,nz2
      double precision alpha



C     output:
C     This subroutine overwrites the arrays lx, ly, d, e
C     and the array r from the common block workspace.
C     lx(i,j,k), ly(i,j,k), d(i,j,k), e(i,j,k) and r(i,j,k)
C     are overwritten for i=0,..,nx+1, j=0,..,ny+1 and k=0,..,nz/2.


C     description:
C     This subroutine computes the incomplete Choleski preconditioner
C     for the ICCG method. The iteration itself is performed in "iccg".
C     After a fourier-transform in the third direction the system of
C     equations, say Ap=r, reads
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
C
C     Here, the righthand-side is complex for k=1,..,nz/2-1,
C     and real for k=0 and k=nz/2.
C     All arrays are extended to the range i=0,..,nx+1,
C     j=0,..,ny+1  Virtual elements, i.e.,
C     elements that fall outside the range i=1,2..,nx,
C     j=1,2,...,ny are set equal to zero (except for d).
C     The incomplete Choleski decomposition is used to compute
C     the preconditioner. The diagonal of the preconditioner
C     is modified according to Gustafson.
C     In general, the preconditioner can be represented as
C       T    -1
C     (U +D)D  (D+U),
C
C     in which D is a diagonal matrix that is determined in such
C     a way that the rowsums of the preconditioner are equal to
C     the rowsums of the matrix A (for alfa=0).
C     Here, the precondioning is made more efficient by eliminating
C     the diagonal D. To that end, the original system Ap=r is
C     re-scaled to
C
C      -0.5  -0.5  0.5      -0.5
C     D    AD    (D   p) = D    r
C
C     The preconditioner of the scaled system reads
C              T
C     (L+I)(I+L )
C
C     where the lower triangular matrix L is given by
C
C     -0.5  -0.5              -0.5  -0.5     T
C     D    AD    = L + diag (D    AD    ) + L
C
C     The arrays lx and ly are overwritten by the non-zero
C     diagonals of L             -0.5
C     Array d is overwritten by D
C     It may be emphasized that this matrix is required for the
C     scaling of the system (see "iccg").
C     For the Eisenstat implementation of CG the matrix
C
C                     -0.5  -0.5
C     E = 2I - diag (D    AD    )
C
C     is computed.
C     Explicit diagonal ordening of the unknowns is applied.
C     Yet, since d is required for the scaling only, the array d
C     is not re-ordened.
C
C     More details can be found in:
C     J.J. Dongarra, I.S. Duff, D.C. Sorensen and H.A. v.d. Vorst
C     Linear system solving on vector and shared memory computers
C     SIAM, 1991.



      nz2=nz/2


C     "computation" of the elements of d, lx, ly  and e
C     outside the range i=1,..,nx, j=1,..,ny.


      do 10 k=0,nz2

      do 11 j=0,ny+1,ny+1
      do 12 i=0,nx+1
       d(i,j,k) = 1.0D0
       e(i,j,k) = 0.0D0
      lx(i,j,k) = 0.0D0
      ly(i,j,k) = 0.0D0
 12   continue
 11   continue

      do 21 i=0,nx+1,nx+1
      do 22 j=0,ny+1
       d(i,j,k) = 1.0D0
       e(i,j,k) = 0.0D0
      lx(i,j,k) = 0.0D0
      ly(i,j,k) = 0.0D0
 22   continue
 21   continue




C     part of the diagonal of the matrix is temporally saved in r


      do 31 j=1,ny
      do 32 i=1,nx
      r(i,j,k) = d(i,j,k)
 32   continue
 31   continue





C     computation of the diagonal (d) of the preconditioner

      alpha=1.0D0+alfa/(ny*nx)



      do 41 j=1,ny
      do 42 i=1,nx
      d(i,j,k) = alpha*r(i,j,k) + e(i,j,k)
     &          -lx(i,j,k)*( lx(i,j,k) + ly(i-1,j+1,k) )/d(i-1,j  ,k)
     &          -ly(i,j,k)*( ly(i,j,k) + lx(i+1,j-1,k) )/d(i  ,j-1,k)
 42   continue
 41   continue




C     the diagonal d is overwritten by 1.0/sqrt(d)


      do 51 j=1,ny
      do 52 i=1,nx
      d(i,j,k) = 1.0D0/sqrt(d(i,j,k))
 52   continue
 51   continue




C     computation of the matrix e that is required for
C     the Eisenstat implementation of CG


      do 61 j=1,ny
      do 62 i=1,nx
      e(i,j,k) = 2.0D0 - (r(i,j,k) + e(i,j,k))*d(i,j,k)*d(i,j,k)
 62   continue
 61   continue



C     computation of the lower trangular matrix L


      do 81 j=1,ny
      do 82 i=1,nx
      lx(i,j,k) = lx(i,j,k)*d(i,j,k)*d(i-1,j  ,k)
      ly(i,j,k) = ly(i,j,k)*d(i,j,k)*d(i  ,j-1,k)
82    continue
81    continue



C     diagonal ordening of lx, ly and e




C     diagonal ordening of lx

      do 91 j=0,nx+1
      do 92 i=0,j
      r(i,j,k) = lx(i,j-i,k)
 92   continue
 91   continue

      do 101 j=nx+2,ny+1
      do 102 i=0,nx+1
      r(i,j,k) = lx(i,j-i,k)
 102  continue
 101  continue

      do 111 j=0,nx
      do 112 i=j+1,nx+1
      r(i,j,k) = lx(i,j-i+2+ny,k)
 112  continue
 111  continue

      do 121 j=0,ny+1
      do 122 i=0,nx+1
      lx(i,j,k) = r(i,j,k)
 122  continue
 121  continue




C     diagonal ordening of ly

      do 131 j=0,nx+1
      do 132 i=0,j
      r(i,j,k) = ly(i,j-i,k)
 132  continue
 131  continue

      do 141 j=nx+2,ny+1
      do 142 i=0,nx+1
      r(i,j,k) = ly(i,j-i,k)
 142  continue
 141  continue

      do 151 j=0,nx
      do 152 i=j+1,nx+1
      r(i,j,k) = ly(i,j-i+2+ny,k)
 152  continue
 151  continue

      do 161 j=0,ny+1
      do 162 i=0,nx+1
      ly(i,j,k) = r(i,j,k)
 162  continue
 161  continue




C     diagonal ordening of e

      do 171 j=0,nx+1
      do 172 i=0,j
      r(i,j,k) = e(i,j-i,k)
 172  continue
 171  continue

      do 181 j=nx+2,ny+1
      do 182 i=0,nx+1
      r(i,j,k) = e(i,j-i,k)
 182  continue
 181  continue

      do 191 j=0,nx
      do 192 i=j+1,nx+1
      r(i,j,k) = e(i,j-i+2+ny,k)
 192  continue
 191  continue

      do 201 j=0,ny+1
      do 202 i=0,nx+1
      e(i,j,k) = r(i,j,k)
 202  continue
 201  continue



 10   continue


      return
      end
