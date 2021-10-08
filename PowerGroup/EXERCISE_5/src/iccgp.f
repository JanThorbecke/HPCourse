      subroutine iccgp(nx,ny,nz,tol)
C
C
      integer nx,ny,nz
C
C     nx is the number of cells in x-direction
C     ny is the number of cells in y-direction
C     nz is the number of cells in z-direction
C
C
      double precision tol
C
C     tol is the tolerance used in the iterative Poisson solver
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
C
C
C     local variables:
      double precision alfa(0:nnz-1),beta(0:nnz-1)
      double precision normresidu(0:nnz-1),normnewresidu(0:nnz-1)
      integer i,j,m,n(0:nnz-1),nz2
      logical contin
C
C
C     This routine overwrites the array p, and the arrays
C     r, q, a and b for the common block workspace.
C     p(i,j,k), r(i,j,k), q(i,j,k), a(i,j,k and b(i,j,k)
C     are overwritten for i=0,..,nx+1,  j=0,..,ny+1 and
C     k=0,...,nz+1.
C
C
C     Before preconditioning and fourier transform the system of
C     equations for the pressure p was given by
C
C     lx(i,j,k)*p(i-1,j,k) +
C
C     ly(i,j,k)*p(i,j-1,k) +
C
C     lz(i,j,k)*p(i,j,k-1) +
C
C     d(i,j,k)*p(i,j,k) +
C
C     lx(i+1,j,k)*p(i+1,j,k) +
C
C     ly(i,j+1,k)*p(i,j+1,k) +
C
C     lz(i,j,k+1)*p(i,j,k+1) =
C
C     r(i,j,k)
C
C     where i=1,2,..,nx, j=1,2,..,ny and k=0,...,nz-1.
C     First, a fast fourier transformation in the third direction is
C     performed. The poisson equation is then solved in the spectral
C     space, and the solution is transformed back into the physical
C     space.
C     The poisson equation in the spectral space reads
C
C     lx(i,j,m)*r(i-1,j,m) +
C
C     ly(i,j,m)*r(i,j-1,m) +
C
C     (d(i,j,m)+e(i,j,m)*r(i,j,m) +
C
C     lx(i+1,j,m)*r(i+1,j,m) +
C
C     ly(i,j+1,m)*r(i,j+1,m)  =
C
C     q(i,j,m)
C
C     where i=1,2,..,nx, j=1,2,..,ny and m=0,...,nz/2.
C     Here, the right-hand-side q is complex.
C     The complex values of q are stored in a real array with the
C     same name:
C
C     q(i,j,m) = the real part of q(i,j,m)
C
C     where i=1,..,nx, j=1,..,ny, and m=0,...,nz/2.
C     and
C
C     q(i,j,m+nz/2) = the imaginary part of q(i,j,m)
C
C     where i=1,..,nx, j=1,..,ny, and m=1,...,nz/2-1.
C
C     The preconditioning overwrites the arrays lx, ly, e,
C     and d (see "preconditioner").
C
C     The explicitly preconditioned system reads
C
C     lx(i,j,k)*p(i-1,j,m) +
C
C     ly(i,j,k)*p(i,j-1,m) +
C
C     p(i,j,m) +
C
C     lx(i+1,j,k)*p(i+1,j,m) +
C
C     ly(i,j+1,k)*p(i,j+1,m) =
C
C     a(i,j,m)
C
C
C     where i=1,2,..,nx, j=1,2,..,ny and k=0,...,nz/2.
C     The right-handside a and the unknown q are related to
C     unscaled rigth-handside r and solution p by
C
C     a = q*d     and    p = r/d
C
C     respectively.
C     All arrays in the preconditioned system, except d, are
C     diagonally ordened (see "preconditioner").
C     This subroutine applies the conjugate gradient method
C     to explicitly preconditioned system.
C     It requires that nx<=ny:
C
      if (nx .gt. ny) stop
C
C
C     For convenience, the arrays p, r, lx, ly, lz, d, a, q and b
C     are extended. That is, elements of these arrays bith indices
C     that lie in the range i=0,..,nx+1, j=0,..,ny+1 are used during
C     the computations.
C     It may be emphasized that elements outside the range
C     i=1,..,nx, j=1,...,ny do not have a physical meaning.
C     The initiation of the virtual elements of lx, ly, e
C     and d has been performed in the subroutine "preconditioner"
C



      nz2 = nz/2


      do 10 m=0,nz-1

C     Initiation of the virtual elements of q and a:

      do 11 j=0,ny+1,ny+1
      do 12 i=0,nx+1
      r(i,j,m) = 0.0D0
      a(i,j,m) = 0.0D0
 12   continue
 11   continue


      do 21 i=0,nx+1,nx+1
      do 22 j=0,ny+1
      r(i,j,m) = 0.0D0
      a(i,j,m) = 0.0D0
 22   continue
 21   continue

 10   continue

      do 20 m=0,nz2


C     Diagonal ordening and scaling of the initial "solution" q
C     and the right-handside a


      do 31 j=0,nx+1
      do 32 i=0,j
      p(i,j,m) = r(i,j-i,m)/d(i,j-i,m)
      a(i,j,m) = d(i,j-i,m)*q(i,j-i,m)
 32   continue
 31   continue

      do 41 j=nx+2,ny+1
      do 42 i=0,nx+1
      p(i,j,m) = r(i,j-i,m)/d(i,j-i,m)
      a(i,j,m) = d(i,j-i,m)*q(i,j-i,m)
 42   continue
 41   continue

      do 51 j=0,nx
      do 52 i=j+1,nx+1
      p(i,j,m) = r(i,j-i+2+ny,m)/d(i,j-i+2+ny,m)
      a(i,j,m) = d(i,j-i+2+ny,m)*q(i,j-i+2+ny,m)
 52   continue
 51   continue


 20   continue

c
cpm  Added for the sake of example
c

      do m=nz2+1,nz-1
         do j = 0,ny+1
            do i = 0,nx+1
               q(i,j,m) = q(i,j,m-nz2)
            enddo
         enddo
      enddo


      do 30 m=nz2+1,nz-1

      do 61 j=0,nx+1
      do 62 i=0,j
      p(i,j,m) = r(i,j-i,m)/d(i,j-i,m-nz2)
      a(i,j,m) = q(i,j-i,m)*d(i,j-i,m-nz2)
 62   continue
 61   continue

      do 71 j=nx+2,ny+1
      do 72 i=0,nx+1
      p(i,j,m) = r(i,j-i,m)/d(i,j-i,m-nz2)
      a(i,j,m) = q(i,j-i,m)*d(i,j-i,m-nz2)
 72   continue
 71   continue

      do 81 j=0,nx
      do 82 i=j+1,nx+1
      p(i,j,m) = r(i,j-i+2+ny,m)/d(i,j-i+2+ny,m-nz2)
      a(i,j,m) = q(i,j-i+2+ny,m)*d(i,j-i+2+ny,m-nz2)
 82   continue
 81   continue


 30   continue



C     computation of the first residu (b)


      do 40 m=0,nz-1

C     first, the virtual elements of b are set equal to zero


      do 91 i=0,nx+1,nx+1
      do 92 j=0,ny+1
      b(i,j,m) = 0.0D0
 92   continue
 91   continue

 40   continue


C     the other elements of the first residue are
C     given by
C                    T
C     b := a - (L + L + 2I - E) p


      do 50 m=0,nz2

      do 101 j=1,ny
      do 102 i=1,nx
      b(i,j,m) = a(i,j,m)
     &          -ly(i,j,m)*p(i  ,j-1,m)
     &          -lx(i,j,m)*p(i-1,j-1,m)
     &          -(2.0D0-e(i,j,m))*p(i,j,m)
     &          -ly(i  ,j+1,m)*p(i,j+1,m)
     &          -lx(i+1,j+1,m)*p(i+1,j+1,m)
 102  continue
 101  continue

      j=0
      do 111 i=1,nx
      b(i,j,m) = a(i,j,m)
     &          -ly(i,j,m)*p(i  ,ny+1,m)
     &          -lx(i,j,m)*p(i-1,ny+1,m)
     &          -(2.0D0-e(i,j,m))*p(i,j,m)
     &          -ly(i  ,j+1,m)*p(i  ,j+1,m)
     &          -lx(i+1,j+1,m)*p(i+1,j+1,m)
 111  continue

      j=ny+1
      do 121 i=1,nx
      b(i,j,m) = a(i,j,m)
     &          -ly(i,j,m)*p(i  ,j-1,m)
     &          -lx(i,j,m)*p(i-1,j-1,m)
     &          -(2.0D0-e(i,j,m))*p(i,j,m)
     &          -ly(i  ,0,m)*p(i  ,0,m)
     &          -lx(i+1,0,m)*p(i+1,0,m)
 121  continue


 50   continue

      do 60 m=nz2+1,nz-1


      do 131 j=1,ny
      do 132 i=1,nx
      b(i,j,m) = a(i,j,m)
     &          -ly(i,j,m-nz2)*p(i  ,j-1,m)
     &          -lx(i,j,m-nz2)*p(i-1,j-1,m)
     &          -(2.0D0-e(i,j,m-nz2))*p(i,j,m)
     &          -ly(i  ,j+1,m-nz2)*p(i,j+1,m)
     &          -lx(i+1,j+1,m-nz2)*p(i+1,j+1,m)
 132  continue
 131  continue

      j=0
      do 141 i=1,nx
      b(i,j,m) = a(i,j,m)
     &          -ly(i,j,m-nz2)*p(i  ,ny+1,m)
     &          -lx(i,j,m-nz2)*p(i-1,ny+1,m)
     &          -(2.0D0-e(i,j,m-nz2))*p(i,j,m)
     &          -ly(i  ,j+1,m-nz2)*p(i  ,j+1,m)
     &          -lx(i+1,j+1,m-nz2)*p(i+1,j+1,m)
 141  continue

      j=ny+1
      do 151 i=1,nx
      b(i,j,m) = a(i,j,m)
     &          -ly(i,j,m-nz2)*p(i  ,j-1,m)
     &          -lx(i,j,m-nz2)*p(i-1,j-1,m)
     &          -(2.0D0-e(i,j,m-nz2))*p(i,j,m)
     &          -ly(i  ,0,m-nz2)*p(i  ,0,m)
     &          -lx(i+1,0,m-nz2)*p(i+1,0,m)
 151  continue


 60   continue

C     transformation of the first residu

C     the virtual elements of the diagonally ordened
C     residu are set equal to zero

      do 70 m=0,nz-1

      do 161 i=1,nx
      q(i,i  ,m) = 0.0D0
      q(i,i-1,m) = 0.0D0
 161  continue


      do 171 i=0,nx+1,nx+1
      do 172 j=0,ny+1
      q(i,j,m) = 0.0D0
 172  continue
 171  continue

 70   continue

C     and the other elements are given by
C
C              -1
C     q := (L+I) b

      do 80 m=0,nz2

      do 181 j=2,nx
      do 182 i=1,j-1
      q(i,j,m) = b(i,j,m)
     &          -ly(i,j,m)*q(i  ,j-1,m)
     &          -lx(i,j,m)*q(i-1,j-1,m)
 182  continue
 181  continue


      do 191 j=nx+1,ny+1
      do 192 i=1,nx
      q(i,j,m) = b(i,j,m)
     &          -ly(i,j,m)*q(i  ,j-1,m)
     &          -lx(i,j,m)*q(i-1,j-1,m)
 192  continue
 191  continue


      j=0
      do 193 i=j+2,nx
      q(i,j,m) = b(i,j,m)
     &          -ly(i,j,m)*q(i  ,ny+1,m)
     &          -lx(i,j,m)*q(i-1,ny+1,m)
 193  continue


      do 201 j=1,nx-2
      do 202 i=j+2,nx
      q(i,j,m) = b(i,j,m)
     &          -ly(i,j,m)*q(i  ,j-1,m)
     &          -lx(i,j,m)*q(i-1,j-1,m)
 202  continue
 201  continue


 80   continue


      do 90 m=nz2+1,nz-1


      do 211 j=2,nx
      do 212 i=1,j-1
      q(i,j,m) = b(i,j,m)
     &          -ly(i,j,m-nz2)*q(i  ,j-1,m)
     &          -lx(i,j,m-nz2)*q(i-1,j-1,m)
 212  continue
 211  continue


      do 221 j=nx+1,ny+1
      do 222 i=1,nx
      q(i,j,m) = b(i,j,m)
     &          -ly(i,j,m-nz2)*q(i  ,j-1,m)
     &          -lx(i,j,m-nz2)*q(i-1,j-1,m)
 222  continue
 221  continue


      j=0
      do 231 i=j+2,nx
      q(i,j,m) = b(i,j,m)
     &          -ly(i,j,m-nz2)*q(i  ,ny+1,m)
     &          -lx(i,j,m-nz2)*q(i-1,ny+1,m)
 231  continue


      do 241 j=1,nx-2
      do 242 i=j+2,nx
      q(i,j,m) = b(i,j,m)
     &          -ly(i,j,m-nz2)*q(i  ,j-1,m)
     &          -lx(i,j,m-nz2)*q(i-1,j-1,m)
 242  continue
 241  continue


 90   continue



      do 100 m=0,nz-1


C     normresidue(m) := inproduct q(.,.,m) and q(.,.,m)

      normresidu(m) = 0.0D0

      bnorm = 0.0
      anorm = 0.0
      pnorm = 0.0

      do 251 j=0,ny+1
      do 252 i=0,nx+1
      normresidu(m) = normresidu(m) + q(i,j,m)*q(i,j,m)
 252  continue
 251  continue



C     first conjugate direction r := q


      do 261 j=0,ny+1
      do 262 i=0,nx+1
      r(i,j,m) = q(i,j,m)
 262  continue
 261  continue


 100  continue


      do 116 m=0,nz-1

C     n counts the number of iterations:
      n(m)=0
 116  continue



C ***   START CG ITERATION  ****



      do 110 m=0,nz-1

      IF (m .le. nz2) THEN

 1111 continue


C     check convergence:

      if (normresidu(m) .gt. tol/((m+1.0E0)**6)) then

C     next:

      n(m) = n(m) + 1

C                     -T
C     compute a = (L+I) r

      do 271 j=nx-2,0,-1
      do 272 i=j+2,nx
      a(i,j,m) = r(i,j,m)
     &          -ly(i  ,j+1,m)*a(i  ,j+1,m)
     &          -lx(i+1,j+1,m)*a(i+1,j+1,m)
 272  continue
 271  continue

      j=ny+1
      do 281 i=1,nx
      a(i,j,m) = r(i,j,m)
     &          -ly(i  ,0,m)*a(i  ,0,m)
     &          -lx(i+1,0,m)*a(i+1,0,m)
 281  continue

      do 291 j=ny,nx+1,-1
      do 292 i=1,nx
      a(i,j,m) = r(i,j,m)
     &          -ly(i  ,j+1,m)*a(i  ,j+1,m)
     &          -lx(i+1,j+1,m)*a(i+1,j+1,m)
 292  continue
 291  continue

      do 301 j=nx,2,-1
      do 302 i=1,j-1
      a(i,j,m) = r(i,j,m)
     &          -ly(i  ,j+1,m)*a(i  ,j+1,m)
     &          -lx(i+1,j+1,m)*a(i+1,j+1,m)
 302  continue
 301  continue

C                     -1
C     compute b = (L+I) (r - Ea)


      do 311 j=2,nx
      do 312 i=1,j-1
      b(i,j,m) =  r(i,j,m)
     &          - e(i,j,m)*a(i  ,j  ,m)
     &          -ly(i,j,m)*b(i  ,j-1,m)
     &          -lx(i,j,m)*b(i-1,j-1,m)
 312  continue
 311  continue

      do 321 j=nx+1,ny+1
      do 322 i=1,nx
      b(i,j,m) =  r(i,j,m)
     &          - e(i,j,m)*a(i  ,j  ,m)
     &          -ly(i,j,m)*b(i  ,j-1,m)
     &          -lx(i,j,m)*b(i-1,j-1,m)
 322  continue
 321  continue

      j=0
      do 331 i=j+2,nx
      b(i,j,m) =  r(i,j,m)
     &          - e(i,j,m)*a(i  ,j   ,m)
     &          -ly(i,j,m)*b(i  ,ny+1,m)
     &          -lx(i,j,m)*b(i-1,ny+1,m)
 331  continue

      do 341 j=1,nx-2
      do 342 i=j+2,nx
      b(i,j,m) =  r(i,j,m)
     &          - e(i,j,m)*a(i  ,j  ,m)
     &          -ly(i,j,m)*b(i  ,j-1,m)
     &          -lx(i,j,m)*b(i-1,j-1,m)
 342  continue
 341  continue


C     compute alfa(m) = -norm residu/inproduct of r and (b + a)

      alfa(m) = 0.0D0

      do 431 j=0,ny+1
      do 432 i=0,nx+1
      alfa(m) =  alfa(m) + r(i,j,m)*(b(i,j,m)+a(i,j,m))
 432  continue
 431  continue

      alfa(m) = -normresidu(m)/alfa(m)


C     update iterate  p := p - alfa*a
C     update residu   q := q + alfa*(b+a)

      do 441 j=0,ny+1
      do 442 i=0,nx+1
      p(i,j,m) = p(i,j,m) - alfa(m)*a(i,j,m)
      q(i,j,m) = q(i,j,m) + alfa(m)*(b(i,j,m)+a(i,j,m))
 442  continue
 441  continue


C     normnewresidu(m) := inproduct of q and q

      normnewresidu(m) = 0.0D0

      do 451 j=0,ny+1
      do 452 i=0,nx+1
      normnewresidu(m) = normnewresidu(m) + q(i,j,m)*q(i,j,m)
 452  continue
 451  continue


      beta(m) = normnewresidu(m)/normresidu(m)


C     next conjugate direction  r := q + beta*r

      do 461 j=0,ny+1
      do 462 i=0,nx+1
      r(i,j,m) = q(i,j,m) + beta(m)*r(i,j,m)
 462  continue
 461  continue


C     update residu:

      normresidu(m) = normnewresidu(m)

C     test convergence:

      if (normresidu(m) .gt. tol/((m+1.0D0)**6)) then

C     stop if CG diverges:

      if (n(m) .eq. 1000) then
       stop "ERROR: ICCG diverges"
      else
       goto 1111
      endif

      endif


C     CG iteration has converged
C     scaling and re-ordening of the solution


      do 471 j=0,nx+1
      do 472 i=0,j
      r(i,j-i,m) = d(i,j-i,m)*p(i,j,m)
 472  continue
 471  continue


      do 481 j=nx+2,ny+1
      do 482 i=0,nx+1
      r(i,j-i,m) = d(i,j-i,m)*p(i,j,m)
 482  continue
 481  continue

      do 491 j=0,nx
      do 492 i=j+1,nx+1
      r(i,j-i+2+ny,m) = d(i,j-i+2+ny,m)*p(i,j,m)
 492  continue
 491  continue

      endif

      ENDIF

cpm  end first part, start second part

      IF (m .gt. nz2) THEN

 2222 continue

C     check convergence:

      if (normresidu(m) .gt. tol/(((m-nz2)+1.0D0)**6)) then

C     next:

      n(m) = n(m) + 1

C                     -T
C     compute a = (L+I) r

      do 351 j=nx-2,0,-1
      do 352 i=j+2,nx
      a(i,j,m) = r(i,j,m)
     &          -ly(i  ,j+1,m-nz2)*a(i  ,j+1,m)
     &          -lx(i+1,j+1,m-nz2)*a(i+1,j+1,m)
 352  continue
 351  continue

      j=ny+1
      do 361 i=1,nx
      a(i,j,m) = r(i,j,m)
     &          -ly(i  ,0,m-nz2)*a(i  ,0,m)
     &          -lx(i+1,0,m-nz2)*a(i+1,0,m)
 361  continue

      do 371 j=ny,nx+1,-1
      do 372 i=1,nx
      a(i,j,m) = r(i,j,m)
     &          -ly(i  ,j+1,m-nz2)*a(i  ,j+1,m)
     &          -lx(i+1,j+1,m-nz2)*a(i+1,j+1,m)
 372  continue
 371  continue

      do 381 j=nx,2,-1
      do 382 i=1,j-1
      a(i,j,m) = r(i,j,m)
     &          -ly(i  ,j+1,m-nz2)*a(i  ,j+1,m)
     &          -lx(i+1,j+1,m-nz2)*a(i+1,j+1,m)
 382  continue
 381  continue

C                     -1
C     compute b = (L+I) (r - Ea)

      do 391 j=2,nx
      do 392 i=1,j-1
      b(i,j,m) =  r(i,j,m)
     &          - e(i,j,m-nz2)*a(i  ,j  ,m)
     &          -ly(i,j,m-nz2)*b(i  ,j-1,m)
     &          -lx(i,j,m-nz2)*b(i-1,j-1,m)
 392  continue
 391  continue

      do 401 j=nx+1,ny+1
      do 402 i=1,nx
      b(i,j,m) =  r(i,j,m)
     &          - e(i,j,m-nz2)*a(i  ,j  ,m)
     &          -ly(i,j,m-nz2)*b(i  ,j-1,m)
     &          -lx(i,j,m-nz2)*b(i-1,j-1,m)
 402  continue
 401  continue

      j=0
      do 411 i=j+2,nx
      b(i,j,m) =  r(i,j,m)
     &          - e(i,j,m-nz2)*a(i  ,j   ,m)
     &          -ly(i,j,m-nz2)*b(i  ,ny+1,m)
     &          -lx(i,j,m-nz2)*b(i-1,ny+1,m)
 411  continue

      do 421 j=1,nx-2
      do 422 i=j+2,nx
      b(i,j,m) =  r(i,j,m)
     &          - e(i,j,m-nz2)*a(i  ,j  ,m)
     &          -ly(i,j,m-nz2)*b(i  ,j-1,m)
     &          -lx(i,j,m-nz2)*b(i-1,j-1,m)
 422  continue
 421  continue

C     compute alfa(m) = -norm residu/inproduct of r and (b + a)

      alfa(m) = 0.0D0

      do 433 j=0,ny+1
      do 434 i=0,nx+1
      alfa(m) =  alfa(m) + r(i,j,m)*(b(i,j,m)+a(i,j,m))
 434  continue
 433  continue

      alfa(m) = -normresidu(m)/alfa(m)

C     update iterate  p := p - alfa*a
C     update residu   q := q + alfa*(b+a)


      do 443 j=0,ny+1
      do 444 i=0,nx+1
      p(i,j,m) = p(i,j,m) - alfa(m)*a(i,j,m)
      q(i,j,m) = q(i,j,m) + alfa(m)*(b(i,j,m)+a(i,j,m))
 444  continue
 443  continue


C     normnewresidu(m) := inproduct of q and q

      normnewresidu(m) = 0.0D0

      do 453 j=0,ny+1
      do 454 i=0,nx+1
      normnewresidu(m) = normnewresidu(m) + q(i,j,m)*q(i,j,m)
 454  continue
 453  continue

      beta(m) = normnewresidu(m)/normresidu(m)

C     next conjugate direction  r := q + beta*r

      do 463 j=0,ny+1
      do 464 i=0,nx+1
      r(i,j,m) = q(i,j,m) + beta(m)*r(i,j,m)
 464  continue
 463  continue


C     update residu:

      normresidu(m) = normnewresidu(m)

C     test convergence:

      if (normresidu(m) .gt. tol/(((m-nz2)+1.0D0)**6))  then

C     stop if CG diverges:

      if (n(m) .eq. 1000) then
       stop "ERROR: ICCG diverges"
      else
       goto 2222
      endif

      endif

      endif

      ENDIF

110   continue

C ***   END CG ITERATION  ***


      do 125 m=nz2+1,nz-1

C     CG iteration has converged
C     scaling and re-ordening of the solution


      do 501 j=0,nx+1
      do 502 i=0,j
      r(i,j-i,m) = d(i,j-i,m-nz2)*p(i,j,m)
 502  continue
 501  continue


      do 511 j=nx+2,ny+1
      do 512 i=0,nx+1
      r(i,j-i,m) = d(i,j-i,m-nz2)*p(i,j,m)
 512  continue
 511  continue

      do 521 j=0,nx
      do 522 i=j+1,nx+1
      r(i,j-i+2+ny,m) = d(i,j-i+2+ny,m-nz2)*p(i,j,m)
 522  continue
 521  continue


 125  continue

      do 600 m=0,nz-1
         write(*,*) m,n(m),normresidu(m)
 600  continue

      print *, ' '
      print *, ' '


      return
      end
