      subroutine fft(nx,ny,nz)

      integer nx,ny,nz

C     last modification:
C     3 May 1996


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




C     local variables:
      double precision vr,vi
      integer i,j,k,m,gamma,ny2,nz2,im,immax,iz
      double precision cc,ss



C     description:
C     this routine computes the fast fourier transform in the third
C     direction of the real array r.
C     the fourier-transform fr is computed according to
C
C     fr(i,j,m) = sum over k=0,..,nz-1 of r(i,j,k)*exp(2*pi*I*k*m/nz)
C
C     where i=1,..,nx, j=1,..,ny, and m=0,...,nz-1.
C     note: I*I=-1, and pi=3.1415...
C     since r(i,j,k) is real, it holds that
C
C     fr(i,j,m) = complex conjugate of fr(i,j,nz-m)
C
C     thus values for m > nz/2 need not to be computed.
C     note that fr(i,j,0) and fr(i,j,nz/2) are real.
C     the result fr(i,j,m) is stored for m=0,..,nz/2 in the array q
C     that is:
C
C     q(i,j,m) = the real part of fr(i,j,m)
C
C     where i=1,..,nx, j=1,..,ny, and m=0,...,nz/2.
C     and
C
C     q(i,j,m+nz/2) = the imaginary part of fr(i,j,m)
C
C     where i=1,..,nx, j=1,..,ny, and m=1,...,nz/2-1.
C     the subroutine inifft must be runned before fft, in order
C     to initialize the variables in the common block fftdat.


      ny2 = ny/2


C     in this routine ny must be even
C     this is checked for first

      if ( ny2 + ny2 .ne. ny ) then
      stop "ny must be even"
      endif


C     start of the fast fourier transform


      gamma = nint(dlog(1.0d0*nz)/dlog(2.0d0))



      nz2 = nz/2
      m=0
      immax=1


      do 10 k=1,gamma


C     to save work, two real fft's are computed simultaneously
C     therefore, r(i,j,m) is packed into the real part of the
C     series that is to be transformed and
C     r(i,j+ny2,m) is packed into the imaginary part


      do 11 im=1,immax


      do 12 iz=1,nz2

      cc = c(m,k)
      ss = s(m,k)

CMIC@    DO ALL SHARED(NZ2,NY2,NX,M,K,CC,R,SS)PRIVATE(IZ,J,I,VR,VI)
      do 13 j=1,ny2
      do 14 i=1,nx

      vr  = cc*r(i,j    ,m+nz2) + ss*r(i,j+ny2,m+nz2)
      vi  = cc*r(i,j+ny2,m+nz2) - ss*r(i,j    ,m+nz2)

      r(i,j    ,m+nz2) = r(i,j    ,m) - vr
      r(i,j+ny2,m+nz2) = r(i,j+ny2,m) - vi
      r(i,j    ,m    ) = r(i,j,    m) + vr
      r(i,j+ny2,m    ) = r(i,j+ny2,m) + vi

 14   continue
 13   continue

      m = m + 1

 12   continue



      m = m + nz2

 11   continue


      nz2 = nz2/2
      m = 0
      immax=immax*2

 10   continue



      do 20 m=0,nz-1

      if (ib(m) .gt. m) then

      do 21 j=1,ny2

      do 22 i=1,nx

      vr               = r(i,j    ,m)
      vi               = r(i,j+ny2,m)
      r(i,j    ,m)     = r(i,j    ,ib(m))
      r(i,j+ny2,m)     = r(i,j+ny2,ib(m))
      r(i,j    ,ib(m)) = vr
      r(i,j+ny2,ib(m)) = vi

 22   continue

 21   continue

      endif

 20   continue


      nz2=nz/2



      do 30 j=1,ny2
      do 31 m=0,nz2-1

      if ( m .gt. 0) then

      do 32 i=1,nx


C     the real part of the transformed data pair
C     ( r(i,j,m), r(i,j+ny2,m) ) is stored into
C     ( q(i,j,m), q(i,j+ny2,m) );
C     the imaginary part is stored into
C     ( q(i,j,m+nz2), q(i,j+ny2,m+mz2) )


      q(i,j    ,m    ) = 0.5D0*( r(i,j    ,nz-m) + r(i,j    ,m) )
      q(i,j    ,m+nz2) = 0.5D0*( r(i,j+ny2,nz-m) - r(i,j+ny2,m) )
      q(i,j+ny2,m    ) = 0.5D0*( r(i,j+ny2,nz-m) + r(i,j+ny2,m) )
      q(i,j+ny2,m+nz2) = 0.5D0*(-r(i,j    ,nz-m) + r(i,j    ,m) )


 32   continue

      else

      do 33 i=1,nx


C     the real part of the transformed data pair
C     ( r(i,j,0), r(i,j+ny2,0) ) is stored into
C     ( q(i,j,0), q(i,j+ny2,0) );
C     it's imaginary part is zero
C     the real part of the transformed data pair
C     ( r(i,j,nz2), r(i,j+ny2,nz2) ) is stored into
C     ( q(i,j,nz2), q(i,j+ny2,nz2) );
C     it's imaginary part is zero

      q(i,j    ,m    ) = r(i,j    ,m)
      q(i,j    ,m+nz2) = r(i,j    ,m+nz2)
      q(i,j+ny2,m    ) = r(i,j+ny2,m)
      q(i,j+ny2,m+nz2) = r(i,j+ny2,m+nz2)


 33   continue

      endif

 31   continue
 30   continue






      return
      end
