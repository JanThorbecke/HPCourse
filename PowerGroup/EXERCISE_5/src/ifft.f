      subroutine ifft(nx,ny,nz)

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
      double precision cc,ss,qnz



C     description:
C     this routine computes the inverse fast fourier transform in the
C     third direction of the real array q.
C     the inverse fourier-transform fq is computed according to
C
C     fq(i,j,m) = sum over k=0,..,nz-1 of
C                                 q(i,j,k)*exp(-2*pi*I*k*m/nz)/nz
C
C     where i=1,..,nx, j=1,..,ny, and m=0,...,nz-1.
C     note: I*I=-1, and pi=3.1415...
C     since here q(i,j,k) must be conjugate even, i.e.
C
C     q(i,j,m) = complex conjugate of q(i,j,nz-m),
C
C     it holds that the inverse fourier transform fq(i,j,m) is real.
C     to save storage the real part of q(i,j,m) is stored in
C     the real array q(i,j,m) where i=1,..,nx, j=1,..,ny, and
C     m=0,...,nz/2; the imaginary part of q(i,j,m) is stored in
C     the real array q(i,j,m+nz/2) where i=1,..,nx, j=1,..,ny,
C     and m=1,...,nz/2-1.
C
C     the results of the inverse transform are stored in q(i,j,m)
C     where i=1,..,nx, j=1,..,ny, and m=0,...,nz-1.
C     the subroutine inifft must be runned before ifft, in order
C     to initialize the variables in the common block fftdat.



      ny2=ny/2
      nz2=nz/2


C     in this routine ny must be even
C     this is checked for first

      if ( ny2 + ny2 .ne. ny ) then
      stop "ny must be even"
      endif


C     the complex input for the inverse fft is constructed as the
C     sum of q(i,j,m) + I*q(i,j+ny2,m)
C     the real part of this sum is stored in the real array r(i,j,m)
C     where j=1,..,ny2; it's imaginary part is stored in r(i,j+ny2,m)
C     where j=1,..,ny2.


      qnz =1.0D0/nz

      do 10 j=1,ny2
      do 11 m=0,nz2-1

      if ( m .gt. 0) then

      do 12 i=1,nx

      r(i,j    ,m   ) = qnz*( q(i,j,m    ) - q(i,j+ny2,m+nz2) )
      r(i,j+ny2,m   ) = qnz*( q(i,j,m+nz2) + q(i,j+ny2,m    ) )

      r(i,j    ,nz-m) = qnz*( q(i,j,m    ) + q(i,j+ny2,m+nz2) )
      r(i,j+ny2,nz-m) = qnz*(-q(i,j,m+nz2) + q(i,j+ny2,m    ) )

 12   continue

      else

      do 13 i=1,nx

      r(i,j    ,m    ) = qnz*( q(i,j    ,m) )
      r(i,j+ny2,m    ) = qnz*( q(i,j+ny2,m) )

      r(i,j    ,nz2+m) = qnz*( q(i,j    ,nz2+m) )
      r(i,j+ny2,nz2+m) = qnz*( q(i,j+ny2,nz2+m) )

 13   continue

      endif

 11   continue
 10   continue




C     start of the fast fourier transform


      gamma = nint(dlog(1.0d0*nz)/dlog(2.0d0))



      nz2 = nz/2
      m=0
      immax=1


      do 20 k=1,gamma

      do 21 im=1,immax

      do 22 iz=1,nz2

      cc =  c(m,k)
      ss = -s(m,k)


CMIC@    DO ALL SHARED(NZ2,NY2,NX,M,K,CC,R,SS)PRIVATE(IZ,J,I,VR,VI)
      do 23 j=1,ny2

      do 24 i=1,nx

      vr       = cc*r(i,j    ,m+nz2) + ss*r(i,j+ny2,m+nz2)
      vi       = cc*r(i,j+ny2,m+nz2) - ss*r(i,j    ,m+nz2)

      r(i,j    ,m+nz2) = r(i,j    ,m) - vr
      r(i,j+ny2,m+nz2) = r(i,j+ny2,m) - vi
      r(i,j    ,m    ) = r(i,j,    m) + vr
      r(i,j+ny2,m    ) = r(i,j+ny2,m) + vi

 24   continue

 23   continue

      m = m + 1

 22   continue

      m = m + nz2

 21   continue

      nz2 = nz2/2
      m = 0
      immax=immax*2

 20   continue



      do 30 j=1,ny2
      do 31 m=0,nz-1

      if (ib(m) .gt. m) then

      do 32 i=1,nx

      vr               = r(i,j    ,m)
      vi               = r(i,j+ny2,m)
      r(i,j    ,m)     = r(i,j    ,ib(m))
      r(i,j+ny2,m)     = r(i,j+ny2,ib(m))
      r(i,j    ,ib(m)) = vr
      r(i,j+ny2,ib(m)) = vi

 32   continue

      endif

 31   continue
 30   continue


      nz2=nz/2

      do 40 j=1,ny2
      do 41 m=1,nz2-1
      do 42 i=1,nx

      vr              =  r(i,j    ,m)
      vi              =  r(i,j+ny2,m)
      r(i,j    ,m)    =  r(i,j    ,nz-m)
      r(i,j+ny2,m)    =  r(i,j+ny2,nz-m)
      r(i,j    ,nz-m) =  vr
      r(i,j+ny2,nz-m) =  vi

 42   continue
 41   continue
 40   continue


      return
      end
