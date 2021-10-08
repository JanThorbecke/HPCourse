      subroutine inifft(nz)

      integer nz


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



C     local variables:
      integer alfa,k,m,gamma,nz2,bitr,im,immax,iz


C     description:
C     this routine initiates the arrays of the common block fftdat.
C     these arrays are needed to compute fast fourier transforms in
C     the z-direction.
C     see for instance the subroutine fft.
C     this routine calles the integer function bitr.

      pi = 4.0D0*atan(1.0D0)

      gamma = nint(dlog(1.0d0*nz)/dlog(2.0d0))
      nz2 = nz/2
      m = 0
      immax=1


      do 10 k=1,gamma


      do 11 im=1,immax

      do 12 iz=1,nz2

      alfa   = bitr(m/nz2,gamma)
      c(m,k) = cos(2.0D0*alfa*pi/nz)
      s(m,k) = sin(2.0D0*alfa*pi/nz)

      m = m + 1

 12   continue

      m = m + nz2

 11   continue

      nz2 = nz2/2
      m  = 0
      immax=immax*2

 10   continue


      do 20 m=0,nz-1
      ib(m) = bitr(m,gamma)
 20   continue


      return
      end
