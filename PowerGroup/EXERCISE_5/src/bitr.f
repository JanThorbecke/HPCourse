      integer function bitr(n,gamma)

      integer n, gamma


C     last modification:
C     31 January 1996


C     local variables
      integer i,j,b,m

C     description:
C     this integer function is called by inifft.

      i = 0
      b = 0
      m = n
 1111 continue
      if (i .lt. gamma) then
         j = m/2
         b = 2*b + (m-2*j)
         m = j
         i = i + 1
      else
         bitr = b
         return
      endif
      goto 1111

      end
