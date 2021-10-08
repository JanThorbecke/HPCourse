      program main

      integer N
      parameter (N=500)
      real  A(N,N,N)
      real  B(N,N,N)
      real  C(N,N,N)
      double precision wallclock_time, t0, t1

      integer i, j, k

      DO i=1, N
         DO j=1, N
            DO k=1, N
               A(i,j,k) = 0.0
               B(i,j,k) = real(i*0.0025+j*0.000125+k*0.00001)
               C(i,j,k) = real(i*0.00125+j*0.00025+k*0.0001)
            END DO
         END DO
      END DO


      t0 = wallclock_time()
      DO i=1, N
         DO j=1, N
            DO k=1, N
               A(i,j,k) = A(i,j,k)+B(i,j,k)*C(i,j,k)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1,1), a(1,N,1), a(N,1,N)
      write(*,*)'time i j k = ', t1-t0

      call setzero(N,A)
      t0 = wallclock_time()
      DO i=1, N
            DO k=1, N
         DO j=1, N
               A(i,j,k) = A(i,j,k)+B(i,j,k)*C(i,j,k)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1,1), a(1,N,1), a(N,1,N)
      write(*,*)'time i k j = ', t1-t0

      call setzero(N,A)
      t0 = wallclock_time()
         DO j=1, N
      DO i=1, N
            DO k=1, N
               A(i,j,k) = A(i,j,k)+B(i,j,k)*C(i,j,k)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1,1), a(1,N,1), a(N,1,N)
      write(*,*)'time j i k = ', t1-t0


      call setzero(N,A)
      t0 = wallclock_time()
         DO j=1, N
            DO k=1, N
      DO i=1, N
               A(i,j,k) = A(i,j,k)+B(i,j,k)*C(i,j,k)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1,1), a(1,N,1), a(N,1,N)
      write(*,*)'time j k i = ', t1-t0

      call setzero(N,A)
      t0 = wallclock_time()
            DO k=1, N
      DO i=1, N
         DO j=1, N
               A(i,j,k) = A(i,j,k)+B(i,j,k)*C(i,j,k)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1,1), a(1,N,1), a(N,1,N)
      write(*,*)'time k i j = ', t1-t0

      call setzero(N,A)
      t0 = wallclock_time()
            DO k=1, N
         DO j=1, N
      DO i=1, N
               A(i,j,k) = A(i,j,k)+B(i,j,k)*C(i,j,k)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1,1), a(1,N,1), a(N,1,N)
      write(*,*)'time k j i = ', t1-t0

      end

      subroutine setzero(N, A)
      integer N,i,j,k
      real A(N,N,N)

       DO k=1, N
         DO j=1, N
            DO i=1, N
               A(i,j,k) = 0.0
            END DO
         END DO
      END DO

      return
      end 
