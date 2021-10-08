      program main

      integer N
      parameter (N=1000)
      real  A(N,N)
      real  B(N,N)
      real  C(N,N)
      double precision wallclock_time, t0, t1

      integer i, j, k

      DO i=1, N
         DO j=1, N
               A(i,j) = 0.0
               B(i,j) = real(i*0.0025+j*0.000125)
               C(i,j) = real(i*0.00125+j*0.00025)
         END DO
      END DO


      t0 = wallclock_time()
      DO i=1, N
         DO j=1, N
            DO k=1, N
               A(i,j) = A(i,j)+B(i,k)*C(k,j)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1), a(1,N), a(N,1)
      write(*,*)'time i j k = ', t1-t0
      DO i=1, N
         DO j=1, N
               A(i,j) = 0.0
         END DO
      END DO


      t0 = wallclock_time()
      DO i=1, N
            DO k=1, N
         DO j=1, N
               A(i,j) = A(i,j)+B(i,k)*C(k,j)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1), a(1,N), a(N,1)
      write(*,*)'time i k j = ', t1-t0
      DO i=1, N
         DO j=1, N
               A(i,j) = 0.0
         END DO
      END DO


      t0 = wallclock_time()
         DO j=1, N
      DO i=1, N
            DO k=1, N
               A(i,j) = A(i,j)+B(i,k)*C(k,j)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1), a(1,N), a(N,1)
      write(*,*)'time j i k = ', t1-t0
      DO i=1, N
         DO j=1, N
               A(i,j) = 0.0
         END DO
      END DO


      t0 = wallclock_time()
         DO j=1, N
            DO k=1, N
      DO i=1, N
               A(i,j) = A(i,j)+B(i,k)*C(k,j)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1), a(1,N), a(N,1)
      write(*,*)'time j k i = ', t1-t0
      DO i=1, N
         DO j=1, N
               A(i,j) = 0.0
         END DO
      END DO

      t0 = wallclock_time()
            DO k=1, N
      DO i=1, N
         DO j=1, N
               A(i,j) = A(i,j)+B(i,k)*C(k,j)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1), a(1,N), a(N,1)
      write(*,*)'time k i j = ', t1-t0
      DO i=1, N
         DO j=1, N
               A(i,j) = 0.0
         END DO
      END DO

      t0 = wallclock_time()
            DO k=1, N
         DO j=1, N
      DO i=1, N
               A(i,j) = A(i,j)+B(i,k)*C(k,j)
            END DO
         END DO
      END DO
      t1 = wallclock_time()

      write(*,*)'A=',a(1,1), a(1,N), a(N,1)
      write(*,*)'time k j i = ', t1-t0

      end


