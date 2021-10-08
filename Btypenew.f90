!  Btypenew.f90 
!
!  FUNCTIONS:
!  Btypenew      - Entry point of console application.
!
!  Example of displaying 'Hello World' at execution time.
!

!****************************************************************************
!
!  PROGRAM: Btypenew
!
!  PURPOSE:  Entry point for 'Hello World' sample console application.
!
!****************************************************************************

    program Btypenew

    implicit none
    
    real t1,t2,t3
    
    real*8,dimension(:,:),allocatable :: a
    real*8,dimension(:,:),pointer :: b
    integer i,k,inner,outer

    inner = 100000
    outer = 10000
    print *, 'Hello World'
    
    call cpu_time(t1)

    allocate( a(inner,3) )
    do k=1,outer
      do i=1,inner
        a(i,1)=5.0
        a(i,2)=3.0
        a(i,3)=1.0
      enddo
    enddo

    write(*,*) " a(538,inner) =", a(538,1), a(538,2), a(538,3) ! , a(538,4)

    call cpu_time(t2)

    do k=1,outer
      a(:,1) = 1.0
      a(:,2) = 2.0
      a(:,3) = 7.0
      ! a(:,4) = 7.0
    enddo

    call cpu_time(t3)

    write(*,*) " a(538,inner) =", a(538,1), a(538,2), a(538,3) ! , a(538,4)

    write(*,*) "CPU array expl.loop/impl.loop:", t2-t1, t3-t2
    
    call cpu_time(t1)

    allocate( b(inner,3) )
    do k=1,outer
      do i=1,inner
        b(i,1)=5.0
        b(i,2)=3.0
        b(i,3)=1.0
      enddo
    enddo

    call cpu_time(t2)
    write(*,*) " b(538,inner) =", b(538,1), b(538,2), b(538,3) 

    do k=1,outer
      b(:,1) = 1.0
      b(:,2) = 2.0
      b(:,3) = 7.0
      ! b(:,4) = 7.0
    enddo

    call cpu_time(t3)

    write(*,*) " b(538,inner) =", b(538,1), b(538,2), b(538,3) ! , b(538,4)

    write(*,*) "CPU p.array expl.loop/impl.loop:", t2-t1, t3-t2
    
    call cpu_time(t1)

    do k=1,outer
      do i=1,inner,5
        b(i,1)=5.0
        b(i+1,1)=5.0
        b(i+2,1)=5.0
        b(i+3,1)=5.0
        b(i+4,1)=5.0
      enddo
      do i=1,inner,5
        b(i,2)=3.0
        b(i+1,2)=3.0
        b(i+2,2)=3.0
        b(i+3,2)=3.0
        b(i+4,2)=3.0
      enddo
      do i=1,inner,5
        b(i,3)=1.0
        b(i+1,3)=1.0
        b(i+2,3)=1.0
        b(i+3,3)=1.0
        b(i+4,3)=1.0
      enddo
    enddo

    call cpu_time(t2)
    write(*,*) " b(538,inner) =", b(538,1), b(538,2), b(538,3) ! , b(538,4)

    do k=1,outer
      b(:,1) = 1.0
      b(:,2) = 2.0
      b(:,3) = 7.0
    enddo

    call cpu_time(t3)
    write(*,*) " b(538,inner) =", b(538,1), b(538,2), b(538,3) ! , b(538,4)

    write(*,*) "CPU p.array unrolled expl.loop/impl.loop:", t2-t1, t3-t2

    end program Btypenew
