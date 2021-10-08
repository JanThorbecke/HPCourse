
program convolution
implicit none
integer, parameter :: wp=kind(1.0d0)  !< wavefunction-type precision
integer l, i1, j, k, i, i2
integer, parameter :: lowfil=-8,lupfil=7, n1=1024, ndat=10000
integer :: mod_arr(lowfil:n1+lupfil)
real(wp), dimension(0:n1,ndat) :: x
real(wp), dimension(ndat,-lupfil:n1-lowfil) :: y
real(wp) :: tt
real(wp) :: fil(lowfil:lupfil)
integer :: cnt, clk_rate, cnt_max, t1, t2
real :: rnumber
DATA fil / &
       8.4334247333529341094733325815816e-7_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.9940415697834003993178616713_wp,&
       -0.604895289196983516002834636e-1_wp, &
       -0.2103025160930381434955489412839065067e-1_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       2.72734492911979659657715313017228e-6_wp /

  do i2=1,ndat
     do i1=0,n1
        call random_number(rnumber)
        x(i1,i2) = rnumber
     enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!
!periodic boundaries

  call system_clock(t1, clk_rate, cnt_max)
  do i2=1,ndat
     do i1=0,n1
        tt = 0.0
        do l=lowfil,lupfil
           k=modulo(i1+l,n1+1)
!           write(*,*)'i1=',i1,' l=',l,'k=',k
           tt=tt+x(k,i2)*fil(l)
        enddo
        y(i2,i1)=tt
     enddo
  enddo
  call system_clock(t2, clk_rate, cnt_max)
  write(*,*)'periodic modulus: y(1,5)=',y(1,5),'y(1,',n1/2,')=',y(1,n1/2),'Time=',real(t2-t1,8)/real(clk_rate,8)

! Implement your faster version of the code in this part
  call system_clock(t1, clk_rate, cnt_max)
  do i2=1,ndat
     do i1=0,n1
        tt = 0.0
        do l=lowfil,lupfil
           k=modulo(i1+l,n1+1)
           tt=tt+x(k,i2)*fil(l)
        enddo
        y(i2,i1)=tt
     enddo
  enddo
  call system_clock(t2, clk_rate, cnt_max)
  write(*,*)'periodic faster: y(1,5)=',y(1,5),'y(1,',n1/2,')=',y(1,n1/2),'Time=',real(t2-t1,8)/real(clk_rate,8)

!!!!!!!!!!!!!!!!!!!!!!!!!
!non-periodic boundaries
  call system_clock(t1, clk_rate, cnt_max)
  do i2=1,ndat
     do i1=0,n1
        tt = 0.0
        do l=max(lowfil,-i1),min(lupfil,n1-i1)
           k=i1+l
           tt=tt+x(k,i2)*fil(l)
        enddo
        y(i2,i1)=tt
     enddo
  enddo
  call system_clock(t2, clk_rate, cnt_max)
  write(*,*)'non-periodic min-max: y(1,5)=',y(1,5),'y(1,',n1/2,')=',y(1,n1/2),'Time=',real(t2-t1,8)/real(clk_rate,8)

! Implement your faster version of the code in this part
  call system_clock(t1, clk_rate, cnt_max)
  do i2=1,ndat
     do i1=0,n1
        tt=0.e0_wp
        do l=max(lowfil,-i1),min(lupfil,n1-i1)
           k=i1+l
           tt=tt+x(k,i2)*fil(l)
        enddo
        y(i2,i1)=tt
     enddo
  enddo
  call system_clock(t2, clk_rate, cnt_max)
  write(*,*)'non-periodic faster: y(1,5)=',y(1,5),'y(1,',n1/2,')=',y(1,n1/2),'Time=',real(t2-t1,8)/real(clk_rate,8)
   
end

