PROGRAM COOLING

    REAL*8 nr,nl,w,Efl,EfR,Pi,Eb,EbR,T,KT1,KT2
    REAL*8 VL,VR,FERMIL,FERMIR,KB,KT,b
    REAL*8 RL,RR
    REAL*8 UP1,UP2,UP3,UP4,UP
    REAL*8 DN1,DN2,DN3,DN4,DN,ELB
    REAL*8 ratio,VBIAS
    REAL*8,ALLOCATABLE,DIMENSION(:) :: ETEMP,WETEMP,EVAL,WGTE
    REAL*8,ALLOCATABLE,DIMENSION(:) :: ETEMP2,WETEMP2,EVAL2,WGTE2
    REAL*8,ALLOCATABLE,DIMENSION(:) :: ETEMP3,WETEMP3,EVAL3,WGTE3
    REAL*8,ALLOCATABLE,DIMENSION(:) :: ETEMP4,WETEMP4,EVAL4,WGTE4
    REAL*8 W1H,W1C,NW,TW,W2H,W2C,W3H,W3C,W4H,W4C
    INTEGER ENUM,J,K
    CHARACTER(len=20) path1,path3
    CHARACTER(len=30) path2
    COMMON /VL/VL/VR/VR/KT1/KT1/KT2/KT2/Eb/Eb/EbR/EbR

!====================================================
!Parameter
!==================================================== 

!    nl=1d-5
    nl=5.97284*1d-2
    nr=5.97284*1d-2 !number divided by 10^24 (10^8)**3, and this is rs=3
    w=0.1d0
    b=1.0d0
    ENUM=500
    KB=8.6170d0*1d-5
    Pi=ACos(-1.0d0)
    WR=13.6d0*(3.d0*Pi**(2.d0)*(nr*(0.529d0)**3.d0))**(2.d0/3.d0)    
    WL=13.6d0*(3.d0*Pi**(2.d0)*(nl*(0.529d0)**3.d0))**(2.d0/3.d0)
    RR=(3.d0/(4*Pi*nr*(0.529d0)**3.d0))**(1.d0/3.d0)
    RL=(3.d0/(4*Pi*nl*(0.529d0)**3.d0))**(1.d0/3.d0)

    OPEN(15, FILE="PARAM.txt")
    WRITE(15,*) "BANDWIDTH(eV)" 
    WRITE(15,*) WL, WR
    WRITE(15,*) "WIGNER-SEITZ RADIUS(in A0)"
    WRITE(15,*) RL, RR
    OPEN(16, FILE="UP_DN.txt")
    OPEN(18, FILE="WORK.txt")
    OPEN(19, FILE="W1H.txt")
    OPEN(20, FILE="W2C.txt")
    OPEN(21, FILE="HC.txt")
    WRITE(16,*) "UP1, UP2, UP3, UP4, DN1, DN2, DN3, DN4"
    T=10.d0
    KT=KB*T
    KT1=KT
    KT2=KT

    WRITE(path1,*) SNGL(b)
    WRITE(path3,*) INT(T)
    path1=ADJUSTL(path1)
    path3=ADJUSTL(path3)
    path2="b="//trim(path1)//"_T="//trim(path3)
          open(17,file=path2,form="formatted",access="sequential", &
          status="replace")
!          WRITE(17,*) "VBIAS, RATIO"
!==========================================================
!BEGIN VOLTAGE LOOP
!===========================================================
 
    DO K=1,81
      VBIAS=-0.05+0.05/40.d0*(K-1)  
      VR=WR
      VL=VR-VBIAS
      EB=VL-WL
      EBR=0.D0
 
      ALLOCATE(EVAL(ENUM),WGTE(ENUM),ETEMP(ENUM),WETEMP(ENUM))
      ALLOCATE(EVAL2(ENUM),WGTE2(ENUM),ETEMP2(ENUM),WETEMP2(ENUM))
      ALLOCATE(EVAL3(ENUM),WGTE3(ENUM),ETEMP3(ENUM),WETEMP3(ENUM))
      ALLOCATE(EVAL4(ENUM),WGTE4(ENUM),ETEMP4(ENUM),WETEMP4(ENUM))

!===================================
!LEFT TO RIGHT (WITH LEFT BEING THE LOWER PART)
!===================================

!==================================
!A CONTROL LOGIC TO SET THE LOWER BOUND OF INTEGRATION
!THE Eb CHANGES WITH RESPECT TO THE BIAS AND IT CAN BE
!LARGER THAN EbR-W. IN SUCH A CASE, THE LOWER BOUND HAS TO BE
!Eb, RATHER THAN EbR-W.
!=================================

    IF(EbR-w.LT.Eb) THEN
     ELB=Eb
    ELSE 
     ELB=EbR-w
    END IF


    CALL GAUSSLEG(ELB,VR+20*KB*T,ENUM,ETEMP,WETEMP)
      DO I=1,ENUM
       EVAL(I)=ETEMP(I)
       WGTE(I)=WETEMP(I)
      END DO

!===================================
!RIGHT TO LEFT (WITH RIGHT BEING THE LOWER PART)
!==================================
    IF(Eb-w.LT.EbR) THEN
     ELB=EbR
    ELSE
     ELB=Eb-w
    END IF

!   WRITE(*,*) "ELB2",ELB, "VR+20*KB*T", VR+20*KB*T

    CALL GAUSSLEG(ELB,VR+20*KB*T,ENUM,ETEMP2,WETEMP2)
      DO I=1,ENUM
       EVAL2(I)=ETEMP2(I)
       WGTE2(I)=WETEMP2(I)
      END DO
!=====================================
!RIGHT ELECTRODE BACKSCATTERING
!=====================================


    CALL GAUSSLEG(EbR,VR+20*KB*T,ENUM,ETEMP3,WETEMP3)
      DO I=1,ENUM
       EVAL3(I)=ETEMP3(I)
       WGTE3(I)=WETEMP3(I)

      END DO
!====================================
!LEFT ELECTRODE BACKSCATTERING
!====================================
  
    CALL GAUSSLEG(Eb,VL+20*KB*T,ENUM,ETEMP4,WETEMP4)
      DO I=1,ENUM
       EVAL4(I)=ETEMP4(I)
       WGTE4(I)=WETEMP4(I)
      END DO

!===================================
!INITIALIZE VARIABLES
!===================================

      UP1=0.D0
      UP2=0.D0
      UP3=0.D0
      UP4=0.D0
      UP=0.D0

      DN1=0.D0
      DN2=0.D0
      DN3=0.D0
      DN4=0.D0
      DN=0.D0
 
      DO I=1,ENUM
         UP1=UP1+WGTE(I)*(1-FERMIL(EVAL(I)))*FERMIR(EVAL(I)+w)        
         UP2=UP2+WGTE2(I)*(1-FERMIR(EVAL2(I)))*FERMIL(EVAL2(I)+w)
         UP3=UP3+WGTE3(I)*(1-FERMIR(EVAL3(I)))*FERMIR(EVAL3(I)+w)*b
         UP4=UP4+WGTE4(I)*(1-FERMIL(EVAL4(I)))*FERMIL(EVAL4(I)+w)*b   
 
         DN1=DN1+WGTE(I)*(FERMIL(EVAL(I))-FERMIR(EVAL(I)+w))
         DN2=DN2+WGTE2(I)*(FERMIR(EVAL2(I))-FERMIL(EVAL2(I)+w))
         DN3=DN3+WGTE3(I)*(FERMIR(EVAL3(I))-FERMIR(EVAL3(I)+w))*b
         DN4=DN4+WGTE4(I)*(FERMIL(EVAL4(I))-FERMIL(EVAL4(I)+w))*b   
       END DO 
         UP=UP1+UP2+UP3+UP4
         DN=DN1+DN2+DN3+DN4
         RATIO=w/(KB*T)*(DLOG(1+DN/UP))**(-1.d0)

         TW=RATIO*T
         NW=1.D0/(DEXP(w/(KB*TW))-1)

         W1C=0.D0
         W1H=0.D0
         W2C=0.D0
         W2H=0.D0
         W3C=0.D0
         W3H=0.D0
         W4C=0.D0
         W4H=0.D0

         DO I=1,ENUM
      W1H=W1H+2.d0*Pi*(1.d0+NW)*(1-FERMIL(EVAL(I)))*FERMIR(EVAL(I)+w)*WGTE(I)
      W1C=W1C+2.d0*Pi*NW*(1-FERMIR(EVAL(I)+w))*FERMIL(EVAL(I))*WGTE(I)
      W2H=W2H+2.d0*Pi*(1.d0+NW)*(1-FERMIR(EVAL2(I)))*FERMIL(EVAL2(I)+w)*WGTE2(I)
      W2C=W2C+2.d0*Pi*NW*(1-FERMIL(EVAL2(I)+w))*FERMIR(EVAL2(I))*WGTE2(I)
    W3H=W3H+b*2.d0*Pi*(1.d0+NW)*(1-FERMIR(EVAL3(I)))*FERMIR(EVAL3(I)+w)*WGTE3(I)
    W3C=W3C+b*2.d0*Pi*NW*(1-FERMIR(EVAL3(I)+w))*FERMIR(EVAL3(I))*WGTE3(I)
    W4H=W4H+b*2.d0*Pi*(1.d0+NW)*(1-FERMIL(EVAL4(I)))*FERMIL(EVAL4(I)+w)*WGTE4(I)
    W4C=W4C+b*2.d0*Pi*NW*(1-FERMIL(EVAL4(I)+w))*FERMIL(EVAL4(I))*WGTE4(I)
         END DO

      WRITE(16,"(8F9.5)")  UP1,UP2,UP3,UP4,DN1,DN2,DN3,DN4
      WRITE(17,"(3F9.6)") VBIAS,WL,ratio 
      WRITE(18,*) WL,W1H,W2C
      WRITE(21,"(4F9.5)") WR,WL,W1H+W2H,W1C+W2C
      DEALLOCATE(EVAL(ENUM),WGTE(ENUM),ETEMP(ENUM),WETEMP(ENUM))
      DEALLOCATE(EVAL2(ENUM),WGTE2(ENUM),ETEMP2(ENUM),WETEMP2(ENUM))
      DEALLOCATE(EVAL3(ENUM),WGTE3(ENUM),ETEMP3(ENUM),WETEMP3(ENUM))
      DEALLOCATE(EVAL4(ENUM),WGTE4(ENUM),ETEMP4(ENUM),WETEMP4(ENUM))
      END DO ! END DO VBIAS LOOP
!    END DO ! END DO TEMP LOOP


END PROGRAM

!=============================================================

      FUNCTION FERMIL(E)
      IMPLICIT NONE
      REAL*8 FERMIL,E,Eb
      real*8 VL,VR,KT1,KT2
      COMMON /VL/VL/VR/VR/KT1/KT1/KT2/KT2/Eb/Eb

       IF(KT1.EQ.0.0d0) then
        IF(E.LT.VL) then
        FERMIL=1.0d0
        else
        FERMIL=0.0d0
        endif

        else if(E.LT.VL.AND.E.GT.Eb) THEN
        FERMIL=1.0d0/(1.0d0+dexp((E-VL)/KT1))

        else if(E.GT.VL) THEN
        FERMIL=dexp(-(E-VL)/KT1)/(1.0d0+dexp(-(E-VL)/KT1))

        else IF(E.LT.Eb) THEN
        FERMIL=0.d0

       ENDIF
       RETURN
      END

      FUNCTION FERMIR(E)
        IMPLICIT NONE
        REAL*8 FERMIR,E,Eb,EbR
        real*8 VL,VR,KT1,KT2
        common /VL/VL/VR/VR/KT1/KT1/KT2/KT2/Eb/Eb/EbR/EbR

        IF(KT2.EQ.0.0d0) then
        IF(E.LT.VR) then
        FERMIR=1.0d0
        else
        FERMIR=0.0d0
        end if
 
 
        else if (E.LT.VR.AND.E.GT.EbR) then
        FERMIR=1.0d0/(1.0d0+dexp((E-VR)/KT2))
        
        else if (E.GT.VR) THEN
        FERMIR=dexp(-(E-VR)/KT2)/(1.0d0+dexp(-(E-VR)/KT2))

        else IF(E.LT.EbR) THEN
        FERMIR=0.d0

        ENDIF

        RETURN

        END

SUBROUTINE GaussLeg(x1,x2,n,x,w)
  Integer n
 Real*8 x1,x2,x(n),w(n)
 Real*8, PARAMETER :: EPS=3.d-14  !** EPS is the relative precision.

 Real*8 ::  p1,p2,p3,pp,xl,xm,z,z1,pi
  Integer :: i,j,m

  pi=acos(-1.d0)
                     !** High precision is a good idea for this routine.
  m=(n+1)/2          !** The roots are symmetric in the interval, so we
  xl=0.5*(x2-x1)     !** only have to do half of them.
  xm=0.5d0*(x2+x1)
  do i=1,m           !** Loop over the desired roots.
     z=cos(pi*(i-.25d0)/(n+.5d0))
     !** Starting with the above approximation to the ith root,
     !** we enter the main loop of refinement by Newton's method.
     z1 = z+10*EPS               !** Cheat to enter loop

     Do While (ABS(z-z1) > EPS)
        p1=1.0d0
        p2=0.0d0

        do j=1,n        !** Loop up the recurrence relation to get the Legendre
           p3=p2        !** polynomial evaluated at z.
           p2=p1
           p1=((2.0*j-1.d0)*z*p2-(j-1.0)*p3)/j
        end do

        !** p1 is now the desired Legendre polynomial. We next compute pp,
        !** its derivative, by a standard relation involving also p2,
        !** the polynomial of one lower order.
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp  !**Newton's method.
     End Do
     x(i)=xm-xl*z                   !** Scale the root to the desired interval,
     x(n+1-i)=xm+xl*z               !** and put in its symmetric counterpart.
     w(i)=2.0*xl/((1.0-z*z)*pp*pp)  !** Compute the weight
     w(n+1-i)=w(i)                  !**and its symmetric counterpart.
  end do
  return
End Subroutine


        SUBROUTINE SPLINE(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=2500)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/&
             (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
       end do
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      return
      END


      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
     return
      END











