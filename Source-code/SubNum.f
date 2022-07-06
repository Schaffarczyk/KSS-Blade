c
C--------------------------------------------------
C--------------------------------------------------
C     NUMERICAL RECIPES ROUTINEN
C--------------------------------------------------
c
c--------------------------------------------------
c    spline interpolation
c--------------------------------------------------
C
       SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
c
c      numerical receipes chapter 3.3, p109
c	N = number of "tabulated data"
c
       INTEGER N,NMAX
       REAL YP1,YPN,X(N),Y(N),Y2(N)
       PARAMETER (NMAX=500)
C
       INTEGER I,K
	REAL P,QN,SIG,UN,U(NMAX)
	IF (YP1.GT..99E30)THEN
            Y2(1) = 0.
            U(1)  = 0.
	ELSE
	    Y2(1)=-0.5
	    U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
	ENDIF
C
	DO I=2,N-1
	   SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
           P=SIG*Y2(I-1)+2.
	   Y2(I)=(SIG-1.)/P
	   U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
C
     +    	 /(X(I)-X(I-1)))/(X(I+1)-x(i-1))-sig*u(i-1))/p
	end do
C
	IF (YPN.GT..99E+30)THEN
	   QN=0.
	   UN=0.
	ELSE
	   QN=0.5
	   UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
	ENDIF
C
	Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
C
	DO K=N-1,1,-1
		Y2(K)=Y2(K)*Y2(K+1)+U(K)
	END DO
C
	RETURN
	END
C
C-----------------------------------------------------
	SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y,YS)
C------------------------------------------------------
c
c      Y2A output from SPLINE -> N must be the same
c
	INTEGER N, K, KHI, KLO
c
	REAL XA(N),Y2A(N),YA(N)
	REAL A,B,H,X,Y,YS, eps
c
	eps = 1.e-5
c
	KLO = 1
	KHI = N
1	IF (KHI-KLO.GT.1) THEN
		K=(KHI+KLO)/2
		IF(XA(K).GT.X)THEN
			KHI = K
		ELSE
			KLO = K
		ENDIF
	GOTO 1
	ENDIF
	H=XA(KHI)-XA(KLO)
	IF(abs(H).lt.eps)then
	    XA(KHI) = XA(KLO) + 2.*eps
            H=XA(KHI)-XA(KLO)
c            write(*,*) 'BAD XA INPUT IN splint',H
	end if
	A=(XA(KHI)-X)/H
	B=(X-XA(KLO))/H
	y=a*ya(klo)+b*ya(khi)+
c
     +  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
c
        YS=(YA(KHI)-YA(KLO))/(XA(KHI)-XA(KLO)) - 
     +     ((3.*A*A-1.)/6.)*(xa(khi)-xa(klo))*y2a(klo) +
     +     ((3.*b*b-1.)/6.)*(xa(khi)-xa(klo))*y2a(khi)
C
c	write(125,*)'splint',y
c
	RETURN
	END
	
