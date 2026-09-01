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
c------------------------------------------------------------------
c------------------------------------------------------------------
c     Elliptic Integral of the First Kind K(k^2)
c------------------------------------------------------------------
      REAL FUNCTION Elliptic_K(m)
      REAL m, m1
      REAL a0, a1, a2, a3, a4
      REAL b0, b1, b2, b3, b4
      ! A&S 17.3.34 Polynomial + Logarithmic Fit
      a0 = 1.38629436112
      a1 = 0.09666344259
      a2 = 0.03590092383
      a3 = 0.03742563713
      a4 = 0.01451196212
      
      b0 = 0.50000000000
      b1 = 0.12498593597
      b2 = 0.06880248576
      b3 = 0.03328355346
      b4 = 0.00441787012

      IF (m .GE. 1.0) THEN
         m1 = 0.000001
      ELSE IF (m .LE. 0.0) THEN
         m1 = 1.0
      ELSE
         m1 = 1.0 - m
      ENDIF

      Elliptic_K = (a0 + m1*(a1 + m1*(a2 + m1*(a3 + m1*a4))))
     +           + (b0 + m1*(b1 + m1*(b2 + m1*(b3 + m1*b4))))
     +             * ALOG(1.0/m1)
      RETURN
      END

c------------------------------------------------------------------
c     Elliptic Integral of the Second Kind E(k^2)
c------------------------------------------------------------------
      REAL FUNCTION Elliptic_E(m)
      REAL m, m1
      REAL c1, c2, c3, c4
      REAL d1, d2, d3, d4
      ! A&S 17.3.36 Polynomial + Logarithmic Fit
      c1 = 0.44325141463
      c2 = 0.06260601220
      c3 = 0.04757383546
      c4 = 0.01736506451

      d1 = 0.24998368310
      d2 = 0.09200180037
      d3 = 0.04069697526
      d4 = 0.00526449639

      IF (m .GE. 1.0) THEN
         m1 = 0.000001
      ELSE IF (m .LE. 0.0) THEN
         m1 = 1.0
      ELSE
         m1 = 1.0 - m
      ENDIF

      Elliptic_E = 1.0 + m1*(c1 + m1*(c2 + m1*(c3 + m1*c4)))
     +           + m1*(d1 + m1*(d2 + m1*(d3 + m1*d4)))
     +             * ALOG(1.0/m1)
      RETURN
      END