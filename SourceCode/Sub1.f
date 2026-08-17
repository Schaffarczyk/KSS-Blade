c
c-----------------------------------------------------------------------------------------
c      Mai, June 2019 First coding at ETS, Montreal, 
c      2021 - 2022   added more details 
c      2023  July: try propeller mode
c      June 2026  Bhima: Integrated radial aero, added rational ur/gamma_t, and split BEM loop.
c      June 2026  Bhima: Added Calc_Induction_Convergence and BEMoutput subroutines.
c      June 24    APS    gam, etc
c      June 2026  Bhima: Unified Physics into Solver & Elliptic Integrals.
c                        Fixed 'kappa' pass-back bug for structural bending.
c                        Added RadialModel toggle (Hansen vs. Vortex Cylinder).
c----------------------------------------------------------------------
       subroutine BEM(cP,cT,errb)
c-----------------------------------------------------------------------
      use mem
      real clint(2),cdint(2),dct(0:1)
      real zero(100),aoabs(100)
      real slope, theta, dFx, dFz, kappa
      real ur, clintth, cdintth
      real cthr, dFzinc, stgp, ctinc, cdl, dz
      real oopdefi, oopdefim1, dum
      real term1, term2, term3, x, z
      real thick1, thick2, kg
      real gamma_t, gam, phi, w2, cn, ct, TL
      real aset, da, a_prev, a1, a2
      real a, dd, BB, cp_loc, ctloc, ctaer
      integer ibs, j, nprs, na, izero
      character*2 iters
      character*5 pstring
      character*20 nout2, nout3, nout4, nout5
      logical bsearch
      bsearch   = .FALSE.

      errb = 0.	
      ac = 0.4
      om    = rpm*pi/30.
      vtip  = om*rtip
      TSR   = vtip/vwind
      write(*,*)'BEM: vwind TSR=',vwind, tsr
      write(*,*)

      IOOUT  = 11
      nout  = './Bem.out'
      OPEN(UNIT=ioout,FILE=Nout,Form='formatted',status='unknown')
      IO2  = 2
      nout2  = './wt-perf.out'
      OPEN(UNIT=io2,FILE=Nout2,Form='formatted',status='unknown')
      IO3  = 3
      nout3  = './FAST.out'
      OPEN(UNIT=io3,FILE=Nout3,Form='formatted',status='unknown')
      io4  = 4        
      nout4 = 'BS.out'
      OPEN(UNIT=io4,FILE=nout4,Form='formatted',status='unknown')
      io5  = 17 
      nout5 = 'bisecAOA.out'
      OPEN(UNIT=io5,FILE=nout5,Form='formatted',status='unknown')

      write(*,101)'r ',' a ','ap ','F','w','chord','twist','phi ',
     +            'aoa','cL','cD','L2D','cNo','cTa','dFn','dFt','Ga',
     +            'iter','err','th','Prof','It','dFx','dFz','kappa',
     +            'cTloc','a_r'
      write(ioout,101)'r ',' a ','ap ','F','w','chord',
     +            'twist','phi ',
     +            'aoa','cL','cD','L2D','cNo','cTa','dFn','dFt','Ga',
     +            'iter','err','th','Prof','It','dFx','dFz','kappa',
     +            'cTloc','a_r'

      thick1 = 1.
      thick2 = 0.	
      thrust = 0.
      torque = 0.
      aold   = .3
      anew   = .3
      apold  = 0.
      apnew  = eps

c     BEGIN section loop
      dr = (rtip-rroot)/float(nsec)
      do i=1,nsec 
         iter = 1
         bsearch = .true.
         do kk = 1,100
            zero(kk) = 0.
            aoabs(kk) = 0.
         end do
         r = rroot + 0.5*dr + (i-1)*dr

         if (DesMode)then
            nspl = Nsec
         else
            nspl = ndes
         end if

         call SPLINT(rsecsp,chsp   ,chS,   Nspl,r,chord,dummy)
         call SPLINT(rsecsp,twistsp,twistS,Nspl,r,twist,dummy)
         call SPLINT(rsecsp,oopdefsp,oopdefS,Nspl,r,oopdefi,dummy)
         call SPLINT(rsecsp,oopdefsp,oopdefS,Nspl,r-dr,oopdefim1,dum)

         if (twistb)then
            call twistbend(vwind,r,eltwist)
            twist = twist - eltwist
         endif

         if(chord.lt.0)chord = 0.001
         rred = r/rtip
         th   = thicksp(rred)

         thick2 = th
         if(thick2.gt.thick1)then
             th = thick1
         else
             thick1 = thick2
         endif
         if(th.lt.minthick) th = minthick

         call Get_Aero_Coeffs(th, xp, clintth, cdintth)
         tsrloc = om*r/vwind

        if (bsearch)then
           dct(0) = 0.
           dct(1) = 0.
           dd     = 1.
           izero = 0
           nprs = indpro
           na = 100
           da = 1.0/float(na)
           do ibs=1,na
              do j = 0,1
                 a    = (ibs+j-1)*da
                 dct(j) = funcBS(a,r,chord,twist,nprs,xp)
              end do
              a1 = (ibs-1)*da
              a2 =     ibs*da
              dd = dct(0)*dct(1)
              if(dd.lt.0.)then
               izero = izero + 1
               zero(izero) = rtbis(r,chord,twist,a1,a2,1.e-8)
               call ctsec(zero(izero),r,chord,twist,nprs,ctloc,ctaer,
     +             aoabs(izero),xp)
              end if
            end do
            write(io4,211)r,izero,(zero(j),j=1,8)
            write(io5,211)r,izero,(aoabs(j),j=1,8)
211         format(f10.2,i6,8f8.3)
        end if

        If (izero.gt.1) then
           da = 0.25*abs(zero(2)-zero(1))
           aset = zero(1)
        else
           da = 0.05
           aset = zero(1)
        endif
        bsearch = .false.	

 10      iter  = iter + 1

         if(iter.gt.1)then
            aold  = anew
            apold = apnew
         end if   

         if(izero.eq.1)then
            iiters = 1
            iters  = ' F'
         else
            iiters = 3
            iters  = ' B'
         end if
         
         iiters = 1
         iters  = ' F'
         if(iter.gt.maxiter-10)then
            iiters=2
            iters = ' N'
         endif

         if (i .eq. 1) then
            a_prev = 0.0
         else
            a_prev = abem(i-1)
         endif

c       -----------------------------------------------------------------------
c       Unified Physics and Numerical Convergence Solver
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         call Calc_Induction_Convergence(iiters, r, chord, twist, dr,
     +           aold, apold, a_prev, th, xp, oopdefi, 
     +           oopdefim1, aset, da, anew, apnew, erri, phi, w2, 
     +           cn, ct, cthr, ur, gam, gamma_t, TL, clintth, cdintth, 
     +           kappa)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(erri.gt.eps.and.iter.lt.maxiter) goto 10

        if(erri.gt.errb)errb=erri

        dT  =  0.5*dens*B*chord*w2*cn*dr
        dFt =  0.5*dens*B*chord*w2*ct*dr
        dTa = dT/ (B*dr)
        dFta= dFt/(B*dr)
        dQ  =  dFt*r
        dT  = dT/1000.
        dFt = dFt/1000.
        dQ  = dQ/1000.
        clo = clintth
        cdo = cdintth

        abem(i) = anew
        apbem(i)= apnew

        dFx = dTa * cos(kappa)
        dFz = dTa * sin(kappa)

        call BEMoutput(r, anew, apnew, TL, w2, chord, twist, phi, 
     +                 clo, cdo, cn, ct, dTa, dFta, gam, iter, 
     +                 erri, th, iters, dFx, dFz, kappa, cthr, ur)

        write(io2,105)chord,th*chord,twist
        pstring = 'PRINT'
        drnodes = dr
        write(io3,106)r,twist,drnodes,chord,indpro,Pstring

         thrust = thrust + dT
         torque = torque + dQ
      end do

      pow = om*torque
      write(*,*)
      write(*,103)'Thrust= ',thrust,' kN'
      write(*,103)'Torque= ',torque,' kNm'
      write(*,103)'Power = ',pow   ,' kW'
      write(ioout,103)'Thrust= ',thrust,' kN'
      write(ioout,103)'Torque= ',torque,' kNm'
      write(ioout,103)'Power = ',pow   ,' kW'

      Pref = 0.5*dens*ar*vwind**3
      Tref = 0.5*dens*ar*vwind**2
      cp   = 1000.*pow/Pref
      ct   = 1000.*thrust/Tref
      write(*,103)'v-wind = ',vwind
      write(*,103)'TSR = ',TSR
      write(*,103)'pitch = ',glopitch
      write(*,103)'cP = ',cp
      write(*,103)'cT = ',ct
      write(*,103)'cP/cT = ',cp/ct
      write(*,103)'cQ = ',cp/tsr
      write(*,103)'errb  = ',errb
      write(ioout,103)'cP = ',cp
      write(ioout,103)'cT = ',ct
      write(ioout,103)'cQ = ',ct/tsr
      write(ioout,103)'TSR = ',TSR
      write(ioout,103)'pitch = ',glopitch
      write(ioout,103)'errb  = ',errb
      write(ioout,103)'cQ = ',ct/tsr
      write(ioout,103)'TSR = ',TSR
      write(ioout,103)'pitch = ',glopitch
      write(ioout,103)'errb  = ',errb
      Close(UNIT=ioout)
      Close(UNIT=io2)
      Close(UNIT=io3)
      Close(UNIT=io4)
      Close(UNIT=io5)

101   format(14a8,             4a9,       a8,    a9,    a5,  a3, 5a9)
102   format(a10,f8.3,a12)
103   format(a20,f12.4,a10)
104   format(2(a15,f12.4))
105   format(3f12.6)
106   format(4f12.6,i5,a10)
500   format(3a12)   
504   format(2f12.6)
888   format(a10,2f12.6)
      end

	function rtbis(r,chord,twist,x1,x2,xacc)
	use mem
	integer jmax, j
	real rtbis, x1, x2, xacc, funcBS
	real dx, f0, fmid, xmid
	external funcBS
	parameter (jmax=20)
        nprs = indpro
	fmid = funcBS(x2,r,chord,twist,nprs,xp)
        f0   = funcBS(x1,r,chord,twist,nprs,xp)
	if (f0*fmid.gt.0.) write(778,*)'root must be inside rtbis'
	if (f0.lt.0.) then
   	   rtbis = x1
	   dx    = x2 - x1
	else
	   rtbis = x2
	   dx    = x1 - x2
	end if
	do j = 1,jmax
	   dx   = 0.5*dx
           xmid = rtbis + dx
	   fmid = funcBS(xmid,r,chord,twist,nprs,xp)
 	   if(fmid.le.0.)rtbis = xmid
	   if(abs(dx).lt.xacc.or.fmid.eq.0.)return
   	enddo
	write(*,*)'too many bisections'
888     format(a10,i4,f12.6)
        end

	function funcBS(a,r,chord,twist,nprs,xp)
	use mem
	nprs = indpro
        call ctsec(a,r,chord,twist,nprs,ctloc,ctaer,aoas,xp)
	funcBS = ctloc - ctaer
        return
	end

        subroutine ctsec(a,r,chord,twist,nprs,ctloc,ctaer,aoas,xp)
	use mem
        real clints(2), cdints(2)
        integer nprs
        aa   = a*(1.-a)/tsrloc
        ap   = sqrt(aa+0.25) - 0.5
        phi  = atan2(vwind*(1.-a), (r*om*(1.+ap)))
        w2   =     (vwind*(1.-a))**2 + (r*om*(1.+ap))**2
        phid = 180.*phi/pi
        aoas = phid-twist
        FP     = FPR(x,phi)
	 do k =0,1
            do j=1,300
               aoasp(j) = 0.
               clsp (j) = 0.
               clss (j) = 0.
               cdsp (j) = 0.
               cdss (j) = 0.
            enddo
            k2 = np(nprs)
	    do j=1,k2
               aoasp(j) = aoain(nprs,j)
               clsp (j) = clin (nprs,j)
               clss (j) = cls  (nprs,j)
               cdsp (j) = cdin (nprs,j)
               cdss (j) = cds  (nprs,j)
            enddo	 
            call SPLINT
     +      (aoasp,clsp,clss,Np(nprs),aoas,clints(k+1),dummy)
            call SPLINT
     +      (aoasp,cdsp,cdss,Np(nprs),aoas,cdints(k+1),dummy)
	 enddo
        clint = xp*clints(1)+(1.-xp)*clints(2)
        cdint = xp*cdints(1)+(1.-xp)*cdints(2)
        cn   = clint*cos(phi) + cdint*sin(phi)  
        cta  = clint*sin(phi) - cdint*cos(phi)  
        ctaer  = cn
        cp = 4.*a*(1.-a)**2
        if(a.gt.ac)then
           ctloc = 4.*a*(1.-0.25*(5.-3.*a)*a)*FP
        else
           ctloc = 4.*a*(1.-a)
        endif
	return
101     format(7f12.4)
	end
	

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Calc_Induction_Convergence(iiters, r, chord, twist, dr,
     +           aold, apold, a_prev, th, xp, oopdefi, 
     +           oopdefim1, aset, da, anew, apnew, erri, phi, w2, 
     +           cn, ct, cthr, ur, gam, gamma_t, TL, clintth, cdintth,
     +           kappa)
c---------------------------------------------------------------------------
      use mem
      integer iiters
      real r, chord, twist, dr, aold, apold, a_prev, th, xp
      real oopdefi, oopdefim1, aset, da
      real anew, apnew, erri, phi, w2, cn, ct, cthr
      real ur, gam, gamma_t, TL, clintth, cdintth, kappa
      
      real dz, phid, cdl
      real R_cyl, k2, k, term1, term2, term3
      real AA, kG, sigH, ahansen, a_guess, a1, a2, denom
      
      real Elliptic_K, Elliptic_E
      external Elliptic_K, Elliptic_E

      a_guess = aold

c     --- Part 1: Kinematics ---
      if (r .le. rroot) then
         kappa = 0.0
      else
         dz    = oopdefi - oopdefim1
         kappa = atan(dz / dr)
      endif 

      phi  = atan2(vwind*(1.-a_guess)*cos(kappa), (r*om*(1.+apold)))
      if (phi.lt.0.) phi = 0.0001
      w2   = (vwind*(1.-a_guess)*cos(kappa))**2 + (r*om*(1.+apold))**2
      phid = 180.*phi/pi
      aoab = phid - twist

      call Get_Aero_Coeffs(th, xp, clintth, cdintth)
      cn   = clintth*cos(phi) + cdintth*sin(phi)         
      ct   = clintth*sin(phi) - cdintth*cos(phi)
      cthr = cn*cos(phi) - ct*sin(phi)
      if (cthr.gt.2.0) cthr = 2.0

c     ==================================================================
c     RADIAL PHYSICS TOGGLE (1 = Hansen, 2 = Vortex Cylinder)
c     ==================================================================
      if (RadialModel .eq. 2) then
c        --- Part 3: Tangential Vorticity (Eq. 20) ---
         gamma_t = 2.0 * vwind * (a_guess - a_prev)

c        --- Part 4: Radial Induction via Elliptic Integrals (Eq. 10 & 11) ---
         R_cyl = r + (dr / 2.0)
         k2    = (4.0 * r * R_cyl) / ((R_cyl + r)**2)
         if (k2 .ge. 1.0) k2 = 0.999999 
         k     = sqrt(k2)
         
         term1 = -(gamma_t) / (2.0 * pi) * sqrt(R_cyl / r)
         term2 = ((2.0 - k2) / k) * Elliptic_K(k2)
         term3 = (2.0 / k) * Elliptic_E(k2)
         ur    = term1 * (term2 - term3)
      else
c        --- Classical Hansen Model (No spanwise wake expansion) ---
         gamma_t = 0.0
         ur      = 0.0
      endif

c     --- Part 2: Circulation ---
      w2  = w2 + ur*ur
      cdl = sqrt(clintth**2 + cdintth**2)
      gam = 0.5 * cdl * sqrt(w2) * chord

c     NUMERICAL ROOT FINDING
      TL = FPR(r/rtip, phi) * RL(r/rtip, phi)
      AA = B * chord * cn / (8. * pi * TL * r * (sin(phi))**2)

      select case(iiters)
      case (1)
         anew  = AA / (1. + AA)
         if (anew.gt.0.4) then
            sigH = chord * B / (2. * pi * r)
            kG   = 4. * TL * sin(phi)*sin(phi) / (sigH * cn)
            anew = ahansen(kG, 0.4)
         end if
         if (anew.ge.0.99) anew = 0.99

      case (2)
         anew  = AA / (1. + AA) 

      case (3)
         a1    = max(0., aset - da)
         a2    = min(aset + da, 1.)
         anew  = rtbis(r, chord, twist, a1, a2, 1.e-6)

      case default
      end select 

      erri  = abs(anew - aold) / abs(anew)
      denom = 8. * TL * pi * r * sin(phi) * cos(phi)
      apnew = (B * chord * ct) / denom
      apnew = apnew / (1. - apnew)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c=======================================================================
      subroutine BEMoutput(r, anew, apnew, TL, w2, chord, twist, phi, 
     +           clo, cdo, cn, ct, dTa, dFta, gam, iter, 
     +           erri, th, iters, dFx, dFz, kappa, cthr, ur)
c=======================================================================
      use mem

      real r, anew, apnew, TL, w2, chord, twist, phi, clo, cdo, cn
      real ct, dTa, dFta, gam, erri, th, dFx, dFz, kappa, cthr, ur
      integer iter
      character*2 iters
      real r2d, a_r

      r2d = 180.0 / pi
c     --- Task 4: Calculate radial induction factor ---
      a_r = ur / vwind

      write(*,100)r,anew,apnew,TL,sqrt(w2),chord,twist,
     +            phi*57.3,aoab,clo,cdo,clo/cdo,cn,ct,dTa,dFta,gam,
     +            iter,erri,th,nameprout,iters,dFx,dFz,r2d*kappa,
     +            cthr,a_r
      write(ioout,100)r,anew,apnew,TL,sqrt(w2),chord,twist,
     +            phi*57.3,aoab,clo,cdo,clo/cdo,cn,ct,dTa,dFta,gam,
     +            iter,erri,th,nameprout,iters,dFx,dFz,r2d*kappa,
     +            cthr,a_r
  100 format(7f8.3,2f8.1,5f8.3,3f9.1,i8,x,e8.1,x,f8.6,x,a4,x,a2, 5f9.3)

      return
      end
c=======================================================================
c     Task 1: Subroutine Get_Aero_Coeffs
c=======================================================================
      subroutine Get_Aero_Coeffs(th, xp, clintth, cdintth)
      use mem
      real th, xp, clintth, cdintth
      real clint(2), cdint(2), dummy, dth
      integer k, kk, k2, j

      do k =1,nopr
         if (th.lt.prothick(k))then 
             continue
         else
             indpro = k
             goto 123
         endif
      end do     

  123 if(k.eq.1) indpro=2      
      nameprout = namepr(k)

      do k =0,1
         do j=1,300
            aoasp(j) = 0.
            clsp (j) = 0.
            clss (j) = 0.
            cdsp (j) = 0.
            cdss (j) = 0.
         enddo
         kk = indpro-k
         k2 = np(kk)
         do j=1,k2
            aoasp(j) = aoain(kk,j)
            clsp (j) = clin (kk,j)
            clss (j) = cls  (kk,j)
            cdsp (j) = cdin (kk,j)
            cdss (j) = cds  (kk,j)
         enddo	 
         call SPLINT
     +   (aoasp,clsp,clss,Np(kk),aoab,clint(k+1),dummy)
         call SPLINT
     +   (aoasp,cdsp,cdss,Np(kk),aoab,cdint(k+1),dummy)
      enddo

      dth =  prothick(indpro)-prothick(indpro-1)
      xp  = (prothick(indpro)-th)/dth
      clintth = xp*clint(2)+(1.-xp)*clint(1)
      if (clintth.gt. 3.)clintth= 3.
      if (clintth.lt.-2.)clintth=-2.
      cdintth = xp*cdint(2)+(1.-xp)*cdint(1)

      return
      end