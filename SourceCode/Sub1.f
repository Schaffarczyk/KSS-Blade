c
c-----------------------------------------------------------------------------------------
c      Mai, June 2019 First coding at ETS, Montreal, 
c      2021 - 2022   added more details 
c
c      2022 binary search if fixed point iteration and Newton does not converge
c
c      2023  July: try propeller mode
c 
c      June 2026  Bhima: Integrated radial aero, added rational ur/gamma_t, and split BEM loop.
c      June 2026 Bhima: Added new subroutine Calc_Induction_Convergence
c      June 2026 Bhima: Added new Subroutine BEMoutput.
c      June 24 APS gam, etc
c----------------------------------------------------------------------
       subroutine BEM(cP,cT,errb)
c-----------------------------------------------------------------------
c
      use mem
c
      real clint(2),cdint(2),dct(0:1)
      real zero(100),aoabs(100)
      real slope, theta, dFx, dFz, kappa
      real ur
      real cthr, dFzinc, ctinc, dz
      real thrsum, ctav, uarg, gam
      real oopdefi, oopdefim1, dum

      real thick1, thick2, kg
c
      character*2 iters
      character*5 pstring
      character*20 nout2, nout3, nout4, nout5
c
      logical bsearch
      bsearch   = .FALSE.
c
c------------------------------------------------------------------
c
      write(*,*)
      write(*,*)'************************'
      write(*,*)'routine BEM ************'
c
      errb = 0.	
c
c     critical a value to for cT(a) correction (Glauert/Hansen)
c
      ac = 0.4
c
c     some auxilary values
c
      om    = rpm*pi/30.
      vtip  = om*rtip
      TSR   = vtip/vwind
      write(*,*)'BEM: vwind TSR=',vwind, tsr
      write(*,*)
c
c     open output files
c
      IOOUT  = 11
      nout  = './Bem.out'
      OPEN(UNIT=ioout,FILE=Nout,Form='formatted',status='unknown')
c
      IO2  = 2
      nout2  = './wt-perf.out'
      OPEN(UNIT=io2,FILE=Nout2,Form='formatted',status='unknown')
c
      IO3  = 3
      nout3  = './FAST.out'
      OPEN(UNIT=io3,FILE=Nout3,Form='formatted',status='unknown')
c          
      io4  = 4        
      nout4 = 'BS.out'
      OPEN(UNIT=io4,FILE=nout4,Form='formatted',status='unknown')

      io5  = 17 
      nout5 = 'bisecAOA.out'
      OPEN(UNIT=io5,FILE=nout5,Form='formatted',status='unknown')
c
      write(*,101)'r ',' a ','ap ','F','w','chord','twist','phi ',
     +            'aoa','cL','cD','L2D','cNo','cTa','dFn','dFt','Ga',
     +            'iter','err','th','Prof','It','dFx','dFz','kappa',
     +            'cTloc','ur_M'
c
      write(ioout,101)'r ',' a ','ap ','F','w','chord',
     +            'twist','phi ',
     +            'aoa','cL','cD','L2D','cNo','cTa','dFn','dFt','Ga',
     +            'iter','err','th','Prof','It','dFx','dFz','kappa',
     +            'cTloc','ur_M'
c----------------------------------------------------------------------------------
c     initialize some data
c
      thick1 = 1.
      thick2 = 0.	
c
      thrust = 0.
      torque = 0.
c
c     radial induction (Madsen): initialize accumulated thrust for
c     the averaged thrust coefficient (Eq. 7) and the radial
c     induced velocity
c
      thrsum = 0.
      ur     = 0.
      kappa  = 0.
      gam    = 0.
c
      aold   = .3
      anew   = .3
c
      apold  = 0.
      apnew  = eps
c
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
c     BEGIN section loop
c
      dr = (rtip-rroot)/float(nsec)
      do i=1,nsec 
         iter = 1
c
	 bsearch = .true.
         do kk = 1,100
            zero(kk) = 0.
            aoabs(kk) = 0.
         end do
c
	 r = rroot + 0.5*dr + (i-1)*dr
c
c        interpolate for chord and twist
c
	 if (DesMode)then
	    nspl = Nsec
         else
            nspl = ndes
         end if
c
         call SPLINT(rsecsp,chsp   ,chS,   Nspl,r,chord,dummy)
         call SPLINT(rsecsp,twistsp,twistS,Nspl,r,twist,dummy)
c        Interpolate out-of-plane deflection for local r and r-dr
         call SPLINT(rsecsp,oopdefsp,oopdefS,Nspl,r,oopdefi,dummy)
         call SPLINT(rsecsp,oopdefsp,oopdefS,Nspl,r-dr,oopdefim1,dum)
c
c        add pitch due to bend twist coupling 
c
 	 if (twistb)then
            call twistbend(vwind,r,eltwist)
c
c           deformation should change twist to more positive angles
c           d. as given always is negative:->  change sign here
c
            twist = twist - eltwist
         endif
c
c        be sure that chord is always a positve number
c
	 if(chord.lt.0)chord = 0.001
c
c        get local thickness
c
	 rred = r/rtip
	 th   = thicksp(rred)
c
c        ensure monotonic decrease of thickness
c
	 thick2 = th
c
	 if(thick2.gt.thick1)then
             th = thick1
	 else
	     thick1 = thick2
         endif
c
c        ensure thickness to be larger than thinnest
c
	 if(th.lt.minthick)then
	    th = minthick
         end if
c
cBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c        BEGIN BEM loop    
c
 10      iter  = iter + 1
         iters = ' F'
c
         if(iter.gt.1)then
            aold  = anew
            apold = apnew
         end if   
c
c     (1) start
c     --- Calculate Kappa (Local Dihedral Angle) ----(Start Code - updated by Bhima on 05.05.2026, 26.05.2026)
c      2026 07: kappa only if Radial = .T. ;
c               comparison now against rroot (not 2.8667)
c
      if (Radial .and. r.gt.rroot) then
         dz    = oopdefi - oopdefim1
         kappa = atan(dz / dr)
      else
         kappa = 0.0
      endif
c ----(End Code - updated by Bhima on 05.05.2026, 26.05.2026)
c
c
c       start with initial values for a and aprime (1/3 and 0.)
c (Line 169 - updated by Bhima on 28.04.2026, 05.05.2026)
c       for Radial = .F. is kappa = 0 and the classical
c       expressions are recovered (cos(kappa) = 1)
c
        phi = atan2(vwind*(1.-aold)*cos(kappa), (r*om*(1.+apold)))
c         phi  = atan2(vwind*(1.-aold),(r*om*(1.+apold)))
c         phi  = atan2((r*om*(1.+apold)),vwind*(1.-aold))
         if (phi.lt.0.)phi = 0.0001
         w2 = (vwind*(1.-aold)*cos(kappa))**2 + (r*om*(1.+apold))**2
         phid = 180.*phi/pi
         aoab = phid-twist
c (Line 182 - updated by Bhima on 05.05.2026)
c      (1) end
c-----------------------------------------------------------------------------------
c        thickness interpolation for cL and cD
c
c        find first profile thinner then local thickness
c
	 do k =1,nopr
	    if (th.lt.prothick(k))then 
                continue
            else
                indpro = k
                goto 123
            endif
         end do     
c
c       1st profile: cylinder with th = 1.
c
123     if(k.eq.1)indpro=2      
	nameprout = namepr(k)
c
c       cL/cD interpolation if thickness is not close (within eps) to a given one
c
	 do k =0,1
            do j=1,300
               aoasp(j) = 0.
               clsp (j) = 0.
               clss (j) = 0.
               cdsp (j) = 0.
               cdss (j) = 0.
            enddo
c
            kk = indpro-k
            k2 = np(kk)
	    do j=1,k2
               aoasp(j) = aoain(kk,j)
               clsp (j) = clin (kk,j)
               clss (j) = cls  (kk,j)
               cdsp (j) = cdin (kk,j)
               cdss (j) = cds  (kk,j)
            enddo	 
c
            call SPLINT
     +      (aoasp,clsp,clss,Np(kk),aoab,clint(k+1),dummy)
            call SPLINT
     +      (aoasp,cdsp,cdss,Np(kk),aoab,cdint(k+1),dummy)
	 enddo
c
         dth =  prothick(indpro)-prothick(indpro-1)
 	 xp   = (prothick(indpro)-th)/dth
c
c        store weight for rtbis/funcBS (was undefined there before)
c
         xpint = xp
c
         clintth = xp*clint(2)+(1.-xp)*clint(1)
c
c        just for safety - may be not appropriate if very special airfoils are used
c        (2026 07: removed duplicated block)
c
	 if (clintth.gt. 3.)clintth= 3.
	 if (clintth.lt.-2.)clintth=-2.
c
         cdintth = xp*cdint(2)+(1.-xp)*cdint(1)
c
c---------------------------------------------------------------------------
c        calculate normal and tangential forces
c        note: ct here is "tangential" not in thrust direction
c    
         cn   = clintth*cos(phi) + cdintth*sin(phi)         
         ct   = clintth*sin(phi) - cdintth*cos(phi)
c        Calculate local thrust component and limit it to 2.0
         cthr = cn*cos(phi) - ct*sin(phi)
         if (cthr.gt.2.0) cthr = 2.0
c
c        tip loss 
c----------------------------------------------------------------------------
 	 x      = r/rtip
         tsrloc = om*r/vwind
c        Prandtl-Glauert
         TL = FPR(x,phi)
c         Burton (does not work at the moment)
c         Tl     = FBu(x,tsrloc,aold)
c
c        root loss
c----------------------------------------------------------------------------
         TL = TL*RL(x,phi)
c----------------------------------------------------------------------------
c
c       scan each section for multiple solutions
c
        if (bsearch)then
           dct(0) = 0.
           dct(1) = 0.
           dd     = 1.
c
	   izero = 0
           nprs = indpro
c
c          a loop
c          scan 0 < a < 1 in steps of 1.e-2
c
           na = 100
           da = 1.0/float(na)
c
	   do ibs=1,na
              do j = 0,1
                 a    = (ibs+j-1)*da
                 dct(j) = funcBS(a,r,chord,twist,nprs,xp)
	      end do
c
              a1 = (ibs-1)*da
              a2 =     ibs*da
              dd = dct(0)*dct(1)
              if(dd.lt.0.)then
               izero = izero + 1
               zero(izero) = rtbis(r,chord,twist,a1,a2,1.e-8)
              call ctsec(zero(izero),r,chord,twist,nprs,ctloc,ctaer,
     +             aoabs(izero),xp)
              end if
c
            end do
            write(io4,211)r,izero,(zero(j),j=1,8)
            write(io5,211)r,izero,(aoaBS(j),j=1,8)
211         format(f10.2,i6,8f8.3)
        end if
c
c       store some values for binary search
c
        If (izero.gt.1) then
           da = 0.25*abs(zero(2)-zero(1))
           aset = zero(1)
        else
           da = 0.05
           aset = zero(1)
        endif
	bsearch = .false.	
c
c
        AA = B*chord*cn/(8.*pi*TL*r*(sin(phi))**2)
c
c      (2) 2026 07:
c      the former in-loop radial correction was removed here:
c      it used ur, clo and cdo before they were defined (undefined
c      variables in the first iteration / first section).
c      The radial contribution (Madsen) is now applied once per
c      section AFTER convergence of the induction iteration,
c      see block (5) below.
c
c       a' indpt from iteration schema
c        
        BB    = B*chord*ct/(8.*TL*pi*r*sin(phi)*cos(phi))
        apnew = BB/(1.-BB)
c
        if(apnew.lt.-.5)apnew = -.5
c
cDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
c---------------------------------------------------------------------------
c       use BINARY SEARCH if izero ne 1 and fixed point else
c       use Newton if only 10 iters left
c---------------------------------------------------------------------------
c
        if(izero.eq.1)then
            iiters = 1
        else
            iiters = 3
        end if
c
c       here seems to be a bug (2023 Feb 06)
c
 	iiters = 1
c
        if(iter.gt.maxiter-10)iiters=2
c
        call Calc_Induction_Convergence(iiters, AA, ac, chord, r,
     +           phi, cn, TL, twist, aset, da, aold, anew, erri, iters)
c       -----------------------------------------------------------------------
        if(erri.gt.eps.and.iter.lt.maxiter) goto 10
c
c       END BEM loop  
cBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c
        if(erri.gt.errb)errb=erri
c
c       (3) start
c--------(Begin Code - added & updated by Bhima on 09.06.2026)
c       -----------------------------------------------------
c       Vortex Cylinder: Tangential Vorticity (Eq. 20 & 38)
c       -----------------------------------------------------
c       Determine the axial induction of the previous section
        if (i .eq. 1) then
c          Ghost section at the root has zero induction
           a_prev = 0.0
        else
c          Pull the converged induction from the previous loop
           a_prev = abem(i-1)
        endif
        
c       Calculate gamma_t using U0 (vwind) and the induction delta
        gamma_t = 2.0 * vwind * (anew - a_prev)
c--
c--------(End Code - added & updated by Bhima on 09.06.2026)
c
c       (3) end
c
c       section thrust and tangential force, all B blades, in N
c
        dT  =  0.5*dens*b*chord*w2*cn*dr
        dFt =  0.5*dens*b*chord*w2*ct*dr
c
c       bound circulation from lift (Kutta-Joukowski)
c
        gam = 0.5*sqrt(w2)*chord*clintth
c
c       (5) start  2026 07
c-----------------------------------------------------------------------
c       simplest algorithm to include the radial velocity component:
c       Madsen correction, Li et al., WES 7, 75-104 (2022)
c
c       Eq. (7): averaged thrust coefficient up to radius r
c                cTav(r) = int_0^r cT 2 pi r* dr* / (pi r**2)
c                        = (accumulated thrust) / (0.5 dens U0**2 pi r**2)
c       Eq. (6): radial induced velocity
c                ur = U0/2.24 * cTav/(4 pi) *
c                     ln[ (0.04**2+(x+1)**2) / (0.04**2+(x-1)**2) ]
c                with x = r/rtip
c       Eq. (5): additional tangential driving force of the
c                non-planar (prebent) blade section
c                dFzinc = dens * gam * ur * tan(kappa) * dr  (per blade)
c-----------------------------------------------------------------------
c
        ur     = 0.
        dFzinc = 0.
        ctinc  = 0.
c
c       V3: store bound circulation and dihedral angle of this
c       section for the vortex cylinder sweep after the loop
c
        gamsec(i) = gam
        kapsec(i) = kappa
c
c       RadMode = 1: Madsen correction, applied per section;
c       RadMode = 2: vortex cylinder model, applied AFTER the
c                    section loop (needs a(i) of ALL sections),
c                    see block (6) below
c
        if (Radial .and. RadMode.eq.1) then
c
c          Eq. (7): accumulate thrust from hub to current section
c
           thrsum = thrsum + dT
           ctav   = thrsum/(0.5*dens*vwind**2*pi*r**2)
c
c          Eq. (6): radial induced velocity (exact log form; the
c          former series expansion 2z+2/3 z**3+2/5 z**5 with
c          z = 2x/(0.04**2+x**2+1) is its 3-term approximation)
c
           x    = r/rtip
           uarg = (0.04**2+(x+1.)**2)/(0.04**2+(x-1.)**2)
           ur   = vwind/2.24 * ctav/(4.*pi) * log(uarg)
c
c          Eq. (5): extra tangential driving force (per blade)
c
           dFzinc = dens*gam*ur*tan(kappa)*dr
c
c          add to tangential force of all B blades and
c          to the tangential force coefficient (diagnostic)
c
           dFt    = dFt + B*dFzinc
           ctinc  = dFzinc/(0.5*dens*w2*chord*dr)
           ct     = ct + ctinc
        endif
c
c       (5) end
c
        dTa = dT/ (B*dr)
        dFta= dFt/(B*dr)
c
        dQ  =  dFt*r
        dT  = dT/1000.
        dFt = dFt/1000.
        dQ  = dQ/1000.
c
        clo = clintth
        cdo = cdintth
c
	abem(i) = anew
        apbem(i)= apnew
c
        dFx = dTa * cos(kappa)
        dFz = dTa * sin(kappa)

c       Write section results to output files
        call BEMoutput(r, anew, apnew, TL, w2, chord, twist, phi, 
     +                 clo, cdo, cn, ct, dTa, dFta, gam, iter, 
     +                 erri, th, iters, dFx, dFz, kappa, cthr, ur)
c--------(End Code - updated by Bhima on 28.04.2026, 22.06.2026)     
c---------------------------------------------------------------------------------
c
c       to compare with wt-perf
c
        write(io2,105)chord,th*chord,twist
c
c       to compare with FAST
c
        pstring = 'PRINT'
        drnodes = dr
        write(io3,106)r,twist,drnodes,chord,indpro,Pstring
c
c----------------------------------------------------------------------------------
c
c        sum up for blade's torque and thrust
c
         thrust = thrust + dT
         torque = torque + dQ
c
      end do
c
c     END section loop
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
c
c     (6) start   V3 2026 07
c-----------------------------------------------------------------------
c     vortex cylinder model (Section 3 of Li et al., WES 7, 2022)
c
c     needs the converged a(i) of ALL sections -> after the loop.
c     Cylinder strengths gam_t(k) = 2 U0 (a(k)-a(k-1)) (Eq. 20/38),
c     radial induced velocity by superposition of the elliptic-
c     integral solution (Eqs. 10+11), see Sub3.f.
c     The extra tangential driving force per blade and section
c     (Eq. 5),  dFzinc = dens*gam*ur*tan(kappa)*dr,  is added to
c     the rotor torque; section details are written to VC.out.
c-----------------------------------------------------------------------
c
      if (Radial .and. RadMode.eq.2) then
c
         call VCradial(nsec,dr,urvc)
c
         io6 = 18
         OPEN(UNIT=io6,FILE='./VC.out',Form='formatted',
     +        status='unknown')
         write(io6,301)'r','a','kappa','ur_VC','dFzinc'
c
         dtorvc = 0.
         do i = 1,nsec
            r      = rroot + 0.5*dr + float(i-1)*dr
            dFzinc = dens*gamsec(i)*urvc(i)*tan(kapsec(i))*dr
            dtorvc = dtorvc + B*dFzinc*r
            write(io6,302)r,abem(i),180.*kapsec(i)/pi,
     +                    urvc(i),dFzinc
         end do
         close(unit=io6)
c
c        add to rotor torque (kNm)
c
         torque = torque + dtorvc/1000.
c
         write(*,103)'VC dTorque= ',dtorvc/1000.,' kNm'
      endif
c
301   format(5a12)
302   format(3f12.4,2f12.5)
c
c     (6) end
c
c     GLOBAL output
c
      pow = om*torque
c
      write(*,*)
      write(*,103)'Thrust= ',thrust,' kN'
      write(*,103)'Torque= ',torque,' kNm'
      write(*,103)'Power = ',pow   ,' kW'
c
      write(ioout,103)'Thrust= ',thrust,' kN'
      write(ioout,103)'Torque= ',torque,' kNm'
      write(ioout,103)'Power = ',pow   ,' kW'
c
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
c
      write(ioout,103)'cP = ',cp
      write(ioout,103)'cT = ',ct
      write(ioout,103)'cQ = ',ct/tsr
      write(ioout,103)'TSR = ',TSR
      write(ioout,103)'pitch = ',glopitch
      write(ioout,103)'errb  = ',errb
c
      Close(UNIT=ioout)
      Close(UNIT=io2)
      Close(UNIT=io3)
      Close(UNIT=io4)
      Close(UNIT=io5)
c ----(Line 523 & 524 - updated by Bhima on 05.05.2026)
cccc--100   format(7f8.3,2f8.1,5f8.3,3f9.1,i8,x,e8.1,x,f8.6,x,a4,x,a2, 5f9.3)----(Deleted by Bhima on 22.06.2026)
101   format(14a8,             4a9,       a8,    a9,    a5,  a3, 5a9)
102   format(a10,f8.3,a12)
103   format(a20,f12.4,a10)
104   format(2(a15,f12.4))
105   format(3f12.6)
106   format(4f12.6,i5,a10)
500   format(3a12)   
504   format(2f12.6)
888   format(a10,2f12.6)
c
      end
c
c------------------------------------------------------------------------
c       2022 11 29 rtbis from NUMERICAL RECIPES pp 347
c-------------------------------------------------------------------------
c
	function rtbis(r,chord,twist,x1,x2,xacc)
c
	use mem
c
	integer jmax, j
	real rtbis, x1, x2, xacc, funcBS
	real dx, f0, fmid, xmid
	external funcBS
	parameter (jmax=20)
c
        nprs = indpro
c
c       2026 07: use xpint from module mem - the former local xp
c       was never set in this routine (undefined variable)
c
	fmid = funcBS(x2,r,chord,twist,nprs,xpint)
        f0   = funcBS(x1,r,chord,twist,nprs,xpint)
c        write(*,*)'x1 x2 f0 fmid',x1,x2,f0,fmid
	if (f0*fmid.gt.0.) write(778,*)'root must be inside rtbis'
	if (f0.lt.0.) then
   	   rtbis = x1
	   dx    = x2 - x1
	else
	   rtbis = x2
	   dx    = x1 - x2
	end if
c
	do j = 1,jmax
	   dx   = 0.5*dx
           xmid = rtbis + dx
	   fmid = funcBS(xmid,r,chord,twist,nprs,xpint)
 	   if(fmid.le.0.)rtbis = xmid
	   if(abs(dx).lt.xacc.or.fmid.eq.0.)return
   	enddo
	write(*,*)'too many bisections'
c
888     format(a10,i4,f12.6)
        end
c
c------------------------------------------------------------------------------------------
	function funcBS(a,r,chord,twist,nprs,xp)
c------------------------------------------------------------------------------------------
	use mem
c
	nprs = indpro
c
        call ctsec(a,r,chord,twist,nprs,ctloc,ctaer,aoas,xp)
	funcBS = ctloc - ctaer
c
        return
	end
c
c------------------------------------------------------------------------------------------
        subroutine ctsec(a,r,chord,twist,nprs,ctloc,ctaer,aoas,xp)
c------------------------------------------------------------------------------------------
	use mem
c
        real clints(2), cdints(2)
        integer nprs
c
        aa   = a*(1.-a)/tsrloc
        ap   = sqrt(aa+0.25) - 0.5
c (Line 590 - updated by Bhima on 28.04.2026)
        phi  = atan2(vwind*(1.-a), (r*om*(1.+ap)))
        w2   =     (vwind*(1.-a))**2 + (r*om*(1.+ap))**2
        phid = 180.*phi/pi
        aoas = phid-twist
c
c       2026 07: x was undefined here before (used by FPR)
c
        x      = r/rtip
        FP     = FPR(x,phi)
c
c       2026 07: interpolate between the two neighbouring profiles
c       (nprs and nprs-1) as in routine BEM - before, both spline
c       evaluations used the same profile nprs and the weights
c       were swapped with respect to BEM
c
	 do k =0,1
            do j=1,300
               aoasp(j) = 0.
               clsp (j) = 0.
               clss (j) = 0.
               cdsp (j) = 0.
               cdss (j) = 0.
            enddo
c
            kk = nprs - k
            if (kk.lt.1) kk = 1
            k2 = np(kk)
	    do j=1,k2
               aoasp(j) = aoain(kk,j)
               clsp (j) = clin (kk,j)
               clss (j) = cls  (kk,j)
               cdsp (j) = cdin (kk,j)
               cdss (j) = cds  (kk,j)
            enddo
c
            call SPLINT
     +      (aoasp,clsp,clss,Np(kk),aoas,clints(k+1),dummy)
            call SPLINT
     +      (aoasp,cdsp,cdss,Np(kk),aoas,cdints(k+1),dummy)
	 enddo
c
c	write(*,*)'clint cdint',clint, cdint
c
        clint = xp*clints(2)+(1.-xp)*clints(1)
        cdint = xp*cdints(2)+(1.-xp)*cdints(1)
c
        cn   = clint*cos(phi) + cdint*sin(phi)  
        cta  = clint*sin(phi) - cdint*cos(phi)  
c
        ctaer  = cn
c
        cp = 4.*a*(1.-a)**2
c
c       2026 07: ac was undefined here (local variable) -
c       set to the same critical value as in routine BEM
c
        ac = 0.4
c
        if(a.gt.ac)then
           ctloc = 4.*a*(1.-0.25*(5.-3.*a)*a)*FP
        else
           ctloc = 4.*a*(1.-a)
        endif
c
c	write(*,*)'ctrs a ctloc ctaer',a,ctloc,ctaer
c
	return
c
101     format(7f12.4)
	end
c	
c
c--------(Start Code - added & updated by Bhima on 09.06.2026, 22.06.2026)
c---------------------------------------------------------------------------
      subroutine Calc_Induction_Convergence(iiters, AA, ac, chord, r,
     +           phi, cn, TL, twist, aset, da, aold, anew, erri, iters)
c---------------------------------------------------------------------------
      use mem

      integer iiters
      real AA, ac, chord, r, phi, cn, TL, twist, aset, da
      real aold, anew, erri
      character*2 iters

      real sigH, kG, dphi, phi2, AA1, AA2, AAp, a1, a2
      real ahansen, rtbis
c
c     (4) 2026 07: the calculation of the radial induction ur was
c     moved to routine BEM, block (5): it now uses the AVERAGED
c     thrust coefficient of Eq. (7) (before: the local one) and the
c     exact logarithmic expression of Eq. (6) (before: 3-term
c     series expansion), and it is switched by the flag "Radial"
c     from Machine.in.
c
      select case(iiters)
c     ----------------------------------------------------------------------
      case (1)
c        (1) Fixed point iteration
c        EQs from WWL '76 pp 23/39 Eqs (2.4.4), (2.4.6), (2.6.1) & (2.6.2)
         anew  = AA/(1.+AA)
c        use Glauert extension for cT(a) if a > ac
         if (anew.gt.ac)then
            sigH = chord*B/(2.*pi*r)
            kG   = 4.*TL*sin(phi)*sin(phi)/(sigH*cn)
            anew = ahansen(Kg,ac)
         end if
c
c        write(*,*)'anew ',anew
c 
c        Assure valid WT operation:  0 < a < 1
         if (anew.ge..99) anew = 0.99

c     ----------------------------------------------------------------------
      case (2)
c        (2) NEWTON's method 
         iters = ' N'
         aold  = aset
         dphi  = 0.01*phi
         phi2  = phi + dphi
         AA1   = AA
         AA2   = B*chord*cn/(8.*TL*pi*r*sin(phi2)*sin(phi2))
         AAp   = (AA2-AA1)/dphi
         anew  = aold - AA/AAp
         anew  = anew/(1.+anew)
         if(anew.ge..99) anew = 0.99

c     ----------------------------------------------------------------------
      case (3)
c        (3) Binary search
         iters = ' B'
         a1    = max(0.,aset - da)
         a2    = min(aset + da,1.)
         aold  = anew
         erri  = 1.e-6
         anew  = rtbis(r,chord,twist,a1,a2,1.e-6)
c
      case default
      end select 

c     Calculate final iteration error
      erri = abs(anew-aold)/abs(anew)
c      write(*,*)'anew aold erri',anew,aold,erri
      return
      end
c------------(subroutine added & updated by Bhima on 22.06.2026)
c
      subroutine BEMoutput(r, anew, apnew, TL, w2, chord, twist, phi, 
     +           clo, cdo, cn, ct, dTa, dFta, gam, iter, 
     +           erri, th, iters, dFx, dFz, kappa, cthr, ur)
c=======================================================================
      use mem

      real r, anew, apnew, TL, w2, chord, twist, phi, clo, cdo, cn
      real ct, dTa, dFta, gam, erri, th, dFx, dFz, kappa, cthr, ur
      integer iter
      character*2 iters
      real r2d

c     Convert kappa to degrees for the output file
      r2d = 180.0 / pi

      write(*,100)r,anew,apnew,TL,sqrt(w2),chord,twist,
     +            phi*57.3,aoab,clo,cdo,clo/cdo,cn,ct,dTa,dFta,gam,
     +            iter,erri,th,nameprout,iters,dFx,dFz,r2d*kappa,
     +            cthr,ur
c
c     added by APS 2026 06 25
c
      write(ioout,100)r,anew,apnew,TL,sqrt(w2),chord,twist,
     +            phi*57.3,aoab,clo,cdo,clo/cdo,cn,ct,dTa,dFta,gam,
     +            iter,erri,th,nameprout,iters,dFx,dFz,r2d*kappa,
     +            cthr,ur
c
  100 format(7f8.3,2f8.1,5f8.3,3f9.1,i8,x,e8.1,x,f8.6,x,a4,x,a2, 5f9.3)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c--------(End Code - added & updated by Bhima on 09.06.2026, 22.06.2026)
