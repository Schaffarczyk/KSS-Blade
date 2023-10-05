c
c-----------------------------------------------------------------------------------------
c      Mai, June 2019 First coding at ETS, Montreal, 
c      2021 - 2022   added more details 
c
c      2022 binary search if fixed point iteration and Newton does not converge
c
c      2023  July: try propeller mode
c
c----------------------------------------------------------------------
      subroutine BEM(cP,cT)
c-----------------------------------------------------------------------
c
      use mem
c
      real clint(2),cdint(2),dct(0:1)
      real zero(100),aoabs(100)
c
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
     +            'iter','err','th','Prof','IterS'
c
      write(ioout,101)'r ',' a ','ap ','F','w','chord',
     +            'twist','phi ',
     +            'aoa','cL','cD','L2D','cNo','cTa','dFn','dFt','Ga',
     +            'iter','err','th','Prof','IterS'
c----------------------------------------------------------------------------------
c     initialize some data
c
      thick1 = 1.
      thick2 = 0.	
c
      thrust = 0.
      torque = 0.
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
c       start with initial values for a and aprime (1/3 and 0.)
c
         phi  = atan(vwind*(1.-aold)/(r*om*(1.+apold)))
         if (phi.lt.0.)phi = 0.0001
         w2   =     (vwind*(1.-aold))**2 + (r*om*(1.+apold))**2
         phid = 180.*phi/pi
         aoab = phid-twist
c
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
         clintth = xp*clint(2)+(1.-xp)*clint(1)
c
c        just for safety - may be not appropriate if very special airfoils are used
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
c
	 cthr = cn*cos(phi) - ct*sin(phi)
	 if (cthr.gt.2) cthr=2.
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
c       aprime indept. of iteration schema
c
        BB    = B*chord*ct/(8.*TL*pi*r*sin(phi)*cos(phi))
        apnew = BB/(1.-BB)
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
c       here's seem to be a bug (2023 Feb 06
c
 	iiters = 1
c
        if(iter.gt.maxiter-10)iiters=2
c
cDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
        AA     = B*chord*cn/(8.*TL*pi*r*sin(phi)*sin(phi))
c
c---------------------------------------------------------------------------
        select case(iiters)
c---------------------------------------------------------------------------
        case (1)
c---------------------------------------------------------------------------
c       (1)  Fixed point iteration
c---------------------------------------------------------------------------
c        EQs from WWL '76 pp 23/39 Eqs (2.4.4), (2.4.6), (2.6.1) and (2.6.2)
c---------------------------------------------------------------------------
c
         anew  = AA/(1.+AA)
c
c        use Glauert extension for cT(a) if a > ac
c        From M. Hansen's book
c
         if (anew.gt.ac)then
            sigH = chord*B/(2.*pi*r)
            kG   = 4.*TL*sin(phi)*sin(phi)/(sigH*cn)
            anew = ahansen(Kg,ac)
         end if 
c
c        NREL's schema
c 	  cthnrel = 8./9.+(4.*tl-40./9.)*anew
c         cthnrel = cthnrel +(50./9.-4.*TL)*anew*anew 
c
c        Assure valid WT operation:  0 < a < 1
c
	 if (anew.ge.1.0)anew = 0.9999
c
c----------------------------------------------------------------------------------
         case (2)
c----------------------------------------------------------------------------------
c        (2) NEWTON's method 
c----------------------------------------------------------------------------------
c
            iters = ' N'
c
c           as usual: initial value is very important !
c
            aold  = aset
c
            dphi  = 0.01*phi
            phi2  = phi + dphi
            AA1   = AA
            AA2   = B*chord*cn/(8.*TL*pi*r*sin(phi2)*sin(phi2))
            AAp   = (AA2-AA1)/dphi
            anew  = aold - AA/AAp
	    anew  = anew/(1.+anew)
c            write(777,*)'r iter anew',r,iter, anew
c
          case (3)
c----------------------------------------------------------------------------------
c        (3) Binary search
c----------------------------------------------------------------------------------
C
            iters  = ' B'
c
            a1   = max(0.,aset - da)
            a2   = min(aset + da,1.)
            aold = anew
            err  = 1.e-6
            anew = rtbis(r,chord,twist,a1,a2,1.e-6)
c
        case default
        end select 
c-----------------------------------------------------------------------------------   
        err = abs(anew-aold)/abs(anew)
        if(err.gt.eps.and.iter.lt.maxiter) goto 10
c
c       END BEM loop  
cBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c
        dT  =  0.5*dens*b*chord*w2*cn*dr
        dFt =  0.5*dens*b*chord*w2*ct*dr
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
        gam = .5*clo*sqrt(w2)*chord
c
        write(*,100)r,anew,apnew,TL,sqrt(w2),chord,twist,
     +              phi*57.3,aoab,clo,cdo,clo/cdo,cn,ct,dTa,dFta,gam,
     +               iter,err,th,nameprout,iters
c
        write(ioout,100)r,anew,apnew,TL,sqrt(w2),chord,twist,
     +              phi*57.3,aoab,clo,cdo,clo/cdo,cn,ct,dTa,dFta,gam,
     +              iter,err,th,nameprout,iters
c
c-------------------------------------------------------------------------------
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
c------------------------------------------------------------
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
c     output
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
      write(*,103)'cQ = ',cp/tsr
c
      write(ioout,103)'cP = ',cp
      write(ioout,103)'cT = ',ct
      write(ioout,103)'cQ = ',ct/tsr
      write(ioout,103)'TSR = ',TSR
      write(ioout,103)'pitch = ',glopitch
c
      Close(UNIT=ioout)
      Close(UNIT=io2)
      Close(UNIT=io3)
      Close(UNIT=io4)
      Close(UNIT=io5)
c
100   format(7f8.3,2f8.1,5f8.3,3f9.1,i8,x,e8.1,x,f8.6,x,a4,x,a2)
101   format(17a8                         ,a8,a8,a4,a6,a6)
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
	fmid = funcBS(x2,r,chord,twist,nprs,xp)
        f0   = funcBS(x1,r,chord,twist,nprs,xp)
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
	   fmid = funcBS(xmid,r,chord,twist,nprs,xp)
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
        phi  = atan(vwind*(1.-a)/(r*om*(1.+ap)))
        w2   =     (vwind*(1.-a))**2 + (r*om*(1.+ap))**2
        phid = 180.*phi/pi
        aoas = phid-twist
c
        FP     = FPR(x,phi)
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
c	    write(*,*)'nprs ',nprs
c
            k2 = np(nprs)
	    do j=1,k2
               aoasp(j) = aoain(nprs,j)
               clsp (j) = clin (nprs,j)
               clss (j) = cls  (nprs,j)
               cdsp (j) = cdin (nprs,j)
               cdss (j) = cds  (nprs,j)
            enddo	 
c
            call SPLINT
     +      (aoasp,clsp,clss,Np(nprs),aoas,clints(k+1),dummy)
            call SPLINT
     +      (aoasp,cdsp,cdss,Np(nprs),aoas,cdints(k+1),dummy)
	 enddo
c
c	write(*,*)'clint cdint',clint, cdint
c
        clint = xp*clints(1)+(1.-xp)*clints(2)
        cdint = xp*cdints(1)+(1.-xp)*cdints(2)
c
        cn   = clint*cos(phi) + cdint*sin(phi)  
        cta  = clint*sin(phi) - cdint*cos(phi)  
c
        ctaer  = cn
c
        cp = 4.*a*(1.-a)**2
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

