c
c      SUBROUTINES
c
c-------------------------------------------------------------------
c      BEM from 2019 coding at ETS, Montreal, 2019
c-------------------------------------------------------------------
c--------------------------------
      subroutine BEM(cP,cT)
c--------------------------------
c
      use mem
c
      real clint(2),cdint(2)
c
      real thick1, thick2, kg
c
      character*20 nout2, nout3

      character*5 pstring
c
      logical deb
      deb = .FALSE.
c
      write(*,*)'routine BEM ************'
c      write(*,104)'vwind',vwind,'pitch',glopitch 
c
      thick1 = 1.
      thick2 = 0.	
c
c     critical a value to for cT(a) correction (Glauert/Hansen)
c
      ac = 0.4
c
      om    = rpm*pi/30.
      vtip  = om*rtip
      TSR   = vtip/vwind
      write(*,*)'BEM: vwind TSR=',vwind, tsr
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
c
      thrust = 0.
      torque = 0.
      aold   = 0.3
      apold  = 0.0000
c
c     loop over sections
c
      write(*,101)'r ',' a ','ap ','F','w','chord','twist','phi ',
     +            'aoa','cL','cD','L2D','cNo','cTa','dFn','dFt','Ga',
     +            'iter','err','th','Prof'
c
      write(ioout,101)'r ',' a ','ap ','F','w','chord',
     +            'twist','phi ',
     +            'aoa','cL','cD','L2D','cNo','cTa','dFn','dFt','Ga',
     +            'iter','err','th','Prof'
c
c--------------------------------------------------------------------------
c     BEGIN section loop
c
      dr = (rtip-rroot)/float(nsec)
      do i=1,nsec 
         iter = 1
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
c        get lokal thickness
c
	 rred = r/rtip
	 th   = thicksp(rred)
c
c        ensure monotonic decrease of thickness
c
	 thick2 = th
c	 write(*,'(a20,4f12.5)')'r th1 th2 thick',r,thick1,thick2,th
	 if(thick2.gt.thick1)then
             th = thick1
	 else
	     thick1 = thick2
         endif
c
c        ensure thickness to be larger than thinnest
c
	 if(th.lt.minthick)then
c            write(*,*)'minthick used',th
	    th = minthick
         end if
c
         anew  = .3
         apnew = eps
c
cBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c        BEGIN main BEM loop  
cBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB    
c
 10      iter = iter + 1
c
c        write(*,888)'a ap old',aold,apold
888      format(a10,2f12.6)
c
         if(iter.gt.1)then
            aold  = anew
            apold = apnew
         end if   
c
c       debug
c	
	if (i.eq.21)then
           write(124,*)anew
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
c       aps 2022 apr 26: 1st profile: cylinder with th = 1.
c
123     if(k.eq.1)indpro=2      
	nameprout = namepr(k)
c
c       two values for interpolation if thickness is not close (within eps) to a given one
c
	 do k =0,1
c
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
     +      (aoasp,clsp,clss,Np(indpro),aoab,clint(k+1),dummy)
            call SPLINT
     +      (aoasp,cdsp,cdss,Np(indpro),aoab,cdint(k+1),dummy)
	 enddo
c
c*************** BEM ******************************************************
c
         dth =  prothick(indpro)-prothick(indpro-1)
 	 x   = (prothick(indpro)-th)/dth
c
         clintth = x*clint(2)+(1.-x)*clint(1)
c
	 if (clintth.gt. 3.)clintth= 3.
	 if (clintth.lt.-2.)clintth=-2.
c
         cdintth = x*cdint(2)+(1.-x)*cdint(1)
c
c-------------------------------------------------------------------------
c        begin debug
	 deb = .False.
         if (deb)then
	    if(th.gt.0.46.and.th.lt..47)then
               write(*,125)'x','th','aoab','clint(1)',
     +         'clint(2)','n1','n2','np1','np2'
               write(*,126) x,th,aoab,clint(1),clint(2),
     +         namepr(indpro-1),namepr(indpro),
     +         np(indpro-1),np(indpro)
   
125            format(9a10)
126            format(5f10.3,2a10,2i10)
            end if
         end if
c        end debug
c------------------------------------------------------------------------
c        note: ct here is "tangential" not in thrust direction
c
         cn   = clintth*cos(phi) + cdintth*sin(phi)         
         ct   = clintth*sin(phi) - cdintth*cos(phi)
c
	 cthr = cn*cos(phi) - ct*sin(phi)
	 if (cthr.gt.2) cthr=2.
c
c        calculate induce axial and tangential velocities from
c        equating momentum and aerodynamic forces
c
c        tip loss 
c
c        no tip loss
c         TL = 1.
c
c        Prandtl / Glauert
 	 x  = r/rtip
         TL = FPR(x,phi)
c
c        Burton
c         tsrloc = om*r/vwind
c         Tl = FBu(x,tsrloc,aold)
c
c        root loss
c
         TL = TL*RL(x,phi)
c	
c        from WWL '76 pp 23/39 Eqs (2.4.4), (2.4.6), (2.6.1) and (2.6.2)
c
         AA = B*chord*cn/(8.*TL*pi*r*sin(phi)*sin(phi))
         BB = B*chord*ct/(8.*TL*pi*r*sin(phi)*cos(phi))
         anew  = AA/(1.+AA)
c
c        iteration schema
c
c        some kind of Glauert extension might be appropriate:
c        From Hansen or from NREL
c
	 cthnrel = 8./9.+(4.*tl-40./9.)*anew
         cthnrel = cthnrel +(50./9.-4.*TL)*anew*anew 
c
c         if (anew.gt.ac.or.cthr.gt.cthnrel)then
         if (anew.gt.ac)then
            sigH = chord*B/(2.*pi*r)
            kG = 4.*TL*sin(phi)*sin(phi)/(sigH*cn)
            anew = ahansen(Kg,ac)
c            anew = anrel(cthr,TL)
c	    write(*,*)'ahansen', anew
         else
	    anew  = AA/(1.+AA)
         endif
c
c----------------------------------------------------------------------------------
c        some kind of Newton if fixed point does not converge after 100 iterations
c----------------------------------------------------------------------------------
c
         if (iter.gt.100)then
            dphi  = 0.001*phi
            phi2  = phi + dphi
            AA1   = AA
            AA2   = B*chord*cn/(8.*TL*pi*r*sin(phi2)*sin(phi2))
            AAp   = (AA2-AA1)/dphi
            anew  = aold - AA/AAp
	    anew  = anew/(1.+anew)
        end if
c
c-----------------------------------------------------------------
c        aprime  
c
         apnew = BB/(1.-BB)
c
c        check limits of BEM
c
c        Assume valid WT operation only for 0 < a < 1
c
	 if (anew.ge.1.0)anew = 0.9999
c
c        negative ap indicates instability from the root vortex 
c
	 if (apnew.lt.0.0)then
            anew = eps
            apnew= eps
         endif
c
c--------------------------------------------------------------
c        DEBUG
c        write(*,888)'new: a ap', anew, apnew
c
         err = abs(anew-aold)/abs(anew)
         if(iter.lt.2)                      goto 10
         if(err.gt.eps.and.iter.lt.maxiter) goto 10
c
c        END main BEM loop
c
         dT  =  0.5*dens*b*chord*w2*cn*dr
         dFt =  0.5*dens*b*chord*w2*ct*dr
         dQ  =  dFt*r
         dT   = dT/1000.
         dFt  = dFt/1000.
         dQ   = dQ/1000.
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
     +               phi*57.3,aoab,clo,cdo,clo/cdo,cn,ct,dT,dFt,gam,
     +               iter,err,th,nameprout
c
         write(ioout,100)r,anew,apnew,TL,sqrt(w2),chord,twist,
     +               phi*57.3,aoab,clo,cdo,clo/cdo,cn,ct,dT,dFt,gam,
     +               iter,err,th,nameprout
c
c       to compare with wt-perf
c
        write(io2,105)chord,th*chord,twist
c
c
c       to compare with FAST
c
        pstring = 'PRINT'
        drnodes = dr
c       +(-1)**i*0.6
        write(io3,106)r,twist,drnodes,chord,indpro,Pstring
c
         thrust = thrust + dT
         torque = torque + dQ
c
      end do
c
c     END section loop
c--------------------------------------------------------------------------
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
c
 100  format(7f8.3,2f8.1,5f8.3,3f8.2,i8,x,e8.1,x,f8.6,x,a4)
 101  format(17a8                         , x,a8,x,a8,x,a4)
 102  format(a10,f8.3,a12)
 103  format(a20,f12.4,a10)
 104  format(2(a15,f12.4))
 105  format(3f12.6)
 106  format(4f12.6,i5,a10)
c
	end
c
c------------------------------------------------------------------------
c       analyse given blade from BlaDes.in
c
c------------------------------------------------------------------------
	subroutine Analysis
c------------------------------------------------------------------------
c      
      use mem
c
      write(*,*) 'routine BLADE ANALYSIS ****************'
c
      NIN   = './BlaDes.in'
      IOIN   = 10
      OPEN(UNIT=IOIN, FILE=NIN, FORM='FORMATTED',STATUS='OLD')
C
      write(*,787)'Blade description from ',NIN
      write(*,786)'r-sec','twist','chord','airfoil name'
786   format(3a14,a14)
787   format(2a25)
c
C     Blade Description if present (BlaDes.in)
C 
      READ(IOIN,*)
c
c     read in nd lines from BladeDes
c
      Nd = 0 
 20   Nd = Nd+1
c     
c     2021 Nov 25:   iprno must be consitent with names
c     -> BlaDes and Prothick have to be "synchronous"
c
c
      READ(IOIN,*,END=21,ERR=21)rsecsp(nd),twistsp(nd),chsp(nd),
c23456
     +                           namedesin(nd)
      Write(*,'(3f14.4,2x,a14)')rsecsp(nd),twistsp(nd),chsp(nd),
     +                           namedesin(nd)
      goto 20
C
 21   close(ioin)
      nd = nd-1
      ndes = nd
      write(*,*)nd,' lines of blade description read'
c
c     add glopitch to twist
c
      do i =1,nd
         twistsp(i) =twistsp(i) + glopitch
      end do	     
c
c     check that all profile names are specified 
c
      do i = 1, nd
         do j = 1,nopr
            if (namedesin(i).eq.namepr(j)) goto 31
         enddo
31       iprno(i) = j
      end do
c
C     prepare SPLINE INTERPOLATION of blade data 
c
c     chord
c
      AB1 =  (chsp( 2) -chsp(1   ))/(rsecsp( 2)-rsecsp(   1))
      AB2 =  (chsp(Nd)- chsp(Nd-1))/(rsecsp(Nd)-rsecsp(Nd-1))
      CALL SPLINE(rsecsp,chsp,Nd,AB1,AB2,chS)
      write(*,*)'spline for chord done'
c
c     twist
c
      AB1 =  (twistsp(2)-twistsp(1  ))/(rsecsp(2)-rsecsp(1))
      AB2 =  (twistsp(Nd)-twistsp(Nd-1))/(rsecsp(Nd)-rsecsp(Nd-1))
      CALL SPLINE(rsecsp,twistsp,Nd,AB1,AB2,twistS)
      write(*,*)'spline for twist done'

      write(*,*)'Blade BLADE ANALYSIS done'
c
	return
	end
c
c---------------------------------------------------------------
c       design Mode: Betz chord and twist distribution 
c---------------------------------------------------------------
c
c---------------------------------------------------------------
	subroutine Design
c---------------------------------------------------------------
c
	use mem
	character*4 prdname
c
        write(*,*)'Routine DESIGN   *******************'
c
        IOUT  = 11
        ioutd = 12
        nout  = "des.out"
        noutd = "BlaDes.out"
C
        OPEN(UNIT=Iout, FILE=Nout, FORM='FORMATTED',STATUS='unknown')
        OPEN(UNIT=Ioutd,FILE=Noutd,FORM='FORMATTED',STATUS='unknown')
c
	write(*,'(6a12)')
     +  'rsec','twist','chord','thick','cldes','aoades'
	write(iout,'(6a12)')
     +  'rsec','twist','chord','thick','cldes','aoades'
	write(ioutd,'(4a12)')
     +  'r','Twist','Chord','Pr-Name'
c
c       define blade (chord, twist) for nsec equal sections
c
        dr   = rtip/float(nsec)
c
	do i = 1,nsec
c
c          radius of section
c
	   rloc = 0.5*dr +(i-1)*dr
	   rsecsp(i)  = rloc
c
           rr = rloc/rtip
c
c          look for local thickness
c
	   thloc = thicksp(rr)
c
c          find profile of closest thickness 
c
	   do k =1,nopr
	      if (thloc.lt.prothick(k))then
                continue
              else
                indpro = k
                goto 321
              endif
           end do
c
321        dth       = prothick(indpro)-prothick(indpro-1)
    	   x         = (prothick(indpro)-thloc)/dth
c
c          angle of attack
c
           aoades    = x*optaoa(indpro-1)+(1.-x)*optaoa(indpro)
c
c----------------------------------------------------------------------
	   aoadesmin = 0.
c----------------------------------------------------------------------
	   aoades    = max(aoades,aoadesmin)
c
c          define twist
c
	   tsrloc     = om*rloc/vwinda
c
c          get aoa = flow angle - twist
c
c----------------------------------------------------------------------
	   twmin = -5.
c--------------------------------------------------------------------
           twistloc = twist(tsrloc) - aoades
c
c           write(*,*)'tsrloc'
c           
           twistsp(i) = max(min(twmax,twistloc),twmin)
c
c          cldes
c
           cldes  = x*optcl(indpro-1)+(1.-x)*optcl(indpro)
c
c	   write(*,*)'cls ',optcl(indpro-1),cldes,optcl(indpro)
c
c	   chord
c
           chloc   = chDes(rloc,cldes) 
           chsp(i) = min(chmax,chloc)  
c
	   rout = rsecsp(i)
           rr   = rout/rtip
c
	   write(*,'(6f12.6)')
     +     rout,twistsp(i),chsp(i),thicksp(rr),cldes,aoades
	   write(iout,'(6f12.6)')
     +     rout,twistsp(i),chsp(i),thicksp(rr),cldes,aoades
c
c          find name of profile of next smallest thicknet
c
c
           thickdessec = thicksp(rr)
           prdname = ""
	   do k =1,nopr
              if(thickdessec.ge.prothick(k))then
                 prdname =namepr(k)
                 goto 888
              end if
           end do
c
888	   write(ioutd,'(3f12.3,a12)')
     +     rout,twistsp(i),chsp(i),prdname
	end do
        close(unit=iout)
        close(unit=ioutd)
c
c     do i =1,nsec
c
c
c--------------------------------------------------------------------------------
C     prepare SPLINE INTERPOLATION of blade data 
c
c     chord
c
      AB1=(chsp(2)-chsp(1  ))/(rsecsp(2)-rsecsp(1))
      AB2=(chsp(Nsec)-chsp(Nsec-1))/(rsecsp(Nsec)-rsecsp(Nsec-1))
      CALL SPLINE(rsecsp,chsp,Nsec,AB1,AB2,chS)
      write(*,*)'spline for chord done'
c
c     twist
c
      AB1=(twistsp(2)-twistsp(1  ))/(rsecsp(2)-rsecsp(1))
      AB2=(twistsp(Nsec)-twistsp(Nsec-1))/(rsecsp(Nsec)-rsecsp(Nsec-1))
      CALL SPLINE(rsecsp,twistsp,Nsec,AB1,AB2,twistS)
c
      write(*,*)'spline for twist done'
      write(*,*)'END design'
c
	return
	end
c
c-------------------------------------------------------------------
c      airfoil analysis
c-------------------------------------------------------------------
c
c-------------------------------------------------------------------
	subroutine Airfoils
c-------------------------------------------------------------------
c
	use mem
c
	real rthpr(20)
	character*20 nout2
c
	write(*,*)'Routine AIRFOILS ******************'
c
	minthick = 1.1
c
c------------------------------------------------------------------
c     read names and thick thickness of aerodynamic profiles
c------------------------------------------------------------------
c
      IOIN = 10
      nin  = "./ProThick.in"
      OPEN(UNIT=IOIN, FILE=NIN, FORM='FORMATTED',STATUS='OLD')
C
      write(*,*)'Profile name and thickness'
      write(*,781)'name','thickness','rel loc','abs loc'
781   format(4a15)  
C 
c     read one line of header !!
c
      READ(IOIN,*)
c----------------------------------------------------------
c     read in nd lines from ProThick.in
c
      IOut = 13
      nout = "th.out"
      OPEN(UNIT=IOut,FILE=Nout,FORM='FORMATTED',STATUS='unknown')
c
      Nd = 1 
 50   READ(IOIN,*,END=51,ERR=51)namepr(nd),prothick(nd)
c      write(*,*)'nd ',nd, namepr(nd),prothick(nd)
c
c     set minthick
c
      if (prothick(nd).lt.minthick)minthick=prothick(nd)
c
c     find location of profile from thickness distribution
c
      dr = 0.001
      nn = int(1./dr)
      do i=0,nn
        rth = dr*i
        thcomp=thicksp(rth)
        if(prothick(nd).gt.thcomp)goto 777
      end do
c
777   rthpr(nd) = rth
      Write(*   ,'(a15,3f15.3)') namepr(nd),prothick(nd),rth,rth*rtip
      Write(iout,'(a15,3f15.3)') namepr(nd),prothick(nd),rth,rth*rtip
c
      Nd = Nd+1
      goto 50
C
 51   nd   = nd-1
      nopr = nd
      write(*,*)nopr,' lines of profile description read'
      Write(*,*)
      close(unit=iout)
c
c------------------------------------------------------------------------------
c     analysis and interpolaton of aerodynamic data 
c------------------------------------------------------------------------------
c
      Iout2 = 15
      Nout2 = 'ProProp.out'
      OPEN(UNIT=IOut2,FILE=Nout2,FORM='FORMATTED',STATUS='unknown')
c
      write(*,*)'profile parameters:'
      write(*    ,'(a8,9a12)')
     +     'profile','location','L2D-max','aoa-des','cL-des',
     +     'AOA-max','cL-max','aoaLzero','slopeL0','No Po'
      write(iout2,'(a8,9a12)')
     +     'profile','location','L2D-max','aoa-des','cL-des',
     +     'AOA-max','cL-max','aoaLzero','slopeL0','No Po' 
c
c     nopr = number of profiles
      do i =1,nopr
         IOIN   = 10
         NIN   = namepr(i) // ".aer"
         OPEN(UNIT=IOIN,FILE=NIN,FORM='FORMATTED',STATUS='OLD')
         Np(i) = 0 
 30      Np(i) = Np(i)+1
         READ(IOIN,*,END=31,ERR=30)aoain(i,np(i)),clin(i,np(i)),
     +                             cdin(i,np(i)), dummy
c
c	 write(*,99)aoain(i,np(i)),clin(i,np(i)),cdin(i,np(i))
99       format(3f10.4)
c
         goto 30
c
31       np(i) = np(i)-1
c	 write(*,*)'airfoil: ',namepr(i),np(i),' lines'
c        np(i) number of lines in profile "i" 
c
c---------------------------------------------------------------
c        search for aoa at which cL/cD = max 
c--------------------------------------------------------------
c
	 gzmax = -1000
	 do j = 1,np(i)
	    gz = clin(i,j)/cdin(i,j)
	    if (gz.gt.gzmax)then
               gzmax = gz
               indopt = j
            end if
	 end do
c
	optcl(i) =  clin(i,indopt)
        optaoa(i)= aoain(i,indopt)
c
	if(abs(optcl(i)).lt.1.e-3)optaoa(i)=0.
c
c---------------------------------------------------------------
c        search for aoa at which cL = max 
c--------------------------------------------------------------
c
        clmax = -1.
	 do j = 1,np(i)
	    if (clin(i,j).gt.clmax)then
               clmax = clin(i,j)
               indmax = j
            end if
	 end do
c
	maxcl(i) =  clin(i,indmax)
        maxaoa(i)= aoain(i,indmax)
c
c---------------------------------------------------------------
c        search for aoa at which cL = 0. 
c--------------------------------------------------------------
c
	 do j = 1,np(i)
            if(aoain(i,j).gt.-20.and.aoain(i,j).lt.20.)then
	       if (clin(i,j).lt.0.)then
                  indzero = j
               end if
             endif
	 end do
c
	if(indzero.lt.1)indzero = 1
        zeroaoa(i)= aoain(i,indzero)
        zeroslope(i) =
     +  (clin(i,indzero+1)-clin(i,indzero))/
     +  (aoain(i,indzero+1)-aoain(i,indzero))
c
c---------------------------------------------------------------    
c       output
c-----------------------------------------------------------------------
c
        optprrad = rtip*rthpr(i)
        write(*    ,123)nin,optprrad,gzmax,optaoa(i),optcl(i),
     +  maxaoa(i),maxcl(i),zeroaoa(i),zeroslope(i),np(i)
        write(iout2,123)nin,optprrad,gzmax,optaoa(i),optcl(i),
     +  maxaoa(i),maxcl(i),zeroaoa(i),zeroslope(i),np(i)
123	format(a8,8f12.4,i12)
c
C        SPLINE INTERPOLATION 
C
c	 write(*,*) 'prepare SPLINE data for Polars '
c
c     cl 
c
         AB1 =  (clin(i,2)-clin(i,1  ))/(aoain(i,2)-aoain(i,1))
         AB2 =  (clin(i,Np(i))-clin(i,np(i)-1))/
     +          (aoain(i,Np(i))-aoain(i,Np(i)-1))
	 do j=1,np(i)
            aoasp(j)=aoain(i,j)
            clsp(j) =clin(i,j)
         enddo
c
         CALL SPLINE(aoasp,clsp,Np(i),AB1,AB2,clSs)
	 do j=1,np(i)
            cls(i,j) = clss(j)
         enddo
C
c     cd
c
         AB1 =  (Cdin(i,2)-CDin(i,1  ))/(aoain(i,2)-aoain(i,1))
         AB2 =  (Cdin(i,Np(i)) -Cdin(i,Np(i)-1))/
     +          (aoain(i,Np(i))-aoain(i,Np(i)-1))
	 do j=1,np(i)
            aoasp(j)=aoain(i,j)
            cdsp(j) =cdin(i,j)
         enddo
         CALL SPLINE(aoasp,Cdsp,Np(i),AB1,AB2,CdSs)
	 do j=1,np(i)
            cds(i,j) = cdss(j)
         enddo
c
      enddo
c     do i=1,nopr
c
      write(*,*)'minthick = ',minthick
c
      write(*,*)'PROFILE Analysis done'
c
      return
      end
c---------------------------------------------
c---------------------------------------------
c   xa = rtip - tiplen
c   ya = chord(xa)
c   yas= slope at xa
c   x (input) : location of chord
c   y (output): chord (x)
c   (1) elliptic shape = infinite slope at rtip
c
c---------------------------------------------
	subroutine TipShapeEl(xa,ya,yas,x,y)
c---------------------------------------------
	use mem
c
	real aa,bb,xm
c
        bb = sqrt(ya**2 +(yas*ya)**2)
        xm = (bb*xa+rtip*ya*yas)/(ya*yas+bb)
        aa = rtip-xm
c
        y = bb*sqrt(1.-((x-xm)/aa)**2)
c
	return
	end	
c---------------------------------------------------    
c   (2) parabolic shape from chord and smooth slope at xa = rtip - tiplen
c       chord at rtip = 0.  
c----------------------------------------------------
        subroutine TipShapePara(xa,ya,yas,x,y)
c----------------------------------------------------
	use mem
c
	real aa,bb,cc
c
	aa = ya-yas*(xa-rtip)
	aa = aa/((xa**2-rtip**2)-2.*xa*(xa-rtip))
c
	bb = yas - 2.*aa*xa
	cc = ya  - aa*xa**2-bb*xa
c
        y = aa*x*x+bb*x+cc
c
	return
	end
c
c--------------------------------------------------
	subroutine Bezier(bb,xt,npb)
c--------------------------------------------------
c
c npb number of sections (points)
c bb: input matrix(2,npb)
c xt: smoothed output same dimensions
c
c--------------------------------------------------
	use mem
c
	real bb(2,0:npb),xt(2,0:npb)
c
	dt =1./float(npb)
        do i =0,npb
           t = float(i)*dt
           do k=1,2
              xt(k,i) = 0.
 	      do j=0,npb
	         xt(k,i)= 
     +           xt(k,i) + fnoverk(npb,j)*t**j*(1.-t)**(npb-j)*bb(k,j)
              end do
           end do
c          write(*,100)   t,(xt(k,i),k=1,2)
c          write(iout,100)(xt(k),k=1,2)
        enddo
c
100	FORMAT(3F12.6)
c
	return
        end
c
c--------------------------------------------------
c
      subroutine smoothCHORD
c	
      use mem
c
      real ab1,ab2
c
      write(*,*)'routine SMOOTH CHORD *********'
c
c     smooth chord
c
      call Bezier(chsp,chspsm,nsec)
c
      do i =1,nsec
         chsp(i)=chspsm(i)
      end do
c
C     prepare SPLINE INTERPOLATION of smoothed chord
c
      AB1=(chsp(2)   -chsp(1     ))/(rsecsp(2)   -rsecsp(1))
      AB2=(chsp(Nsec)-chsp(Nsec-1))/(rsecsp(Nsec)-rsecsp(Nsec-1))
      CALL SPLINE(rsecsp,chsp,Nsec,AB1,AB2,chS)
c
      return
      end
c
c--------------------------------------------------
c
      subroutine smoothTWIST
c	
      use mem
c
      real ab1,ab2
c
      write(*,*)'routine SMOOTH TWIST *********'
c
      nt = nsec
c
c     smooth chord
c
      call Bezier(twistsp,twistspsm,nt)
c
      do i =1,nt
         twistsp(i)=twistspsm(i)
      end do
c
C     prepare SPLINE INTERPOLATION of smoothed chord
c
      AB1=(twistsp(2) -twistsp(1   ))/(rsecsp(2) -rsecsp(1))
      AB2=(twistsp(Nt)-twistsp(Nt-1))/(rsecsp(Nt)-rsecsp(Nt-1))
      CALL SPLINE(rsecsp,twistsp,Nt,AB1,AB2,twistS)
c
      return
      end
c
c--------------------------------------------------
c
      subroutine improveTWIST
c	
c     2022 May  23: twist limiter changed
c     2022 June 14: implementation of Jens' 
c
c--------------------------------------------------
      use mem
c
      real ab1,ab2
c
      write(*,*)'START improve TWIST *********'
c
      nt = nsec
c
c     use from BEM: d phi = -1/(tsr_loc (1+ a')) d a
c
      write(*,100)'r','a-old','da','dphi'
c
      do i =1,nt
         rr = rsecsp(i)
         x  = rr/rtip
         da     =  aopt(DesSchema,x)- abem(i)
         tsrloc = tsr*rr/rtip
         dphi   =-da/(tsrloc*(1.+apbem(i)))
         dphi   = 180.*dphi/pi
c
c        2022 06 17: no negative twist change
c
         if(dphi.lt.0.)dphi=0.
         twistsp(i)= twistsp(i) + dphi
c
c        25 05 2022: mintwist arbitrarily set to -7°
c
	 if(twistsp(i).lt.-7.)twistsp(i) = -7.
	 write(*,101)rsecsp(i),abem(i),da,dphi
      end do
      write(*,*)
c
C     update SPLINE INTERPOLATION data 
c
      AB1=(twistsp(2) -twistsp(1   ))/(rsecsp(2) -rsecsp(1))
      AB2=(twistsp(Nt)-twistsp(Nt-1))/(rsecsp(Nt)-rsecsp(Nt-1))
      CALL SPLINE(rsecsp,twistsp,Nt,AB1,AB2,twistS)
c
      write(*,*)'END improve TWIST *********'
c
100	format(4a12)
101	format(4f12.6)
c
      return
      end
c
c----------------------------------------------------------------------
	subroutine ThickPrep(filename)
c----------------------------------------------------------------------
c
	use mem
	character*30 filename
c
      Intemp  = 19
      OPEN(UNIT=intemp,FILE=filename,Form='formatted',status='unknown')
      irth = 1
      write(*,'(2a10)')'R','Thick'
99    read(intemp,*,END=98,ERR=98)dthr(irth),dthth(irth)
	  write(*,'(2f10.4)')dthr(irth),dthth(irth)
         irth = irth+1
         goto 99
98    Close(UNIT=intemp)
      irth = irth-1
      write(*,*)'From ',nread,' : ',irth,' data read'
c
      AB1=(dthth(2)   - dthth(1))     /(dthr(2)   -dthr(1))
      AB2=(dthth(irth)- dthth(irth-1))/(dthr(irth)-dthr(irth-1))
      CALL SPLINE(dthr,dthth,irth,AB1,AB2,dthsp)
c
      write(*,*)'thickness spline GENERATED'
c
	return
	end
c
c-----------------------------------------------------
c-----------------------------------------------------
c     
c     FUNCTIONS
c
c-----------------------------------------------------
c     Prandtl's F
c     from (for example) Hansen (3rd Ed) p 44  
c-----------------------------------------------------
c
c--------------------------
      function FPR(x,phi)    
c--------------------------
c
      use mem
c
c     Eq. (6.34) and WWL76 3rd Eq on page 35
      f   = 0.5*B*(1.-x)/(x*sin(phi))
c     Eq. (6.33)
      FPR = 2./pi*acos(exp(-f))
c
      return
      end
c
c-----------------------------------------------------
c     Burton's F
c     from Ramdin (MSc Thesis, 2017) p 28  
c-----------------------------------------------------
c
c--------------------------------
      function FBu(x,tsrl,a)    
c--------------------------------
c
      use mem
c
      sq = sqrt(1.+(tsrl*x/(1.-a))**2)
      f   = 0.5*B*(1.-x)*sq/x
c     
      FBu = 2./pi*acos(exp(-f))
c
      return
      end
c
c-----------------------------------------------------
c     Burton's  root loss
c     page 84 Eq (3.80) 
c     
c-----------------------------------------------------
      function RL(x,phi)
c
      use mem
c
      xr  = rroot/rtip
      if (x.gt.xr)then
         f   = 0.5*B*(x-xr)/(x*sin(phi))
         RL  = 2./pi*acos(exp(-f))
      else
         RL = 1.e-4
      endif
c
      return
      end
c
c-----------------------------------------------------
c     thickness distribution of blade
c-----------------------------------------------------
c
c---------------------------------
      function thickPoly(rred)    
c---------------------------------
c
      real k(0:3)
c
c     initial cubic regression polynom for thickness from 4 points:
c     r = 1    th =  .17
c     r = 2/3  th =  .21
c     r = 1/3  th =  .4 (thickest meaningful profile)
c     r = 0    th = 1.  (cylinder)
c
c     0 \le rred \le 1 as well as thick 
c
      k(0) =  1. 
      k(1) = -2.6782
      k(2) =  3.023 
      k(3) = -1.1748
c
      thick = 0.
      do i=0,3
        thick = thick + k(i)*rred**i
      end do
c
      return
      end
c
c---------------------------------
      function thicksp(rred)    
c---------------------------------
c
      use mem
c
c	write(*,*)'AB1 und 2 ',AB1,AB2
c	write(*,*)'r 12 th 12 ',rth(1),rth(2),th(1),th(2)
c
      call SPLINT(dthr,dthth,dthSp,irth,rred,thint,dummy)
      thicksp = thint
      return
      end
c
c-----------------------------------------------------
c     twist from Betz
c-----------------------------------------------------
c
      function twist(tsrloc)    
c
	use mem
c
      twist = atan(2./(3.*tsrloc))      
      twist = 180.*twist/pi
c
      return
      end
c
c-------------------------------------------------------------
c     chord from Betz or Schmitz and geometric tip reduction
c--------------------------------------------------------------
c
      function chDes(rloc,cldes)    
c
      use mem
c
      real chloc1, chloc2, xa, ya, dr, yas
      real chtip
c
      dr = 1.e-3
c
      chloc1 = chBe(rloc,cldes)
      xa     = rtip-tiplen
      ya     = chloc1
      dr     = eps
      chloc2 = chBe(rloc+dr,cldes)
      yas    = (chloc2 - chloc1)/dr
c
      chDes = chloc1
c
      itip = 0
      if(tipshape.eq.'E')itip = 1
      if(tipshape.eq.'P')itip = 2
c
      if(rloc.gt.xa) then
         SELECT CASE(ITIP)
            CASE(1)
               call TipShapeEl(xa,ya,yas,rloc,chtip)
               chDes = chtip
            CASE(2)
               call TipShapePara(xa,ya,yas,rloc,chtip)
               chDes = chtip
            CASE DEFAULT
               chDEs = chloc1
         end SELECT
      endif
c
      return
      end
c
c------------------------------------------------------------
c    see Gasch, 8. Ed, p 206 Eq. 5.64 (Betz)
c------------------------------------------------------------
      function chBe(rloc,cldes)    
c------------------------------------------------------------
c
      use mem
c
      chBe = 16.*pi*rtip/(9.*B*cldes*tsr)
c      write(*,888)'B cldes tsr', B,cldes,tsr
      chBe = chBe/(sqrt(4./9. +(TSR*rloc/rtip)**2))
c  
c      write(*,888)'om TSR chbe ',om, tsr, chbe
c
888   format(a30,3f8.3)
      return
      end		
c
c-----------------------------------------------------------
c    see Gasch, 8. Ed, p 206 Eq. 5.65 (Schmitz)
c-----------------------------------------------------------
      function chSc(rloc,cldes)    
c-----------------------------------------------------------
c
      use mem
c
      al1  = atan(rtip/(tsr*rloc))
      chSc = 16.*pi*rloc/(B*cldes)
      chSc = chSc*(sin(al1/3.))**2
c
      write(*,*)'rloc, cldes chsc',rloc,cldes,chSc
c
      return
      end		
c
c---------------------------------------
c       numbers grow very fast therefore shift to FLOAT
c--------------------------------------
	real function fnoverk(n,k)
c-------------------------------------
c
	integer n,k
c
	if(2*k.gt.n)then k=n-k
        fnoverk = 1.
c
	do i=1,k
           t = Float(n+1-i)/Float(i)
           fnoverk = fnoverk*t
	end do
c
c	write(*,*)'n over k',fnoverk
c
	return
	end
c
c----------------------------------------------------------
c       glauert correction from Hansen Eqs (8.18/44/45) 
c       for a > ac
c----------------------------------------------------------
	real function ahansen(Kg,ac)
c----------------------------------------------------------
c
	real kg, ah
c
        ah = kg*(1.-2.*ac)+2.
        ah = ah*ah+4.*(kg*ac*ac-1.)
        if(ah.lt.0.)ah=0.
        ah = sqrt(ah)
        ahansen = 0.5*(2.+kg*(1.-2.*ac)-ah)
c
	return
        end
c
c
c----------------------------------------------------------
c       glauert correction from Nrel TP 36881 Eq 15 p 7  
c       for a > ac
c----------------------------------------------------------
	real function anrel(ct,F)
c----------------------------------------------------------
c
	real ct, F
c
c        anrel = ct*(50.-36*F)+12.*F*(3.*F-4.)
c        anrel = (18.*F-20.-3.*sqrt(anrel))/(36.*F-50.)
        anrel = 25.*ct-24.*F-18.*cT*F+18.*F*F
        if(anrel.lt.0)anrel=0.
        anrel = -20.+18.*F-3.*sqrt(2.)*sqrt(anrel)
        anrel = anrel/(36.*F-50.)
c
c	write(*,*)'ct F anrel ',ct, f, anrel
c
	return
        end
c
c----------------------------------------------------------
c       twist bend
c----------------------------------------------------------
	subroutine twistbend(vw,rr,adpi)
c----------------------------------------------------------
c
	use mem
c
	real  mr, mwind
c
	adpi = 0.
c
c	find index of nearest location
c
	do i =0,12
          rtw1 = rstw(i)  +rroot
          rtw2 = rstw(i+1)+rroot
          if(rr.gt.rtw1.and.rr.lt.rtw2)goto 111
        end do
c
111     nrend = i+1
        nranf = i
        ranf  = rstw(nranf)
        rend  = rstw(nrend)
c
c       find index of nearest windspeed
c
        nwanf = aint(vw)
        nwend = nwanf + 1
        vanf  = float(nwanf)
        vend  = float(nwend)
c
c       2D linear interpolation according to wiki
c
        fq11 = eltw(nranf,nwanf)
        fq12 = eltw(nranf,nwend)
        fq21 = eltw(nrend,nwanf)
        fq22 = eltw(nrend,nwend)
c
	adpi =        fq22*(rr  -ranf)*(vw  -vanf)
c	write(*,*)'1',adpi
        adpi = adpi + fq12*(rend-rr  )*(vw  -vanf)
c	write(*,*)'2',adpi
        adpi = adpi + fq21*(rr  -ranf)*(vend-vw)
c	write(*,*)'3',adpi
        adpi = adpi + fq11*(rend-rr  )*(vend-vw)
c	write(*,*)'4',adpi
        adpi = adpi/((rend-ranf)*(vend-vanf))
c
c	write(*,*)'twist bend'
c	write(*,*)'ra re va ve',nranf, nrend, nwanf, nwend
c	write(*,*)'11 12 21 22',fq11,fq12,fq21,fq22
c	write(*,*)'r1 r2 v1 v2',ranf,rend,vanf,vend
c	write(*,*)'rr vw adpi (5)',rr,vw,adpi
c
	return
        end
c
c----------------------------------------------------------
c       a opt
c----------------------------------------------------------
	real function aopt(is,x)
c----------------------------------------------------------
c
	use mem
c
	real k(0:4)
c
	k(0) =  0.2502
	k(1) =  0.5297
	k(2) = -1.0601
	k(3) =  0.6539    
	k(4) =  0.026
c
        select case(is)
          case(1)
c
c         Betz: a = 1 / 3
c
             aopt = 1./3.
          case(2)
c
c         Jens (2022)
c
             aopt = 0.
             do i=0,4
                aopt = aopt+k(i)*x**i
             end do   
             if (aopt.lt.0.25)aopt =0.25
             if (aopt.gt.0.40)aopt =0.40
          case(3)
c
c         trial for Gamma = const 
c
c              write(321,*)'v rt pi b tsr ',vwind,rtip,pi,b,tsr
c              normalized Gamma by vin * rtip
c
              gaopt = 8.*pi/(9.*b*tsr)
c              write(321,*)'gaopt = ',gaopt
              tsrl = tsr*x
c
c             this is a very crude approximation
c             valid only if ap << 1
c             a limiter may be appropriate
c
              ap = 2./(9.*tsrl*tsrl)
              if(ap.gt.1.)ap=1.
              cc = tsr*B*gaopt*(1.+ap)/(4.*pi)
c              write(321,*)'tsrl ap cc = ',tsrl,ap,cc
              if(cc.lt..25)then
                 aopt = .5-sqrt(0.25-cc)
              else
c                 write(321,*)'cc failed'
                 aopt = 0.5
              endif
              if (aopt.lt.0.  )aopt =0.0
              if (aopt.gt.0.50)aopt =0.5
c              write(321,*)'aopt = ',aopt
          case default   
             aopt = 1./3.
        end select
	return
        end
c




