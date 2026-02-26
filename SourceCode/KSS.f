C
C     ÉTS 02 June 2019
C     BEM Code for NREl Baseline
c     
c     extended to KSS (aerodynamic design)
c
c     guided by lecture notes from M. Korjahn
c     Senvion / SGRE / bewind
c
c
c     12 Nov 2021 start coding: more modern features
c     16 Nov      thickness interpolation of cL, cD
c                 1. design loop from Betz
c     21 Jan 2022 Enforce strict monotonic decreasing thickess; hub -> tip
c                 (output AOA,cL)(max)
c
c     01 Feb 2022 limit lower value of thickness to thinnest profile
c
c     07 Feb 2022 twist as function of location and wind
c
c     17 Feb 2022 ap < 0 near hub
c                 multiple solutions 
c
c     24 Feb 2022 V1 incl. report
c
c     06 Jul 2022 V2 incl. report
c                 features have been added, fore example
c                 BEM interation shifts to Newton-Raphson after 100 not succesful 
c                 fixed-point iteration
c                 tested for: SWT (r-tip = 5.1m), NREL 5MW, DTU 10MW, NREL 15 MW,
c                             r-tip= 110, 124 and 133 m)
c
c         Nov/Dec implementation of binary search
c
c      20 Nov 2025 V4 optimus syria and minor changes
c	report at overleaf
c      02 Dec 2025
c       root chord
c
c------------------------------------------------------------------------------------------
	program KSS
c
	use mem
c
	character*10 inpstring
	character*30 filename
        character*20 nout2,nout3,nout4,nout5,nout6
c
        real cpmax,tsrmax,pitchmax
c
c-----------------------------------------------------------------
c      allocate memory
c----------------------------------------------------------------
c
	nsp = 300
        allocate(rsecsp(1:nsp)   ,stat=status)
	allocate(twistsp(1:nsp)  ,stat=status)
        allocate(chsp(1:nsp)     ,stat=status)
        allocate(chs(1:nsp)      ,stat=status)
        allocate(twists(1:nsp)   ,stat=status)
        allocate(chspsm(1:nsp)   ,stat=status)
        allocate(twistspsm(1:nsp) ,stat=status)
c
	npol=300
        allocate(aoasp(1:npol)  ,stat=status)
        allocate(clsp(1:npol)   ,stat=status)
        allocate(cdsp(1:npol)   ,stat=status)

        allocate(clss(1:npol)   ,stat=status) 
        allocate(cdss(1:npol)   ,stat=status)

        allocate(abem(1:npol)   ,stat=status) 
        allocate(apbem(1:npol)  ,stat=status)	
c
	call date_and_time(date,timea,zone,values)
c
	write(*,*)
	write(*,*)'================================================ '
        write(*,*)'KSS  START ***** V5 *****                        '
	write(*,'(a6,a4,a1,a2,a1,a2,a6,a10)')
     +  'date ',date(1:4),' ',date(5:6),' ',date(7:8),' time ',timea

	write(*,*)'================================================ '
	write(*,*)
c
c----------------------------------------------------------------------
c     basic numerical data
c----------------------------------------------------------------------
c
      pi      = 4.*atan(1.)
      eps     = 1.e-5
      maxiter = 200
c
      cpmax = -1.
      ctmax = -1.
c
c     number of BEM-iterations performed
c
      ibem = 0	
c
c        Assure valid WT operation: only for 0 < a < 1
c
	 if (anew.ge.1.0)anew = 0.9999
c
c----------------------------------------------------------------------
c     read in machine and job data
c----------------------------------------------------------------------
c
      Intemp  = 19
      nread   = './Machine.in'
      OPEN(UNIT=intemp,FILE=Nread,Form='formatted',status='unknown')
c
c----------------  analyse/desgin parameter -----------------------------------
c
         read(intemp,*)inpstring, Desmode
         read(intemp,*)inpstring, B
         read(intemp,*)inpstring, dens
         read(intemp,*)inpstring, rroot
         read(intemp,*)inpstring, rtip
         read(intemp,*)inpstring, vwinda, vwinde, nwind
         read(intemp,*)inpstring, rpm
         read(intemp,*)inpstring, glopitcha,glopitche, npitch
         read(intemp,*)inpstring, nsec
c
c----------------  purely design parameter -----------------------------------
c
         read(intemp,*)inpstring, tipshape
         read(intemp,*)inpstring, tiplen
         read(intemp,*)inpstring, twmax
         read(intemp,*)inpstring, chmax
         read(intemp,*)inpstring, ImpChord
         read(intemp,*)inpstring, ImpTwist, noimptw
         read(intemp,*)inpstring, ImpThick
c
c----------------  twist bend flag -----------------------------------
c
         read(intemp,*)inpstring, TwistB
c
         read(intemp,*)inpstring, DesSchema

c
c----------------  root chord -----------------------------------
c
         read(intemp,*)inpstring, chroot
c
c
      Close(UNIT=intemp)
c
c     derive some more data
c
      ar   = pi*rtip**2
c
      write(*,*)'General Data from Machine.in'
      write(*,*)  'DesMode     ',DesMode
      write(*,*)  'TwistBend  ',TwistB
      write(*,108)'windrange =', vwinda, vwinde, nwind
      write(*,*)vwinda, vwinde, nwind
      write(*,102)'RPM       = ',rpm
      write(*,102)'B         =', B
      write(*,102)'Rtip      =', Rtip
      write(*,102)'Tilength  =', tiplen
      write(*,108)'Pitch     = ',glopitcha,glopitche,npitch
      write(*,106)'Tipshape  =', tipshape
      write(*,102)'max twist = ',twmax
      write(*,102)'max chord = ',chmax
      write(*,105)'NSEC      = ',Nsec
      write(*,105)'Des Schema= ',DesSchema
c 
      write(*,*)  'ImpTwist    ',ImpTwist
      write(*,*)  'ImpThick    ',ImpThick
      write(*,*)  'chroot      ',chroot
c
      write(*,*)
c
c
c-----------------------------------------------------------------------
c      read first pre-definded thickness distribution
c-----------------------------------------------------------------------
c
      filename   = './ThickDis.in'
c      write(*,*)'data from ThickDis.in '
      call ThickPrep(filename)
c
c----------------------------------------------------------------------
c       analyse a given set of airfoils
c-----------------------------------------------------------------------
c
	call airfoils
c
c-----------------------------------------------------------------------
c       chose mode:
c
c       ANALYSIS or DESIGN 
c
c-----------------------------------------------------------------------
c
c     ANALYZE a GIVEN BLADE
c
      if(.not.DesMode)then
c
c       read twist bend data
c
	if(twistb) then
c
	write(*,*)
	write(*,*)'***** enter twistb routine *****'
c
c       sectional radius (m) at which torsional deformation is given
c
c          Type 1: 12 sections
c          Type 2: 13
c          Type 3: 11
c          Type 4: 40
c-----------------------------------------------------------
c          Type 4 IEA 22MW RWT: 40 sections
c
	   nseltw =40
c
	rstw( 1 )  =   0.000
	rstw( 2 )  =   2.831
	rstw( 3 )  =   7.080
	rstw( 4 )  =  10.835
 	rstw( 5 )  =  14.588
 	rstw( 6 )  =  18.340
 	rstw( 7 )  =  22.090
 	rstw( 8 )  =  25.838
 	rstw( 9 )  =  29.586
	rstw(10 )  =  33.333
 	rstw(11 )  =  37.080
 	rstw(12 )  =  40.827
 	rstw(13 )  =  44.573
 	rstw(14 )  =  48.320
 	rstw(15 )  =  52.066
 	rstw(16 )  =  55.813
 	rstw(17 )  =  59.560
 	rstw(18 )  =  63.308
 	rstw(19 )  =  67.055
 	rstw(20 )  =  70.803
 	rstw(21 )  =  74.786
 	rstw(22 )  =  78.769
 	rstw(23 )  =  82.753
 	rstw(24 )  =  86.738
 	rstw(25 )  =  90.724
 	rstw(26 )  =  94.712
 	rstw(27 )  =  98.702
 	rstw(28 )  = 102.694
 	rstw(29 )  = 106.689
 	rstw(30 )  = 110.687
 	rstw(31 )  = 114.690
 	rstw(32 )  = 118.697
	rstw(33 )  = 122.710
 	rstw(34 )  = 126.729
 	rstw(35 )  = 130.756
 	rstw(36 )  = 134.793
 	rstw(37 )  = 139.110
 	rstw(38 )  = 140.553
 	rstw(39 )  = 141.276
 	rstw(40 )  = 142.000
c
c          Type 3 GW110 : 11 sections
c
c	   rstw(0) =   0.
c	   rstw(1) =  10.5
c	   rstw(2) =  20.0
c	   rstw(3) =  30.0
c	   rstw(4) =  40.
c	   rstw(5) =  50.
c	   rstw(6) =  60.
c	   rstw(7) =  70.
c	   rstw(8) =  80.
c	   rstw(9) =  90.
c	   rstw(10)=  100.
c	   rstw(11)= 110.5
c
c          hub radius added later
c
c          Type 2 Sinoma 124: 13 sections
c
c	   rstw(1) =  10.5049
c	   rstw(2) =  20.514
c	   rstw(3) =  30.5191
c	   rstw(4) =  40.5222
c	   rstw(5) =  50.5245
c	   rstw(6) =  60.5258
c	   rstw(7) =  70.5265
c	   rstw(8) =  78.0273
c	   rstw(9) =  82.0276
c	   rstw(10)=  90.0286
c	   rstw(11)= 100.036
c	   rstw(12)= 110.107
c	   rstw(13)= 122.806
c	   rstw(14)= 125.0
c
c          Type 1 Goldwind 124: 12 sections
c
c          rstw(0)  =  0.0
c  	   rstw(1)  = 10.5022
c	   rstw(2)  = 20.5038
c	   rstw(3)  = 30.5071
c	   rstw(4)  = 40.5143
c	   rstw(5)  = 50.525
c	   rstw(6)  = 60.5403
c	   rstw(7)  = 70.5587
c	   rstw(8)  = 80.5793
c	   rstw(9)  = 90.5995
c	   rstw(10) =100.625
c	   rstw(11) =110.659
c	   rstw(12) =122.231
c	   rstw(13) =125.0
c
           Intemp  = 19
           nread   = './El-Twist.dat'
           OPEN(UNIT=intemp,FILE=Nread,Form='formatted',
     +     status='unknown')
c
c          2022 April 26:  k = int(winspeed)
c                          icutout = 23
c
	   icutout = 25
	   do k =3,icutout
              eltw(0,k) = 0.0
           end do
	   write(*,*)'eltw r = 0'
c
c         i # sections    = nseltw = 40
c         j # wind-speeds = icutout= 23 
c
	   do i =1,nseltw
               read(intemp,*) (eltw(i,j),j=3,icutout)
c	       write(*,'(i3,23f6.2)')i,(eltw(i,j),j=3,icutout)
           end do           
           close (unit=intemp)
c
           do k =3,icutout
              del = (eltw(nseltw-1,k)-eltw(nseltw-2,k))/
     +              (rstw(nseltw-1)-rstw(nseltw-2))
              del = (rstw(nseltw)-rstw(nseltw-1))*del
              eltw(nseltw+1,k) = eltw(nseltw,k)+del
           end do
	   write(*,*)'eltw ',nseltw+1,' data set'
c
           write(*,*)'*****************'
           write(*,*)'eltw ',nseltw,' data read'
c
	   do i =0,nseltw+1
               write(*,888)(eltw(i,j),j=3,icutout)
           end do    
c
	write(*,*)'***** twistb routine left *****'
	write(*,*)
c
        end if
c
c       end if Twist Bend
c
888     format(23f6.2)  
c
         nout     = './cP-tsr-pi.gnu'
         nout3    = './cT-tsr-pi.gnu'
         nout6    = './cQ-tsr-pi.gnu'
         iouttemp  = 22
         iouttemp3 = 23
         iouttemp6 = 26
         OPEN(UNIT=iouttemp,FILE=Nout,
     +   Form='formatted',status='unknown')
         OPEN(UNIT=iouttemp3,FILE=Nout3,
     +   Form='formatted',status='unknown')
         OPEN(UNIT=iouttemp6,FILE=Nout6,
     +   Form='formatted',status='unknown')
c
c        loops for wind- and pitch
c
         dvwind = (vwinde-vwinda)/float(nwind)
         dpitch = (glopitche-glopitcha)/float(npitch)
c          write(*,*)glopitche, glopitcha,npitch,dpitch
c          write(iouttemp,109)'tsr','pitch','cP'
         do j = 1,npitch
            glopitch = glopitcha +(j-1)*dpitch
            call Analysis
c	     write(*,*)'KKS gp ',glopitch
            do i = 1,nwind
               vwind = vwinda +(i-1)*dvwind
               om   = rpm*pi/30.
               vtip = om*rtip
               TSR  = vtip/vwind
               write(*,*)'TSR ',TSR
c
cBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
               call bem(cP,cT,errb)
               ibem = ibem +1
cBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c
c              sort out some strange values
c
	       if(cp.lt.0)  cp = 0.
	       if(cp.gt.0.6)cp = 0.
	       if(ct.lt.0.) ct = 0.
	       if(ct.gt.2.) ct = 2.
c
c               if(errb.gt.10.)cP = 0.
c
c             search only "around" maximum
c             avoid turbulent wake state etc
c
	      tsrmin =  1.
              tsrmax = 15.
c
               if (tsr.gt.tsrmin.and.tsr.lt.tsrmax)then
	          if(cp.gt.cpmax.and.cT.gt.0.and.cp.lt..55)then
                     cpmax      = cp
                     tsrmaxcp   = tsr
                     pitchmaxcp = glopitch
                  end if
	          if(ct.gt.ctmax)then
                     ctmax      = ct
                     tsrmaxct   = tsr
                     pitchmaxct = glopitch
                  end if
               endif
c
	       write(iouttemp ,107)tsr,glopitch,cp
	       write(iouttemp3,107)tsr,glopitch,ct
	       write(iouttemp6,107)tsr,glopitch,cp/tsr
            end do
         write(iouttemp ,107)
         write(iouttemp3,107)
         end do
         close(Unit=iouttemp)
         close(Unit=iouttemp3)
         close(Unit=iouttemp6)
       endif
c
c      END not.design
c
         nout2     = './cP-max.dat'
         iouttemp2 = 23
         OPEN(UNIT=iouttemp2,FILE=Nout2,
     +   Form='formatted',status='unknown')
c
         nout4     = './cT-max.dat'
         iouttemp4 = 24
         OPEN(UNIT=iouttemp4,FILE=Nout4,
     +   Form='formatted',status='unknown')
c
       write(*,*)
       write(*,*)'MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM'
       write(iouttemp2,112)tsrmaxcp,pitchmaxcp,cpmax
       write(iouttemp4,112)tsrmaxct,pitchmaxct,ctmax
111    format(a10,f10.4,a20,2f8.3)
112    format(3f10.4)
       close(Unit=iouttemp2)
       close(Unit=iouttemp4)
c
c---------------------------------------------------------------------
c
c       Design mode 
c
c----------------------------------------------------------------------
c
        if(DesMode)then
        nout5     = './cP-KSS.dat'
        iout5 = 25
        OPEN(UNIT=iout5,FILE=Nout5,Form='formatted',status='unknown')
c
	vwind = vwinda
	om    = rpm*pi/30.
        vtip  = om*rtip
        tsr   = vtip/vwind
c       write(*,*)'KSS des ',TSR
c
c       initialize chord and twist from Betz/Schmitz
c
        call Design
c
c---------------------------------------------------------------------
c
c       smooth chord and/or twist using Bezier polynoms 
c
c----------------------------------------------------------------------
c
        call smoothCHORD
c
        call smoothTWIST
c
c       initial BEM
c     
        write(*,*)'Initial BEM'
	call bem(cP,cT,errb)
        ibem = ibem + 1
        write(iout5,201)ibem,'. iter cP/cT = ',cP,cT
c
c-----------------------------------------------------------------------
c      perform noimptw times twist improvement loops towards a = 1/3
c-----------------------------------------------------------------------
c
	   if(ImpTwist)then
             write(*,*)'IMPROVE twist'
             do I =1,noimptw
	       write(*,*)'twist loop no: ',i
               call improveTWIST
	       call bem(cP,cT,errb)
               ibem = ibem + 1
               write(iout5,201)ibem,'. iter cP/cT = ',cP,cT
             end do
           end if
c
	else
        end if 
c
c-----------------------------------------------------------------------
c      improve THICKNESS: read improved thickness distribution
c-----------------------------------------------------------------------
c
	   if(impThick)then
              write(*,*)'IMPROVE thickness'
              filename   ='./ThickDis.new'
               call ThickPrep(filename)
	       call bem(cP,cT,errb)
               ibem = ibem + 1
               write(iout5,201)ibem,'. iter cP/cT = ',cP,cT
 	   end if
c
c-----------------------------------------------------------------------
c      perform additional twist improvement loops because of new thickness
c-----------------------------------------------------------------------
c
c------------------------------------------------------------------------
c       finished
c------------------------------------------------------------------------
	call date_and_time(date,timee,zone,values)
c
	write(*,*)
	write(*,*)'================================================ '
	write(*,'(a10,a4,a1,a2,2(a1,a2))')'date ',date(1:4),' ',
     +           date(5:6),' ',date(7:8)
        write(*,'(2a10)')'Start: ', timea
        write(*,'(2a10)')'End:   ', timee
        Write(*,*)'KSS END                                          '
	write(*,*)'================================================ '
	write(*,*)
c
100     format(7f8.3,2f8.1,4f8.3,2f8.2,i8,x,e8.1,x,a4)
101     format(16a8,x,a8,x,a4)
102     format(a10,f8.3,a12)
103     format(a20,f12.4,a10)
105     format(a10,i5)
106     format(a10,a1)
107     format(2f12.3,f12.6)
108     format(a20,2f12.3,i4)
109     format(3a12)
110     format(a120,f12.6,a30,2f10.3)
c
200     format(a18,2f10.4)
201     format(i2,a16,2f10.4)
c
999     format(a6,a8,a6)
998     format(34f12.4)
997     format(f6.3,f8.3,f6.1)
996     format(4a12,a6)    
c
      end
c
