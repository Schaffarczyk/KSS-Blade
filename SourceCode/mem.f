c
c	module mem instead of "evil" common
c
	module mem
c	
      integer status,ioout, maxiter,nsec,ndes
      integer indpro
c
      logical DesMode,ImpChord,ImpTwist,ImpThick
c
      real b, dens
      real rtip, tiplen, ar, rroot 
      real glopitch, glopitcha, glopitche
      integer npitch
      real aoab
      real vtip, vwinda, vwinde, vwind
      real RPM, om, tsr,tsrloc
      integer nwind
      integer DesSchema
 
      real twmax, chmax
      real pi, eps

      CHARACTER*100 nin, nameprout,namedesin(100)
      Character*20  nout,noutd,nread
      character*10  timea,timee
      character*8   date
      character*5   zone
      CHARACTER*4   namepr(1:20)
      CHARACTER*1   tipshape
c     
      REAL, allocatable :: rsecsp(:),aoasp(:)
      REAL, allocatable :: chsp(:)  ,twistsp(:)
      REAL, allocatable :: chspsm(:),twistspsm(:)
      REAL, allocatable :: clsp(:)  ,cdsp(:)
      REAL, allocatable :: chs(:),twists(:),clss(:),cdss(:)
      REAL, allocatable :: abem(:),apbem(:)
c
c     profile data
c     1st index: max number of profiles
c     2nd index: max number of aoas
c
      REAL aoain(10,300),clin(10,300), cdin(10,300)
      real cls  (10,300),cds (10,300)
c
      real prothick(100)
      real optcl(100),optaoa(100), maxcl(100),maxaoa(100)
      real zeroaoa(100),zeroslope(100)
c
c     nopr: number of profiles, nd = lines in BlaDes
c
      integer nopr,nd
      integer values(8)
c
c     np: number of data in profile file, iprno: profile number 
c
      integer np(0:100),iprno(100)
c
c     profile thickness
c
      real dthr(100),dthth(100),dthsp(100)
      real minthick
      integer irth
c
c     elastic pitch due to twist bend
c
      real eltw(0:15,25),rstw(0:15)
      logical twistb
c
	end module mem
