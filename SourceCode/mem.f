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
c     flag to switch inclusion of the radial velocity component
c     (Madsen correction, Li et al. WES 7, 75-104, 2022, Eqs. 5-7)
c     on (.T.) or off (.F.); read from Machine.in
c
      logical Radial
c
c     V3: RadMode selects the radial-induction model (if Radial=.T.)
c         RadMode = 1 : Madsen correction        (Eqs. 6+7, V2)
c         RadMode = 2 : vortex cylinder model    (Section 3, V3)
c     read from Machine.in, default 1
c
      integer RadMode
c
c     V3: per-section storage for the vortex cylinder sweep
c     (bound circulation, dihedral angle, VC radial induction)
c
      real, allocatable :: gamsec(:), kapsec(:), urvc(:)
c
c     interpolation weight for binary-search routines (rtbis/funcBS);
c     set in BEM (was an undefined local in rtbis before)
c
      real xpint
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

      CHARACTER*100 nin, nameprout,namedesin(200)
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
      Real, allocatable :: oopdefsp(:),oopdefS(:)
c
c     profile data
c     1st index: max number of profiles
c     2nd index: max number of aoas
c
      REAL aoain(20,300),clin(20,300), cdin(20,300)
      real cls  (20,300),cds (20,300)
c
      real prothick(200)
      real optcl(200),optaoa(200), maxcl(200),maxaoa(200)
      real zeroaoa(200),zeroslope(200)
c
c     nopr: number of profiles, nd = lines in BlaDes
c
      integer nopr,nd
      integer values(8)
c
c     np: number of data in profile file, iprno: profile number 
c
      integer np(0:100),iprno(200)
c
c     profile thickness
c
      real dthr(200),dthth(200),dthsp(200)
      real minthick
      integer irth
c
c     elastic pitch due to twist bend
c
      real eltw(0:50,30),rstw(0:50)
      logical twistb
c
	end module mem
