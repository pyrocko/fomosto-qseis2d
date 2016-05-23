c     GLOBAL INDEX PARAMETERS FOR DEFINING ARRAYS
c     ===========================================
c     nzmax: max. interface index;
c     lmax: max. number of total homogeneous layers (lmax <= nzmax-2);
c     nfmax: max. number of frequency samples
c     nkmax: max. number of slowness samples
c
      integer nzmax,lmax,nfmax,nkmax
      parameter(lmax=2500)
      parameter(nzmax=lmax+2)
      parameter(nfmax=2048)
	  parameter(nkmax=100000)
c
c     EARTH RADIUS IN METER
c     =====================
c
      double precision rr0,km2m
      parameter(rr0=6.371d+06,km2m=1.0d+03)
c
c     ATMOSPHERIC PARAMETERS
c     ======================
c     double precision roair,vpair
c     parameter(roair=0.1300d+01,vpair=0.3318d+03)
      double precision roair,vpair
      parameter(roair=0.d0,vpair=0.d0)
c
c     THE MININUM VS/VP RATIO: VSPMIN
c     ===============================
c     (if vs/vp < vspmin, then fluid medium is assumed)
c
      double precision vspmin
      parameter(vspmin=0.05d0)
c
c     FOR FLAT-EARTH-TRANSFORMATION
c     =============================
c
      integer ndens
      parameter(ndens=1)
c
      double complex accair,cvpair,kpair,comega
      common /airpara/ accair,cvpair,kpair,comega
c
c     zr: receiver depth
c     lzr: sublayer no of receiver
c
      integer lzr
      double precision zr,dismax
	  common /ireceiver/ lzr
      common /dreceiver/ zr,dismax
c
      integer lp,nno(nzmax)
      double precision hp(nzmax)
	  common /isublayer/ lp,nno
      common /dsublayer/ hp
c
c     original model parameters
c
      integer l0
      double precision z1(lmax),z2(lmax),ro1(lmax),ro2(lmax)
      double precision vp1(lmax),vp2(lmax),vs1(lmax),vs2(lmax)
      double precision qp1(lmax),qp2(lmax),qs1(lmax),qs2(lmax)
	  common /imodel0/l0
      common /dmodel0/z1,z2,ro1,ro2,vp1,vp2,vs1,vs2,qp1,qp2,qs1,qs2
c       
c     layered model parameter:
c     n0: number of homogeneous layers
c
      integer n0
      double precision h(lmax),ro(lmax),vp(lmax),vs(lmax)
      double precision qp(lmax),qs(lmax)
	  common /imodel/n0
      common /dmodel/h,ro,vp,vs,qp,qs
c
      double complex acc(lmax),kp(lmax),ks(lmax),cla(lmax),cmu(lmax)
      double complex cvp(lmax),cvs(lmax),wa(lmax),wb(lmax)
      common /dpara/ acc,kp,ks,cla,cmu,cvp,cvs,wa,wb
c
c     for partial solution
c
      integer ipartial,npar,ipsv(nzmax)
      double precision zup(nzmax),zlow(nzmax)
	  common /ipartial/ ipartial,npar,ipsv
      common /dpartial/ zup,zlow
c
      logical pup(nzmax),pdw(nzmax),svup(nzmax),svdw(nzmax),sh(nzmax)
      common /psvfilter/ pup,pdw,svup,svdw,sh
c
c     slowness cut-offs
c
      double precision slw(4)
      common /slwcutoffs/ slw
c
c     source parameters
c
      integer ls
      double precision zs,zs0
      double precision sfct0(6,4),sfct1(6,4)
	  common /isource/ ls
      common /dsource/ zs,zs0,sfct0,sfct1
c
c     path filtering
c
      integer iflat,ipath,lpath,isurf
      double precision pathdepth
	  common /ipathfilter/iflat,ipath,lpath,isurf
      common /dpathfilter/pathdepth
c
      integer nt,nf
      double precision dt,df,fi
	  common /isampling/nt,nf
      common /dsampling/dt,df,fi
c
c     input and output data files
c
      character*800 inputfile
      common /inputdata/ inputfile
      character*800 greeninfo,greenfile
      common /outdata/ greeninfo,greenfile
