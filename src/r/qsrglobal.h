c     GLOBAL INDEX PARAMETERS FOR DEFINING ARRAYS
c     ===========================================
c     nzmax: max. interface index;
c     lmax: max. number of total homogeneous layers (lmax <= nzmax-2);
c     nfmax: max. number of frequency samples
c
      integer nzmax,lmax,nfmax,ndtransmax
      parameter(lmax=2500)
      parameter(nzmax=lmax+2)
      parameter(nfmax=8192)
      parameter(ndtransmax=4)
c
c     INDEX PARAMETERS FOR BESSEL FUNCTION TABLES
c     ===========================================
c
      integer nbsjmax
      parameter(nbsjmax=500000)
c       
c     layered model parameter:
c     n0: number of homogeneous layers
c
      integer n0
      double precision h(lmax),ro(lmax),vp(lmax),vs(lmax)
      double precision qp(lmax),qs(lmax)
	  common /imodel/ n0
      common /dmodel/ h,ro,vp,vs,qp,qs
c
      integer lp,lzr,nno(nzmax)
      double precision hp(nzmax)
	  common /isublayer/ lp,lzr,nno
      common /dsublayer/ hp
c
c     original model parameters
c
      integer l0
      double precision z1(lmax),z2(lmax),ro1(lmax),ro2(lmax)
      double precision vp1(lmax),vp2(lmax),vs1(lmax),vs2(lmax)
      double precision qp1(lmax),qp2(lmax),qs1(lmax),qs2(lmax)
      common /model0/z1,z2,ro1,ro2,vp1,vp2,vs1,vs2,qp1,qp2,qs1,qs2,l0
c
      double complex acc(lmax),kp(lmax),ks(lmax),cla(lmax),cmu(lmax)
      double complex cvp(lmax),cvs(lmax),wa(lmax),wb(lmax)
      common /cpara/ acc,kp,ks,cla,cmu,cvp,cvs,wa,wb
c
      integer nt,nf,ndtrans
      double precision dt,df,fi,twindow
	  double complex comega
	  common /isampling/ nt,nf,ndtrans
      common /dsampling/ dt,df,fi,twindow,comega
c
      double precision stdis,stazi,stbazi
      common /eqstparas/ stdis,stazi,stbazi
c
c     table of J_n(x), n = -1, 0, 1, 2, 3
c
      double precision bsj(nbsjmax,-1:3)
      common /bessels/ bsj
c
      double precision seist(2*nfmax,3)
      double complex seisf(nfmax,3)
	  common /outputdata/ seist,seisf
c
c     constants
c
      double precision PI,PI2
      parameter(PI=3.14159265358979d0,PI2=6.28318530717959d0)
      double precision DEG2RAD,KM2M
      parameter(DEG2RAD=1.745329251994328d-02,KM2M=1.0d+03)
      double precision REARTH
      parameter(REARTH=6.371d+06)
