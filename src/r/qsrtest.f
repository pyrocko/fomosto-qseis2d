      program qsrtest
      implicit none
c
c     SOURCE PARAMETERS
c     =================
c     (lats, lons) = geographic coordinates of epicenter
c     mxx,myy,mzz,mxy,myz,mzx = moment tensor
c     greenfk = f-k spectra of incident waves
c
      double precision lats,lons
      double precision mxx,myy,mzz,mxy,myz,mzx
      character*800 greenfk
c
c     iwavelet = selection of wavelet forms for source time function
c     nwvl = number of time samples of user-defined source wavelet (<= 1024)
c     tau = wavelet duration
c     wvl = time series of user-defined source wavelet
c
      integer nwvlmax
      parameter(nwvlmax=1024)
      integer iwavelet,nwvl
      double precision tau
      double precision wvl(nwvlmax)
c
c     RECEIVER PARAMETERS
c     ===================
c
c     iout = selection of output format
c     ibndp = order of bandpass filter
c     (latr, lonr, depr) = 3D location coordinates of station
c     flw, fup = lower and upper cutoff frequencies
c     outfile = output file
c
      integer iout,ibndp
      double precision latr,lonr,depr,tstart,flw,fup
      character*800 outfile
c
c     nl = number of data lines representing the model
c     dep = depth
c     vp = P wave velocity, vs = S wave velocity, ro = density
c     qp = Q of P wave, qs = Q of S wave
c
      integer nlmax
      parameter(nlmax=100)
      integer nl0
      double precision z0(nlmax),vp0(nlmax),vs0(nlmax),ro0(nlmax)
      double precision qp0(nlmax),qs0(nlmax)
c
c     work space
c
      integer i,j,runtime
      character*800 inputfile
      character*1800 dataline
c
      integer time
c
      write(*,'(a,$)')' the input data file is '
      read(*,'(a)')inputfile
c
c     read input file
c
      open(10,file=inputfile,status='old')
c
      call getdata(10,dataline)
      read(dataline,*)lats,lons
      call getdata(10,dataline)
      read(dataline,*)mxx,myy,mzz,mxy,myz,mzx
      call getdata(10,dataline)
      read(dataline,*)greenfk
c
      call getdata(10,dataline)
      read(dataline,*)tau,iwavelet
      if(iwavelet.lt.0.or.iwavelet.gt.2)then
        stop ' Error in qsrmain: wrong wavelet selection!'
      else if(iwavelet.eq.0)then
        call getdata(10,dataline)
        read(dataline,*)nwvl
        if(nwvl.gt.nwvlmax)then
          stop ' too long time series of source wavelet!'
        endif
        read(10,*)(wvl(i),i=1,nwvl)
      endif
c
      call getdata(10,dataline)
      read(dataline,*)latr,lonr,depr
      depr=depr*1000.d0
      call getdata(10,dataline)
      read(dataline,*)tstart
      call getdata(10,dataline)
      read(dataline,*)ibndp,flw,fup
      call getdata(10,dataline)
      read(dataline,*)iout
      call getdata(10,dataline)
      read(dataline,*)outfile
c
      call getdata(10,dataline)
      read(dataline,*)nl0
      if(nl0.gt.nlmax)then
        stop ' too large number of datalines for structure model!'
      endif
      do i=1,nl0
        call getdata(10,dataline)
        read(dataline,*)j,z0(i),vp0(i),vs0(i),ro0(i),qp0(i),qs0(i)
c
c       input units:    -,km,  km/s, km/s, g/cm^3,-,-
c
        z0(i)=1000.d0*z0(i)
        vp0(i)=1000.d0*vp0(i)
        vs0(i)=1000.d0*vs0(i)
        ro0(i)=1000.d0*ro0(i)
      enddo
      close(10)
c
      call qsrmain(lats,lons,mxx,myy,mzz,mxy,myz,mzx,greenfk,
     &             iwavelet,tau,nwvl,wvl,
     &             latr,lonr,depr,tstart,iout,outfile,
     &             ibndp,flw,fup,
     &             nl0,z0,vp0,vs0,ro0,qp0,qs0)
c
      stop
      end
