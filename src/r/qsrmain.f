      subroutine qsrmain(lats,lons,mxx,myy,mzz,mxy,myz,mzx,greenfk,
     &                   iwavelet,tau,nwvl,wvl,
     &                   latr,lonr,depr,tstart,iout,outfile,
     &                   ibndp,flw,fup,
     &                   nl0,z0,vp0,vs0,ro0,qp0,qs0)
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
      character*80 greenfk
c
c     iwavelet = selection of wavelet forms for source time function
c     nwvl = number of time samples of user-defined source wavelet (<= 1024)
c     tau = wavelet duration
c     wvl = time series of user-defined source wavelet
c
      integer iwavelet,nwvl
      double precision tau
      double precision wvl(nwvl)
c
c     RECEIVER PARAMETERS
c     ===================
c     iout = selection of output format
c     ibndp = order of bandpass filter
c     (latr, lonr, depr) = 3D location coordinates of station
c     flw, fup = lower and upper cutoff frequencies
c     outfile = output file
c
      integer iout,ibndp
      double precision latr,lonr,depr,tstart,flw,fup
      character*80 outfile
c
c     nl0 = number of data lines representing the model
c     isurf = selection of surface condition
c     z0 = depth
c     vp0 = P wave velocity, vs0 = S wave velocity, ro0 = density
c     qp0 = Q of P wave, qs0 = Q of S wave
c
      integer nl0
      double precision z0(nl0),vp0(nl0),vs0(nl0),ro0(nl0),
     &                 qp0(nl0),qs0(nl0)
c===================================================================
      include 'qsrglobal.h'
c
c     work space
c
      integer runtime,ierr
      double precision dise,disn
      character*180 dataline
c
      integer time
c
      runtime=time()
c
      print *,'#######################################################'
      print *,'#                                                     #'
      print *,'#               Welcome to the program                #'
      print *,'#                                                     #'
      print *,'#                                                     #'
      print *,'#      QQQ     SSSS   EEEEE   III    SSSS   RRRR      #'
      print *,'#     Q   Q   S       E        I    S       R    R    #'
      print *,'#     Q Q Q    SSS    EEEE     I     SSS    RRRR      #'
      print *,'#     Q  QQ       S   E        I        S   R  R      #'
      print *,'#      QQQQ   SSSS    EEEEE   III   SSSS    R   R     #'
      print *,'#                                                     #'
      print *,'#                  (Version 2014)                     #'
      print *,'#                                                     #'
      print *,'#                                                     #'
      print *,'#                      by                             #'
      print *,'#                 Rongjiang Wang                      #'
      print *,'#              (wang@gfz-potsdam.de)                  #'
      print *,'#                                                     #'
      print *,'#           GeoForschungsZentrum Potsdam              #'
      print *,'#           Last modified: December 2014              #'
      print *,'#######################################################'
c
      call disazi(REARTH,lats,lons,latr,lonr,disn,dise)
      stdis=dsqrt(disn**2+dise**2)
      stazi=datan2(dise,disn)
      call disazi(REARTH,latr,lonr,lats,lons,disn,dise)
      stbazi=datan2(-dise,-disn)
      call qsrsublay(nl0,z0,vp0,vs0,ro0,qp0,qs0,depr)
      call qsrwvint(mxx,myy,mzz,mxy,myz,mzx,greenfk,ibndp,flw,fup)
      call qsrfftinv(iwavelet,tau,nwvl,wvl,tstart,iout,outfile)
c
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #      End of computations with QseisR      #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
      stop
      end
