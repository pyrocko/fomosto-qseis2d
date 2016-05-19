      program main
      implicit none
c
      include 'qssglobal.h'
c
c     work space
c
      integer i,runtime
      double precision srate
      integer time
c
c     read input file file
c
      print *,'#######################################################'
      print *,'#                                                     #'
      print *,'#               Welcome to the program                #'
      print *,'#                                                     #'
      print *,'#                                                     #'
      print *,'#       QQQ     SSSS   EEEEE   III    SSSS    SSSS    #'
      print *,'#      Q   Q   S       E        I    S       S        #'
      print *,'#      Q Q Q    SSS    EEEE     I     SSS     SSS     #'
      print *,'#      Q  QQ       S   E        I        S       S    #'
      print *,'#       QQQQ   SSSS    EEEEE   III   SSSS    SSSS     #'
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
      print *,'                          '
      write(*,'(a,$)')' the input data file is '
      read(*,'(a)')inputfile
      runtime=time()
c
      open(10,file=inputfile,status='old')
      call qssgetinp(10,srate)
      close(10)
c
      call qsswvint(srate)
c
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #      End of computations with QseisS      #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
      stop
      end
