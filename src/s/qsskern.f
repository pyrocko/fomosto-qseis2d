      subroutine qsskern(yc,f,k)
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     yc(3,4): solution vector (complex)
c     k: wave number
c     f: frequency. 2*pi*f = circlar frequency
c     fi: imaginary part of the frequency
c
      double precision f,k
      double complex yc(3,4)
c
      include 'qssglobal.h'
c
      integer i,istp,lup,llw
c
      call qsswaveno(f,k)
c
      do istp=1,4
        do i=1,3
          yc(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
c     determination of starting upper sublayer
c
      lup=1
      llw=lp
      call qsspsv(yc,k,lup,llw)
      call qsssh(yc,k,lup,llw)
c
      return
      end