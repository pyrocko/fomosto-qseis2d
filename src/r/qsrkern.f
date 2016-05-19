      subroutine qsrkern(f,k,ak,yk)
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ak(3,4), yk(6,4): solution vector (complex)
c     k: wave number
c     f: frequency. 2*pi*f = circular frequency
c     fi: imaginary part of the frequency
c
      double precision f,k
      double complex ak(3,4),yk(6,4)
c
      include 'qsrglobal.h'
c
      integer i,istp
c
      call qsrwaveno(f,k)
c
      do istp=1,4
        do i=1,6
          yk(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
c     determination of starting upper sublayer
c
      call qsrpsv(yk,ak,f,k)
      call qsrsh(yk,ak,f,k)
c
      return
      end