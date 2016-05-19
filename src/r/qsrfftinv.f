      subroutine qsrfftinv(iwavelet,tau,nwvl,wvl,tstart,iout,outfile)
      implicit none
c
      integer iwavelet,nwvl,iout
      double precision tau,tstart,flw,fup
      double precision wvl(nwvl)
      character*80 outfile
c
      include 'qsrglobal.h'
c
      integer nn
      parameter (nn=2*nfmax)
      integer lf,mf,i,j,it
      double precision t,slw0,omi
      double precision f(nn)
      double complex s,bcs,bss
      double complex yb(3)
      double precision dy(2*nn,3)
      double complex cy(nn,3)
      equivalence(dy,cy)
c
      double complex wvf(nfmax)
c
      if(iout.ne.1)then
c
c       convert cylindrical to cartesian system
c
        bcs=dcmplx(dcos(stbazi),0.d0)
        bss=dcmplx(dsin(stbazi),0.d0)
        do lf=1,nf
          do i=1,3
            yb(i)=seisf(lf,i)
          enddo
          seisf(lf,1)=yb(2)*bss-yb(3)*bcs           !east
          seisf(lf,2)=yb(2)*bcs+yb(3)*bss           !north
          seisf(lf,3)=-yb(1)                        !up
        enddo
      endif
      do lf=1,nf
        f(lf)=dble(lf-1)*df
        do i=1,3
          cy(lf,i)=seisf(lf,i)
        enddo
      enddo
c
c     for time reduction
c
      do lf=1,nf
        s=cdexp(dcmplx(-fi,f(lf))*dcmplx(PI2*tstart,0.d0))
        do i=1,3
          cy(lf,i)=cy(lf,i)*s
        enddo
      enddo
c
c     muliplication with wavelet spectrum
c
      call qsrwavelet(iwavelet,tau,nwvl,wvl,wvf,nf)
c
      do lf=1,nf
        do i=1,3
          cy(lf,i)=cy(lf,i)*wvf(lf)
        enddo
      enddo
c
      mf=1
      do lf=2*nf,nf+2,-1
        mf=mf+1
        do i=1,3
          cy(lf,i)=dconjg(cy(mf,i))
        enddo
      enddo
      do i=1,3
        cy(nf+1,i)=(0.d0,0.d0)
      enddo
c
c     convention for Fourier transform:
c     f(t)=\int F(f) exp(i2\pi f t) df
c
      do i=1,3
        call four1(dy(1,i),2*nf,+1)
      enddo
c
      open(20,file=outfile,status='unknown')
      if(iout.ne.1)then
        write(20,'(a)')'       Time[sec]       '
     &               //'East[m/s]      North[m/s]         Up[m/s]'
      else
        write(20,'(a)')'       Time[sec]       '
     &               //'Down[m/s]     Radial[m/s]  Azimuthal[m/s]'
      endif
c
      omi=PI2*fi
      do it=1,nt
        t=dble(it-1)*dt
        do i=1,3
          seist(it,i)=df*dreal(cy(it,i))*dexp(-omi*t)
        enddo
        write(20,'(f16.6,3E16.8)')t+tstart,
     &                            seist(it,1),seist(it,2),seist(it,3)
      enddo
      close(20)
      return
      end
