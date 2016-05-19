      subroutine qsrwavelet(iwavelet,tau,nwvl,wvl,wvf,mm2)
      implicit none
      integer iwavelet,nwvl,mm2
      double precision tau
      double precision wvl(nwvl)
      double complex wvf(mm2)
c
      include 'qsrglobal.h'
c
      integer l,n
      double precision f,omi,x,dt0
      double complex alfa,beta,gamma,eta
c
      double precision eps
      data eps/1.0d-04/
c
      if(iwavelet.ne.0)then
c
c       for wavelet: normalized square half-sinus
c
        do l=1,mm2
          f=df*dble(l-1)
          x=f*tau
          if(x.eq.0.d0)then
            wvf(l)=(1.d0,0.d0)
          else if(x.ge.1.d0-eps.and.x.le.1.d0+eps)then
            wvf(l)=dcmplx(-1.d0/x/(1+x),0.d0)
          else
            wvf(l)=dcmplx(0.d0,1.d0/(PI2*x*(1.d0+x)*(1.d0-x)))
     &               *(cdexp(dcmplx(0.d0,-PI2*x))-(1.d0,0.d0))
          endif
        enddo
      else
c
c       user's own wavelet function
c
        dt0=tau/dble(nwvl-1)
c
        wvf(1)=dcmplx(0.5d0*(wvl(1)+wvl(nwvl)),0.d0)
        do n=2,nwvl-1
          wvf(1)=wvf(1)+dcmplx(wvl(n),0.d0)
        enddo
        wvf(1)=wvf(1)*dcmplx(dt0,0.d0)
c
        do l=2,mm2
          wvf(l)=(0.d0,0.d0)
          omi=PI2*df*dble(l-1)
          alfa=cdexp(dcmplx(0.d0,-omi*dt0))
          beta=(alfa-(1.d0,0.d0))*dcmplx(0.d0,1.d0/omi)
          gamma=alfa*dcmplx(0.d0,1.d0/omi)
     &           -beta*dcmplx(0.d0,1.d0/omi/dt0)
          eta=(1.d0,0.d0)
          do n=1,nwvl-1
            wvf(l)=wvf(l)+eta*(dcmplx(wvl(n),0.d0)*(beta-gamma)
     &            +dcmplx(wvl(n+1),0.d0)*gamma)
            eta=eta*alfa
          enddo
        enddo
      endif
      return
      end
