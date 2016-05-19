      subroutine qsrhksh(hk,z,n)
      implicit none
c
      integer n
      double precision z
      double complex hk(2,2)
c
      include 'qsrglobal.h'
c
      double complex cx,cem,cch,csh
c
      cx=ks(n)*dcmplx(2.d0*z,0.d0)
      if(z.gt.0.d0)then
        cem=cdexp(-cx)
        cch=(0.5d0,0.d0)*((1.d0,0.d0)+cem)
        csh=(0.5d0,0.d0)*((1.d0,0.d0)-cem)
      else
        cem=cdexp(cx)
        cch=(0.5d0,0.d0)*((1.d0,0.d0)+cem)
        csh=-(0.5d0,0.d0)*((1.d0,0.d0)-cem)
      endif
c
c     propagator matrix for SH waves
c
      hk(1,1)=cch
      hk(1,2)=csh/(cmu(n)*ks(n))
      hk(2,1)=csh*cmu(n)*ks(n)
      hk(2,2)=cch
c
      return
      end
