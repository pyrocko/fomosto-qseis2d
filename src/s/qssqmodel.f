      subroutine qssqmodel(f)
      implicit none
c
c     calculate q based on the constant q model
c
c     f = frequency
c
      double precision f
c
      include 'qssglobal.h'
c
      integer n
      double complex cmp,cms
c
      double precision pi2
      data pi2/6.28318530717959d0/
c
      comega=dcmplx(pi2*f,pi2*fi)
      cvpair=dcmplx(vpair,0.d0)
c
      do n=1,n0
        if(qp(n).le.0.d0)then
          cmp=(1.d0,0.d0)
        else
          cmp=dcmplx(1.d0,0.5d0/qp(n))
        endif
        if(qs(n).le.0.d0.or.vs(n).le.vspmin*vp(n))then
          cms=(1.d0,0.d0)
        else
          cms=dcmplx(1.d0,0.5d0/qs(n))
        endif
        cvp(n)=dcmplx(vp(n),0.d0)*cmp
        cvs(n)=dcmplx(vs(n),0.d0)*cms
        cmu(n)=dcmplx(ro(n)*vs(n)**2,0.d0)*cms**2
        cla(n)=dcmplx(ro(n)*vp(n)**2,0.d0)*cmp**2-(2.d0,0.d0)*cmu(n)
        acc(n)=dcmplx(ro(n),0.d0)*comega**2
      enddo
c
      return
      end
