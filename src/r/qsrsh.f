      subroutine qsrsh(yk,ak,f,k)
      implicit none
c
c     calculation of response to sh source
c     yc(3,4): solution vector (complex)
c     k: wave number
c
      double precision f,k
      double complex ak(3,4),yk(6,4)
c
      include 'qsrglobal.h'
c
c     work space
c
      integer i,istp,l,n,key
      double complex cnorm
      double complex y0(2)
      double complex y1(2),yup(2)
      double complex hk(2,2,nzmax),c0
c
c===============================================================================
c
c     matrix propagation from surface to source
c
      do l=1,lp-1
        n=nno(l)
        call qsrhksh(hk(1,1,l),hp(l),n)
      enddo
c
c     yup: the starting solution vector
c
      yup(1)=(1.d0,0.d0)
      yup(2)=(0.d0,0.d0)
      if(lzr.eq.1)call cmemcpy(yup,y0,2)
c
      do l=2,lp
        n=nno(l-1)
c
c       determination of propagation matrix
c
        call caxcb(hk(1,1,l-1),yup,2,2,1,y1)
        call cmemcpy(y1,yup,2)
        if(l.gt.lzr)then
c
c         normalization
c
          cnorm=cdexp(-ks(n)*dcmplx(hp(l-1),0.d0))
          y0(1)=y0(1)*cnorm
          y0(2)=y0(2)*cnorm
        else if(l.eq.lzr)then
          call cmemcpy(yup,y0,2)
        endif
      enddo
      n=nno(lp)
      c0=(0.5d0,0.d0)*(yup(1)+yup(2)/(cmu(n)*ks(n)))
      do istp=1,4
        yk(5,istp)=y0(1)*ak(3,istp)/c0
        yk(6,istp)=y0(2)*ak(3,istp)/c0
      enddo
      return
      end
