      subroutine qsrpsv(yk,ak,f,k)
      implicit none
c
c     calculation of response to p-sv source
c     ak(3,4),yk(6,4): solution vector (complex)
c     k: wave number
c
      double precision f,k
      double complex ak(3,4),yk(6,4)
c
      include 'qsrglobal.h'
c
c     work space
c
      integer i,istp,j,l,n,key
      double complex cfac,ck,ch0,pwave,swave
      double complex y0(4,2),c0(4,2),c1(4,2),b(2,4)
      double complex y1(4,2),yup(4,2),orth(2,2)
      double complex coef(2,2)
c
      double complex c2
	data c2/(2.d0,0.d0)/
c
      ck=dcmplx(k,0.d0)
c
c===============================================================================
c
c     matrix propagation from surface to source
c
c     determination of starting upper sublayer
c
      do j=1,2
        do i=1,4
          yup(i,j)=(0.d0,0.d0)
        enddo
      enddo
      yup(1,1)=(1.d0,0.d0)
      yup(3,2)=(1.d0,0.d0)
      if(lzr.eq.1)call cmemcpy(yup,y0,8)
c
      do l=2,lp
        ch0=dcmplx(hp(l-1),0.d0)
        n=nno(l-1)
c
c       determination of propagation matrix
c
        call qsrve2am(n,ck,0.d0,yup,c0,2)
        pwave=cdexp(-kp(n)*ch0)
        swave=cdexp(-ks(n)*ch0)
c
c       orthonormalization of the p-sv modes
c
        cfac=(1.d0,0.d0)/(c0(3,2)*c0(1,1)-c0(1,2)*c0(3,1))
        orth(1,1)=c0(3,2)*cfac
        orth(1,2)=-c0(1,2)*cfac
        orth(2,1)=-c0(3,1)*cfac
        orth(2,2)=c0(1,1)*cfac
        call caxcb(c0,orth,4,2,2,c1)
        if(l.gt.lzr)then
c
c         additional normalization to avoid overflow
c
          do i=1,2
            orth(i,1)=orth(i,1)*pwave
            orth(i,2)=orth(i,2)*swave
          enddo
          call caxcb(y0,orth,4,2,2,y1)
          call cmemcpy(y1,y0,8)
        endif
c
c        c1(1,1)=c1(1,1)
        c1(2,1)=c1(2,1)*pwave*pwave
        c1(3,1)=(0.d0,0.d0)
        c1(4,1)=c1(4,1)*pwave*swave
c
        c1(1,2)=(0.d0,0.d0)
        c1(2,2)=c1(2,2)*swave*pwave
c        c1(3,2)=c1(3,2)
        c1(4,2)=c1(4,2)*swave*swave
c
        call qsram2ve(n,ck,hp(l-1),yup,c1,2)
        if(l.eq.lzr)call cmemcpy(yup,y0,8)
      enddo
c
      n=nno(lp)
      call qsrve2am(n,ck,0.d0,yup,c0,2)
      do istp=1,4
        b(1,istp)=ak(1,istp)
        b(2,istp)=ak(2,istp)
      enddo
      coef(1,1)=c0(1,1)
      coef(1,2)=c0(1,2)
      coef(2,1)=c0(3,1)
      coef(2,2)=c0(3,2)
      key=0
      call cdgemp(coef,b,2,4,0.d0,key)
      if(key.eq.0)then
        print *,'warning in qsrpsv: anormal exit from cdgemp!'
        return
      endif
      do istp=1,4
        do i=1,4
          yk(i,istp)=(0.d0,0.d0)
          do j=1,2
            yk(i,istp)=yk(i,istp)+b(j,istp)*y0(i,j)
          enddo
        enddo
      enddo
c
      return
      end
