      subroutine qsssh(yc,k,lup,llw)
      implicit none
c
c     calculation of response to sh source
c     yc(3,4): solution vector (complex)
c     k: wave number
c
      integer lup,llw
      double precision k
      double complex yc(3,4)
c
      include 'qssglobal.h'
c
c     work space
c
      integer i,istp,l,n,key
      double complex cnorm
      double complex y0(2),b(2,4),b0(2,4),y(2,4)
      double complex y1(2),yup(2),ylw(2)
      double complex hk(2,2,nzmax),coef(2,2),c0
c
c===============================================================================
c
c     matrix propagation from surface to source
c
      do l=lup,ls-1
        n=nno(l)
        call qsshksh(hk(1,1,l),hp(l),n)
      enddo
      do l=ls,llw-1
        n=nno(l)
        call qsshksh(hk(1,1,l),-hp(l),n)
      enddo
c
c     yup: the starting solution vector
c
      yup(1)=(1.d0,0.d0)
      if(lup.eq.1.and.isurf.eq.0)then
        yup(2)=(0.d0,0.d0)
      else
        n=nno(lup)
        yup(2)=cmu(n)*ks(n)
      endif
      if(lup.eq.lzr)call cmemcpy(yup,y0,2)
c
      do l=lup+1,ls
        n=nno(l-1)
c
c       determination of propagation matrix
c
        call caxcb(hk(1,1,l-1),yup,2,2,1,y1)
        call cmemcpy(y1,yup,2)
c
c       normalization
c
        if(l.gt.lzr)then
          cnorm=cdexp(-ks(n)*dcmplx(hp(l-1),0.d0))
          y0(1)=y0(1)*cnorm
          y0(2)=y0(2)*cnorm
        else if(l.eq.lzr)then
          call cmemcpy(yup,y0,2)
        endif
      enddo
c
c===============================================================================
c
c     matrix propagation from half-space to source
c
c     ylw: the starting solution vector
c
      n=nno(llw)
      ylw(1)=(1.d0,0.d0)
      if(sh(llw))then
        ylw(2)=-cmu(n)*ks(n)
      else
        ylw(2)=(0.d0,0.d0)
      endif
      if(llw.gt.ls.and.llw.eq.lzr)call cmemcpy(ylw,y0,2)
c
      do l=llw-1,ls,-1
        n=nno(l)
c
c       determination of propagation matrix
c
        call caxcb(hk(1,1,l),ylw,2,2,1,y1)
        call cmemcpy(y1,ylw,2)
        if(l.lt.lzr)then
          cnorm=cdexp(-ks(n)*dcmplx(hp(l),0.d0))
          y0(1)=y0(1)*cnorm
          y0(2)=y0(2)*cnorm
        else if(l.gt.ls.and.l.eq.lzr)then
          call cmemcpy(ylw,y0,2)
        endif
      enddo
c
c===============================================================================
c     source function
c===============================================================================
c
      do istp=1,4
        do i=1,2
          b(i,istp)=dcmplx(sfct0(i+4,istp)+k*sfct1(i+4,istp),0.d0)
        enddo
      enddo
      do i=1,2
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
      enddo
      if(ipath.eq.1)call cmemcpy(b,b0,8)
      key=0
      call cdgemp(coef,b,2,4,0.d0,key)
      if(key.eq.0)then
        print *,'warning in qsssh: anormal exit from cdgemp!'
        return
      endif
      if(lzr.le.ls)then
        do istp=1,4
          do i=1,2
            y(i,istp)=b(1,istp)*y0(i)
          enddo
        enddo
      else
        do istp=1,4
          do i=1,2
            y(i,istp)=b(2,istp)*y0(i)
          enddo
        enddo
      endif
      if(ipath.eq.1)then
        n=nno(lpath)
        ylw(1)=(1.d0,0.d0)
        if(sh(lpath))then
          ylw(2)=-cmu(n)*ks(n)
        else
          ylw(2)=(0.d0,0.d0)
        endif
        if(lpath.gt.ls.and.lpath.eq.lzr)call cmemcpy(ylw,y0,2)
        do l=lpath-1,ls,-1
          n=nno(l)
c
c         determination of propagation matrix
c
          call caxcb(hk(1,1,l),ylw,2,2,1,y1)
          call cmemcpy(y1,ylw,2)
c
          if(l.lt.lzr)then
            cnorm=cdexp(-ks(n)*dcmplx(hp(l),0.d0))
            y0(1)=y0(1)*cnorm
            y0(2)=y0(2)*cnorm
          else if(l.gt.ls.and.l.eq.lzr)then
            call cmemcpy(ylw,y0,2)
          endif
        enddo
        do i=1,2
          coef(i,1)=yup(i)
          coef(i,2)=-ylw(i)
        enddo
        key=0
        call cdgemp(coef,b0,2,4,0.d0,key)
        if(key.eq.0)then
          print *,'warning in qsssh: anormal exit from cdgemp!'
          return
        endif
        if(lzr.le.ls)then
          do istp=1,4
            do i=1,2
              y(i,istp)=y(i,istp)-b0(1,istp)*y0(i)
            enddo
          enddo
        else
          do istp=1,4
            do i=1,2
              y(i,istp)=y(i,istp)-b0(2,istp)*y0(i)
            enddo
          enddo
        endif
      endif
c
      n=nno(lzr)
      do istp=1,4
        yc(3,istp)=(0.5d0,0.d0)*(y(1,istp)+y(2,istp)/(cmu(n)*ks(n)))
      enddo
c
      return
      end
