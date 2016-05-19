      subroutine qsspsv(yc,k,lup,llw)
      implicit none
c
c     calculation of response to p-sv source
c     yc(4,4): solution vector (complex)
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
      integer i,istp,j,l,n,key
      double complex cfac,ck,ch0,pwave,swave
      double complex y0(4,2),c0(4,2),c1(4,2),b(4,4),b0(4,4),y(4,4)
      double complex y1(4,2),yup(4,2),ylw(4,2),orth(2,2)
      double complex coef(4,4),cnorm(2)
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
c
      if(lup.eq.1)then
        do j=1,2
          do i=1,4
            yup(i,j)=(0.d0,0.d0)
          enddo
        enddo
        yup(1,1)=(1.d0,0.d0)
        if(vpair.gt.0.d0)yup(2,1)=-accair/kpair
        yup(3,2)=(1.d0,0.d0)
      else
        n=nno(lup)
c
        yup(1,1)=kp(n)
        yup(2,1)=wa(n)
        yup(3,1)=ck
        yup(4,1)=wb(n)*kp(n)
c
        yup(1,2)=ck
        yup(2,2)=wb(n)*ks(n)
        yup(3,2)=ks(n)
        yup(4,2)=wa(n)
      endif
      if(lup.eq.lzr)call cmemcpy(yup,y0,8)
c
      do l=lup+1,ls
        ch0=dcmplx(hp(l-1),0.d0)
        n=nno(l-1)
c
c       determination of propagation matrix
c
        call qssve2am(n,ck,0.d0,yup,c0,2)
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
c
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
c       for partial solution only!
c
        if(.not.pup(l-1))then
          c1(2,1)=(0.d0,0.d0)
          c1(4,1)=(0.d0,0.d0)
          if(l.gt.lzr)then
            do j=1,4
              y0(j,1)=(0.d0,0.d0)
            enddo
          endif
        endif
        if(.not.svup(l-1))then
          c1(2,2)=(0.d0,0.d0)
          c1(4,2)=(0.d0,0.d0)
          if(l.gt.lzr)then
            do j=1,4
              y0(j,2)=(0.d0,0.d0)
            enddo
          endif
        endif
        if(.not.pdw(l-1))then
          c1(2,1)=(0.d0,0.d0)
          c1(2,2)=(0.d0,0.d0)
          if(l-1.eq.lzr)then
            call qssve2am(n,ck,0.d0,y0,c0,2)
            c0(2,1)=(0.d0,0.d0)
            c0(2,2)=(0.d0,0.d0)
            call qssam2ve(n,ck,0.d0,y0,c0,2)
          endif
        endif
        if(.not.svdw(l-1))then
          c1(4,1)=(0.d0,0.d0)
          c1(4,2)=(0.d0,0.d0)
          if(l-1.eq.lzr)then
            call qssve2am(n,ck,0.d0,y0,c0,2)
            c0(4,1)=(0.d0,0.d0)
            c0(4,2)=(0.d0,0.d0)
            call qssam2ve(n,ck,0.d0,y0,c0,2)
          endif
        endif
c
c       end partial solution procedure!
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
        call qssam2ve(n,ck,hp(l-1),yup,c1,2)
c
        if(l.eq.lzr)call cmemcpy(yup,y0,8)
      enddo
c
c===============================================================================
c
c     matrix propagation from half-space to source
c
      n=nno(llw)
      ylw(1,1)=-kp(n)
      ylw(2,1)=wa(n)
      ylw(3,1)=ck
      ylw(4,1)=-wb(n)*kp(n)
c
      ylw(1,2)=ck
      ylw(2,2)=-wb(n)*ks(n)
      ylw(3,2)=-ks(n)
      ylw(4,2)=wa(n)
c
      if(llw.gt.ls.and.llw.eq.lzr)call cmemcpy(ylw,y0,8)
c
      do l=llw-1,ls,-1
        ch0=dcmplx(hp(l),0.d0)
        n=nno(l)
c
c       determination of propagation matrix
c
        call qssve2am(n,ck,0.d0,ylw,c0,2)
        pwave=cdexp(-kp(n)*ch0)
        swave=cdexp(-ks(n)*ch0)
c
c       orthonormalization of the p-sv modes
c
        cfac=(1.d0,0.d0)/(c0(4,2)*c0(2,1)-c0(2,2)*c0(4,1))
        orth(1,1)=c0(4,2)*cfac
        orth(1,2)=-c0(2,2)*cfac
        orth(2,1)=-c0(4,1)*cfac
        orth(2,2)=c0(2,1)*cfac
        call caxcb(c0,orth,4,2,2,c1)
c
        if(l.lt.lzr)then
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
c       for partial solution only!
c
        if(.not.pup(l))then
          c1(1,1)=(0.d0,0.d0)
          c1(1,2)=(0.d0,0.d0)
          if(l+1.eq.lzr)then
            call qssve2am(n,ck,0.d0,y0,c0,2)
            c0(1,1)=(0.d0,0.d0)
            c0(1,2)=(0.d0,0.d0)
            call qssam2ve(n,ck,0.d0,y0,c0,2)
          endif
        endif
        if(.not.svup(l))then
          c1(3,1)=(0.d0,0.d0)
          c1(3,2)=(0.d0,0.d0)
          if(l+1.eq.lzr)then
            call qssve2am(n,ck,0.d0,y0,c0,2)
            c0(3,1)=(0.d0,0.d0)
            c0(3,2)=(0.d0,0.d0)
            call qssam2ve(n,ck,0.d0,y0,c0,2)
          endif
        endif
        if(.not.pdw(l))then
          c1(1,1)=(0.d0,0.d0)
          c1(3,1)=(0.d0,0.d0)
          if(l.lt.lzr)then
            do j=1,4
              y0(j,1)=(0.d0,0.d0)
            enddo
          endif
        endif
        if(.not.svdw(l))then
          c1(1,2)=(0.d0,0.d0)
          c1(3,2)=(0.d0,0.d0)
          if(l.lt.lzr)then
            do j=1,4
              y0(j,2)=(0.d0,0.d0)
            enddo
          endif
        endif
c
c       end partial solution procedure!
c
        c1(1,1)=c1(1,1)*pwave*pwave
c        c1(2,1)=c1(2,1)
        c1(3,1)=c1(3,1)*pwave*swave
        c1(4,1)=(0.d0,0.d0)
c
        c1(1,2)=c1(1,2)*swave*pwave
        c1(2,2)=(0.d0,0.d0)
        c1(3,2)=c1(3,2)*swave*swave
c        c1(4,2)=c1(4,2)
c
        call qssam2ve(n,ck,-hp(l),ylw,c1,2)
        if(l.gt.ls.and.l.eq.lzr)call cmemcpy(ylw,y0,8)
      enddo
c
c===============================================================================
c     source function
c===============================================================================
c
      do istp=1,4
        do i=1,4
          b(i,istp)=dcmplx(sfct0(i,istp)+k*sfct1(i,istp),0.d0)
        enddo
      enddo
      do i=1,4
        do j=1,2
          coef(i,j)=yup(i,j)
          coef(i,j+2)=-ylw(i,j)
        enddo
      enddo
      if(ipath.eq.1)call cmemcpy(b,b0,16)
      key=0
      call cdgemp(coef,b,4,4,0.d0,key)
      if(key.eq.0)then
        print *,'warning in qsspsv: anormal exit from cdgemp!'
        return
      endif
      if(lzr.le.ls)then
        do istp=1,4
          do i=1,4
            y(i,istp)=(0.d0,0.d0)
            do j=1,2
              y(i,istp)=y(i,istp)+b(j,istp)*y0(i,j)
            enddo
          enddo
        enddo
      else
        do istp=1,4
          do i=1,4
            y(i,istp)=(0.d0,0.d0)
            do j=1,2
              y(i,istp)=y(i,istp)+b(j+2,istp)*y0(i,j)
            enddo
          enddo
        enddo
      endif
c
      if(ipath.eq.1)then
        n=nno(lpath)
        ylw(1,1)=-kp(n)
        ylw(2,1)=wa(n)
        ylw(3,1)=ck
        ylw(4,1)=-wb(n)*kp(n)
c
        ylw(1,2)=ck
        ylw(2,2)=-wb(n)*ks(n)
        ylw(3,2)=-ks(n)
        ylw(4,2)=wa(n)
c
        if(lpath.gt.ls.and.lpath.eq.lzr)call cmemcpy(ylw,y0,8)
        do l=lpath-1,ls,-1
          ch0=dcmplx(hp(l),0.d0)
          n=nno(l)
c
c         determination of propagation matrix
c
          call qssve2am(n,ck,0.d0,ylw,c0,2)
          pwave=cdexp(-kp(n)*ch0)
          swave=cdexp(-ks(n)*ch0)
c
c         orthonormalization of the p-sv modes
c
          cfac=(1.d0,0.d0)/(c0(4,2)*c0(2,1)-c0(2,2)*c0(4,1))
          orth(1,1)=c0(4,2)*cfac
          orth(1,2)=-c0(2,2)*cfac
          orth(2,1)=-c0(4,1)*cfac
          orth(2,2)=c0(2,1)*cfac
c
          call caxcb(c0,orth,4,2,2,c1)
          if(l.lt.lzr)then
c
c           additional normalization to avoid overflow
c
            do i=1,2
              orth(i,1)=orth(i,1)*pwave
              orth(i,2)=orth(i,2)*swave
            enddo
            call caxcb(y0,orth,4,2,2,y1)
            call cmemcpy(y1,y0,8)
          endif
c
c         for partial solution only!
c
          if(.not.pup(l))then
            c1(1,1)=(0.d0,0.d0)
            c1(1,2)=(0.d0,0.d0)
            if(l+1.eq.lzr)then
              call qssve2am(n,ck,0.d0,y0,c0,2)
              c0(1,1)=(0.d0,0.d0)
              c0(1,2)=(0.d0,0.d0)
              call qssam2ve(n,ck,0.d0,y0,c0,2)
            endif
          endif
          if(.not.svup(l))then
            c1(3,1)=(0.d0,0.d0)
            c1(3,2)=(0.d0,0.d0)
            if(l+1.eq.lzr)then
              call qssve2am(n,ck,0.d0,y0,c0,2)
              c0(3,1)=(0.d0,0.d0)
              c0(3,2)=(0.d0,0.d0)
              call qssam2ve(n,ck,0.d0,y0,c0,2)
            endif
          endif
          if(.not.pdw(l))then
            c1(1,1)=(0.d0,0.d0)
            c1(3,1)=(0.d0,0.d0)
            if(l.lt.lzr)then
              do j=1,4
                  y0(j,1)=(0.d0,0.d0)
                enddo
              endif
            endif
          if(.not.svdw(l))then
            c1(1,2)=(0.d0,0.d0)
            c1(3,2)=(0.d0,0.d0)
            if(l.lt.lzr)then
              do j=1,4
                y0(j,2)=(0.d0,0.d0)
             enddo
            endif
          endif
c
c         end partial solution procedure!
c
          c1(1,1)=c1(1,1)*pwave*pwave
c          c1(2,1)=c1(2,1)
          c1(3,1)=c1(3,1)*pwave*swave
          c1(4,1)=(0.d0,0.d0)
c
          c1(1,2)=c1(1,2)*swave*pwave
          c1(2,2)=(0.d0,0.d0)
          c1(3,2)=c1(3,2)*swave*swave
c          c1(4,2)=c1(4,2)
c
          call qssam2ve(n,ck,-hp(l),ylw,c1,2)
          if(l.gt.ls.and.l.eq.lzr)call cmemcpy(ylw,y0,8)
        enddo
        do i=1,4
          do j=1,2
            coef(i,j)=yup(i,j)
            coef(i,j+2)=-ylw(i,j)
          enddo
        enddo
        key=0
        call cdgemp(coef,b0,4,4,0.d0,key)
        if(key.eq.0)then
          print *,'warning in qsspsv: anormal exit from cdgemp!'
          return
        endif
        if(lzr.le.ls)then
          do istp=1,4
            do i=1,4
              do j=1,2
                y(i,istp)=y(i,istp)-b0(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,4
            do i=1,4
              do j=1,2
                y(i,istp)=y(i,istp)-b0(j+2,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
      endif
c
      n=nno(lzr)
      call qssve2am(n,ck,0.d0,y,coef,4)
      do istp=1,4
        yc(1,istp)=coef(1,istp)
        yc(2,istp)=coef(3,istp)
      enddo
c
      return
      end
