      subroutine qsrbsj(dk,nk,geospr)
      implicit none
      integer nk
      double precision dk,geospr
c
      include 'qsrglobal.h'
c
      integer i,ik
      double precision k,x
      double precision bessj0,bessj1,bessj
c
      do ik=1,nk
        k=dble(ik)*dk
        x=k*stdis
        bsj(ik,0)=bessj0(x)
        bsj(ik,1)=bessj1(x)
        if(x.gt.2.d0)then
          bsj(ik,2)=bsj(ik,1)*2.d0/x-bsj(ik,0)
        else
          bsj(ik,2)=bessj(2,x)
        endif
        if(x.gt.3.d0)then
          bsj(ik,3)=bsj(ik,2)*4.d0/x-bsj(ik,1)
        else
          bsj(ik,3)=bessj(3,x)
        endif
        bsj(ik,-1)=-bsj(ik,1)
        do i=-1,3
          bsj(ik,i)=bsj(ik,i)*geospr
        enddo
      enddo
c
      return
      end