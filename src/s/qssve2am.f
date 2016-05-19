      subroutine qssve2am(n,ck,z,y,c,nj)
      implicit none
c
c     converting displacement-stress vectors to wave amplitudes
c     by Langer block-diagonal decomposition technique for com-
c     putational efficiency
c
      integer n,nj
      double precision z
      double complex ck
      double complex y(4,nj),c(4,nj)
c
      include 'qssglobal.h'
c
      integer i,j,m,key
      double complex ca,cadp,cads
      double complex swap(4),a(4,4),b(4,4)
c
      ca=(0.5d0,0.d0)/acc(n)
      cadp=ca/kp(n)
      cads=ca/ks(n)
      do j=1,nj
        swap(1)=(-wa(n)*y(1,j)+ck*y(4,j))*cadp
        swap(2)=(-y(2,j)+wb(n)*y(3,j))*ca
        swap(3)=(-ck*y(2,j)+wa(n)*y(3,j))*cads
        swap(4)=(wb(n)*y(1,j)-y(4,j))*ca
c
        c(1,j)= swap(1)+swap(2)
        c(2,j)=-swap(1)+swap(2)
        c(3,j)=-swap(3)+swap(4)
        c(4,j)= swap(3)+swap(4)
      enddo
c
      return
      end
