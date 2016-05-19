      subroutine qsssource(ros,vps,vss)
      implicit none
c
      double precision ros,vps,vss
c
      include 'qssglobal.h'
c
      integer i,istp
      double precision pi,pi2
c
      do istp=1,4
        do i=1,6
          sfct0(i,istp)=0.d0
          sfct1(i,istp)=0.d0
        enddo
      enddo
c
      pi=4.d0*datan(1.d0)
      pi2=2.d0*pi
c
c     istp = 1
c     explosion source (m11=m22=m33=1)
c
      sfct0(1,1)=-1.d0/(pi2*ros*vps*vps)
      sfct1(4,1)=-(vss/vps)**2/pi
c
c     istype = 2
c     strike-slip (m12=m21=1)
c
      sfct1(4,2)=1.d0/pi2
      sfct1(6,2)=-sfct1(4,2)
c
c     istype = 3
c     dip-slip (m13=m31=1)
c
      sfct0(3,3)=-1.d0/(pi2*ros*vss*vss)
      sfct0(5,3)=sfct0(3,3)
c
c     istp = 4
c     compensated linear vector dipole (CLVD) (m11=m22=-1/2, M33=1)
c
      sfct0(1,4)=-1.d0/(pi2*ros*vps*vps)
      sfct1(4,4)=(3.d0-4.d0*(vss/vps)**2)/(2.d0*pi2)
c
      return
      end
