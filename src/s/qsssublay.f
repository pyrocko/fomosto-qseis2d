      subroutine qsssublay(resolut,fcut)
      implicit none
c
      double precision resolut(3),fcut
c
      include 'qssglobal.h'
c
c     work space
c
      integer i,i0,l,ivp,ivs,iro
      double precision dh,dro,dvp,dvs,dqp,dqs,z,dz,halfwvlen
c
      n0=0
      do l=1,l0-1
        dz=z2(l)-z1(l)
        dvp=2.d0*dabs(vp2(l)-vp1(l))/(vp2(l)+vp1(l))
        if(vs1(l).gt.vspmin*vp1(l))then
          dvs=2.d0*dabs(vs2(l)-vs1(l))/(vs2(l)+vs1(l))
        else
          dvs=0.d0
        endif
        dro=2.d0*dabs(ro2(l)-ro1(l))/(ro2(l)+ro1(l))
c
        if(dvp.le.0.d0)then
          ivp=1
        else
          halfwvlen=0.25d0*(vp2(l)+vp1(l))/fcut
          ivp=1+idint(dz/halfwvlen)
          if(resolut(1).gt.0.d0)ivp=min0(ivp,1+idint(dvp/resolut(1)))
        endif
        if(dvs.le.0.d0)then
          ivs=1
        else
          halfwvlen=0.25d0*(vs2(l)+vs1(l))/fcut
          ivs=1+idint(dz/halfwvlen)
          if(resolut(2).gt.0.d0)ivs=min0(ivs,1+idint(dvs/resolut(2)))
        endif
        iro=1
        if(dro.gt.0.d0.and.resolut(3).gt.0.d0)then
          iro=1+idint(dro/resolut(3))
        endif
c
        i0=max0(ivp,ivs,iro)
        dro=(ro2(l)-ro1(l))/dz
        dvp=(vp2(l)-vp1(l))/dz
        dvs=(vs2(l)-vs1(l))/dz
        dqp=(qp2(l)-qp1(l))/dz
        dqs=(qs2(l)-qs1(l))/dz
        dh=dz/dble(i0)
        do i=1,i0
          n0=n0+1
          if(n0.ge.lmax)then
            stop ' Max. number of layers (lmax) too small defined!'
          endif
          h(n0)=dh
          z=(dble(i)-0.5d0)*dh
          ro(n0)=ro1(l)+dro*z
          vp(n0)=vp1(l)+dvp*z
          vs(n0)=vs1(l)+dvs*z
          qp(n0)=qp1(l)+dqp*z
          qs(n0)=qs1(l)+dqs*z
        enddo
      enddo
c
c     last layer is halfspace
c
      n0=n0+1
      h(n0)=0.d0
      ro(n0)=ro1(l0)
      vp(n0)=vp1(l0)
      vs(n0)=vs1(l0)
      qp(n0)=qp1(l0)
      qs(n0)=qs1(l0)
c
      return
      end
