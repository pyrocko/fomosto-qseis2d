      subroutine qsrsublay(nl0,z0,vp0,vs0,ro0,qp0,qs0,depr)
      implicit none
c
c     nl0 = number of data lines representing the model
c     z0 = depth
c     vp0 = P wave velocity, vs0 = S wave velocity, ro0 = density
c     qp0 = Q of P wave, qs0 = Q of S wave
c     depr = receiver depth
c
      integer nl0
      double precision z0(nl0),vp0(nl0),vs0(nl0),ro0(nl0),
     &                 qp0(nl0),qs0(nl0),depr
c
      include 'qsrglobal.h'
c
c     work space
c
      integer i,i0,l,n,li,lp0,ivp,ivs,iro
      double precision dh,dro,dvp,dvs,dqp,dqs,z,dz,zswap
      double precision dep(nzmax),dep0(nzmax)
c
c     determine upper und lower parameter values of each layer
c
      if(depr.gt.z0(nl0))then
        stop 'qsrsublay: too large receiver depth!'
      endif
      l0=1
      z1(l0)=0.d0
      do i=2,nl0
        if(z0(i).gt.z0(i-1))then
          z1(l0)=z0(i-1)
          vp1(l0)=vp0(i-1)
          vs1(l0)=vs0(i-1)
          ro1(l0)=ro0(i-1)
          qp1(l0)=qp0(i-1)
          qs1(l0)=qs0(i-1)
c
          z2(l0)=z0(i)
          vp2(l0)=vp0(i)
          vs2(l0)=vs0(i)
          ro2(l0)=ro0(i)
          qp2(l0)=qp0(i)
          qs2(l0)=qs0(i)
          l0=l0+1
        else
          z1(l0)=z0(i)
          vp1(l0)=vp0(i)
          vs1(l0)=vs0(i)
          ro1(l0)=ro0(i)
          qp1(l0)=qp0(i)
          qs1(l0)=qs0(i)
        endif
      enddo
      z1(l0)=z0(nl0)
      vp1(l0)=vp0(nl0)
      vs1(l0)=vs0(nl0)
      ro1(l0)=ro0(nl0)
      qp1(l0)=qp0(nl0)
      qs1(l0)=qs0(nl0)
c
      n0=0
      do l=1,l0-1
        dz=z2(l)-z1(l)
        dvp=2.d0*dabs(vp2(l)-vp1(l))/(vp2(l)+vp1(l))
        dvs=2.d0*dabs(vs2(l)-vs1(l))/(vs2(l)+vs1(l))
        dro=2.d0*dabs(ro2(l)-ro1(l))/(ro2(l)+ro1(l))
c
        if(dvp.le.0.d0)then
          ivp=1
        else
          ivp=1+idint(dvp/0.01d0)
        endif
        if(dvs.le.0.d0)then
          ivs=1
        else
          ivs=1+idint(dvs/0.01d0)
        endif
        iro=1
        if(dro.gt.0.d0)then
          iro=1+idint(dro/0.05d0)
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
      lp0=1
      dep0(lp0)=0.d0
      do n=1,n0-1
        lp0=lp0+1
        dep0(lp0)=dep0(lp0-1)+h(n)
      enddo
      lp0=lp0+1
      dep0(lp0)=depr
c
c     sort the z0-profile
c
      do l=1,lp0-1
        do li=l+1,lp0
          if(dep0(li).lt.dep0(l))then
            zswap=dep0(l)
            dep0(l)=dep0(li)
            dep0(li)=zswap
          endif
        enddo
      enddo
c
c     delete duplicates
c
      lp=1
      dep(lp)=0.d0
      do l=2,lp0
        if(dep0(l).gt.dep(lp))then
          hp(lp)=dep0(l)-dep(lp)
          lp=lp+1
          dep(lp)=dep0(l)
        endif
      enddo
      hp(lp)=0.d0
c
c     determine lzr
c
      lzr=1
      do l=1,lp
        if(dep(l).eq.depr)lzr=l
      enddo
c
c     determine layer no of each depth
c
      li=1
      zswap=h(1)
      nno(1)=1
      do l=2,lp
        if(dep(l).ge.zswap.and.li.lt.n0)then
          li=li+1
          zswap=zswap+h(li)
        endif
        nno(l)=li
      enddo
      write(*,*)' Receiver-site layered structure:'
      write(*,'(7a)')'    no ','  z(km)  ',
     &               '  vp(km/s) ','  vs(km/s) ',' ro(g/cm^3)',
     &               '    qp   ','    qs'
      z=0.d0
      do i=1,n0
        write(*,1000)2*i-1,z/KM2M,vp(i)/KM2M,
     &               vs(i)/KM2M,ro(i)/KM2M,qp(i),qs(i)
        z=z+h(i)
        if(i.lt.n0)then
          write(*,1000)2*i,z/KM2M,vp(i)/KM2M,
     &               vs(i)/KM2M,ro(i)/KM2M,qp(i),qs(i)
        endif
      enddo
1000  format(i5,f12.2,3f11.4,2f8.1)
c
      return
      end
