      subroutine qsrwvint(mxx,myy,mzz,mxy,myz,mzx,greenfk,
     &                    ibndp,flw,fup)
      implicit none
c
      integer ibndp
      double precision mxx,myy,mzz,mxy,myz,mzx,flw,fup
      character*800 greenfk
c
      include 'qsrglobal.h'
c
      integer i,ik,istp,lf,nbsj,nd,ik1,ik2,nk1,nk2,idtrans
      integer ics(4),ms(4)
      double precision f,zs,rs,k,dk,delta,fac,ss1,cs1,ss2,cs2,dismax
      double precision srate,geospr,expl,ss00,ss45,ds00,ds90,clvd
      double precision slw(4),kcut(4)
      double complex ck,ck2,cdk,c2dk,cdk2,czs2,cfac,bcs,bss
      double complex cics(4),cm2(4,4),grns(3,4),yb(3)
      double complex ak(3,4),yk(6,4),y(3,4,nbsjmax),swap(nbsjmax)
      double complex bpf(nfmax)
c
      double precision taper
c
      double complex c2
      data c2/(2.d0,0.d0)/
c
c     ics = 1  when the azmuth-factor is cos(ms*theta) for poloidal mode
c             (psv) and sin(ms*theta) for the toroidal mode (sh);
c     ics = -1 otherwise.
c
      ms(1)=0
      ics(1)=1
      ms(2)=2
      ics(2)=-1
      ms(3)=1
      ics(3)=1
      ms(4)=0
      ics(4)=1
c
      do istp=1,4
        cics(istp)=dcmplx(dble(ics(istp)),0.d0)
        cm2(1,istp)=dcmplx(dble(ms(istp)**2),0.d0)
        cm2(2,istp)=dcmplx(dble((ms(istp)-1)**2),0.d0)
        cm2(3,istp)=dcmplx(dble((ms(istp)+1)**2),0.d0)
        cm2(4,istp)=cm2(1,istp)
      enddo
c
      open(30,file=greenfk,
     &     form='unformatted',status='old')
      read(30)zs
      read(30)nf,df,fi
      read(30)nbsj,dk
      read(30)(slw(i),i=1,4),srate
      read(30)nd
c
      nt=2*nf
      dt=1.d0/(dble(2*nf)*df)
      twindow=dble(nt-1)*dt
c
      cdk=dcmplx(dk,0.d0)
      c2dk=dcmplx(2.d0*dk,0.d0)
      cdk2=dcmplx(dk*dk,0.d0)
      czs2=dcmplx(zs*zs,0.d0)
c
      if(stdis.gt.PI2/(dmax1(srate,1.d0)*dk))then
        print *,' Warning in qsrwvint:'
     &        //' too large epicentral distance or'
     &        //' too low slowness sampling rate in f-k spectra!'
      endif
c
      ndtrans=min0(idint(PI2/(dk*stdis*5.d0)),nd,4)
c
      if(nbsj.gt.nbsjmax)then
        stop ' Error in qsrwvint: nbsjmax too small!'
      else
        geospr=1.d0/(zs**2+stdis**2)**ndtrans
        call qsrbsj(dk,nbsj,geospr)
      endif
c
      do lf=1,nf
        do i=1,3
          seisf(lf,i)=(0.d0,0.d0)
        enddo
      enddo
c
      expl=(mxx+myy+mzz)/3.d0
      ss00=mxy
      ds00=mzx
      clvd=mzz-expl
      ss45=0.5d0*(mxx-myy)
      ds90=myz
      cs1=dcos(stazi)
      ss1=dsin(stazi)
      cs2=dcos(2.d0*stazi)
      ss2=dsin(2.d0*stazi)
      do istp=1,4
        do i=1,3
          ak(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
      do lf=2,nf
        read(30)f,nk1,nk2
        call qsrqmodel(f)
        do ik=nk1,nk2
          k=dble(ik)*dk
          ck=dcmplx(k,0.d0)
          read(30)ak(1,1),ak(2,1),
     &            ak(1,2),ak(2,2),ak(3,2),
     &            ak(1,3),ak(2,3),ak(3,3),
     &            ak(1,4),ak(2,4)
          call qsrkern(f,k,ak,yk)
          do istp=1,4
            y(1,istp,ik)=yk(1,istp)
            y(2,istp,ik)=( yk(3,istp)+cics(istp)*yk(5,istp))/c2
            y(3,istp,ik)=(-yk(3,istp)+cics(istp)*yk(5,istp))/c2
          enddo
        enddo
        write(*,'(i6,a,E13.6,a,i7)')lf,'.',f,
     &      'Hz: slowness samples = ',1+nk2-nk1
c
        do idtrans=1,ndtrans
          ik1=max0(1,nk1)+idtrans
          ik2=nk2-idtrans
          do istp=1,4
            do i=1,3
              do ik=ik1-1,ik2+1
                swap(ik)=y(i,istp,ik)
              enddo
              do ik=ik1,ik2
                ck=dcmplx(dble(ik)*dk,0.d0)
                ck2=ck*ck
                y(i,istp,ik)=swap(ik)*(czs2+cm2(i,istp)/ck2)
     &            -(swap(ik+1)-swap(ik-1))/c2dk/ck
     &            -(swap(ik+1)-c2*swap(ik)+swap(ik-1))/cdk2
              enddo
            enddo
          enddo
        enddo
c
        do istp=1,4
          do i=1,3
            grns(i,istp)=(0.d0,0.d0)
          enddo
        enddo
c
        ik1=max0(1,nk1+ndtrans)
        ik2=nk2-ndtrans
        do ik=ik1,ik2
          k=dble(ik)*dk
          do i=1,4
            kcut(i)=PI2*f*slw(i)
          enddo
          fac=k*dk
          fac=fac*taper(k,kcut(1),kcut(2),kcut(3),kcut(4))
          cfac=dcmplx(fac,0.d0)
c
          do istp=1,4
            do i=1,3
              y(i,istp,ik)=y(i,istp,ik)*cfac
            enddo
            yb(1)=y(1,istp,ik)*dcmplx(bsj(ik,ms(istp)),0.d0)
            yb(2)=y(2,istp,ik)*dcmplx(bsj(ik,ms(istp)-1),0.d0)
            yb(3)=y(3,istp,ik)*dcmplx(bsj(ik,ms(istp)+1),0.d0)
c
            grns(1,istp)=grns(1,istp)+yb(1)
            grns(2,istp)=grns(2,istp)+yb(2)+yb(3)
            grns(3,istp)=grns(3,istp)-cics(istp)*(yb(2)-yb(3))
          enddo
        enddo
        seisf(lf,1)=seisf(lf,1)+grns(1,1)*dcmplx(expl,0.d0)
     &              +grns(1,4)*dcmplx(clvd,0.d0)
     &              +grns(1,2)*dcmplx(ss00*ss2+ss45*cs2,0.d0)
     &              +grns(1,3)*dcmplx(ds00*cs1+ds90*ss1,0.d0)
        seisf(lf,2)=seisf(lf,2)+grns(2,1)*dcmplx(expl,0.d0)
     &              +grns(2,4)*dcmplx(clvd,0.d0)
     &              +grns(2,2)*dcmplx(ss00*ss2+ss45*cs2,0.d0)
     &              +grns(2,3)*dcmplx(ds00*cs1+ds90*ss1,0.d0)
        seisf(lf,3)=seisf(lf,3)
     &              +grns(3,2)*dcmplx(ss00*cs2-ss45*ss2,0.d0)
     &              +grns(3,3)*dcmplx(ds00*ss1-ds90*cs1,0.d0)
      enddo
c
      close(30)
c
c     Butterworth bandpass filter
c
      if(ibndp.gt.0)then
        call bandpass(ibndp,flw,fup,df,nf,bpf)
      else
        do lf=1,nf
          bpf(lf)=(1.d0,0.d0)
        enddo
      endif
c
c     amplitude correction when using the flat-earth transform
c     see Mueller (1977) for n = -2
c
      rs=REARTH*dexp(-zs/REARTH)
      delta=stdis/REARTH
      cfac=dcmplx((REARTH/rs)**2*dsqrt(delta/dsin(delta)),0.d0)
      do i=1,3
        do lf=1,nf
          seisf(lf,i)=seisf(lf,i)*cfac*bpf(lf)
        enddo
      enddo
c
      return
      end
