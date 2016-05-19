      subroutine qsswvint(srate)
      implicit none
c
      double precision srate
c
      include 'qssglobal.h'
c
      integer i,istp,lf,ik,nk1,nk2,nnk
      integer unit(4)
      double precision f,k,dk
      double precision pi,pi2
      double precision kcut(4)
      double complex yc(3,4)
c
      integer nd
      data nd/5/
c
      pi=4.d0*datan(1.d0)
      pi2=2.d0*pi
      dk=pi/(srate*dismax)
      nnk=idint(pi2*dble(nf)*df*slw(4)/dk)+6
c
      write(*,'(a)')'  calculate f-k spectra ...'
c
      open(20,file=greeninfo,status='unknown')
      write(20,'(a)')'# input file'
      write(20,'(a)')inputfile
      write(20,'(a)')'# zs[km]'
      write(20,'(a)')zs0/1000.d0
      write(20,'(a)')'# nf, df[Hz], fi[Hz]'
      write(20,'(i6,2E16.8)')nf,df,fi
      write(20,'(a)')'# sw1-[s/km], srate'
      write(20,'(4f10.4)')(slw(i)*1000.d0,i=1,4),srate
      write(20,'(a)')'# nnk, dk[1/m]'
      write(20,'(i8,E16.8)')nnk,dk
      write(20,'(a)')'# nd'
      write(20,'(i4)')nd
      close(20)
c
      open(30,file=greenfile,
     &        form='unformatted',status='unknown')
      write(30)zs0
      write(30)nf,df,fi
      write(30)nnk,dk
      write(30)(slw(i),i=1,4),srate
      write(30)nd
c
      do lf=2,nf
        f=dble(lf-1)*df
c
        call qssqmodel(f)
c
        do i=1,4
          kcut(i)=pi2*f*slw(i)
        enddo
        nk1=max0(1,idint(kcut(1)/dk)-nd)
        nk2=1+idint(kcut(4)/dk)+nd
c
        write(30)f,nk1,nk2
c
        do ik=nk1,nk2
          k=dble(ik)*dk
          call qsskern(yc,f,k)
          write(30)yc(1,1),yc(2,1),
     &             yc(1,2),yc(2,2),yc(3,2),
     &             yc(1,3),yc(2,3),yc(3,3),
     &             yc(1,4),yc(2,4)
        enddo
c
        write(*,'(i6,a,E13.6,a,i7)')lf,'.',f,
     &      'Hz: slowness samples = ',1+nk2-nk1
      enddo
      close(30)
c
      return
      end