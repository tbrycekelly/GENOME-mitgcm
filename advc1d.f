
      subroutine advc1d(i,j,n,nt,dz,vrbos,errcon,deposit)
c
c --- ppm-based 1-dim transport routine, extracted from fct3d.f
c --- originator: rainer bleck
c
c --- kk       - number of layers
c --- fld(kk)  - mixing ratio of dependent variable in each layer
c --- z(kk+1)  - interface depth; z(k+1)-z(k) = thickness of layer k
c --- dz(kk+1) - elements of -fld- travel from z(k)-dz(k) to z(k) during
c ---            present time step
c --- vrbos  - if c.true., print diagnostic messages
c --- errcon - set to .true. if error condition encountered
c
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
      implicit none
      include 'common_blocks.h'
c
      integer, intent(IN) :: i,j,nt,n
      real,    intent(IN) :: dz(kdm+1)
      logical, intent(IN) :: vrbos,deposit
      logical,intent(OUT) :: errcon
c
      real thk(kdm),flux(kdm+1),div(kdm)
      real a(kdm),b(kdm),c(kdm),dx,fcdx,yl,yr
      real amount,bfore,after,scale,slab,dslab
      integer k,kp,kd
      real, parameter :: athird=1./3.
      real, parameter :: small=9.806e-8   !  (for z,dz given in press)
c
      if (vrbos) write (*,100)
     & 'entering advc1d:  old loc''n  new loc''n   variable',
     & (k,p(i,j,k)/onem,(p(i,j,k)+dz(k))/onem,
     & tracer(i,j,k,n,nt),k=1,kk),
     &  kk+1,p(i,j,kk+1)/onem,(p(i,j,kk+1)+dz(kk+1))/onem
 100  format (a/(i15,2f15.7,es15.6))
c
      bfore=0.
      scale=1.e-33
      do k=1,kk
        thk(k)=p(i,j,k+1)-p(i,j,k)
        bfore=bfore+tracer(i,j,k,n,nt)*thk(k)
        scale=max(scale,abs(tracer(i,j,k,n,nt)))
      end do
c
cc --- exit if parcels overtake each other
c      errcon=.false.
      do k=2,kk
        if (p(i,j,k)+dz(k).lt.p(i,j,k-1)+dz(k-1) .or.
     &      p(i,j,k)+dz(k).gt.p(i,j,k+1)+dz(k+1)) then
c          write (*,'(a,i3)') 'error advc1d -- crossover at k =',k
c        write (*,100)
c     & 'entering advc1d:  old loc''n  new loc''n   variable',
c     & (k,p(i,j,kd)/onem,(p(i,j,kd)+dz(kd))/onem,
c     & tracer(i,j,kd,n,nt),kd=1,kk),
c     &  kk+1,p(i,j,kk+1)/onem,(p(i,j,kk+1)+dz(kk+1))/onem

          errcon=.true.
          return
        end if
      end do
c
c --- start by filling massless cells with data from layer(s) above or below
c
      do 17 k=kk-1,1,-1
 17   tracer(i,j,k,n,nt)=(tracer(i,j,k,n,nt)*thk(k)+
     &                    tracer(i,j,k+1,n,nt)*small)
     &      /(       thk(k)+         small)
      do 18 k=2,kk
 18   tracer(i,j,k,n,nt)=(tracer(i,j,k,n,nt)*thk(k)+
     &                    tracer(i,j,k-1,n,nt)*small)
     &      /(       thk(k)+         small)
c
c --- fit 0th, 1st, or 2nd deg. polynomial to -fld- in each cell
      a(1 )=tracer(i,j,1,n,nt)
      b(1 )=0.
      c(1 )=0.
      a(kk)=tracer(i,j,kk,n,nt)
      b(kk)=0.
      c(kk)=0.
c
      do 16 k=2,kk-1
c --- uncomment one of the following 3 options to activate pcm,plm,ppm resp.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise constant method:
ccc      a(k)=fld(k)
ccc      b(k)=0.
ccc      c(k)=0.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise linear method:
c --- fit linear function a+bx to tracer in each cell (-.5 < x < +.5)
ccc      a(k)=fld(k)
ccc      b(k)=0.
ccc      if (fld(k).le.min(fld(k-1),fld(k+1)) .or.
ccc     &    fld(k).ge.max(fld(k-1),fld(k+1))) then
ccc        b(k)=0.
ccc      else if ((fld(k+1)-fld(k-1))*(fld(k-1)+fld(k+1)
ccc     &  -2.*fld(k)).gt.0.) then
ccc        b(k)=fld(k)-fld(k-1)
ccc      else
ccc        b(k)=fld(k+1)-fld(k)
ccc      end if
ccc      c(k)=0.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise parabolic method:
c --- fit parabola a+bx+cx^2 to tracer in each cell (-.5 < x < +.5)
      yl=.5*(tracer(i,j,k-1,n,nt)+tracer(i,j,k,n,nt))
      yr=.5*(tracer(i,j,k+1,n,nt)+tracer(i,j,k,n,nt))
      a(k)=1.5*tracer(i,j,k,n,nt)-.25*(yl+yr)
      b(k)=yr-yl
      c(k)=6.*(.5*(yl+yr)-tracer(i,j,k,n,nt))
      if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-tracer(i,j,k,n,nt))) then
c --- apex of parabola lies inside interval [-.5,+.5], implying an over-
c --- or undershoot situation. change curve to prevent over/undershoots.
        if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-tracer(i,j,k,n,nt))) then
c --- put apex of parabola on edge of interval [-.5,+.5]
          if ((yr-yl)*(.5*(yl+yr)-tracer(i,j,k,n,nt)) .gt. 0.) then
c --- apex at x=-.5
            a(k)=.25*(3.*tracer(i,j,k,n,nt)+yl)
            c(k)=3.*(tracer(i,j,k,n,nt)-yl)
            b(k)=c(k)
          else
c --- apex at x=+.5
            a(k)=.25*(3.*tracer(i,j,k,n,nt)+yr)
            c(k)=3.*(tracer(i,j,k,n,nt)-yr)
            b(k)=-c(k)
          end if
        else                    !  -1/6 < x < +1/6
c --- moving apex won't help. replace parabola by constant.
          a(k)=tracer(i,j,k,n,nt)
          b(k)=0.
          c(k)=0.
        end if
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 16   continue
c
c
c --- get flux by summing -fld- over upstream slab of thickness -dz-
c
      do 23 k=1,kk+1
      slab=small
      if (dz(k).lt.0.) then                     ! -fld- moves up
        if (k.eq.kk+1) then
          amount=0.                             ! influx from outside domain
        else
          amount=slab*tracer(i,j,k,n,nt)
          kp=k-1
 24       kp=kp+1
          if (kp.le.kk) then
            if (thk(kp).gt.0.) then
              dslab=min(slab+thk(kp),-dz(k),p(i,j,kk+1)-p(i,j,k))
     &             -min(slab        ,-dz(k),p(i,j,kk+1)-p(i,j,k))
              dx=dslab/thk(kp)
              fcdx=a(kp)
     &          +b(kp)*.5*(dx-1.)               !  not needed in pcm
     &          +c(kp)*(.25-dx*(.5-dx*athird))  !  not needed in pcm,plm
              amount=amount+fcdx*dslab
              slab=slab+dslab
            end if
          else                                  ! reaching outside domain
            slab=-dz(k)
          end if
          if (slab.lt.-dz(k)) go to 24
        end if
      else if (dz(k).gt.0.) then                ! -fld- moves down
        if (k.eq.1) then
          amount=0.                             ! influx from outside domain
        else
          amount=slab*tracer(i,j,k-1,n,nt)
          kp=k
 25       kp=kp-1
          if (kp.ge.1) then
            if (thk(kp).gt.0.) then
              dslab=min(slab+thk(kp), dz(k),p(i,j,k)-p(i,j,1))
     &             -min(slab        , dz(k),p(i,j,k)-p(i,j,1))
              dx=dslab/thk(kp)
              fcdx=a(kp)
     &          +b(kp)*.5*(1.-dx)               !  not needed in pcm
     &          +c(kp)*(.25-dx*(.5-dx*athird))  !  not needed in pcm,plm
              amount=amount+fcdx*dslab
              slab=slab+dslab
            end if
          else                                  ! reaching outside domain
            slab=dz(k)
          end if
          if (slab.lt.dz(k)) go to 25
        end if
      else                                      !  dz = 0
        amount=0.
      end if
 23   flux(k)=dz(k)*amount/slab
c
      do 26 k=1,kk
 26   div(k)=flux(k+1)-flux(k)
c
      do 4 k=1,kk
      amount=max(0.,tracer(i,j,k,n,nt)*thk(k)-div(k))
 4    tracer(i,j,k,n,nt)=(tracer(i,j,k,n,nt)*small+amount)/
     .                    (small+thk(k))
c
      after=flux(kk+1)-flux(1)                  !  account for outflow loss
      do k=1,kk
        after=after+tracer(i,j,k,n,nt)*thk(k)
      end do
c      if (abs(bfore-after)*kk.gt.1.e-9*scale*p(i,j,kk+1)) then
c        write (*,104) '  advc1d - bad column intgl.:',i,j,bfore,after
c        write (*,101)
c     & 'exiting  advc1d:  old loc''n  new loc''n   variable',
c     & (k,p(i,j,k)/onem,(p(i,j,k)+dz(k))/onem,
c     &  tracer(i,j,k,n,nt),dz(k)/onem,flux(k),k=1,kk),
c     &  kk+1,p(i,j,kk+1)/onem,(p(i,j,kk+1)+dz(kk+1))/onem,0,
c     &  dz(kk+1)/onem,flux(kk+1)
c
c        errcon=.true.
c      end if
 104  format (a,2i4,2es19.10)
c
      if (vrbos) write (*,101)
     & 'exiting  advc1d:  old loc''n  new loc''n   variable',
     & (k,p(i,j,k)/onem,(p(i,j,k)+dz(k))/onem,
     &  tracer(i,j,k,n,nt),dz(k)/onem,flux(k),k=1,kk),
     &  kk+1,p(i,j,kk+1)/onem,(p(i,j,kk+1)+dz(kk+1))/onem,0,
     &  dz(kk+1)/onem,flux(kk+1)
 101  format (a/(i15,2f15.7,3es15.6))
c
      return
    end
c end subroutine advc1d
