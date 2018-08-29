subroutine initrc(mnth)
use mod_xc    ! HYCOM communication interface
use mod_pipe  ! HYCOM debugging interface
implicit none
c
include 'common_blocks.h'
c
integer mnth
c
c --- --------------------------
c --- initializatize all tracers
c --- --------------------------
c
logical    lpipe_initrc
parameter (lpipe_initrc=.false.)
c
character ptxt*12,cformat*99
integer   i,ibio,nbio,j,k,ktr,l,bmodel
real      bio_n,bio_p,zk
real      pwij(kk+1),trwij(kk,ntracr),
&          prij(kk+1),trcij(kk,ntracr)
cvic
integer   mxtracr
c
if (ntracr.eq.0) then
  return  ! no tracer
endif
c
c --- expand trcflg to allow for number of biology fields.
c
nbio = 0
ibio = 0
do ktr= 1,ntracr+1
  if     (ktr.ne.ntracr+1 .and.
&          trcflg(min(ktr,ntracr)).eq.9) then
    if     (ibio.eq.0) then !start biology
      ibio = ktr
    endif
  elseif (ibio.ne.0) then !end biology
    nbio = ktr-ibio
    if     (nbio.eq.22) then
c ---       Stukel 5-Phyto Diazotroph Model
      trcflg(ibio)   =  923
      trcflg(ibio+1) = -923
      trcflg(ibio+2) = -923
      trcflg(ibio+3) = -923
      trcflg(ibio+4) = -923
      trcflg(ibio+5) = -923
      trcflg(ibio+6) = -923
      trcflg(ibio+7) = -923
      trcflg(ibio+8) = -923
      trcflg(ibio+9) = -923
      trcflg(ibio+10) = -923
      trcflg(ibio+11) = -923
      trcflg(ibio+12) = -923
      trcflg(ibio+13) = -923
      trcflg(ibio+14) = -923
      trcflg(ibio+15) = -923
      trcflg(ibio+16) = -923
      trcflg(ibio+17) = -923
      trcflg(ibio+18) = -923
      trcflg(ibio+19) = -923
      trcflg(ibio+20) = -923
      trcflg(ibio+21) = -923
      trcflg(ibio+22) = -923
      ibio = 0
elseif (nbio .gt.22) then
c ---     ROCA model
      trcflg(ibio)   =  100
      do k=ibio+1,ibio+nbio-1
  trcflg(k)   =  -100
enddo
    else
c ---       unknown standard biology.
      if (mnproc.eq.1) then
      write(lp,'(/ 2a,i3 /)')
&        'error - trcflg=9 (standard biology) expects',
&        ' 3/4/6/7/8 consecutive tracers but have',nbio
*    &        ' 3/4/6/7/8/9 consecutive tracers but have',nbio
      call flush(lp)
      endif !1st tile
      call xcstop('(trcini)')
             stop '(trcini)'
    endif
  endif
enddo
c
if (mnproc.eq.1) then
write(lp,*)
do k= 1,ntracr
  write(lp,'(a,i3,i6)') 'initrc: k,trcflg =',k,trcflg(k)
enddo
write(lp,*)
endif !1st tile
c
if     (nbio.gt.0) then
c
c ---   input bio-tracer parameters.
c ---   note that multiple sets of bio-tracers are allowed,
c ---   each is read from tracer.input in tracer order.
c
  open(unit=uoff+99,file=trim(flnminp)//'tracer.input')
  do ktr= 1,ntracr
    if     (trcflg(ktr).eq.923) then
c ---       Stukel Diazo Model
c            call trcupd_923(1,2,-ktr)
elseif (trcflg(ktr).eq.100) then
c ---       ROCA Model
      call trcupd_100(1,2,-ktr)
    endif
  enddo
  close(unit=uoff+99)
endif
c
if     (trcrin) then
  return  ! tracer from restart
endif
c
margin = 0
c
if     (iniflg.eq.2) then  ! use climatology
  call rdrlax(mnth,1)
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,ktr,pwij,trwij,prij,trcij)
!$OMP&           SCHEDULE(STATIC,jblk)
  do j=1-margin,jj+margin
    do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
        prij(1)=0.0
        do k=1,kk
          prij(k+1)=prij(k)+dp(i,j,k,1)
          pwij(k)  =pwall(i,j,k,1)
          do ktr= 1,ntracr
            trwij(k,ktr)=trwall(i,j,k,1,ktr)
          enddo !ktr
        enddo !k
        pwij(kk+1)=prij(kk+1)
*             call plctrc(trwij,pwij,kk,ntracr,
*    &                    trcij,prij,kk        )
        call plmtrc(trwij,pwij,kk,ntracr,
&                    trcij,prij,kk        )
        do k=1,kk
          do ktr= 1,ntracr
            tracer(i,j,k,1,ktr)=trcij(k,ktr)
            tracer(i,j,k,2,ktr)=trcij(k,ktr)
          enddo !ktr
        enddo !k
      enddo !i
    enddo !l
  enddo !j
!$OMP   END PARALLEL DO
else ! analytic inititalization
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,ktr)
!$OMP&           SCHEDULE(STATIC,jblk)
  do j=1-margin,jj+margin
    do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
        p(i,j,1)=0.0
        do k=1,kk
          p(i,j,k+1)=p(i,j,k)+dp(i,j,k,1)
          do ktr= 1,ntracr
            if     (trcflg(ktr).eq.0) then !100% in the mixed layer
              if     (p(i,j,k).le.dpmixl(i,j,1)) then
                tracer(i,j,k,1,ktr)=10.0
                tracer(i,j,k,2,ktr)=10.0
              else
                tracer(i,j,k,1,ktr)=0.0
                tracer(i,j,k,2,ktr)=0.0
              endif
            elseif (trcflg(ktr).eq.1) then !20 below euphotic zone
              if     (p(i,j,k)*betabl(jerlv0).lt.4.0) then
                tracer(i,j,k,1,ktr)=0.0
                tracer(i,j,k,2,ktr)=0.0
              else
                tracer(i,j,k,1,ktr)=20.0  ! mg/m^3
                tracer(i,j,k,2,ktr)=20.0  ! mg/m^3
              endif
            elseif (trcflg(ktr).eq.2) then !temperature
              tracer(i,j,k,1,ktr)=temp(i,j,k,1)
              tracer(i,j,k,2,ktr)=temp(i,j,k,1)
            elseif (trcflg(ktr).eq.3) then !fully passive
              tracer(i,j,k,1,ktr)=0.0 !should never get here
              tracer(i,j,k,2,ktr)=0.0 !should never get here
            elseif (trcflg(ktr).eq.904 .or.
&                    trcflg(ktr).eq.903     ) then !NPZD or NPZ
              zk = 0.5*(p(i,j,k+1)+p(i,j,k))*qonem
              if     (zk.le.300.0) then
                ! 0.1 at 300m, 1.0 at 100m, 2.025 at 0m
                bio_p = 0.1 + (300.0-zk)**2 * (0.9/200.0**2)
              elseif (zk.le.900.0) then
                ! 0.1 at 300m, 0.0 at 900m
                bio_p = (900.0-zk) * 0.1/600.0
              else
                bio_p = 0.0
              endif
              if     (temp(i,j,k,1).lt. 6.0) then
                bio_n = 37.0
              elseif (temp(i,j,k,1).gt.27.0) then
                bio_n =  0.0
              else
*                     bio_n = (27.0-temp(i,j,k,1)) * 37.0/21.0
                bio_n = 39.3116-1.335*temp(i,j,k,1)
              endif
              tracer(i,j,k,1,ktr  )=bio_n  !N
              tracer(i,j,k,2,ktr  )=bio_n
              tracer(i,j,k,1,ktr+1)=bio_p  !P
              tracer(i,j,k,2,ktr+1)=bio_p
              tracer(i,j,k,1,ktr+2)=bio_p  !Z=P
              tracer(i,j,k,2,ktr+2)=bio_p
              if     (trcflg(ktr).eq.904) then
                tracer(i,j,k,1,ktr+3)=bio_p + 1.0  !D=P+1
                tracer(i,j,k,2,ktr+3)=bio_p + 1.0
              endif
            endif !trcflg
          enddo !ktr
        enddo !k
      enddo !i
    enddo !l
  enddo !j
!$OMP   END PARALLEL DO
endif !iniflg.eq.2:else
c
if     (lpipe .and. lpipe_initrc) then
   do ktr= 1,ntracr
     do k= 1,kk
       write (ptxt,'(a4,i2.2,a3,i3)') 'trc.',ktr,' k=',k
       call pipe_compare_sym1(tracer(1-nbdy,1-nbdy,k,1,ktr),
&                              ip,ptxt)
     enddo !k
   enddo !ktr
 endif !lpipe.and.lpipe_initrc
c
cvic
mxtracr=min(ntracr,23)
if     (itest.gt.0 .and. jtest.gt.0) then
   write(cformat,'(a,i2,a,i2,a)')
   &     'a / (23x,i3,2f8.2,', mxtracr,'f8.4))'
&     '(i9,2i5,a,',mxtracr,
   write (lp,cformat)
&     nstep,i0+itest,j0+jtest,
&     '  istate:  thkns    dpth',
&     ('  tracer',ktr=1,mxtracr),
&     (k,
&      dp(itest,jtest,k,1)*qonem,
&      (p(itest,jtest,k+1)+p(itest,jtest,k))*0.5*qonem,
&      (tracer(itest,jtest,k,1,ktr),ktr=1,mxtracr),
&      k=1,kk)
   write(lp,'(23x,a,8x,f8.2)') 'bot',depths(itest,jtest)
endif !test tile
cvic
call xcsync(flush_lp)
c
return
end
c end subroutine initrc
