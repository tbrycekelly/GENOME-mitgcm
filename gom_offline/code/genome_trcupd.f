      subroutine trcupd(m,n)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'

      integer m
      integer n

! --- -----------------------------------------------------------
! --- tracer-specific operations (side-wall relaxation in thermf)
! --- -----------------------------------------------------------

      integer i
      integer j
      integer k
      integer ktr
      integer l
      real    beta_b
      real    pijk
      real    pijkp
      real    q

      margin = 0  ! no horizontal derivatives

      !!!
      !!! Section 1 - Relaxation boundary conditions
      !!!
      ! Loop through each tracer and apply the type of relaxation
      ! condition applicable to them (based on trcflg):
      !     trcflg == 0   ??based on bc??
      !     trcflg == 1   30 day euphotic half life
      !
      do ktr = 1,ntracr
        !!!
        !!! Section 1.1 - Relaxation based on bc
        !!!
        if (trcflg(ktr).eq.0) then
          if (trcrlx) then ! If trcrlx is on then do the relaxation
            ! tracer always trwall, when non-zero, at surface
!$OMP       PARALLEL DO PRIVATE(j,k,l,i,ktr,q)
!$OMP&      SCHEDULE(STATIC,jblk)
            do j=1-margin,jj+margin
              do l=1,isp(j)
                do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                  q = trwall(i,j,1,lc0,ktr) * wc0
     &               + trwall(i,j,1,lc1,ktr) * wc1
     &               + trwall(i,j,1,lc2,ktr) * wc2
     &               + trwall(i,j,1,lc3,ktr) * wc3
                  if (q.gt.0.0) then
                    tracer(i,j,1,n,ktr) = q
                  endif
                enddo !i
              enddo !l
            enddo !j
!$OMP       END PARALLEL DO

          elseif (.not. trcrlx) then
            ! tracer always 10.0 at surface
!$OMP       PARALLEL DO PRIVATE(j,k,l,i,ktr)
!$OMP&      SCHEDULE(STATIC,jblk)
            do j=1-margin,jj+margin
              do l=1,isp(j)
                do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                  tracer(i,j,1,n,ktr) = 10.0
                enddo !i
              enddo !l
            enddo !j
!$OMP       END PARALLEL DO
          endif !trcrlx:else

        !!!
        !!! Section 1.2 Relaxation based on time
        !!!
        elseif (trcflg(ktr).eq.1) then
! ---     psudo-silicate, half-life of 30 days in euphotic zone
          q = 1.0 - delt1 / (30.0 * 86400.0)

!$OMP     PARALLEL DO PRIVATE(j,k,l,i,ktr,pijk,pijkp,beta_b)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do l=1,isp(j)
              ! Look at stencil
              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                if (jerlv0.eq.0) then
                  beta_b = qonem * ( akpar(i,j,lk0) * wk0
     &                              +akpar(i,j,lk1) * wk1
     &                              +akpar(i,j,lk2) * wk2
     &                              +akpar(i,j,lk3) * wk3 )
                else
                  beta_b = betabl(jerlov(i,j))
                endif

                pijkp = 0.0

                do k = 1,kk
                  pijk = pijkp
                  pijkp = pijk + dp(i,j,k,n)

                  ! Apply value to euphotic cells
                  if (0.5*(pijk+pijkp)*beta_b .lt. 4.0) then
                    tracer(i,j,k,n,ktr) = q * tracer(i,j,k,n,ktr)
                  else
                    exit  !too deep
                  endif
                enddo ! k
              enddo !i
            enddo !l
          enddo !j
!$OMP     END PARALLEL DO


          !!! Section 1.3 - Run Roca Model
          call trcupd_100(m,n,ktr)
        enddo !ktr
        !!! End Section 1
        return
      end subroutine trcupd
