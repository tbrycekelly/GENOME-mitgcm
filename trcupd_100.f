        subroutine trcupd_100(m,n, ibio)
        use mod_xc  ! HYCOM communication interface
        implicit none
        include 'common_blocks.h'
        integer m
        integer n
        integer ibio
c
c --- -------------------------------------------------
c --- tracer-specific operations for ROCA
c --- -------------------------------------------------
c
      integer, parameter :: itermax=3             ! number of iterations
      integer, parameter :: mxsubs=8              ! max number of substrates - real number must be less than sbio
      integer, parameter :: mxproc=20             ! more than mxsubs
      integer, parameter :: mxgene=25	            !max space alloted for genes
      real, save :: chldel               ! number of days to lag chlorophyll
      real, save :: rivdin               ! river DIN to input
      real, save :: rivtss               ! river tss to input
      real, save :: rivdetrs             ! river small detritus to input
      real, save :: rivdetrl             ! river large detritus to input
      real, save :: rivdon(mxsubs)       ! river organic matter to input
      real, save :: alpha(mxsubs)        ! uptake sent to O and IO nuts
      real, save :: alphao(mxsubs)       ! nonassimilated to each kind of DOM
      real, save :: alphanh4			 ! uptake sent to IO and O nuts in nitrification
      real, save :: ae(mxsubs)           ! assim efficiency for each substrate
      real, save :: aei(mxsubs)          ! 1.- assim efficiency for each substrate
      real, save :: respi(mxsubs)        ! growth efficiency for each substrate
      real, save :: ga(mxsubs)           ! excretion partitioning to NH4
      real, save :: gai(mxsubs)          ! excretion partitioning to DON
      real, save :: betab                ! partitioning of dead bio to detritus
      real, save :: betap                ! partitioning of dead bio to NH4
      real, save :: betai                ! partitioning of dead bio to don
      real, save :: bigdet               ! cutoff for bug size to go to big det
      real, save :: mxldet               ! max large det conc for sinking fraction
      real, save :: betasl(mxtrcr)       ! fraction of dead bug going to small det
      real, save :: betarl(mxsubs)       ! fraction of dead bug going to each dom sum to 1
      real, save :: domtocdom(mxsubs)    ! fraction of dom that is colored
      real, save :: bp(mxproc,mxtrcr)    ! array of parameters read in for each bio
      real, save :: trans(mxgene)		 ! transcription number for each gene
      real, save :: aggreg			     ! aggregation rate for small detritus to large
      real, save :: gmax			     ! fraction of community biomass below which an organism is removed

      real, parameter :: rs  = 0.6            ! fraction of incident Q that is longwave
      real, parameter :: rs2 = 0.4            ! 580 nm spectral division in pslight
!  &  wkr = 1.6666,                           ! red Ks for water - bad!
      real, parameter :: wkr2= 0.35           ! red Ks for water in pslight
      real, parameter :: wkb = 0.03           ! blue Ks for water
      real, parameter :: dsec = 86400.0       ! seconds in one day
      real, parameter :: par = 0.45           ! fraction of incident sw that is par
      real, parameter :: wk = 0.022           ! blue Ks for chl
      real, parameter :: ntochl = 1.59        !nitrogen to chl conversion
      real, parameter :: basalmet = 1.1e-7    ! .01 per day basal metabolism

      integer, save :: ichl                 ! place for chlorophyl
      integer, save :: icdm                 ! place for cdom estimate
      integer, save :: inh4e                ! place for nitrification
      integer, save :: ilight               ! place for light
      integer, save :: ihet                 ! growth dure to het
      integer, save :: iauto                ! growth due to autotrophy
      integer, save :: ino3,					        ! growth due to nitrate
      integer, save :: isflxb,					      ! sinking flux due to n2-fix
      integer, save :: ibug,				        ! bug ID location
      integer, save :: imot,					    ! gene for motility/particle attach
      integer, save :: irhod,					      ! gene for rhodopsin
      integer, save :: ngene,					      ! number of genes
      integer, save :: gene(mxgene,mxtrcr)		! gene array

c - process places
      integer, parameter :: jbug = 1					       ! bug number
      integer, parameter :: jsink = 2               ! place for sinking rate
      integer, parameter :: jmaxgrow = 3            ! place for max growth
      integer, parameter :: jdin = 4					       ! DIN uptake
      integer, parameter :: jnh4 = 5                ! ammonium uptake
      integer, parameter :: jdetrs = 6              ! small detritus
      integer, parameter :: jdetrl = 7              ! large detritus
      integer, parameter :: jn2 = 8                 ! nitrogen fixation
      integer, parameter :: jnh4e = 9               ! nitrification
      integer, parameter :: jsize = 10              ! size
      integer, parameter :: jdeath = 11             ! individual death rate (virus)
      integer, parameter :: jddeath = 12            ! density dependent death
      integer, parameter :: jorgs = 13              ! place to start organic subs
      integer, parameter :: jorge = 14              ! place to end organics
      integer, parameter :: jphoto = 15             ! place for photosynthesis
      integer, parameter :: jphotonh4e = 16				 ! photoinhibition of nitrification
      integer, parameter :: jsubf = 17              ! penalty for substrate uptake to account for extra N need
      integer, parameter :: jbasal = 18				       ! basal metabolism

      integer, save :: sbio                 ! begin live biology
      integer, save :: dtb                  ! bio timestep
      integer, save :: irein 		 	         ! how often to reinitialize
      integer, save :: itss			           ! index for sediment
      integer, save :: idin			           ! index for DIN
      integer, save :: inh4			           ! index for Ammonium
      integer, save :: idetrs			         ! index for small detritus
      integer, save :: idetrl			         ! index for large detritus
      integer, save :: sorg			           ! index for start of dissolved organic components
      integer, save :: eorg			           ! index for end of dissolved organic components
      integer, save :: in2			             ! index for nitrogen fixation
      integer, save :: nbproc               ! number of "processes" with coefficients
      integer, save :: bpi(mxproc,mxtrcr)   ! array of binary on offs for each bio process
      integer, save :: itt(30)              !array indices for outputting transcription
      integer, save :: jtt(30)                 !array indices for outputting transcription
      integer, save :: ntt                     !array indices for outputting transcription
      integer, save :: nbb                     !array indices for outputting transcription
      integer, save :: idiag                    !array indices for outputting transcription

c - model params
      integer i
      integer j
      integer k
      integer l
      integer inb
      integer iter
      integer ktr
      integer ktra
      integer io
      integer jo
      integer itbio
      integer iyear
      integer iday
      integer ihour

      real bm(mxtrcr)
      real bn(mxtrcr)
      real bi(mxtrcr)         ! temp bio arrays
      real, save :: growth(mxtrcr,mxsubs+3)        ! temp substrate sums
      real sink(kdm+1),bf(kdm+1),dps(kdm+1),rf,bsl  ! arrays for sinking
      real bsum(idm,jdm,kdm) 				                ! sum of all bio

c - light vars
      real radfl
      real swfl
      real rior
      real riob
      real z1
      real z2 ! radiation, shortwave, red and blue light, depths of top and bottom of layer
      real rlt(kdm)
      real ll(mxtrcr)	!avg light in layer, light limitation
      real, save :: pkr
      real, save :: pkb
      real, save :: kr
      real, save :: kb !red and blue light params
      real, save :: cdkr
      real, save :: sedkr
      real, save :: sedkb !cdom K, sediment K
      real, save :: cdkb
      real, save :: a_CDM
      real, save :: a_sed   ! absorption for sed and cdom
      real, save :: airtinc
      real, save :: precipinc ! increment for air t and precip

c - misc
      logical llw
      logical errcon1
      logical ldiag
      logical lbioloc !output diagnostics
      real t1
      real t2
      real t3
      real t4
      real t5
      real t6
      real t7
      real t8
      real t(mxgene)
      real sum     !placeholder terms
      real ao
      character*48 bionm

c--------------------------------------------------------------------
      if     (ibio.lt.0) then !initialize only
c
c ---   read from tracer_NN.input:
c ---   'biotyp' = type (90X=std.bio,X=3,4,9) must be 100
c
c        write(6,*)'entered trcupd_100 - Initializing'
        i = -ibio
        dtb = delt1

        call blkini(k, 'biotyp')

        if (k.ne.100) then
          if (mnproc.eq.1) then
            write(lp,'(/ a /)') 'error - biotyp must be 100'
            call flush(lp)
          endif ! 1st tile

          call xcstop('(trcini)')

          stop '(trcini)'
        endif ! k.ne.100

        call blkini(addtracr , 'addtra')
        call blkini(ntt ,      'ntt   ')

        do io = 1,ntt
          call blkini(itt(io) ,  'itt   ')
          call blkini(jtt(io) ,  'jtt   ')
        enddo !io

        if (mnproc.eq.1) then
          write(lp,*) 'dtime ', dtime
          call forday(dtime,yrflag, iyear,iday,ihour)
          nbb = ntracr
          bionm = "bioloc.out"
          open(unit=uoff+109,file=trim(flnminp)//bionm,status='new')
        endif !mnproc

  ! - read
        ichl   = ntracr+1       ! place for chlorophyl
        icdm   = ntracr+2       ! place for cdom estimate
        ilight = ntracr+3       ! place for light
        iauto  = ntracr+4       ! growth due to photosynthesis
        ihet   = ntracr+5		! growth due to heterotrophy
        ino3   = ntracr+6       ! growth due to nitrate
        isflxb = ntracr+7		! sinking flux
        ibug   = ntracr+7		! place to write transcription

        call blkini(itss     , 'itss  ')
        call blkini(idin     , 'idin  ')
        call blkini(inh4     , 'inh4  ')
        call blkini(idetrs   , 'idetrs')
        call blkini(idetrl   , 'idetrl')
        call blkini(sorg     , 'sorg  ')
        call blkini(eorg     , 'eorg  ')
        call blkini(in2      , 'in2   ')
        call blkini(imot     , 'imot  ')
        call blkini(irhod    , 'irhod ')
        call blkini(sbio     , 'sbio  ')

        inh4e = sbio                !place for nitrification

        call blkini(nbproc   , 'nproc ')
        call blkinr(chldel   , 'chldel','(a6," =",f10.4," ")')

        chldel = 1.0 / (chldel*dtb)

        call blkinr(a_CDM    , 'a_CDM ','(a6," =",f10.4," ")')
        call blkinr(a_sed    , 'a_CDM ','(a6," =",f10.4," ")')
        call blkinr(rivtss   , 'rivtss','(a6," =",f10.4," ")')
        call blkinr(rivdin   , 'rivDIN','(a6," =",f10.4," ")')
        call blkinr(rivdetrs , 'rivdts','(a6," =",f10.4," ")')
        call blkinr(rivdetrl , 'rivdtl','(a6," =",f10.4," ")')

        do io = sorg,eorg
          call blkinr(ao      , 'rivdon','(a6," =",f10.4," ")')
          rivdon(io) = ao
        enddo !io

        do io = 1,sbio
          call blkinr(ao    , 'alpha ','(a6," =",f10.4," ")')
          alpha(io) = ao
        enddo !io

        do io = sorg,eorg
          call blkinr(ao     , 'alphao','(a6," =",f10.4," ")')
          alphao(io) = ao
        enddo !io

        call blkinr(alphanh4 , 'alpha4','(a6," =",f10.4," ")')

        do io = 1,sbio
          call blkinr(ao     , 'aeo   ','(a6," =",f10.4," ")')
          ae(io) = ao
          aei(io) = 1.-ao
        enddo !io

        do io = 1,sbio
          ! in2 is 1, inh4e is sbio
          call blkinr(ao     , 'resp  ','(a6," =",f10.4," ")')
          respi(io) = ao
          bioid(io) = io
        enddo !io

        do io = 1,sbio
          call blkinr(ao     , 'ga    ','(a6," =",f10.4," ")')
          ga(io) = ao
          gai(io) = 1.0-ao
        enddo !io

        call blkinr(betab     , 'betab ','(a6," =",f10.4," ")')
        call blkinr(betap     , 'betap ','(a6," =",f10.4," ")')

        betai = 1.0 - betap - betab

        do io = sorg,eorg
          call blkinr(ao     , 'betarl','(a6," =",f10.4," ")')
          betarl(io) = ao
        enddo !io

        call blkinr(bigdet   , 'bigdet','(a6," =",f10.4," ")')
        call blkinr(mxldet   , 'mxldet','(a6," =",f10.4," ")')
        call blkinr(bsl       , 'betasl','(a6," =",f10.4," ")')

        do io = sorg,eorg
          call blkinr(ao     , 'cdom  ','(a6," =",f10.4," ")')
          domtocdom(io) = ao
        enddo !io

        do io = 1,sbio-1
          call blkinr(ao     , 'sink  ','(a6," =",f10.4," ")')
          bp(jsink,io)=ao/dsec
        enddo !io

        call blkinr(aggreg   , 'aggreg','(a6," =",f10.4," ")')

        aggreg = aggreg/dsec

        call blkinr(gmax     , 'gmax  ','(a6," =",f10.4," ")')
        call blkini(irein    , 'irein ') !timesteps to reinitialize
        call blkini(idiag    , 'idiag ') !timesteps to reinitialize
        call blkinr(airtinc     , 'airt  ','(a6," =",f10.4," ")')
        call blkinr(precipinc     ,    'pinc  ','(a6," =",f10.4," ")')

        airtinc = airtinc / (100.d0*365.d0*86400.d0)
        precipinc = precipinc / (100.d0*365.d0*86400.d0)

        close(unit=uoff+99)

        open(unit=uoff+99,file=trim(flnminp)//'bio100.input')

        call blkini(ngene     , 'ngene ') !gene transcript type

        do io = 2,sbio-1
          betasl(io)=1.
        enddo !io

        betasl(idetrl)=.5

        do inb = sbio,ntracr
          do io = 1,ngene
            call blkini(l     , 'gene  ') !gene transcript type
            gene(io,inb)=l
          enddo !io
          do io = 1,nbproc
            call blkinr(ao    ,'bproc ','(a6," =",e14.5," ")')
            bp(io,inb)=ao
            call blkini(l     ,'bproci')
            bpi(io,inb)=l
          enddo   !io
          nbb = max(nbb,int(bp(1,inb)))
  !        betasl(inb)=max(.3,((bigdet*bp(jsize,inb)+(.9-bigdet))))
          betasl(inb) = 1.
          bioid(inb)=bp(1,inb)
          if(bp(jsize,inb).gt.bigdet) then
            betasl(inb)=bsl
          endif
        enddo	!inb
        close(unit=uoff+99)
        return
      endif !ibio.lt.0

! ============================================
! =============== start main code ============
! ============================================

      ! - vic add false climate forcing increment at each time step
      teinc = teinc + airtinc
      pinc = pinc + precipinc
      margin = 0  ! no horizontal derivatives
	    itbio = 18
      ldiag = .false.
      lbioloc = .false.
      if (mod(nstep,idiag) .lt. 1) then
        ldiag = .true.
      endif

      ! - get biomass at each location
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            do k = 1,kk
              bsum(i,j,k)=0.
              do inb = sbio, ntracr
                bsum(i,j,k) = bsum(i,j,k)+tracer(i,j,k,n,inb)
              enddo ! inb
            enddo !k
		  enddo !i
		enddo !l
      enddo !j

      ! - test for bio reinitialization
      
      if (mod(nstep,irein) .lt. 1 .and. gmax .gt. 0 ) then
      ! - fast forward to start of new bugs...
        open(unit=uoff+99,file=trim(flnminp)//'bio100.rein')
        call blkini(ngene     , 'ngene ') !gene transcript type

        do inb = sbio,nbb
          do io = 1,ngene
            call blkini(l     , 'gene  ') !gene transcript type
          enddo !io
          do io = 1,nbproc
            call blkinr(ao    ,'bproc ','(a6," =",e14.5," ")')
            call blkini(l     ,'bproci')
          enddo   !io
        enddo !inb

      ! - Find tracers that are a small fraction gmax of the community
        do inb = sbio, ntracr
!         for organism inb
          t1=0
          do j=1-margin,jj+margin
            do l=1,isp(j)
              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                do k = 1,kk
!                 for each layer
                  if (bsum(i,j,k).gt.1.e-8) then
!                   if there is real biomass, then
                    t2 = tracer(i,j,k,n,inb)/bsum(i,j,k)
                    t1 = max(t2,t1)
                    if (t1 .eq. t2) then
!                     save index of maximum contribution of organism inb
                      t4 = i
                      t5 = j
                      t6 = k
                    endif !t1=t2
                  endif ! bsum
                enddo !k
		          enddo !i
		        enddo !l
          enddo !j

          if (t1.le.gmax) then
      ! - reinitialize bio - add new organism
            t2 = bp(1,inb)
            do io = 1,ngene
              call blkini(l     , 'gene  ') !gene transcript type
              gene(io,inb)=l
            enddo !io
            do io = 1,nbproc
              call blkinr(ao    ,'bproc ','(a6," =",e14.5," ")')
              bp(io,inb)=ao
              call blkini(l     ,'bproci')
              bpi(io,inb)=l
            enddo   !io
      ! - Adjust biomass for new biomass addition
            do j=1-margin,jj+margin
              do l=1,isp(j)
                do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                  do k = 1,kk
                    t3 = tracer(i,j,k,n,inb)
                    tracer(i,j,k,n,inb) = (bsum(i,j,k)*1.5*gmax)
                    tracer(i,j,k,m,inb) = tracer(i,j,k,n,inb)
!                   t3 = fractional change to biomass
                    t3 = (tracer(i,j,k,n,inb)-t3)/bsum(i,j,k)
! - loop through inb and take away mass from other groups proportional to their biomass here...
!  ie. we know tracer(inb), and we know bsum, so, each bug/bsum times traacer(inb) would conserve mass.
                    do jo = sbio,ntracr
                       tracer(i,j,k,n,jo)=tracer(i,j,k,n,jo)-
     &                    (t3*tracer(i,j,k,n,jo))
                    enddo !jo
                  enddo !k
		            enddo !i
		          enddo !l
            enddo !j

            betasl(inb)=1.
            if(bp(jsize,inb).gt.bigdet) then
              betasl(inb)=bsl
            endif
            nbb=nbb+1
            bioid(inb) = bp(1,inb)

            if (mnproc .eq. 1 ) then
              write(lp,*) nstep,' swapped old bug',int(t2),
     &        ' for',int(bp(1,inb)),t4,t5,t6
!              write(lp,*)'genes',(gene(io,inb), io=1,ngene)
            endif   !mnproc

          endif !t1<=gmax
        enddo ! inb

        close(unit=uoff+99)
        
! write over bio100.input file
        if (mnproc.eq.1) then
          open(unit=uoff+108,file=trim(flnminp)//'bio100.input',
     &       status='unknown')
          write(uoff+108,'(i11,a23)') ngene, " 'ngene ' = gene params"
          do inb = sbio,ntracr
            do io = 1,ngene
              write(uoff+108,'(i11,a,i4)') gene(io,inb),
     &           " 'gene ' = gene params for ",inb
            enddo !io
            do io = 1,nbproc
              write(uoff+108,'(e14.5,a33,i5,a5,f5.0)') bp(io,inb),
     &        " 'bproc ' = bio params  for slot ",inb,
     &        " bug ",bp(1,inb)
              write(uoff+108,'(i11,a23,i5)') bpi(io,inb),
     &        " 'bproci' = bio onoff  for bug ",inb
            enddo   !io
          enddo   !inb
          close(uoff+108)
        endif ! mnproc
      endif ! reinitialize mod nstep,rein

!      write(lp,*)'Entered trcupd_100 - Doing Biology'
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            if (jj .gt. 1) then
              radfl = radflx(i,j,l0)*w0 + radflx(i,j,l1)*w1
     &           + radflx(i,j,l2)*w2 + radflx(i,j,l3)*w3

              if (flxflg.ne.3) then
                radfl = sswflx(i,j)
              endif  !flxflg

              rior = par*radfl*rs2           !Red fraction
              riob = par*radfl*(1.0-rs2)      !Blue fraction
              z1  = 0

              do k=1,kk
                if (riob.gt.0.001) then
                    z2 = dp(i,j,k,n)/onem
                    pkr = 0.
                    pkb =  wk    * ntochl*max(0.,tracer(i,j,k,n,ichl))
                    cdkr = a_CDM * max(0.,tracer(i,j,k,n,icdm))
                    cdkb = a_CDM * max(0.,tracer(i,j,k,n,icdm))
                    sedkr = a_sed * max(0.,tracer(i,j,k,n,itss))
                    sedkb = a_sed * max(0.,tracer(i,j,k,n,itss))
                    kr = wkr2 + pkr + cdkr + sedkr   !Calc. Kt (water + phyto + CDOM) for red.
                    kb = wkb + pkb + cdkb + sedkb   !Calc. Kt (water + phyto + CDOM) for blue.
              rlt(k) = max(0.,-rior*(exp(-kr*z2) - 1)
    &                      /(kr*(z2))
    &                          -riob*(exp(-kb*z2) - 1)
    &                      /(kb*(z2)))
                      rior = rior*exp(-kr*(z2))
                      riob = riob*exp(-kb*(z2))
                  else
                      rlt(k)=0.
                  endif  !riob.gt.0.001
              enddo  !k

! apply nutrients to river inflows.
! rivers are m/sec
              rf = max(0.,abs(rivers(i,j,lc0)*wc0+rivers(i,j,lc1)*wc1
     &        + rivers(i,j,lc2)*wc2+rivers(i,j,lc3)*wc3)*
     &        (dtb*onem / dp(i,j,1,n)))

              tracer(i,j,1,n,idin) = (tracer(i,j,1,n,idin) +
     &        rf*rivdin) / (1.0+rf)

              tracer(i,j,1,n,itss) = (tracer(i,j,1,n,itss) +
     &        rf * rivtss) / (1.0+rf)

              do inb=sorg,eorg
                tracer(i,j,1,n,inb) = (tracer(i,j,1,n,inb) +
     &            rf*rivdon(inb))/(1.0+rf)
              enddo !inb

              tracer(i,j,1,n,idetrs) = (tracer(i,j,1,n,idetrs) +
     &        rf*rivdetrs) / (1.0+rf)
              tracer(i,j,1,n,idetrl) = (tracer(i,j,1,n,idetrl) +
     &        rf*rivdetrl) / (1.0+rf)
            else
              rlt(1) = tracer(i,j,1,n,ilight)
			      endif  !jj>1
! Set up bm_* as the "new" timestep
            do k=1,kk
              if(dp(i,j,k,n) .gt. 0.1*onem) then
                do inb=1,ntracr
                  bm(inb)=max(0.,tracer(i,j,k,n,inb))
                  bn(inb)=bm(inb)
                enddo !inb

! -- calculate growth penalty due to energy availability as sum of all energy sources
! light
                do inb = sbio,ntracr
! -- with photoinhibition
      			      ll(inb)=bpi(jphoto,inb)*
           &            (1.-exp(-rlt(k)/max(.001,bp(jphoto,inb)   )))*
           &            (   exp(-rlt(k)/max(.001,bp(jphoto,inb)**2)))
! -- without photoinhibition
!                ll(inb)=bpi(jphoto,inb)*
!     &            (1.-exp(-rlt(k)/max(.001,bp(jphoto,inb)   )))
!			  if (mnproc.eq.1 .and.itest.eq.i.and.jtest.eq.j
!     &           .and. mod(nstep,idiag) .lt. 1.) then
!       write(lp,'(a,i10,2i4,4f12.5)') 'lightlim',nstep,k,inb,
!     & bp(jphoto,inb), rlt(k), ll(inb),  bpi(jphoto,inb)*
!     &            (1.-exp(-rlt(k)/max(.001,bp(jphoto,inb)   )))
!              endif
			             enddo !inb
!                 if (mnproc.eq.1 .and.itest.eq.i.and.jtest.eq.j
!     &          .and. k.eq.1. .and. mod(nstep,idiag) .lt. 1.) then
!               sum=0.
!               do inb = 2,ntracr
!                 sum = sum+bn(inb)
!               enddo
!               write(lp,'(a,e15.8)') "before iter sum",sum
!               endif
                do iter = 1,itermax
                  do inb=1,ntracr
                    bi(inb)=bm(inb)
                  enddo !inb
! -- calculate terms in bio equations before applying them
                  do inb=sbio,ntracr
! -- din uptake
                    t(jdin) = bpi(jdin,inb)*respi(idin)
! -- ammonium uptake
                    t(jnh4) = bpi(jnh4,inb)*respi(inh4)
! -- nitrogen fixation
                    t(jn2)  = bpi(jn2,inb) *respi(in2)
! -- nitrification
! -- t1 = 1 for no light inhibition, = 0 until light gets to 1 then 1
! -- for light inhibition
                    t1 = max(0.,
     & (1.-(bpi(jphotonh4e,inb)*bp(jphotonh4e,inb)*rlt(k))))
                    t(jnh4e) = bpi(jnh4e,inb) * respi(inh4e) *  t1
                    jo = jorgs
! -- organic substrate uptake
                    do io = sorg,eorg
                      t(jo) = bpi(jo,inb)*respi(io)
                      jo = jo+1
                    enddo !io
! -- particulate uptake
                    t(jdetrs) = bpi(jdetrs,inb)*respi(idetrs)
                    t(jdetrl) = bpi(jdetrl,inb)*respi(idetrl)
! -- calculate growth
                    t7 = bp(jmaxgrow,inb)*bn(inb)
                    growth(inb,idin) = t7*ll(inb)
     & *t(jdin)/(bn(idin)+(bp(jsubf,inb)*bp(jdin,inb)))
                    t1 = growth(inb,idin)*bn(idin)
!			  if (mnproc.eq.1 .and.itest.eq.i.and.jtest.eq.j
!     &          .and. k.eq.1. .and. mod(nstep,idiag) .lt. 1.) then
!                write(lp,*)'din',t1
!              endif
                    growth(inb,inh4) = t7*ll(inb)
     & *t(jnh4)/(bn(inh4)+(bp(jsubf,inb)*bp(jnh4,inb)))
				            t1 = max(t1,growth(inb,inh4)*bn(inh4))
!			  if (mnproc.eq.1 .and.itest.eq.i.and.jtest.eq.j
!     &          .and. k.eq.1. .and. mod(nstep,idiag) .lt. 1.) then
!                write(lp,*)'nh4',growth(inb,inh4)*bn(inh4),t1
!              endif
                    growth(inb,in2) = t7*ll(inb)*t(jn2)
                    t1 = max(t1,growth(inb,in2))
!			  if (mnproc.eq.1 .and.itest.eq.i.and.jtest.eq.j
!     &          .and. k.eq.1. .and. mod(nstep,idiag) .lt. 1.) then
!                write(lp,*)'n2',growth(inb,in2),t1
!              endif
                    growth(inb,inh4e)=t7*t(jnh4e)/
     & (bn(inh4)+(bp(jsubf,inb)*bp(jnh4e,inb)))
				            t1 = max(t1,growth(inb,inh4e)*bn(inh4))
!			  if (mnproc.eq.1 .and.itest.eq.i.and.jtest.eq.j
!     &          .and. k.eq.1. .and. mod(nstep,idiag) .lt. 1.) then
!                write(lp,*)'nh4e',growth(inb,inh4e)*bn(inh4),t1
!              endif
				            jo = jorgs
                    do io=sorg,eorg
                      growth(inb,io) = t7*t(jo)/
     & (bn(io)+(bp(jsubf,inb)*bp(jo,inb)))
	                    t1 = max(t1,growth(inb,io)*bn(io))
!			  if (mnproc.eq.1 .and.itest.eq.i.and.jtest.eq.j
!     &          .and. k.eq.1. .and. mod(nstep,idiag) .lt. 1.) then
!                write(lp,*)'org',io,growth(inb,io)*bn(io),t1
!              endif
                      jo=jo+1
                    enddo !io
                    growth(inb,idetrs) = t7*t(jdetrs)/
     & (bn(idetrs)+(bp(jsubf,inb)*bp(jdetrs,inb)))
				            t1 = max(t1,growth(inb,idetrs)*bn(idetrs))
!			  if (mnproc.eq.1 .and.itest.eq.i.and.jtest.eq.j
!     &          .and. k.eq.1. .and. mod(nstep,idiag) .lt. 1.) then
!                write(lp,*)'dets',growth(inb,idetrs)*bn(idetrs),t1
!              endif
                    growth(inb,idetrl) = t7*t(jdetrl)/
     & (bn(idetrl)+(bp(jsubf,inb)*bp(jdetrl,inb)))
				            t1 = max(t1,growth(inb,idetrl)*bn(idetrl))
!			  if (mnproc.eq.1 .and.itest.eq.i.and.jtest.eq.j
!     &          .and. k.eq.1. .and. mod(nstep,idiag) .lt. 1.) then
!                write(lp,*)'detl',growth(inb,idetrl)*bn(idetrl),t1
!              endif
!-- Set growth for each organism to the growth term which yields maximum growth
!-- Could also keep all terms to express all the genes all the time
!			  if (mnproc.eq.1 .and.itest.eq.i.and.jtest.eq.j
!     &          .and. k.eq.1. .and. mod(nstep,idiag) .lt. 1. ) then
!       write(lp,'(a,i4,9e10.3)') 'growth',inb,
!     & growth(inb,idin)*bn(idin)*86400.,
!     & growth(inb,inh4)*bn(inh4)*86400.,
!     & growth(inb,in2)*86400.,
!     & growth(inb,inh4e)*bn(inh4)*86400.,
!     & growth(inb,sorg)*bn(sorg)*86400.,
!     & growth(inb,eorg)*bn(eorg)*86400.,
!     & growth(inb,idetrs)*bn(idetrs)*86400.,
!     & growth(inb,idetrl)*bn(idetrl)*86400.,
!     & t1*86400.
!       write(lp,'(a,9e10.3)') 'sub    ',bn(inb),
!     & bn(idin),
!     & bn(inh4),
!     & 0.,
!     & bn(inh4),
!     & bn(sorg),
!     & bn(eorg),
!     & bn(idetrs),
!     & bn(idetrl)
!              endif
                    if (growth(inb,idin)*bn(idin).lt.t1) then
                      growth(inb,idin) = 0.
                    endif
                    if (growth(inb,inh4)*bn(inh4).lt.t1) then
                      growth(inb,inh4) = 0.
                    endif
                    if (growth(inb,in2).lt.t1) then
                      growth(inb,in2) = 0.
                    endif
                    if (growth(inb,inh4e)*bn(inh4).lt.t1) then
                      growth(inb,inh4e) = 0.
                    endif
                    do io=sorg,eorg
                      if (growth(inb,io)*bn(io).lt.t1) then
                        growth(inb,io) = 0.
                      endif
                    enddo !io
                    if (growth(inb,idetrs)*bn(idetrs).lt.t1) then
                      growth(inb,idetrs) = 0.
                    endif
                    if (growth(inb,idetrl)*bn(idetrl).lt.t1) then
                      growth(inb,idetrl) = 0.
                    endif
                  enddo !inb

! -- bio uptake of DIN
			            t7 = 0.
                  do inb=sbio,ntracr
			              t7 = t7 + growth(inb,idin)*dtb
			            enddo !inb
			            bn(idin)=bi(idin)/(1.+t7)
                  bi(idin)=bn(idin)
                  t2 = bn(idin)

                  do inb=sbio,ntracr
			              t1 = growth(inb,idin)*t2*dtb
                    growth(inb,idin) = (1.-ae(idin)-alpha(idin))*t1
				            bn(inb) = bi(inb)+growth(inb,idin)
                    bi(inb) = bn(inb)

				            do io=sorg,eorg
                      bn(io)=bi(io)+(1.-ga(idin))*alphao(io)*alpha(idin)*t1
                      bi(io)=bn(io)
                    enddo !io

                    bn(inh4)=bi(inh4)+ga(idin)*alpha(idin)*t1
                    bi(inh4)=bn(inh4)
                    bn(idetrs)=bi(idetrs)+ae(idin)*betasl(idin)*t1
                    bi(idetrs)=bn(idetrs)
                    bn(idetrl)=bi(idetrl)+ae(idin)*(1.-betasl(idin))*t1
                    bi(idetrl)=bn(idetrl)
			            enddo !inb

! -- bio uptake of NH4 for nuts
			            t7 = 0.

                  do inb=sbio,ntracr
			              t7 = t7+growth(inb,inh4)*dtb
                  enddo

                  bn(inh4)=bi(inh4)/(1.+t7)
                  bi(inh4)=bn(inh4)
                  t2 = bn(inh4)

                  do inb=sbio,ntracr
                    t1 = growth(inb,inh4)*t2*dtb
                    growth(inb,inh4) = (1.-ae(inh4)-alpha(inh4))*t1
				            bn(inb) =bi(inb)+growth(inb,inh4)
                    bi(inb) =bn(inb)

				            do io=sorg,eorg
                      bn(io)=bi(io)+(1.-ga(inh4))*alphao(io)*alpha(inh4)*t1
                      bi(io)=bn(io)
                    enddo ! io

                    bn(inh4)=bi(inh4)+ga(inh4)*alpha(inh4)*t1
                    bi(inh4)=bn(inh4)
                    bn(idetrs)=bi(idetrs)+ae(inh4)*betasl(inh4)*t1
                    bi(idetrs)=bn(idetrs)
                    bn(idetrl)=bi(idetrl)+ae(inh4)*(1.-betasl(inh4))*t1
                    bi(idetrl)=bn(idetrl)
			            enddo !inb

! -- bio nitrogen fixation
                  t7 = 0.
                  do inb=sbio,ntracr
			              t7 = t7+growth(inb,in2)*dtb
                  enddo

! increment N2 gas here :)
                  do inb=sbio,ntracr
! limit growth by N2 gas availability here...
                    t1=growth(inb,in2)*dtb
                    growth(inb,in2) = (1.-ae(in2)-alpha(in2))*t1
                    bn(inb) =bi(inb)+ growth(inb,in2)
                    bi(inb) =bn(inb)
                    do io=sorg,eorg
                      bn(io)=bi(io)+(1.-ga(in2))*(alphao(io))*alpha(in2)*t1
                      bi(io)=bn(io)
                    enddo !io
                    bn(inh4)=bi(inh4)+ga(in2)*alpha(in2)*t1
                    bi(inh4)=bn(inh4)
                    bn(idetrs)=bi(idetrs)+(ae(in2))*betasl(in2)*t1
                    bi(idetrs)=bn(idetrs)
                    bn(idetrl)=bi(idetrl)+(ae(in2))*(1.-betasl(in2))*t1
                    bi(idetrl)=bn(idetrl)
                  enddo !inb

! Nitrification
!			   if (mnproc.eq.1 .and.itest.eq.i.and.jtest.eq.j
!     &          .and. k.eq.1. .and. mod(nstep,idiag) .lt. 1.) then
!       write(lp,'(a,3e10.3)') 'before nh4e',bn(inh4),bn(idin),
!     &     bn(eorg)
!               endif
			             t7 = 0.
!			  t6 = 0.
!			  t2 = 0.
                   do inb=sbio,ntracr
			               t7 = t7+growth(inb,inh4e)*dtb
                   enddo !inb

                   bn(inh4)=bi(inh4)/(1.+t7)
                   bi(inh4)=bn(inh4)
                   t2 = bn(inh4)

                   do inb=sbio,ntracr
                     t1 = growth(inb,inh4e)*t2*dtb
                     growth(inb,inh4e) = (1.-ae(inh4e)-alpha(inh4e))*t1
				             bn(inb) =bi(inb)+growth(inb,inh4e)
                     bi(inb) =bn(inb)

                     do io=sorg,eorg
                       bn(io)=bi(io)+(1.-alphanh4)*(1.-ga(inh4e))*
     &                          (alphao(io))*alpha(inh4e)*t1
                       bi(io)=bn(io)
                     enddo !io

                     bn(inh4)=bi(inh4)+ga(inh4e)*alpha(inh4e)*t1
                     bi(inh4)=bn(inh4)
                     bn(idetrs)=bi(idetrs)+ae(inh4e)*betasl(inh4e)*t1
                     bi(idetrs)=bn(idetrs)
                     bn(idetrl)=bi(idetrl)+
     &                          ae(inh4e)*(1.-betasl(inh4e))*t1
                     bi(idetrl)=bn(idetrl)
				             bn(idin)=bi(idin)+alphanh4*alpha(inh4e)*t1
				             bi(idin)=bn(idin)
                   enddo !inb

			             bi(inh4)=bn(inh4)

! -- bio uptake of DON
                   do jo = sorg,eorg
				             t7 = 0.
                     do inb=sbio,ntracr
				               t7 = t7+growth(inb,jo)*dtb
                     enddo !inb
                     bn(jo) = bi(jo)/(1.+t7)
                     bi(jo) = bn(jo)
                     t2 = bn(jo)
                     do inb=sbio,ntracr
                       t1 = growth(inb,jo)*t2*dtb
                       growth(inb,jo) = (1.-ae(jo)-alpha(jo))*t1
				               bn(inb) =bi(inb)+growth(inb,jo)
                       bi(inb) =bn(inb)
                       do io=sorg,eorg
                         bn(io)=bi(io)+(1.-ga(jo))*(alphao(io))*alpha(jo)*t1
                         bi(io)=bn(io)
                       enddo !io
                       bn(inh4)=bi(inh4)+ga(jo)*alpha(jo)*t1
                       bi(inh4)=bn(inh4)
                       bn(idetrs)=bi(idetrs)+ae(jo)*betasl(jo)*t1
                       bi(idetrs)=bn(idetrs)
                       bn(idetrl)=bi(idetrl)+ae(jo)*(1.-betasl(jo))*t1
                       bi(idetrl)=bn(idetrl)
                     enddo !inb
                   enddo !jo

! -- bio uptake of Detr - non-assimilated only goes to small.
			             do jo = idetrs,idetrl
				             t7 = 0.
				             do inb=sbio,ntracr
				               t7 = t7+growth(inb,jo)*dtb
                     enddo !inb
                     bn(jo)=bi(jo)/(1.+t7)
                     bi(jo)=bn(jo)
                     t2 = bn(jo)
                     do inb=sbio,ntracr
                       t1 = growth(inb,jo)*t2*dtb
                       growth(inb,jo) = (1.-ae(jo)-alpha(jo))*t1
				               bn(inb) =bi(inb)+growth(inb,jo)
                       bi(inb) =bn(inb)
                       do io=sorg,eorg
                         bn(io)=bi(io)+(1.-ga(jo))*(alphao(io))*alpha(jo)*t1
                         bi(io)=bn(io)
                       enddo !io
                       bn(inh4)=bi(inh4)+ga(jo)*alpha(jo)*t1
                       bi(inh4)=bn(inh4)
                       bn(idetrs)=bi(idetrs)+ae(jo)*betasl(jo)*t1
                       bi(idetrs)=bn(idetrs)
                       bn(idetrl)=bi(idetrl)+ae(jo)*(1.-betasl(jo))*t1
                       bi(idetrl)=bn(idetrl)
                     enddo !inb
                   enddo! jo

! -- bio death
! -- aggregation
                   t2 = aggreg*bi(idetrs)*bi(idetrs)
                   bn(idetrs) = bi(idetrs)-t2
                   bn(idetrl) = bi(idetrl)+t2
                   bi(idetrs)=bn(idetrs)
                   bi(idetrl)=bn(idetrl)
!-- bio squared death - viruses
			             do inb=sbio,ntracr
				             growth(inb,mxsubs+1)= bp(jdeath,inb)*bi(inb)*bi(inb)
				             bn(inb)=bi(inb)/(1.+dtb*bp(jdeath,inb)*bi(inb))
				             t1=bi(inb)-bn(inb)
				             bi(inb)=bn(inb)
! -- partitioning of dying bio
                     do io=sorg,eorg
                       bn(io)=bi(io)+betarl(io)*betai*t1
                       bi(io)=bn(io)
                     enddo!io
                     bn(inh4)=bi(inh4)+betap*t1
                     bi(inh4)=bn(inh4)
                     bn(idetrs)=bi(idetrs)+betasl(inb)*betab*t1
                     bn(idetrl)=bi(idetrl)+(1.-betasl(inb))*betab*t1
                     bi(idetrl)=bn(idetrl)
                     bi(idetrs)=bn(idetrs)
                   enddo !inb

!-- bio death - basal metabolism
			             do inb=sbio,ntracr
                     t2 = bp(jbasal,inb)*basalmet
!-- reduce motility and basal metabolism when gene imot in on and there is det
                     t1 = gene(imot,inb)*(tracer(i,j,k,n,idetrs)+
     &                tracer(i,j,k,n,idetrl)+bsum(i,j,k))
                     if (t1 .gt. mxldet) then
                       t2 = basalmet
                     endif
                     if (gene(irhod,inb).ge.1 .and. ll(inb).le..001) then
                       t2 = basalmet
                     endif
				             bn(inb)=bi(inb)/(1.+dtb*t2)
				             t1=bi(inb)-bn(inb)
				             bi(inb)=bn(inb)
                     growth(inb,mxsubs+1)=t1
! -- partitioning of dying bio
                     do io=sorg,eorg
                       bn(io)=bi(io)+betarl(io)*betai*t1
                       bi(io)=bn(io)
                     enddo!io
                     bn(inh4)=bi(inh4)+betap*t1
                     bi(inh4)=bn(inh4)
                     bn(idetrs)=bi(idetrs)+betasl(inb)*betab*t1
                     bn(idetrl)=bi(idetrl)+(1.-betasl(inb))*betab*t1
                     bi(idetrl)=bn(idetrl)
                     bi(idetrs)=bn(idetrs)
                   enddo !inb

! -- overall density dependent biodeath - t3 is summed biomass
			             do inb=sbio,ntracr
                     growth(inb,mxsubs+2)=bp(jddeath,inb)*bsum(i,j,k)*
     &                                bi(inb)
                     bn(inb)=bi(inb)/(1.+dtb*bp(jddeath,inb)*bsum(i,j,k))
!				t6=bn(inb)-bi(inb)
                     t1=bi(inb)-bn(inb)
                     bi(inb)=bn(inb)
! -- partitioning of dying bio
                     do io=sorg,eorg
                       bn(io)=bi(io)+betarl(io)*betai*t1
!				  t6=t6+(bn(io)-bi(io))
                       bi(io)=bn(io)
                     enddo!io
                     bn(inh4)=bi(inh4)+betap*t1
                     bi(inh4)=bn(inh4)
                     bn(idetrs)=bi(idetrs)+betasl(inb)*betab*t1
!				t6=t6+(bn(idetrs)-bi(idetrs))
                     bn(idetrl)=bi(idetrl)+(1.-betasl(inb))*betab*t1
!				t6=t6+(bn(idetrl)-bi(idetrl))
                     bi(idetrl)=bn(idetrl)
                     bi(idetrs)=bn(idetrs)
                   enddo !inb
                 enddo !iter

c- forward step eqns - keep any negative values to keep conservation
                 do inb=1,ntracr
                   tracer(i,j,k,n,inb)=bn(inb)
                 enddo!inb
c1c-- calculate chlorophyll by estimating fraction of growth due to photo
c1c--synthesis
c1c -- other diagnstics too.
                 t4 = 0.0
			           t1 = 0.0
			           t2 = 0.0
			           t7 = 0.0
                 t5 = 0.0
                 do inb = sbio,ntracr
                   growth(inb,idin)=growth(inb,idin)/dtb
                   growth(inb,inh4)=growth(inb,inh4)/dtb
                   growth(inb,in2)=growth(inb,in2)/dtb
                   growth(inb,inh4e)=growth(inb,inh4e)/dtb
                   t1 = growth(inb,idin)+growth(inb,inh4)
     &           +growth(inb,in2)
                   t5=t5+growth(inb,idin)
                   t3 = growth(inb,mxsubs+1)+growth(inb,mxsubs+2)
                   if ( t1 .gt. 0) then
                     t7 = t7+bn(inb)
                   endif
                   t4 = t4+t1+growth(inb,inh4e)
                   if (t1+growth(inb,inh4e) .gt. 0) then
                     t1 = t1-growth(inb,mxsubs+2)
                   else
                     t2=t2-growth(inb,mxsubs+2)
                   endif
                   do io=idetrs,eorg
                     growth(inb,io)=growth(inb,io)/dtb
                     t2=t2+growth(inb,io)
                   enddo !io
                 enddo !inb
                 tracer(i,j,k,n,ichl)  = t7
			           tracer(i,j,k,n,iauto) = t4
			           tracer(i,j,k,n,ihet)  = t2
			           tracer(i,j,k,n,ino3)  = t5

! - Use extra space in addtrc
                 tracer(i,j,k,n,ilight) = rlt(k)
!-- calculate total dom pool
                 tracer(i,j,k,n,icdm)=0
                 do io = sorg,eorg
                   tracer(i,j,k,n,icdm)=tracer(i,j,k,n,icdm)+
     &          tracer(i,j,k,n,io)*domtocdom(io)
                 enddo ! io
!------------transcription
!            if (mnproc.eq.1
!     &         .and. mod(nstep,idiag) .lt. 1.) then
                 do inb = 1,ngene
                   trans(inb)=0.
                   t(inb) = 0.
                 enddo
                 if (ldiag) then
                   lbioloc = .false.
                   do io = 1,ntt
                     if (itt(io).eq.i .and.jtt(io).eq.j .and. k.lt.6)then
                       lbioloc=.true.
                     endif
                   enddo !io
                 else
                   lbioloc = .false.
                 endif !ldiag
                 do inb = sbio,ntracr
                   t2=0
                   do jo = 1,sbio
                     t2 = t2+growth(inb,jo)
                   enddo
c              t2 = t2+growth(inb,mxsubs+1)
c-- gene 1-3 pcb-hl,pcb-ll, pbs - assume transcription is increased for low light levels to enhance funneling of light
                   t(1) = gene(1,inb)*(1.-ll(inb))*
     &               (growth(inb,idin)+growth(inb,inh4)+growth(inb,in2))
                   t(2) = gene(2,inb)*(1.-ll(inb))*
     &               (growth(inb,idin)+growth(inb,inh4)+growth(inb,in2))
                   t(3) = gene(3,inb)*(1.-ll(inb))*
     &               (growth(inb,idin)+growth(inb,inh4)+growth(inb,in2))
c - gene 4 = vana  - assume transcription is expressed proportional to growth on that substrate
                   t(4) = gene(4,inb)*growth(inb,sorg)
c     &               +gene(4,inb)*tracer(i,j,k,n,inb)*(2.3e-5)
c gene5 = AMA - organic pools of amino acids, carbs, everything non terrestrial - cellulosic
                   t(5) = gene(5,inb)*growth(inb,eorg)
     &               +gene(5,inb)*tracer(i,j,k,n,inb)*(2.3e-5)
c - gene 6 = detrl  - assume transcription is expressed proportional to growth on that substrate
                   t(6) = gene(6,inb)*growth(inb,idetrl)
     &               +gene(6,inb)*tracer(i,j,k,n,inb)*(2.3e-5)+
     &               +growth(inb,idetrs)
     &               +gene(7,inb)*tracer(i,j,k,n,inb)*(2.3e-5)
c gene7 = rhod - low level light
                   t(7) = gene(7,inb)*t2
                   if (ll(inb).lt..001) then
                     t(7)=0
                   endif
c gene 8  pbs-ll
                   t(8) = gene(8,inb)*(1.-ll(inb))*
     &               (growth(inb,idin)+growth(inb,inh4)+growth(inb,in2))
c gene 9  nif
                   t(9) = gene(9,inb)*(growth(inb,in2))
c gene 10 mot-grazing gene
                   t(10) = gene(10,inb)*t2
c gene 11 sil gene shell formation
                   t(11) = gene(11,inb)*t2
c gene 12 nat gene - nitrate uptake
                   t(12) = gene(12,inb)*
     &        ( growth(inb,idin)*(1.1-
     &            (bn(idin)/(bp(jsubf,inb)*bp(jdin,inb)+bn(idin))))
     &         +gene(12,inb)*(5.8e-6)*tracer(i,j,k,n,inb))
c gene 13 amtb ammonium transport - high affinity
                  t(13) = gene(13,inb)*
     &        (growth(inb,inh4)*(1.1-
     &        (bn(inh4)/(bp(jsubf,inb)*bp(jnh4,inb)+bn(inh4))))
     &         +gene(13,inb)*(5.8e-6)*tracer(i,j,k,n,inb))
c gene 14 nirk gene - nitrate uptake - high affinity
                  t(14) = gene(14,inb)*(growth(inb,idin)*
     &        (1.1-(bn(idin)/(bp(jsubf,inb)*bp(jdin,inb)+bn(idin))))
     &         +gene(14,inb)*(5.8e-6)*tracer(i,j,k,n,inb))
c gene 15 amt ammonium transport
                  t(15)=gene(15,inb)*
     &        (growth(inb,inh4)*(1.1-
     &        (bn(inh4)/(bp(jsubf,inb)*bp(jnh4,inb)+bn(inh4))))
     &         +gene(15,inb)*(5.8e-6)*tracer(i,j,k,n,inb))
c gene imot mot gene - for particle attached
                  t3=(tracer(i,j,k,n,idetrs)+tracer(i,j,k,n,idetrl)+
     &            bsum(i,j,k))
                  t4 = basalmet*(bp(jbasal,inb)**3)*86400
                  if (t3 .gt. mxldet) then
                    t4 = basalmet*86400
                  endif
                  t(imot)=gene(imot,inb)*t2*t4
c - gene 17,18 = amoA  - nitrification - set as nitrate return
c transcription = nitrification rate
c transcription matched to the observations
                  t5 = (alphanh4*alpha(inh4e))/(1.-alpha(inh4e))
                  t(17) = gene(17,inb)*(growth(inb,inh4e))*t5
     &                +gene(17,inb)*1.1e-7*tracer(i,j,k,n,inb)
                  t(18) = gene(18,inb)*(growth(inb,inh4e))*t5
     &                +gene(18,inb)*1.1e-7*tracer(i,j,k,n,inb)
c gene 19 cheA/B gene -reduce mortality through toxin
                  t(19) = gene(19,inb)*t2
c gene 20 chisyn gene -reduce sinking
                  t(20) = gene(20,inb)*t2
c - energy acquisition based transcription
c--calculate transcription - hardwired....
                  if (mnproc.eq.1 .and. lbioloc) then
                    write(uoff+109,
     &            '(i10,3i4,f6.0,f7.3,8e10.3,20e10.3)') nstep,i,j,k,
     &            bp(1,inb),saln(i,j,k,n),
     &            (tracer(i,j,k,n,jo),jo=1,sbio-1),tracer(i,j,k,n,inb),
     &            (t(io)*86400.,io=1,ngene)
                    write(uoff+109,
     &            '(i10,3i4,f6.0,f7.3,8e10.3,20i3)') -1*nstep,
     &            i,j,k,bp(1,inb),saln(i,j,k,n),
     &            (tracer(i,j,k,n,jo),jo=1,sbio-1),tracer(i,j,k,n,inb),
     &            (gene(io,inb),io=1,ngene)
                    write(uoff+109,
     &            '(i10,3i4,f6.0,f7.3,8e10.3,20e10.3)') 0,
     &            i,j,k,bp(1,inb),saln(i,j,k,n),
     &            (tracer(i,j,k,n,jo),jo=1,sbio-1),tracer(i,j,k,n,inb),
     &            (t(io)*86400./tracer(i,j,k,n,inb),io=1,ngene)
                    write(uoff+109,
     &            '(i10,3i4,f6.0,f7.3,8e10.3,20e10.3)') 1,
     &            i,j,k,bp(1,inb),saln(i,j,k,n),
     &            (tracer(i,j,k,n,jo),jo=1,sbio-1),tracer(i,j,k,n,inb)
     &           ,gene(1,inb)*(1.-ll(inb))
     &           ,gene(2,inb)*(1.-ll(inb))
     &           ,gene(3,inb)*(1.-ll(inb))
     &           ,gene(4,inb),gene(5,inb),gene(6,inb),gene(7,inb)
     &           ,gene(8,inb)*(1.-ll(inb))
     &           ,gene(9,inb),gene(10,inb),gene(11,inb)
     &           ,gene(12,inb)*(1.1-
     &            (bn(idin)/(bp(jsubf,inb)*bp(jdin,inb)+bn(idin))))
     &           ,gene(13,inb)*(1.1-
     &            (bn(inh4)/(bp(jsubf,inb)*bp(jnh4,inb)+bn(inh4))))
     &           ,gene(14,inb)*(1.1-
     &            (bn(idin)/(bp(jsubf,inb)*bp(jdin,inb)+bn(idin))))
     &           ,gene(15,inb)*(1.1-
     &            (bn(inh4)/(bp(jsubf,inb)*bp(jnh4,inb)+bn(inh4))))
     &           ,gene(imot,inb)*t4
     &           ,gene(17,inb),gene(18,inb),gene(19,inb),gene(20,inb)
                 endif  !transcription output
                 do io = 1,ngene
                   trans(io) = trans(io)+t(io)
                 enddo
               enddo   !inb
               do inb = 1,ngene
                 tracer(i,j,k,n,ibug+inb)=trans(inb)
               enddo
             endif !dp gt 0.1*onem
          enddo !k
        enddo !i
      enddo !l
    enddo !j

c!$OMP END PARALLEL DO
    llw=.false.
    do j=1-margin,jj+margin
      do l=1,isp(j)
        do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
          
          do k = 1,kk
            dps(k) = (max(((p(i,j,k+1)-p(i,j,1))*qonem),.1))/dtb
            tracer(i,j,k,n,isflxb)=0.
          enddo
          
          do inb = 1,ntracr
          
            do k = 1,kk+1
              sink(k)=bp(jsink,inb)*onem*dtb
              bf(k)= tracer(i,j,k,n,inb)
            enddo
            
            errcon1=.false.
            
            if (sink(1).gt.0) then
              call advc1d(i,j,n,inb,sink,llw,errcon1,.false.)
              t2 = 0
              do k = 1,kk
                t1= (bf(k)- tracer(i,j,k,n,inb))*dps(k)+t2
                t2 = t1
                tracer(i,j,k,n,isflxb)=tracer(i,j,k,n,isflxb)+t1
              enddo
            endif
          
          enddo !inb
          
        enddo !i
      enddo !l
      enddo !j
      return
      end subroutine trcupd_100
!====================================================s
! ======== END ROCA =================================
!====================================================
