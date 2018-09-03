      subroutine genome_reseed(m,n, ibio)
      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"

      integer myIter
      integer n
      integer ibio
!
! --- -------------------------------------------------
! --- tracer-specific operations for ROCA
! --- -------------------------------------------------
!
      integer, parameter :: itermax=3    ! number of iterations
      integer, parameter :: mxsubs=8     ! max number of substrates - real number must be less than sbio
      integer, parameter :: mxproc=20    ! more than mxsubs
      integer, parameter :: mxgene=25    !max space alloted for genes
      real chldel               ! number of days to lag chlorophyll
      real rivdin               ! river DIN to input
      real rivtss               ! river tss to input
      real rivdetrs             ! river small detritus to input
      real rivdetrl             ! river large detritus to input
      real rivdon(mxsubs)       ! river organic matter to input
      real alpha(mxsubs)        ! uptake sent to O and IO nuts
      real alphao(mxsubs)       ! nonassimilated to each kind of DOM
      real alphanh4            ! uptake sent to IO and O nuts in nitrification
      real ae(mxsubs)           ! assim efficiency for each substrate
      real aei(mxsubs)          ! 1.- assim efficiency for each substrate
      real respi(mxsubs)        ! growth efficiency for each substrate
      real ga(mxsubs)           ! excretion partitioning to NH4
      real gai(mxsubs)          ! excretion partitioning to DON
      real betab                ! partitioning of dead bio to detritus
      real betap                ! partitioning of dead bio to NH4
      real betai                ! partitioning of dead bio to don
      real bigdet               ! cutoff for bug size to go to big det
      real mxldet               ! max large det conc for sinking fraction
      real betasl(PTRACERS_num)       ! fraction of dead bug going to small det
      real betarl(mxsubs)       ! fraction of dead bug going to each dom sum to 1
      real domtocdom(mxsubs)    ! fraction of dom that is colored
      real bp(mxproc,PTRACERS_num)    ! array of parameters read in for each bio
      real trans(mxgene)         ! transcription number for each gene
      real aggreg              ! aggregation rate for small detritus to large
      real gmax                ! fraction of community biomass below which an organism is removed

      real, parameter :: rs  = 0.6            ! fraction of incident Q that is longwave
      real, parameter :: rs2 = 0.4            ! 580 nm spectral division in pslight
      !real, parameter :: wkr = 1.6666,                           ! red Ks for water - bad!
      real, parameter :: wkr2= 0.35           ! red Ks for water in pslight
      real, parameter :: wkb = 0.03           ! blue Ks for water
      real, parameter :: dsec = 86400.0       ! seconds in one day
      real, parameter :: par = 0.45           ! fraction of incident sw that is par
      real, parameter :: wk = 0.022           ! blue Ks for chl
      real, parameter :: ntochl = 1.59        !nitrogen to chl conversion
      real, parameter :: basalmet = 1.1e-7    ! .01 per day basal metabolism

      integer ichl                 ! place for chlorophyl
      integer icdm                 ! place for cdom estimate
      integer inh4e                ! place for nitrification
      integer ilight               ! place for light
      integer ihet                 ! growth dure to het
      integer iauto                ! growth due to autotrophy
      integer ino3                ! growth due to nitrate
      integer isflxb              ! sinking flux due to n2-fix
      integer ibug                ! bug ID location
      integer imot                ! gene for motility/particle attach
      integer irhod               ! gene for rhodopsin
      integer ngene               ! number of genes
      integer gene(mxgene,PTRACERS_num) ! gene array

      ! - process places
      integer, parameter :: jbug = 1                ! bug number
      integer, parameter :: jsink = 2               ! place for sinking rate
      integer, parameter :: jmaxgrow = 3            ! place for max growth
      integer, parameter :: jdin = 4                ! DIN uptake
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
      integer, parameter :: jphotonh4e = 16         ! photoinhibition of nitrification
      integer, parameter :: jsubf = 17              ! penalty for substrate uptake to account for extra N need
      integer, parameter :: jbasal = 18             ! basal metabolism

      integer sbio                ! begin live biology
      integer dtb                 ! bio timestep
      integer irein                ! how often to reinitialize
      integer itss                 ! index for sediment
      integer idin                 ! index for DIN
      integer inh4                 ! index for Ammonium
      integer idetrs               ! index for small detritus
      integer idetrl               ! index for large detritus
      integer sorg                 ! index for start of dissolved organic components
      integer eorg                 ! index for end of dissolved organic components
      integer in2                ! index for nitrogen fixation
      integer nbproc              ! number of "processes" with coefficients
      integer bpi(mxproc,PTRACERS_num)  ! array of binary on offs for each bio process
      integer itt(30)             !array indices for outputting transcription
      integer jtt(30)             !array indices for outputting transcription
      integer ntt                 !array indices for outputting transcription
      integer nbb                 !array indices for outputting transcription
      integer idiag               !array indices for outputting transcription

      ! - model params
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

      ! temp bio arrays
      real bm(PTRACERS_num)
      real bn(PTRACERS_num)
      real bi(PTRACERS_num)
      real growth(PTRACERS_num,mxsubs+3) ! temp substrate sums

      ! arrays for sinking
      real sink(Nr+1)
      real bf(Nr+1)
      real dps(Nr+1)
      !real rf ! Already defined in MITgcm?
      real bsl

      real bsum(Nx,Ny,Nr) ! sum of all bio

      ! - light vars
      real radfl        ! radiation
      real swfl         ! shortwave
      real rior         ! red
      real riob         ! blue
      real z1           ! depths of top of layer
      real z2           ! depths of bottom of later
      real rlt(Nr)     ! avg light in layer
      real ll(PTRACERS_num)   ! light limitation
      real pkr
      real pkb
      real kr      ! red light param
      real kb      ! blue light params
      real cdkr    ! cdom K
      real sedkr   ! sediment red K
      real sedkb   ! sediment blue K
      real cdkb
      real a_CDM
      real a_sed   ! absorption for sed and cdom

      ! climate change
      real airtinc
      real precipinc  ! increment for air t and precip

      ! - misc
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


      !!!
      !!! Section 2 - Reseed if necessary
      !!!
      ! Count biomass for each organism (bsum), check if it's a small
      ! fraction of total growth everywhere, reseed everywhere if
      ! necessary.

      ! - vic add false climate forcing increment at each time step
      !teinc = teinc + airtinc
      !pinc = pinc + precipinc

      margin = 0  ! no horizontal derivatives
	    itbio = 18
      ldiag = .false.
      lbioloc = .false.
      if (mod(nstep,idiag) .lt. 1) then
        ldiag = .true.
      endif

      ! Section 2.1 - Get biomass at each location
      do j=1,Nx
        do i=1,Ny
          !do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            do k = 1,Nr
              bsum(i,j,k) = 0.0
              do inb = sbio, ntracr
                bsum(i,j,k) = bsum(i,j,k)+pTracer(i,j,k,n,inb)
              enddo ! inb
            enddo !k
		      !enddo !i
		    enddo !l
      enddo !j

      ! Section 2.2 - Test for bio reinitialization
      if (mod(nstep,irein).lt.1 .and. gmax.gt.0) then
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

        ! Section 2.2.1 - Reseeding
        ! Find tracers that are a small fraction gmax of the community
        ! Save index of cell where the maximum relative tracer is
        ! present. Then check if it is less than gmax, reinitialize if
        ! not.
        do inb = sbio, ntracr
          !for organism inb

          ! Section 2.2.1.1 - Max relative contribution
          t1=0
          do j=1-margin,jj+margin
            do l=1,isp(j)
              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                do k = 1,kk
                  !for each layer
                  if (bsum(i,j,k).gt.1.e-8) then
                    !if there is real biomass, then
                    t2 = tracer(i,j,k,n,inb)/bsum(i,j,k)
                    t1 = max(t2,t1)
                    if (t1 .eq. t2) then
                      !save index of maximum contribution of organism
                      t4 = i
                      t5 = j
                      t6 = k
                    endif !t1=t2
                  endif ! bsum
                enddo !k
		          enddo !i
		        enddo !l
          enddo !j

          ! Section 2.2.1.2 - Reseed and adjust
          ! Reseed with new organism if max relative growth is less
          ! than gmax. First sum affected biomass and and reseed at 1.5
          ! times gmax. Rescale all organisms in order to conserve
          ! mass.
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
                    !t3 = fractional change to biomass
                    t3 = (tracer(i,j,k,n,inb)-t3)/bsum(i,j,k)

                    ! Loop through inb and take away mass from other
                    ! groups proportional to their biomass here...
                    ! ie. we know tracer(inb), and we know bsum, so,
                    ! each bug/bsum times tracer(inb) would conserve
                    ! mass.
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
      ! End section 2 - Reseeding


      return
      end subroutine genome_reseed
      !====================================================
      ! ======== END ROCA =================================
      !====================================================
