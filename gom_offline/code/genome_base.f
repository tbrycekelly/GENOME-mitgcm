      subroutine genome_base(myIter)
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
      integer, parameter :: mxgene=25	   !max space alloted for genes
      real chldel               ! number of days to lag chlorophyll
      real rivdin               ! river DIN to input
      real rivtss               ! river tss to input
      real rivdetrs             ! river small detritus to input
      real rivdetrl             ! river large detritus to input
      real rivdon(mxsubs)       ! river organic matter to input
      real alpha(mxsubs)        ! uptake sent to O and IO nuts
      real alphao(mxsubs)       ! nonassimilated to each kind of DOM
      real alphanh4			       ! uptake sent to IO and O nuts in nitrification
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
      real trans(mxgene)		     ! transcription number for each gene
      real aggreg			         ! aggregation rate for small detritus to large
      real gmax			           ! fraction of community biomass below which an organism is removed

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
      integer ino3					      ! growth due to nitrate
      integer isflxb					    ! sinking flux due to n2-fix
      integer ibug				        ! bug ID location
      integer imot					      ! gene for motility/particle attach
      integer irhod					      ! gene for rhodopsin
      integer ngene					      ! number of genes
      integer gene(mxgene,PTRACERS_num)	! gene array

      ! - process places
      integer, parameter :: jbug = 1					      ! bug number
      integer, parameter :: jsink = 2               ! place for sinking rate
      integer, parameter :: jmaxgrow = 3            ! place for max growth
      integer, parameter :: jdin = 4					      ! DIN uptake
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
      integer, parameter :: jphotonh4e = 16				  ! photoinhibition of nitrification
      integer, parameter :: jsubf = 17              ! penalty for substrate uptake to account for extra N need
      integer, parameter :: jbasal = 18				      ! basal metabolism

      integer sbio                ! begin live biology
      integer dtb                 ! bio timestep
      integer irein 		 	         ! how often to reinitialize
      integer itss			           ! index for sediment
      integer idin			           ! index for DIN
      integer inh4			           ! index for Ammonium
      integer idetrs			         ! index for small detritus
      integer idetrl			         ! index for large detritus
      integer sorg			           ! index for start of dissolved organic components
      integer eorg			           ! index for end of dissolved organic components
      integer in2			           ! index for nitrogen fixation
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
      real ll(PTRACERS_num)	  ! light limitation
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
      !!! Section 1 - Initialize
      !!!
      ! Run if first timestep
      if (myIter.eq.1) then
        !call genome_init()
        print*, 'Iter=',myIter,' **Initialization**'
      end if

      !!!
      !!! Section 2 - Reseed if necessary
      !!!
      ! Run every nth time step
      if (mod(myIter,10).eq.0) then
        !call genome_reseed()
        print*, 'Iter=',myIter,' **Reseed**'
      end if


      !!!
      !!! Section 3 - Run standard biology
      !!!
      ! Run for each lat/lon position
      !do j=1-margin,jj+margin
      !  do k=1-margin,kk+margin
      !    ! Water column biology
      !    call genome_biology()
      !    ! sink that watercolumn
      !    call genome_sink()
      !  end do !k
      !end do !j
      do j=1,Nx
        do k=1,Ny
          call genome_biology(j,k)
          call genome_sink(j,k)
        end do !k
      end do !j


      return
      end subroutine genome_base
      !====================================================
      ! ======== END ROCA =================================
      !====================================================
