      subroutine genome_init(m,n, ibio)
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
      !!! Section 1 - Initialization
      !!!
      if (ibio.lt.0) then !initialize only
        i = -ibio
        dtb = delt1

        call blkini(k, 'biotyp')

!        if (k.ne.100) then
!          if (mnproc.eq.1) then
!            write(lp,'(/ a /)') 'error - biotyp must be 100'
!            call flush(lp)
!          endif ! 1st tile

!          call xcstop('(trcini)')
!
!          stop '(trcini)'
!        endif ! k.ne.100

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

        ! - read paramters
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
      !!! End Section 1 - Initialization


      return
      end subroutine genome_init
      !====================================================
      ! ======== END ROCA =================================
      !====================================================
