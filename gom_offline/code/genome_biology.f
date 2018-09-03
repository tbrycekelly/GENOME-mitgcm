      subroutine genome_biology(i,j)
      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"

      integer i
      integer j
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
      real rlt(Nr)      ! avg light in layer
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
      !!! Section 3 - Biology
      !!!
      do l=1,Nr
        do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))

          !!!
          !!! Section 3.1 -
          !!!
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

                  !Calc. Kt (water + phyto + CDOM) for red and blue
                  kr = wkr2 + pkr + cdkr + sedkr
                  kb = wkb + pkb + cdkb + sedkb

                  rlt(k) = max(0.,-rior*(exp(-kr*z2) - 1) / (kr*(z2))
    &                          -riob*(exp(-kb*z2) - 1) / (kb*(z2)))
                      rior = rior*exp(-kr*(z2))
                      riob = riob*exp(-kb*(z2))
                  else
                      rlt(k)=0.
                  endif  !riob.gt.0.001
              enddo  !k

              ! Apply nutrients to river inflows.
              ! NB. Rivers are m/sec
              rf = max(0.,abs(rivers(i,j,lc0)*wc0+rivers(i,j,lc1)*wc1
     &                    + rivers(i,j,lc2)*wc2+rivers(i,j,lc3)*wc3)*
     &                    (dtb*onem / dp(i,j,1,n)))

              tracer(i,j,1,n,idin) = (tracer(i,j,1,n,idin) +
     &                                rf*rivdin) / (1.0+rf)

              tracer(i,j,1,n,itss) = (tracer(i,j,1,n,itss) +
     &                                rf * rivtss) / (1.0+rf)

              do inb=sorg,eorg
                tracer(i,j,1,n,inb) = (tracer(i,j,1,n,inb) +
     &                                rf*rivdon(inb))/(1.0+rf)
              enddo !inb

              tracer(i,j,1,n,idetrs) = (tracer(i,j,1,n,idetrs) +
     &                                rf*rivdetrs) / (1.0+rf)
              tracer(i,j,1,n,idetrl) = (tracer(i,j,1,n,idetrl) +
     &                                rf*rivdetrl) / (1.0+rf)
            else
              rlt(1) = tracer(i,j,1,n,ilight)
			      endif  !jj>1
			      ! End Basic Light



            !!!
            !!! Section 3.2 - Limitations
            !!!
            ! Set up bm_* as the "new" timestep
            do k=1,kk
              if(dp(i,j,k,n) .gt. 0.1*onem) then
                do inb=1,ntracr
                  bm(inb)=max(0.,tracer(i,j,k,n,inb))
                  bn(inb)=bm(inb)
                enddo !inb

                ! Calculate growth penalty due to energy availability
                ! as sum of all energy sources light
                do inb = sbio,ntracr
                  ! With photoinhibition
      			      ll(inb)=bpi(jphoto,inb)*
           &            (1.-exp(-rlt(k)/max(.001,bp(jphoto,inb)   )))*
           &            (   exp(-rlt(k)/max(.001,bp(jphoto,inb)**2)))
                  ! Without photoinhibition
!                 ll(inb)=bpi(jphoto,inb)*
!     &            (1.-exp(-rlt(k)/max(.001,bp(jphoto,inb)   )))
			          enddo !inb



                !!!
                !!! Section 3.2.1 - Backward Euler
                !!!
                do iter = 1,itermax
                  do inb=1,ntracr
                    bi(inb)=bm(inb)
                  enddo !inb
                  ! Calculate terms in bio equations before applying
                  ! them
                  do inb=sbio,ntracr
                    ! Din uptake
                    t(jdin) = bpi(jdin,inb)*respi(idin)
                    ! Ammonium uptake
                    t(jnh4) = bpi(jnh4,inb)*respi(inh4)
                    ! Nitrogen fixation
                    t(jn2)  = bpi(jn2,inb) *respi(in2)
                    ! Nitrification
                    ! t1 = 1 for no light inhibition, = 0 until light
                    ! gets to 1 then 1 for light inhibition
                    t1 = max(0.,
     & (1.-(bpi(jphotonh4e,inb)*bp(jphotonh4e,inb)*rlt(k))))
                    t(jnh4e) = bpi(jnh4e,inb) * respi(inh4e) *  t1
                    jo = jorgs
                    ! Organic substrate uptake
                    do io = sorg,eorg
                      t(jo) = bpi(jo,inb)*respi(io)
                      jo = jo+1
                    enddo !io
                    ! Particulate uptake
                    t(jdetrs) = bpi(jdetrs,inb)*respi(idetrs)
                    t(jdetrl) = bpi(jdetrl,inb)*respi(idetrl)

                    ! Calculate growth
                    t7 = bp(jmaxgrow,inb)*bn(inb)
                    growth(inb,idin) = t7*ll(inb)*t(jdin)/
     &                        (bn(idin)+(bp(jsubf,inb)*bp(jdin,inb)))

                    t1 = growth(inb,idin)*bn(idin)
                    growth(inb,inh4) = t7*ll(inb)*t(jnh4)
     &                       /(bn(inh4)+(bp(jsubf,inb)*bp(jnh4,inb)))

				            t1 = max(t1,growth(inb,inh4)*bn(inh4))
                    growth(inb,in2) = t7*ll(inb)*t(jn2)

                    t1 = max(t1,growth(inb,in2))
                    growth(inb,inh4e)=t7*t(jnh4e)/(bn(inh4)
     &                              +(bp(jsubf,inb)*bp(jnh4e,inb)))

				            t1 = max(t1,growth(inb,inh4e)*bn(inh4))
				            jo = jorgs

                    do io=sorg,eorg
                      growth(inb,io) = t7*t(jo) / (bn(io)
     &                               +(bp(jsubf,inb)*bp(jo,inb)))
	                    t1 = max(t1,growth(inb,io)*bn(io))
                      jo=jo+1
                    enddo !io
                    growth(inb,idetrs) = t7*t(jdetrs)/(bn(idetrs)
     &                              +(bp(jsubf,inb)*bp(jdetrs,inb)))
				            t1 = max(t1,growth(inb,idetrs)*bn(idetrs))
                    growth(inb,idetrl) = t7*t(jdetrl)/(bn(idetrl)
     &                              +(bp(jsubf,inb)*bp(jdetrl,inb)))
				            t1 = max(t1,growth(inb,idetrl)*bn(idetrl))

                    ! Set growth for each organism to the growth term
                    ! which yields maximum growth
                    ! Could also keep all terms to express all the
                    ! genes all the time

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

                  ! Bio uptake of DIN
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

                  ! Bio uptake of NH4 for nuts
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

                  ! Bio nitrogen fixation
                  t7 = 0.
                  do inb=sbio,ntracr
			              t7 = t7+growth(inb,in2)*dtb
                  enddo

                  ! Increment N2 gas here :)
                  do inb=sbio,ntracr
                    ! Limit growth by N2 gas availability here...
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
			             t7 = 0.
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

                   ! Bio uptake of DON
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

                   ! Bio uptake of Detr - non-assimilated only goes to
                   ! small.
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

                   !!! Bio death
                   ! Aggregation
                   t2 = aggreg*bi(idetrs)*bi(idetrs)
                   bn(idetrs) = bi(idetrs)-t2
                   bn(idetrl) = bi(idetrl)+t2
                   bi(idetrs)=bn(idetrs)
                   bi(idetrl)=bn(idetrl)
                   ! Bio squared death - viruses
			             do inb=sbio,ntracr
				             growth(inb,mxsubs+1)= bp(jdeath,inb)*bi(inb)*bi(inb)
				             bn(inb)=bi(inb)/(1.+dtb*bp(jdeath,inb)*bi(inb))
				             t1=bi(inb)-bn(inb)
				             bi(inb)=bn(inb)
                     ! Partitioning of dying bio
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

                   ! Bio death - basal metabolism
			             do inb=sbio,ntracr
                     t2 = bp(jbasal,inb)*basalmet
                     ! Reduce motility and basal metabolism when gene
                     ! imot in on and there is det.
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
                     ! Partitioning of dying bio
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

                   ! Overall density dependent biodeath - t3 is summed biomass
			             do inb=sbio,ntracr
                     growth(inb,mxsubs+2)=bp(jddeath,inb)*bsum(i,j,k)*
     &                                bi(inb)
                     bn(inb)=bi(inb)/(1.+dtb*bp(jddeath,inb)*bsum(i,j,k))
!				              t6=bn(inb)-bi(inb)
                     t1=bi(inb)-bn(inb)
                     bi(inb)=bn(inb)
                     ! Partitioning of dying bio
                     do io=sorg,eorg
                       bn(io)=bi(io)+betarl(io)*betai*t1
!				                t6=t6+(bn(io)-bi(io))
                       bi(io)=bn(io)
                     enddo!io
                     bn(inh4)=bi(inh4)+betap*t1
                     bi(inh4)=bn(inh4)
                     bn(idetrs)=bi(idetrs)+betasl(inb)*betab*t1
!				              t6=t6+(bn(idetrs)-bi(idetrs))
                     bn(idetrl)=bi(idetrl)+(1.-betasl(inb))*betab*t1
!				              t6=t6+(bn(idetrl)-bi(idetrl))
                     bi(idetrl)=bn(idetrl)
                     bi(idetrs)=bn(idetrs)
                   enddo !inb
                 enddo !iter
                ! End 3.2.1 - Backward Euler


                 !!!
                 !!! Section 3.2.2 - Forward step eqns
                 !!!
                 ! Keep any negative values to keep conservation
                 do inb=1,ntracr
                   tracer(i,j,k,n,inb)=bn(inb)
                 enddo!inb
                 ! Calculate chlorophyll by estimating fraction of
                 ! growth due to photosynthesis, other diagnstics too.
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

                 ! Use extra space in addtrc
                 tracer(i,j,k,n,ilight) = rlt(k)
                 ! Calculate total dom pool
                 tracer(i,j,k,n,icdm)=0
                 do io = sorg,eorg
                   tracer(i,j,k,n,icdm)=tracer(i,j,k,n,icdm)+
     &                              tracer(i,j,k,n,io)*domtocdom(io)
                 enddo ! io



                 !!!
                 !!! Section 3.2.2.1 - Transcription
                 !!!
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
!                   t2 = t2+growth(inb,mxsubs+1)

                   ! Gene 1-3 pcb-hl,pcb-ll, pbs - assume transcription
                   ! is increased for low light levels to enhance
                   ! funneling of light
                   t(1) = gene(1,inb)*(1.-ll(inb))*
     &               (growth(inb,idin)+growth(inb,inh4)+growth(inb,in2))
                   t(2) = gene(2,inb)*(1.-ll(inb))*
     &               (growth(inb,idin)+growth(inb,inh4)+growth(inb,in2))
                   t(3) = gene(3,inb)*(1.-ll(inb))*
     &               (growth(inb,idin)+growth(inb,inh4)+growth(inb,in2))

                   ! Gene 4 = vana
                   ! Assume transcription is expressed proportional
                   ! to growth on that substrate
                   t(4) = gene(4,inb)*growth(inb,sorg)
!     &               +gene(4,inb)*tracer(i,j,k,n,inb)*(2.3e-5)

                   ! Gene5 = AMA
                   ! organic pools of amino acids, carbs, everything
                   ! non terrestrial - cellulosic
                   t(5) = gene(5,inb)*growth(inb,eorg)
     &               +gene(5,inb)*tracer(i,j,k,n,inb)*(2.3e-5)

                   ! Gene 6 = detrl
                   ! assume transcription is expressed proportional to
                   ! growth on that substrate
                   t(6) = gene(6,inb)*growth(inb,idetrl)
     &               +gene(6,inb)*tracer(i,j,k,n,inb)*(2.3e-5)+
     &               +growth(inb,idetrs)
     &               +gene(7,inb)*tracer(i,j,k,n,inb)*(2.3e-5)

                   ! Gene7 = rhod - low level light
                   t(7) = gene(7,inb)*t2
                   if (ll(inb).lt..001) then
                     t(7)=0
                   endif

                   ! Gene 8  pbs-ll
                   t(8) = gene(8,inb)*(1.-ll(inb))*
     &               (growth(inb,idin)+growth(inb,inh4)+growth(inb,in2))
                   ! Gene 9  nif
                   t(9) = gene(9,inb)*(growth(inb,in2))

                   ! Gene 10 mot-grazing gene
                   t(10) = gene(10,inb)*t2

                   ! gene 11 sil gene shell formation
                   t(11) = gene(11,inb)*t2

                   ! gene 12 nat gene - nitrate uptake
                   t(12) = gene(12,inb)*
     &        ( growth(inb,idin)*(1.1-
     &            (bn(idin)/(bp(jsubf,inb)*bp(jdin,inb)+bn(idin))))
     &         +gene(12,inb)*(5.8e-6)*tracer(i,j,k,n,inb))

                  ! Gene 13 amtb ammonium transport - high affinity
                  t(13) = gene(13,inb)*
     &        (growth(inb,inh4)*(1.1-
     &        (bn(inh4)/(bp(jsubf,inb)*bp(jnh4,inb)+bn(inh4))))
     &         +gene(13,inb)*(5.8e-6)*tracer(i,j,k,n,inb))

                  ! Gene 14 nirk gene - nitrate uptake - high affinity
                  t(14) = gene(14,inb)*(growth(inb,idin)*
     &        (1.1-(bn(idin)/(bp(jsubf,inb)*bp(jdin,inb)+bn(idin))))
     &         +gene(14,inb)*(5.8e-6)*tracer(i,j,k,n,inb))

                  ! Gene 15 amt ammonium transport
                  t(15)=gene(15,inb)*
     &        (growth(inb,inh4)*(1.1-
     &        (bn(inh4)/(bp(jsubf,inb)*bp(jnh4,inb)+bn(inh4))))
     &         +gene(15,inb)*(5.8e-6)*tracer(i,j,k,n,inb))

                  ! Gene imot mot gene - for particle attached
                  t3=(tracer(i,j,k,n,idetrs)+tracer(i,j,k,n,idetrl)+
     &            bsum(i,j,k))
                  t4 = basalmet*(bp(jbasal,inb)**3)*86400
                  if (t3 .gt. mxldet) then
                    t4 = basalmet*86400
                  endif
                  t(imot)=gene(imot,inb)*t2*t4

                  ! Gene 17,18 = amoA  - nitrification - set as nitrate return
                  ! transcription = nitrification rate
                  ! transcription matched to the observations
                  t5 = (alphanh4*alpha(inh4e))/(1.-alpha(inh4e))
                  t(17) = gene(17,inb)*(growth(inb,inh4e))*t5
     &                +gene(17,inb)*1.1e-7*tracer(i,j,k,n,inb)
                  t(18) = gene(18,inb)*(growth(inb,inh4e))*t5
     &                +gene(18,inb)*1.1e-7*tracer(i,j,k,n,inb)
                  ! Gene 19 cheA/B gene -reduce mortality through toxin
                  t(19) = gene(19,inb)*t2

                  ! Gene 20 chisyn gene -reduce sinking
                  t(20) = gene(20,inb)*t2
                  ! Energy acquisition based transcription
                  ! Calculate transcription - hardwired....
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

      ! End Section 3 - Biology

!$OMP END PARALLEL DO
      return
      end subroutine genome_biology
      !====================================================
      ! ======== END ROCA =================================
      !====================================================
