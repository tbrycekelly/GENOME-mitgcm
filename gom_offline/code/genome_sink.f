      subroutine genome_sink(m,n, ibio)
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
      real mxldet               ! max large det conc for sinking fraction
      real betasl(PTRACERS_num)       ! fraction of dead bug going to small det
      real betarl(mxsubs)       ! fraction of dead bug going to each dom sum to 1
      real domtocdom(mxsubs)    ! fraction of dom that is colored
      real bp(mxproc,PTRACERS_num)    ! array of parameters read in for each bio
      real trans(mxgene)         ! transcription number for each gene
      real aggreg              ! aggregation rate for small detritus to large


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
      !!! Section 4 - Sinking
      !!!
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
      ! End Section 4 - Sinking


      return
      end subroutine genome_sink
      !====================================================
      ! ======== END ROCA =================================
      !====================================================
