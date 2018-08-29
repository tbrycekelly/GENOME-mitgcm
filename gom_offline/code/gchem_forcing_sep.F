C $Header: /u/gcmpack/MITgcm/pkg/gchem/gchem_forcing_sep.F,v 1.41 2015/07/18 21:42:04 jmc Exp $
C $Name: checkpoint65v $

#include "GCHEM_OPTIONS.h"
#ifdef ALLOW_DIC
#include "DIC_OPTIONS.h"
#endif
#ifdef ALLOW_DARWIN
#include "DARWIN_OPTIONS.h"
#endif

CBOP
C !ROUTINE: GCHEM_FORCING_SEP
C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_FORCING_SEP( myTime, myIter, myThid )

C !DESCRIPTION:
C     calls subroutine that will update passive tracers values
C     with a separate timestep. Since GCHEM_FORCING_SEP is now
C     called before DO_FIELDS_BLOCKING_EXCHANGES, the passive
C     tracer values in the halo regions are not up to date and
C     must not be used.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
C#include "KPP.h"
C#include "KPP_PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "GCHEM.h"
#include "GCHEM_FIELDS.h"
#ifdef ALLOW_DIC
#include "DIC_VARS.h"
#endif /* ALLOW_DIC */
#ifdef ALLOW_DARWIN
#include "DARWIN_FLUX.h"
#include "DARWIN_SIZE.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      _RL myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_GCHEM
#ifdef GCHEM_SEPARATE_FORCING

C!LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  bi,bj                :: tile indices
C  k                    :: vertical level

      REAL:: vol_sum
      INTEGER bi,bj,iMin,iMax,jMin,jMax
      INTEGER i, j, k, iter, dd
      PARAMETER( iMin = 1 , iMax = sNx )
      PARAMETER( jMin = 1 , jMax = sNy )
      character(len=46):: part1
      character(len=33):: part2
      character(len=79)::fname_l
C
C *** NEMUOR Variables ***
      REAL:: NO_v1, NO_v2, NH_v1, NH_v2, pass_v1
      REAL:: sp_v1, sp_v2, lp_v1, lp_v2, th_v1
      REAL:: sz_v1, sz_v2, lz_v1, lz_v2, pass_v2
      REAL:: pz_v1, pz_v2, SI_v1, SI_v2, th_v2
      REAL:: DON_v1, DON_v2, PON_v1, PON_v2, OP_v1, OP_v2
      REAL:: delta_NO, delta_NH, delta_SI
      REAL:: delta_PON, delta_DON, delta_OP
      REAL:: delta_sp, delta_lp
      REAL:: delta_sz, delta_lz, delta_pz
      REAL:: NO_lim_sp, NH_lim_sp
      REAL:: NO_lim_lp, NH_lim_lp, SI_lim_lp
      REAL:: l_lim_sp, t_lim_v_sp, t_lim_r_sp, t_lim_m_sp
      REAL:: l_lim_lp, t_lim_v_lp, t_lim_r_lp, t_lim_m_lp
      REAL:: np_frac_sp, np_frac_lp
      REAL:: t_lim_g_sz_sp, g_lim_sz_sp, graz_sz_sp
      REAL:: NS_lim_comp
      REAL:: t_lim_g_lz_sp, g_lim_lz_sp, graz_lz_sp
      REAL:: t_lim_g_lz_lp, g_lim_lz_lp, graz_lz_lp
      REAL:: t_lim_g_pz_lp, g_lim_pz_lp, graz_pz_lp
      REAL:: t_lim_m_sz, t_lim_m_lz, t_lim_m_pz
      REAL:: t_lim_g_lz_sz, g_lim_lz_sz, graz_lz_sz
      REAL:: t_lim_g_pz_sz, g_lim_pz_sz, graz_pz_sz
      REAL:: t_lim_g_pz_lz, g_lim_pz_lz, graz_pz_lz
      REAL:: t_lim_nit, nitr, t_lim_dec_PON_NH
      REAL:: t_lim_dec_DON_NH, t_lim_dec_PON_DON
      REAL:: t_lim_dec_OP_SI, tmp
      REAL:: l_dum, volume1, volume2, sink_frac1, sink_frac2, cell_depth
      REAL,DIMENSION(PTRACERS_numInUse):: var
      _RL temp_PON(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL temp_OP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL light_adj(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C *** END NEMURO Variables ***

C *** Chl:C Model
      REAL:: v_chl_sp, v_chl_lp, sp_chl2c, lp_chl2c
      REAL:: sp_chl, lp_chl

C *** Diagnostics ***
      REAL:: sp_diag, lp_diag, sz_diag, lz_diag, pz_diag
      REAL:: pp_NO, pp_NH
      REAL:: l_lim_sp_diag, l_lim_lp_diag
      REAL:: NO_lim_sp_diag, NO_lim_lp_diag
      REAL:: NH_lim_sp_diag, NH_lim_lp_diag
      REAL:: SI_lim_lp_diag
C *** End Diagnostics ***

#if (defined ALLOW_OBCS) || (defined ALLOW_DIAGNOSTICS)
      INTEGER iTr
#endif
#ifdef ALLOW_DIAGNOSTICS
      CHARACTER*8 diagName
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('GCHEM_FORCING_SEP',myThid)
#endif

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
C--   fill-in tracer diagnostics before any GChem udate
       DO iTr = 1,PTRACERS_numInUse
        diagName = '        '
        WRITE(diagName,'(A5,A2)') 'GC_Tr', PTRACERS_ioLabel(iTr)
        CALL DIAGNOSTICS_FILL( pTracer(1-OLx,1-OLy,1,1,1,iTr), diagName,
     &                         0, Nr, 0, 1, 1, myThid )
       ENDDO
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

ccccccccccccccccccccccccc
c global calculations   c
ccccccccccccccccccccccccc
#ifdef ALLOW_OLD_VIRTUALFLUX
#ifdef ALLOW_DIC
# ifdef ALLOW_AUTODIFF
      IF ( .NOT.useDIC ) STOP 'ABNORMAL END: S/R GCHEM_FORCING_SEP (1)'
# else /* ALLOW_AUTODIFF */
      IF ( useDIC ) THEN
# endif /* ALLOW_AUTODIFF */
c find global surface averages
       gsm_s = 0. _d 0
       gsm_dic = 0. _d 0
       gsm_alk = 0. _d 0
       CALL GCHEM_SURFMEAN(salt,gsm_s,myThid)
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,1), gsm_dic, myThid )
       print*,'mean surface dic', gsm_dic,gsm_s
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,2), gsm_alk, myThid )
# ifndef ALLOW_AUTODIFF
      ENDIF
# endif /* ALLOW_AUTODIFF */
#endif /* ALLOW_DIC */
#ifdef ALLOW_DARWIN
c     IF ( useDARWIN ) THEN
c find global surface averages
       gsm_s = 0. _d 0
       gsm_dic = 0. _d 0
       gsm_alk = 0. _d 0
       CALL GCHEM_SURFMEAN(salt,gsm_s,myThid)
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,iDIC), gsm_dic, myThid )
       print*,'mean surface dic', gsm_dic,gsm_s
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,iALK), gsm_alk, myThid )
c     ENDIF
#endif
ccccccccccccccccccccccccccccccccccccccccccc
#endif /* ALLOW_OLD_VIRTUALFLUX */

#ifdef ALLOW_DARWIN
      IF ( useDARWIN ) THEN
        CALL DARWIN_CONS( myIter, myTime, myThid )
      ENDIF
#endif

ccccccccccccccccccccccccc
c chemical forcing      c
ccccccccccccccccccccccccc
C$taf loop = parallel
      DO bj=myByLo(myThid),myByHi(myThid)
C$taf loop = parallel
       DO bi=myBxLo(myThid),myBxHi(myThid)

ccccccccccccccccccccccccccc DIC cccccccccccccccccccccccccccccccc

#ifdef ALLOW_DIC
# ifdef ALLOW_AUTODIFF
        IF (.NOT.useDIC) STOP 'ABNORMAL END: S/R GCHEM_FORCING_SEP (2)'
# else /* ALLOW_AUTODIFF */
        IF ( useDIC ) THEN
# endif /* ALLOW_AUTODIFF */
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL('DIC_BIOTIC_FORCING',myThid)
#endif
#ifdef ALLOW_FE
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,6),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
#else
#ifdef ALLOW_O2
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
#else
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
#endif
#endif
# ifndef ALLOW_AUTODIFF
        ENDIF
# endif /* ALLOW_AUTODIFF */
#endif /* ALLOW_DIC */
cccccccccccccccccccccccccc END DIC cccccccccccccccccccccccccccccccccc

#ifdef ALLOW_DARWIN
        IF ( useDARWIN ) THEN
#ifdef NUT_SUPPLY
c articficial supply of nutrients
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL('DARWIN_NUT_SUPPLY',myThid)
#endif
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
#endif
ccccccccccccccc
C darwin_forcing operates on bi,bj part only, but needs to get full
C array because of last (iPtr) index
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL('DARWIN_FORCING',myThid)
#endif
          CALL DARWIN_FORCING(  pTracer(1-OLx,1-OLy,1,1,1,1),
     &                          bi, bj, iMin, iMax, jMin, jMax,
     &                          myIter, myTime, myThid )
        ENDIF
#endif /* ALLOW_DARWIN */

#ifdef ALLOW_OBCS
C--   Apply (again) open boundary conditions for each passive tracer
C Note: could skip this 2nd call to OBCS_APPLY if all DIC/DARWIN
C       updates of ptracers were only done in the interior (i.e. with
C       tendency multiplied by maskInC)
        IF ( useOBCS .AND. .NOT.useDIC ) THEN
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL('OBCS_APPLY_PTRACER',myThid)
#endif
          DO iTr = 1,PTRACERS_numInUse
            CALL OBCS_APPLY_PTRACER(
     I                bi, bj, 0, iTr,
     U                pTracer(1-OLx,1-OLy,1,bi,bj,iTr),
     I                myThid )
          ENDDO
        ENDIF
#endif /* ALLOW_OBCS */

       ENDDO
      ENDDO

#ifdef ALLOW_DARWIN
      IF ( useDARWIN ) THEN
         CALL DARWIN_CONS( myIter, myTime, myThid )
#ifdef ALLOW_CARBON
         CALL DIC_ATMOS( 1, myTime, myIter, myThid )
#endif
      ENDIF
#endif /* ALLOW_DARWIN */

#ifdef ALLOW_DIC
# ifdef ALLOW_AUTODIFF
      IF ( .NOT.useDIC ) STOP 'ABNORMAL END: S/R GCHEM_FORCING_SEP (3)'
# else /* ALLOW_AUTODIFF */
      IF ( useDIC ) THEN
# endif /* ALLOW_AUTODIFF */
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('DIC_ATMOS',myThid)
#endif
       CALL DIC_ATMOS( myTime, myIter, myThid )
# ifdef COMPONENT_MODULE
       CALL DIC_STORE_FLUXCO2( myTime, myIter, myThid )
# endif
# ifdef ALLOW_COST
       CALL DIC_COST( myTime, myIter, myThid )
# endif
# ifndef ALLOW_AUTODIFF
      ENDIF
# endif /* ALLOW_AUTODIFF */
#endif /* ALLOW_DIC */

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('GCHEM_FORCING_SEP',myThid)
#endif

#endif /* GCHEM_SEPARATE_FORCING */
#endif /* ALLOW_GCHEM */

C       Apply Correction
        do bj=myByLo(myThid),myByHi(myThid)
        do bi=myBxLo(myThid),myBxHi(myThid)
        do k = 1, Nr
        do j = 1, sNy
        do i = 1, sNx

C       Nitrate Correction
        pTracer(i,j,k,bi,bj,1)=pTracer(i,j,k,bi,bj,1)+
     &  (totSurfCorPTr(1)+totSurfCorPTr(2)+totSurfCorPTr(4)+
     &  totSurfCorPTr(5)+totSurfCorPTr(7)+totSurfCorPTr(8)+
     &  totSurfCorPTr(9)+totSurfCorPTr(10)+totSurfCorPTr(11))*
     &  deltaT/globalVolume*maskC(i,j,k,bi,bj)

C       Silica Correction
        pTracer(i,j,k,bi,bj,3)=pTracer(i,j,k,bi,bj,3)+
     &  (totSurfCorPTr(3)+totSurfCorPTr(6)+totSurfCorPTr(8)*r_SI_N)*
     &  deltaT/globalVolume*maskC(i,j,k,bi,bj)

C NOTE: correction equation is simplified, original equation would look
C like correction=totSurfCorPTr*deltaT*vcell/vtot*1/vcell
C totSurfCorPTr is the total correction in units of mmol N/s for the
C tile not the cell. Since vcel/vtot*1/vcell simplify to 1/vtot the
C equation above is correction = totSurfCorPTr*deltaT*1/vtot

        end do
        end do
        end do
        end do
        end do

C       River Contribution (Surface Cell Only)
        pTracer(:,:,1,:,:,1)=pTracer(:,:,1,:,:,1)+
     &  rivers*deltaT
        pTracer(:,:,1,:,:,3)=pTracer(:,:,1,:,:,3)+
     &  rivers*r_SI_N_riv*deltaT

C       Reset Diagnostic Fields
        if (mod(myIter,floor(dumpFreq/deltaT))==1) then
        print*, '* Reset Diagnostic Arrays (gchem_forcing_sep.f) *'

        phy_chl(:,:,:,:,:)=0.0
        l_lim_sp_field(:,:,:,:,:)=0.0
        l_lim_lp_field(:,:,:,:,:)=0.0

        print*, 'Set Arrays = 0'

        day_count=day_count+1
        day_count_tot=day_count_tot+1

        if (day_count>=366) then
        day_count=1
        year_count=year_count+1
        end if

        end if

C NOTE: if myIter = 1 the bio iteration is 0, so when myIter = 25 there
C have been 24 bio iterations and it is time to reset the diagnostic
C field, hence mod(25,24)=1.

C   ##########    Biological Time Step ############   C
        if (mod(myIter,floor(dt_b_sec/deltaT))==0) then

        var  = totSurfCorPTr*deltaT/globalVolume*maskC(i,j,k,bi,bj)

        print*, '***************************'
        print*, '******** BIO STEP *********'
        print*, '***************************'

C       Stats:
        print*, 'Year=',year_count,'Day = ', day_count
        print*, 'Day Number Total = ', day_count_tot
        print*, 'myIter = ', myIter
        print*, 'Total N Surface Correction = ',
     &  sum(var(1:2))+sum(var(4:5))+sum(var(7:11))
        print*, 'Total Si Surface Correction = ',
     &  var(3)+var(6)+var(8)
        print*, 'Light = ', minval(light) , ':', maxval(light)
        print*, 'Temperatue = ', minval(theta), ':', maxval(theta)
        print*, 'Salinity = ', minval(salt), ':', maxval(salt)

C       Biological Loop
        do bj=myByLo(myThid),myByHi(myThid)
        do bi=myBxLo(myThid),myBxHi(myThid)
        do k = 1, Nr
        do j = 1-OLy, sNy+OLy
        do i = 1-OLx, sNx+OLx

        if (maskC(i,j,k,bi,bj) == 1.0) then

        NO_v1=max(pTracer(i,j,k,bi,bj,1),1.0e-5)
        NH_v1=max(pTracer(i,j,k,bi,bj,2),1.0E-5)
        SI_v1=max(pTracer(i,j,k,bi,bj,3),1.0E-5)
        DON_v1=max(ptracer(i,j,k,bi,bj,4),1.0E-5)
        PON_v1=max(ptracer(i,j,k,bi,bj,5),1.0E-5)
        OP_v1=max(ptracer(i,j,k,bi,bj,6),1.0E-5)
        sp_v1=max(ptracer(i,j,k,bi,bj,7),1.0E-5)
        lp_v1=max(ptracer(i,j,k,bi,bj,8),1.0E-5)
        sz_v1=max(ptracer(i,j,k,bi,bj,9),1.0E-5)
        lz_v1=max(ptracer(i,j,k,bi,bj,10),1.0E-5)
        pz_v1=max(ptracer(i,j,k,bi,bj,11),1.0E-5)
        pass_v1=max(ptracer(i,j,k,bi,bj,12),1.0E-5)
        th_v1=max(ptracer(i,j,k,bi,bj,13),1.0E-5)

        NO_v2=NO_v1
        NH_v2=NH_v1
        SI_v2=SI_v1
        DON_v2=DON_v1
        PON_v2=PON_v1
        OP_v2=OP_v1
        sp_v2=sp_v1
        lp_v2=lp_v1
        sz_v2=sz_v1
        lz_v2=lz_v1
        pz_v2=pz_v1
        th_v2 = th_v1
        pass_v2 = pass_v1
        tmp=theta(i,j,k,bi,bj)

C       Light Attenuation
        if (k==1) then
        light_adj(i,j,bi,bj)=light(i,j,bi,bj)*PAR_frac
     &  *exp(-1.0*(ext_sp*sp_v1+ext_lp*lp_v1+ext_w)
     &  *(depth_uv(k)-0.0))
        else
        light_adj(i,j,bi,bj)=light_adj(i,j,bi,bj)
     &  *exp(-1.0*(ext_sp*sp_v1+ext_lp*lp_v1+ext_w)
     &  *(depth_uv(k)-depth_uv(k-1)))
        end if

        l_dum=light_adj(i,j,bi,bj)

C       Backward Euler Loop
        do iter=1,5

C       %%%%%%%%%%%%%%%%%%%%
C       %%%  N E M U R O %%%
C       %%%%%%%%%%%%%%%%%%%%

C       ### Functional Groups: ###

C       *** Nutrient Uptake Terms ***

C       % --- Small Phytoplankton --- %
C       Calculate Limitations
        l_lim_sp = (1.0-exp(-1.0*photo_chem_sp*l_dum/vmax_sp))
     &  *exp(-1.0*photo_inh_sp*l_dum/vmax_sp)
        t_lim_v_sp = exp(tc_v_sp*tmp)*scale_tmp_phyto
C        NO_lim_sp = NO_v2/(NO_v2+k_NO_sp)*exp(-1.0*inh_NH_NO_sp*NH_v2)
        NO_lim_sp = NO_v2/(NO_v2+k_NO_sp)*(1.0/(1.0+NH_v2/k_NH_sp))
        NH_lim_sp = NH_v2/(NH_v2+k_NH_sp)

C       sp NO Uptake
        delta_NO = NO_v1-(NO_v1/(1.0+dt_b/NO_v2*(vmax_sp*t_lim_v_sp
     &  *l_lim_sp*NO_lim_sp*sp_v2)))

C       sp NH Uptake
        delta_NH = NH_v1-(NH_v1/(1.0+dt_b/NH_v2*(vmax_sp*t_lim_v_sp
     &  *l_lim_sp*NH_lim_sp*sp_v2)))

C       sp Excretion
        delta_sp = ext_excr_sp*(delta_NO+delta_NH)

C       Update
        NO_v2 = NO_v1-delta_NO
        NH_v2 = NH_v1-delta_NH
        DON_v2 = DON_v1+delta_sp
        sp_v2 = sp_v1+delta_NO+delta_NH-delta_sp

C       Other
        np_frac_sp = delta_NO/(max(min_val,(delta_NO+delta_NH)))

C       Diagnostics
C       % ---------------------- %

C       % --- Large Phytoplankton --- %
C       Calculate Limitations
        t_lim_v_lp = exp(tc_v_lp*tmp)*scale_tmp_phyto
        l_lim_lp = (1.0-exp(-1.0*photo_chem_lp*l_dum/vmax_lp))
     &  *exp(-1.0*photo_inh_lp*l_dum/vmax_lp)
C        NO_lim_lp = NO_v2/(NO_v2+k_NO_lp)*exp(-inh_NH_NO_lp*NH_v2)
        NO_lim_lp = NO_v2/(NO_v2+k_NO_lp)*(1.0/(1.0+NH_v2/k_NH_lp))
        NH_lim_lp = NH_v2/(NH_v2+k_NH_lp)
        SI_lim_lp = SI_v2/(SI_v2+k_SI_lp)
        NS_lim_comp = min(1.0,(SI_lim_lp
     &  /(max(min_val,(NO_lim_lp+NH_lim_lp)))))

C       lp NO Uptake
        delta_NO = NO_v2-(NO_v2/(1.0+dt_b/NO_v2*(vmax_lp*t_lim_v_lp
     &  *l_lim_lp*NO_lim_lp*NS_lim_comp*lp_v2)))

C       lp NH Uptake
        delta_NH = NH_v2-(NH_v2/(1.0+dt_b/NH_v2*(vmax_lp*t_lim_v_lp
     &  *l_lim_lp*NH_lim_lp*NS_lim_comp*lp_v2)))

C       lp Excretion
        delta_lp = ext_excr_lp*(delta_NO+delta_NH)

C       Update
        NO_v2 = NO_v2-delta_NO
        NH_v2 = NH_v2-delta_NH
        SI_v2 = SI_v1-r_SI_N*(delta_NO+delta_NH)+r_SI_N*delta_lp
        DON_v2 = DON_v2+delta_lp
        lp_v2 = lp_v1+delta_NO+delta_NH-delta_lp

C       Other
        np_frac_lp = delta_NO/(max(min_val,(delta_NO+delta_NH)))

C       Diagnostics
C       % ---------------------- %

C       *** END Nutrient Uptake Terms *** %

C       *** Phytoplankton Mortality Terms *** %

C       Calculate Limiations
        t_lim_m_sp = exp(tc_m_sp*tmp)*scale_tmp_phyto

C       sp Mortality
        delta_sp = sp_v2-(sp_v2/(1.0+dt_b/sp_v2*(ref_mort_sp*t_lim_m_sp
     &  *sp_v2)))

C       Calculate Limitations
        t_lim_m_lp = exp(tc_m_lp*tmp)*scale_tmp_phyto

C       lp Mortality
        delta_lp = lp_v2-(lp_v2/(1.0+dt_b/lp_v2*(ref_mort_lp*t_lim_m_lp
     &  *lp_v2)))

C       Update
        sp_v2 = sp_v2-delta_sp
        lp_v2 = lp_v2-delta_lp
        PON_v2 = PON_v1+delta_sp+delta_lp
        OP_v2 = OP_v1+delta_lp*r_SI_N

C       *** END of Mortality Terms ***

C       *** Phytoplankton Respiration Terms ***

C       Calculate Limitations
        t_lim_r_sp = exp(tc_r_sp*tmp)*scale_tmp_phyto2

C       sp Respiration
        delta_sp = sp_v2-(sp_v2/(1.0+dt_b/sp_v2*(ref_resp_sp*t_lim_r_sp
     &  *sp_v2)))

C       Calculate Limitations
        t_lim_r_lp = exp(tc_r_lp*tmp)*scale_tmp_phyto2

C       lp Respiration
        delta_lp = lp_v2-(lp_v2/(1.0+dt_b/lp_v2*(ref_resp_lp*t_lim_r_lp
     &  *lp_v2)))

C       Update
        sp_v2 = sp_v2-delta_sp
        lp_v2 = lp_v2-delta_lp
        NO_v2 = NO_v2+delta_sp*np_frac_sp+delta_lp*np_frac_lp
        NH_v2 = NH_v2+delta_sp*(1.0-np_frac_sp)+delta_lp*
     &  (1.0-np_frac_lp)
        SI_v2 = SI_v2+delta_lp*r_SI_N

C       Diagnostics

C       *** END PHytoplankton Respiration Terms ***

C       *** Zooplankton Grazing Terms ***

C       % --- Small Zooplankton --- %
C       Calculate Limitations
        t_lim_g_sz_sp = exp(tc_g_sz_sp*tmp)*scale_tmp_zoo
        g_lim_sz_sp = max(0.0,(1.0-exp(iv_sz_sp*(thresh_sz_sp-sp_v2))))

C       sz Grazing on sp
        delta_sp = sp_v2-(sp_v2/(1.0+dt_b/sp_v2*(gmax_sz_sp
     &  *t_lim_g_sz_sp*g_lim_sz_sp*sz_v2)))

C       Update
        sp_v2 = sp_v2-delta_sp
        sz_v2 = sz_v1+gge_sz*delta_sp
        NH_v2 = NH_v2+(ae_sz-gge_sz)*delta_sp
        PON_v2 = PON_v2+(1.0-ae_sz)*delta_sp

C       Diagnostics
C       % ------------------------ %

C       % --- Large Zooplankton --- %
C       Calculate Limitations
        t_lim_g_lz_sp = exp(tc_g_lz_sp*tmp)*scale_tmp_zoo
        g_lim_lz_sp = max(0.0,(1.0-exp(iv_lz_sp*(thresh_lz_sp-sp_v2))))

C       lz Grazing on sp
        delta_sp = sp_v2-(sp_v2/(1.0+dt_b/sp_v2*(gmax_lz_sp
     &  *t_lim_g_lz_sp*g_lim_lz_sp*lz_v2)))

C       Calculate Limitations
        t_lim_g_lz_lp = exp(tc_g_lz_lp*tmp)*scale_tmp_zoo
        g_lim_lz_lp = max(0.0,(1.0-exp(iv_lz_lp*(thresh_lz_lp-lp_v2))))

C       lz Grazing on lp
        delta_lp = lp_v2-(lp_v2/(1.0+dt_b/lp_v2*(gmax_lz_lp
     &  *t_lim_g_lz_lp*g_lim_lz_lp*lz_v2)))

C       Calculate Limitations
        t_lim_g_lz_sz = exp(tc_g_lz_sz*tmp)*scale_tmp_zoo
        g_lim_lz_sz = max(0.0,(1.0-exp(iv_lz_sz*(thresh_lz_sz-sz_v2))))

C       lz Grazing on sz
        delta_sz = sz_v2-(sz_v2/(1.0+dt_b/sz_v2*(gmax_lz_sz
     &  *t_lim_g_lz_sz*g_lim_lz_sz*lz_v2)))

C       Update
        sp_v2 = sp_v2-delta_sp
        lp_v2 = lp_v2-delta_lp
        sz_v2 = sz_v2-delta_sz
        lz_v2 = lz_v1+gge_lz*(delta_sp+delta_lp+delta_sz)
        NH_v2 = NH_v2+(ae_lz-gge_lz)*(delta_sp+delta_lp+delta_sz)
        PON_v2 = PON_v2+(1.0-ae_lz)*(delta_sp+delta_lp+delta_sz)
        OP_v2 = OP_v2+delta_lp*r_SI_N

C       Diagnostics
C       % -------------------------%

C       % --- Predatroy Zooplankton --- %
C       Calculate Limitations
        t_lim_g_pz_lp = exp(tc_g_pz_lp*tmp)*scale_tmp_zoo
        g_lim_pz_lp = max(0.0,((1.0-exp(iv_pz_lp*(thresh_pz_lp-lp_v2)))
     &  *exp(-inh_szlz_lp*(sz_v2+lz_v2))))

C       pz Grazing on lp
        delta_lp = lp_v2-(lp_v2/(1.0+dt_b/lp_v2*(gmax_pz_lp
     &  *t_lim_g_pz_lp*g_lim_pz_lp*pz_v2)))

C       Calculate Limitations
        t_lim_g_pz_sz = exp(tc_g_pz_sz*tmp)*scale_tmp_zoo
        g_lim_pz_sz = max(0.0,((1.0-exp(iv_pz_sz*(thresh_pz_sz-sz_v2)))
     &  *exp(-inh_lz_sz*lz_v2)))

C       pz Grazing on sz
        delta_sz = sz_v2-(sz_v2/(1.0+dt_b/sz_v2*(gmax_pz_sz
     &  *t_lim_g_pz_sz*g_lim_pz_sz*pz_v2)))

C       Calculate Limitations
        t_lim_g_pz_lz = exp(tc_g_pz_lz*tmp)*scale_tmp_zoo
        g_lim_pz_lz = max(0.0,(1.0-exp(iv_pz_lz*(thresh_pz_lz-lz_v2))))

C       pz Grazing on lz
        delta_lz = lz_v2-(lz_v2/(1.0+dt_b/lz_v2*(gmax_pz_lz
     &  *t_lim_g_pz_lz*g_lim_pz_lz*pz_v2)))

C       Update
        lp_v2 = lp_v2-delta_lp
        sz_v2 = sz_v2-delta_sz
        lz_v2 = lz_v2-delta_lz
        pz_v2 = pz_v1+gge_pz*(delta_lp+delta_sz+delta_lz)
        NH_v2 = NH_v2+(ae_pz-gge_pz)*(delta_lp+delta_sz+delta_lz)
        PON_v2 = PON_v2+(1.0-ae_pz)*(delta_lp+delta_sz+delta_lz)
        OP_v2 = OP_v2+delta_lp*r_SI_N

C       Diagnostics
C       % ----------------------------- %

C       *** END of Grazing Terms ***

C       *** Zooplankton Mortality Terms ***

C       Calculate Limitations
        t_lim_m_sz = exp(tc_m_sz*tmp)*scale_tmp_zoo

C       sz Mortality
        delta_sz = sz_v2-(sz_v2/(1.0+dt_b/sz_v2*(ref_mort_sz*t_lim_m_sz
     &  *sz_v2)))

C       Calculate Limitations
        t_lim_m_lz = exp(tc_m_lz*tmp)*scale_tmp_zoo

C       lz Mortality
        delta_lz = lz_v2-(lz_v2/(1.0+dt_b/lz_v2*(ref_mort_lz*t_lim_m_lz
     &  *lz_v2)))

C       Calculate Limitations
        t_lim_m_pz = exp(tc_m_pz*tmp)*scale_tmp_zoo

C       lz Mortality
        delta_pz = pz_v2-(pz_v2/(1.0+dt_b/pz_v2*(ref_mort_pz*t_lim_m_pz
     &  *pz_v2*pz_v2)))

C       Update
        sz_v2 = sz_v2-delta_sz
        lz_v2 = lz_v2-delta_lz
        pz_v2 = pz_v2-delta_pz
        PON_v2 = PON_v2+delta_sz+delta_lz+delta_pz
C       % ----------------- %

C       *** END of Mortality Terms ***

C       ### END Functional Groups ###

C       ### Chemical Transformations ###

C       Nitrification
        t_lim_nit = exp(tc_nitr*tmp)*scale_tmp_decomp
        delta_NH = NH_v2-(NH_v2/(1.0+dt_b/NH_v2*(ref_nitr*t_lim_nit
     &  *NH_v2)))
        NH_v2 = NH_v2-delta_NH
        NO_v2 = NO_v2+delta_NH

C       Decomp of PON to NH
        t_lim_dec_PON_NH = exp(tc_dec_PON_NH*tmp)*scale_tmp_decomp
        delta_PON = PON_v2-(PON_v2/(1.0+dt_b/PON_v2*(t_lim_dec_PON_NH
     &  *ref_dec_PON_NH*PON_v2)))
        PON_v2 = PON_v2-delta_PON
        NH_v2 = NH_v2+delta_PON

C       Decomp of DON to NH
        t_lim_dec_DON_NH = exp(tc_dec_DON_NH*tmp)*scale_tmp_decomp
        delta_DON = DON_v2-(DON_v2/(1.0+dt_b/DON_v2*(t_lim_dec_DON_NH
     &  *ref_dec_DON_NH*DON_v2)))
        DON_v2 = DON_v2-delta_DON
        NH_v2 = NH_v2+delta_DON

C       Decomp of PON to DON
        t_lim_dec_PON_DON = exp(tc_dec_PON_DON*tmp)*scale_tmp_decomp
        delta_PON = PON_v2-(PON_v2/(1.0+dt_b/PON_v2*(t_lim_dec_PON_DON
     &  *ref_dec_PON_DON*PON_v2)))
        PON_v2 = PON_v2-delta_PON
        DON_v2 = DON_v2+delta_PON

C       Decomp of OP to SI (First OP Loss Term)
        t_lim_dec_OP_SI = exp(tc_dec_OP_SI*tmp)*scale_tmp_decomp
        delta_OP = OP_v2-(OP_v2/(1.0+dt_b/OP_v2*(ref_dec_OP_SI
     &  *t_lim_dec_OP_SI*OP_v2)))
        OP_v2 = OP_v2-delta_OP
        SI_v2 = SI_v2+delta_OP


C       ### END Chemical Transformations ###

       if (NO_v2 /= NO_v2 .or. NO_v2 < 0.0 .or.
     & NH_v2 /= NH_v2 .or. NH_v2 < 0.0 .or.
     & SI_v2 /= SI_v2 .or. SI_v2 < 0.0 .or.
     & DON_v2 /= DON_v2 .or. DON_v2 < 0.0 .or.
     & PON_v2 /= PON_v2 .or. PON_v2 < 0.0 .or.
     & OP_v2 /= OP_v2 .or. OP_v2 < 0.0 .or.
     & sp_v2 /= sp_v2 .or. sp_v2 < 0.0 .or.
     & lp_v2 /= lp_v2 .or. lp_v2 < 0.0 .or.
     & sz_v2 /= sz_v2 .or. sz_v2 < 0.0 .or.
     & lz_v2 /= lz_v2 .or. lz_v2 < 0.0 .or.
     & pz_v2 /= pz_v2 .or. pz_v2 < 0.0 ) then

       print*, 'Index (i,j,k,iter) = ', i, j, k, iter
       print* ,' '
       print*,'NO_v1 = ', NO_v1
       print*,'NH_v1 = ', NH_v1
       print*,'SI_v1 = ', SI_v1
       print*,'DON_v1 = ', DON_v1
       print*,'PON_v1 = ', PON_v1
       print*,'OP_v1 = ', OP_v1
       print*,'sp_v1 = ', sp_v1
       print*,'lp_v1 = ', lp_v1
       print*,'sz_v1 = ', sz_v1
       print*,'lz_v1 = ', lz_v1
       print*,'pz_v1 = ', pz_v1
       print*,'pass_v1 = ', pass_v1
       print*,'th_v1 = ', th_v1
       print*, ' '
       print*,'NO_v2 = ', NO_v2
       print*,'NH_v2 = ', NH_v2
       print*,'SI_v2 = ', SI_v2
       print*,'DON_v2 = ', DON_v2
       print*,'PON_v2 = ', PON_v2
       print*,'OP_v2 = ', OP_v2
       print*,'sp_v2 = ', sp_v2
       print*,'lp_v2 = ', lp_v2
       print*,'sz_v2 = ', sz_v2
       print*,'lz_v2 = ', lz_v2
       print*,'pz_v2 = ', pz_v2
       print*,'pass_v2 = ', pass_v2
       print*,'th_v2 = ', th_v2
       print*, ' '
       pause
       end if

C iter
        end do

C       Thorium
      th_v2=th_v1/(1.0-0.003-0.0003*(sp_v2+lp_v2))+0.5

C       Update State Variables
        pTracer(i,j,k,bi,bj,1)=NO_v2
        pTracer(i,j,k,bi,bj,2)=NH_v2
        pTracer(i,j,k,bi,bj,3)=SI_v2
        pTracer(i,j,k,bi,bj,4)=DON_v2
        pTracer(i,j,k,bi,bj,5)=PON_v2
        pTracer(i,j,k,bi,bj,6)=OP_v2
        pTracer(i,j,k,bi,bj,7)=sp_v2
        pTracer(i,j,k,bi,bj,8)=lp_v2
        pTracer(i,j,k,bi,bj,9)=sz_v2
        pTracer(i,j,k,bi,bj,10)=lz_v2
        pTracer(i,j,k,bi,bj,11)=pz_v2
        pTracer(i,j,k,bi,bj,13)=th_v2

C       Calculate Diagnostics
        l_lim_sp_field(i,j,k,bi,bj)=l_lim_sp_field(i,j,k,bi,bj)+l_lim_sp
        l_lim_lp_field(i,j,k,bi,bj)=l_lim_lp_field(i,j,k,bi,bj)+l_lim_lp

C Note: light limitation and chl diagnostics are still added together so
C I can output a daily average light limitation and chl

C maskC(i,j,k,bi,bj)==1
        end if

C i
        end do
C j
        end do
C k
        end do
C bi
        end do
C bj
        end do

C       %%%%%%%%%%%%%%%%%%%%%%%%%%
C       %%% E N D  N E M U R O %%%
C       %%%%%%%%%%%%%%%%%%%%%%%%%%

C       *** Sinkng Detritus ***
        temp_PON(:,:,:,:,:)=pTracer(:,:,:,:,:,5)
        temp_OP(:,:,:,:,:)=pTracer(:,:,:,:,:,6)

        do bj=myByLo(myThid),myByHi(myThid)
        do bi=myBxLo(myThid),myBxHi(myThid)
        do k = 1, Nr
        do j = 1-OLy, sNy+OLy
        do i = 1-OLx, sNx+OLx

        if (maskC(i,j,k,bi,bj) == 1.0) then
        cell_depth = depth(i,j,bi,bj)
        PON_v1 = pTracer(i,j,k,bi,bj,5)
        OP_v1 = pTracer(i,j,k,bi,bj,6)
        volume1= rA(i,j,bi,bj)*drF(k)*hFacC(i,j,k,bi,bj)
        sink_frac1 = sink*dt_b/drF(k)
        NO_v1=pTracer(i,j,k,bi,bj,1)
        SI_v1=pTracer(i,j,k,bi,bj,3)

C       Surface Cells
        if (k==1) then
        if (cell_depth>depth_w(k+1)) then
        temp_PON(i,j,k,bi,bj)=PON_v1-(PON_v1*sink_frac1)
        temp_OP(i,j,k,bi,bj)=OP_v1-(OP_v1*sink_frac1)
        end if
        if (cell_depth<=depth_w(k+1)) then
        temp_PON(i,j,k,bi,bj)=PON_v1-(PON_v1*sink_frac1)*
     &  (sed_frac+diagen_frac)
        pTracer(i,j,k,bi,bj,1)=NO_v1+PON_v1*sink_frac1*diagen_frac

        temp_OP(i,j,k,bi,bj)=OP_v1-(OP_v1*sink_frac1)*
     &  (sed_frac+diagen_frac)
        pTracer(i,j,k,bi,bj,3)=SI_v1+OP_v1*sink_frac1*diagen_frac
        end if
        end if

C       Water Column Cells
        if (cell_depth>depth_w(k+1) .and. k>1) then
        PON_v2 = pTracer(i,j,k-1,bi,bj,5)
        OP_v2 = pTracer(i,j,k-1,bi,bj,6)
        volume2 = rA(i,j,bi,bj)*drF(k-1)*hFacC(i,j,k-1,bi,bj)
        sink_frac2 = sink*dt_b/drF(k-1)

        temp_PON(i,j,k,bi,bj)=(PON_v1*volume1
     &  +(PON_v2*volume2*sink_frac2)
     &  -(PON_v1*volume1*sink_frac1))/volume1

        temp_OP(i,j,k,bi,bj)=(OP_v1*volume1
     &  +(OP_v2*volume2*sink_frac2)
     &  -(OP_v1*volume1*sink_frac1))/volume1
        end if

C       Bottom Cells
        if (cell_depth==depth_w(k+1) .and. k>1) then
        PON_v2 = pTracer(i,j,k-1,bi,bj,5)
        OP_v2 = pTracer(i,j,k-1,bi,bj,6)
        volume2 = rA(i,j,bi,bj)*drF(k-1)*hFacC(i,j,k-1,bi,bj)
        sink_frac2 = sink*dt_b/drF(k-1)

        temp_PON(i,j,k,bi,bj)=(PON_v1*volume1
     &  +(PON_v2*volume2*sink_frac2)
     &  -(PON_v1*volume1*sink_frac1)*(sed_frac+diagen_frac))/volume1
        pTracer(i,j,k,bi,bj,1)=NO_v1+PON_v1*sink_frac1*diagen_frac

        temp_OP(i,j,k,bi,bj)=(OP_v1*volume1
     &  +(OP_v2*volume2*sink_frac2)
     &  -(OP_v1*volume1*sink_frac1)*(sed_frac+diagen_frac))/volume1
        pTracer(i,j,k,bi,bj,3)=SI_v1+OP_v1*sink_frac1*diagen_frac
        end if

        end if

        end do
        end do
        end do
        end do
        end do

        pTracer(:,:,:,:,:,5)=temp_PON(:,:,:,:,:)
        pTracer(:,:,:,:,:,6)=temp_OP(:,:,:,:,:)

C       *** Chl:C Model Diagnostic ***
        do bj=myByLo(myThid),myByHi(myThid)
        do bi=myBxLo(myThid),myBxHi(myThid)
        do k = 1, Nr
        do j = 1-OLy, sNy+OLy
        do i = 1-OLx, sNx+OLx

        if (maskC(i,j,k,bi,bj) == 1.0) then


        NO_v1=pTracer(i,j,k,bi,bj,1)
        NH_v1=pTracer(i,j,k,bi,bj,2)
        SI_v1=pTracer(i,j,k,bi,bj,3)
        sp_v1=ptracer(i,j,k,bi,bj,7)
        lp_v1=ptracer(i,j,k,bi,bj,8)
        tmp=theta(i,j,k,bi,bj)

C       Light Attenuation
        if (k==1) then
        light_adj(i,j,bi,bj)=light(i,j,bi,bj)*PAR_frac
     &  *exp(-1.0*(ext_sp*sp_v1+ext_lp*lp_v1+ext_w)
     &  *(depth_uv(k)-0.0))
        else
        light_adj(i,j,bi,bj)=light_adj(i,j,bi,bj)
     &  *exp(-1.0*(ext_sp*sp_v1+ext_lp*lp_v1+ext_w)
     &  *(depth_uv(k)-depth_uv(k-1)))
        end if

        l_dum=light_adj(i,j,bi,bj)

C       Small Phytoplankton Chl:C
        t_lim_v_sp = exp(tc_v_sp*tmp)*scale_tmp_phyto
C        NO_lim_sp = NO_v1/(NO_v1+k_NO_sp)*exp(-1.0*inh_NH_NO_sp*NH_v1)
        NO_lim_sp = NO_v1/(NO_v1+k_NO_sp)*(1.0/(1.0+NH_v1/k_NH_sp))
        NH_lim_sp = NH_v1/(NH_v1+k_NH_sp)

        v_chl_sp=vmax_sp*t_lim_v_sp*(NO_lim_sp+NH_lim_sp)*(photo_chem_sp
     &  /(photo_chem_sp+photo_inh_sp))*(photo_inh_sp/(photo_chem_sp
     &  +photo_inh_sp))**(photo_inh_sp/photo_chem_sp)

        sp_chl2c=sp_chl2c_max/(1.0+sp_chl2c_max*alpha_chl*l_dum*0.5
     &  /v_chl_sp)

        if (sp_chl2c<sp_chl2c_min) then
        sp_chl2c=sp_chl2c_min
        end if
        sp_chl=sp_v1/1000*6.625*12.011*sp_chl2c*1000

C       Large Phytoplankton Chl:C
        t_lim_v_lp = exp(tc_v_lp*tmp)*scale_tmp_phyto
C        NO_lim_lp = NO_v1/(NO_v1+k_NO_lp)*exp(-inh_NH_NO_lp*NH_v1)
        NO_lim_lp = NO_v1/(NO_v1+k_NO_lp)*(1.0/(1.0+NH_v1/k_NH_lp))
        NH_lim_lp = NH_v1/(NH_v1+k_NH_lp)
        SI_lim_lp = SI_v1/(SI_v1+k_SI_lp)
        NS_lim_comp = min(1.0,(SI_lim_lp
     &  /(max(min_val,(NO_lim_lp+NH_lim_lp)))))

        v_chl_lp=vmax_lp*t_lim_v_lp*(NO_lim_lp+NH_lim_lp)*NS_lim_comp
     &  *(photo_chem_lp/(photo_chem_lp+photo_inh_lp))*(photo_inh_lp
     &  /(photo_chem_lp+photo_inh_lp))**(photo_inh_lp/photo_chem_lp)

        lp_chl2c=lp_chl2c_max/(1.0+lp_chl2c_max*alpha_chl*l_dum*0.5
     &  /v_chl_lp)

        if (lp_chl2c<lp_chl2c_min) then
        lp_chl2c=lp_chl2c_min
        end if
        lp_chl=lp_v1/1000*6.625*12.011*lp_chl2c*1000

        phy_chl(i,j,k,bi,bj)=phy_chl(i,j,k,bi,bj)+sp_chl+lp_chl

        end if

        end do
        end do
        end do
        end do
        end do

       bio_count=bio_count+1

C       *** Diagnostics ****
        if (mod(myIter,floor(dumpFreq/dt_b_sec))==0) then

        print*, '* Filling Bio Diagnostic Array (gchem_forcing_sep.F) *'
        print*, 'Bio Count = ',bio_count

C       Fill Diagnostics
        phy_chl=phy_chl/real(bio_count)
        l_lim_sp_field=l_lim_sp_field/real(bio_count)
        l_lim_lp_field=l_lim_lp_field/real(bio_count)

        CALL DIAGNOSTICS_FILL(phy_chl,'phy_chl   ',0,Nr,0,1,1,
     &  myThid)
        CALL DIAGNOSTICS_FILL(l_lim_sp_field,'l_lim_sp_field   '
     &  ,0,Nr,0,1,1,myThid)
        CALL DIAGNOSTICS_FILL(l_lim_lp_field,'l_lim_lp_field   '
     &  ,0,Nr,0,1,1,myThid)

        bio_count=0

        end if

        end if
C  END Bio Time Step:  if (mod(myIter,floor(14400/deltaT))==1) then

      RETURN
      END
