C $Header: /u/gcmpack/MITgcm/pkg/ptracers/ptracers_apply_forcing.F,v 1.3 2016/08/23 18:29:49 jahn Exp $
C $Name:  $

#include "PTRACERS_OPTIONS.h"

CBOP
C !ROUTINE: PTRACERS_APPLY_FORCING

C !INTERFACE: ==========================================================
      SUBROUTINE PTRACERS_APPLY_FORCING(
     U                    gPtracer,
     I                    surfForcPtr,
     I                    iMin,iMax,jMin,jMax, k, bi, bj,
     I                    iTracer, myTime, myIter, myThid )

C !DESCRIPTION:
C     Apply passive tracer forcing, i.e., sources and sinks of tracer,
C      by adding forcing terms to the tendency array

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"

C !INPUT PARAMETERS: ===================================================
C  gPtracer             :: the tendency array
C  surfForcPtr          :: surface forcing term
C  iMin iMax jMin jMax  :: working range of tile for applying forcing
C  k                    :: vertical level number
C  bi,bj                :: tile indices
C  iTracer              :: tracer number
C  myIter               :: time-step number
C  myTime               :: model time
C  myThid               :: thread number
      _RL gPtracer   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfForcPtr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin,iMax,jMin,jMax
      INTEGER k, bi,bj, iTracer, kk
      _RL myTime, vol_sum
      INTEGER myIter
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  gPtracer       :: updates tendency array

#ifdef ALLOW_PTRACERS

C !LOCAL VARIABLES: ====================================================
C  i,j            :: loop indices
      INTEGER i,j
C     number of surface interface layer
      INTEGER kSurface
CEOP

      IF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF

C--   Surface forcing term surfForcPtr is previously computed by
C     PTRACERS_FORCING_SURF (and stored as surfaceForcingPTr)
C     because it is needed by KPP_TRANSPORT_PTR.

#ifdef ALLOW_GCHEM
      IF ( useGCHEM ) THEN
        CALL GCHEM_ADD_TENDENCY(
     U                 gPtracer,
     I                 iMin,iMax,jMin,jMax, k, bi, bj,
     I                 iTracer, myTime, myIter, myThid )
      ENDIF
#endif /* ALLOW_GCHEM */

      IF ( k .EQ. kSurface ) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
          gPtracer(i,j) = gPtracer(i,j)
     &                  + surfForcPtr(i,j)
     &                   *recip_drF(k)*recip_hFacC(i,j,k,bi,bj)
         ENDDO
        ENDDO
      ENDIF

C      IF (PTRACERS_linFSConserve(iTracer)) THEN
C        IF (k .EQ. 1) THEN

C        vol_sum=0
C        do bj=myByLo(myThid),myByHi(myThid)
C        do bi=myBxLo(myThid),myBxHi(myThid)
C        do kk = 1, Nr
C        do j = 1-OLy, sNy+OLy
C        do i = 1-OLx, sNx+OLx

C        vol_sum=vol_sum+rA(i,j,bi,bj)*drF(kk)*hFacC(i,j,kk,bi,bj)

C        end do 
C        end do
C        end do  
C        end do
C        end do

C        end if 

C          DO j=0,sNy+1
C           DO i=0,sNx+1
C           gTracer(i,j,k) = gTracer(i,j,k)+
C     &     sum(totSurfCorPTR)/vol_sum
C           ENDDO
C          ENDDO

C       ENDIF

      IF (PTRACERS_linFSConserve(iTracer)) THEN
        IF (k .EQ. kSurface) THEN
           DO j=0,sNy+1
            DO i=0,sNx+1
             gPtracer(i,j) = gPtracer(i,j)
C     &       +meanSurfCorPTr(iTracer)*recip_drF(k)
C     &       *_recip_hFacC(i,j,k,bi,bj)
           ENDDO
          ENDDO
        ENDIF
      ENDIF

#ifdef ALLOW_RBCS
      IF ( useRBCS ) THEN
         CALL RBCS_ADD_TENDENCY(
     U                 gPtracer,
     I                 k, bi, bj, iTracer+2,
     I                 myTime, myIter, myThid )
      ENDIF
#endif /* ALLOW_RBCS */

#endif /* ALLOW_PTRACERS */

      RETURN
      END
