C $Header: /u/gcmpack/MITgcm/pkg/gchem/gchem_init_fixed.F,v 1.12 2016/05/23 13:09:31 jmc Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"

C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_INIT_FIXED(myThid )

C !DESCRIPTION:
C calls subroutines that initialized fixed variables for any tracer
C experiment

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GCHEM.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
      INTEGER        diagNum
      INTEGER        diagMate
      CHARACTER*8    diagName
      CHARACTER*16   diagCode
      CHARACTER*16   diagUnits
      CHARACTER*(80) diagTitle
CEOP

#ifdef ALLOW_GCHEM

C     !LOCAL VARIABLES:

C jmc: simpler to put MASTER thread directive within individual S/R :
C      that leave the choice, for each S/R, to do either single (master)
C      thread tasks or Multi-threaded tasks.
c        _BEGIN_MASTER( myThid )
ccccccccccccccccccccccccccc
c initialize parameters   c
ccccccccccccccccccccccccccc
#ifdef ALLOW_CFC
      IF ( useCFC ) THEN
         print*,'QQ start cfc params'
         CALL CFC_PARAM( myThid )
         print*,'QQ read cfc atmos data'
         CALL CFC_ATMOS( myThid )
      ENDIF
#endif /* ALLOW_CFC */
#ifdef ALLOW_DIC
      IF ( useDIC ) THEN
         print*,'QQ load dic parameters, initial fixed'
         CALL DIC_INIT_FIXED( myThid )
      ENDIF
#endif /* ALLOW_DIC */
#ifdef ALLOW_BLING
      IF ( useBLING ) THEN
         CALL BLING_INIT_FIXED( myThid )
      ENDIF
#endif /* ALLOW_BLING */
c         _END_MASTER( myThid )

#ifdef ALLOW_DARWIN
      IF ( useDARWIN ) THEN
         CALL DARWIN_INIT_FIXED(myThid)
      ENDIF
#endif

#ifdef ALLOW_DIAGNOSTICS
C     Define diagnostics Names :
      IF ( useDiagnostics ) THEN
        CALL GCHEM_DIAGNOSTICS_INIT( myThid )

C     ***  Edited by Taylor Shropshire ***
C     New Diagnostics
C      diagName  = 'sp_npp     '
C      diagTitle = 'Small Phytoplankton Net Primary Production'
C      diagUnits = 'mmol N/m3 per day'
C      diagCode  = 'SMR     MR     '
C      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C      diagName  = 'lp_npp     '
C      diagTitle = 'Large Phytoplankton Net Primary Production'
C      diagUnits = 'mmol N/m3 per day'
C      diagCode  = 'SMR     MR     '
C      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'phy_chl     '
      diagTitle = 'Chlorophyll-a Concentration'
      diagUnits = 'mg Chl-a/m3'
      diagCode  = 'SMR     MR     '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C      diagName  = 'np_frac     '
C      diagTitle = 'New Production Fraction (f-ratio)'
C      diagUnits = 'Non-Dimensional Number'
C      diagCode  = 'SMR     MR     '
C      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C      diagName  = 'sz_graze     '
C      diagTitle = 'Microzooplankton Total Grazing'
C      diagUnits = 'mmol N/m3 per day'
C      diagCode  = 'SMR     MR     '
C      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C      diagName  = 'lz_graze     '
C      diagTitle = 'Mesozooplankton Total Grazing'
C      diagUnits = 'mmol N/m3 per day'
C      diagCode  = 'SMR     MR     '
C      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C      diagName  = 'pz_graze     '
C      diagTitle = 'Predatory Zooplankton Total Grazing'
C      diagUnits = 'mmol N/m3 per day'
C      diagCode  = 'SMR     MR     '
C      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'l_lim_sp_field     '
      diagTitle = 'Small Phytoplankton Light Limitation'
      diagUnits = 'mmol N/m3 per day'
      diagCode  = 'SMR     MR     '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'l_lim_lp_field     '
      diagTitle = 'Large Phytoplankton Light Limitation'
      diagUnits = 'mmol N/m3 per day'
      diagCode  = 'SMR     MR     '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C      diagName  = 'NO_lim_sp_field     '
C      diagTitle = 'Small Phytoplankton Nitrate Limitation'
C      diagUnits = 'mmol N/m3 per day'
C      diagCode  = 'SMR     MR     '
C      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C      diagName  = 'NH_lim_sp_field     '
C      diagTitle = 'Small Phytoplankton Ammonium Limitation'
C      diagUnits = 'mmol N/m3 per day'
C      diagCode  = 'SMR     MR     '
C      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C      diagName  = 'NO_lim_lp_field     '
C      diagTitle = 'Large Phytoplankton Nitrate Limitation'
C      diagUnits = 'mmol N/m3 per day'
C      diagCode  = 'SMR     MR     '
C      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C      diagName  = 'NH_lim_lp_field     '
C      diagTitle = 'Large Phytolankton Ammonium Limitation'
C      diagUnits = 'mmol N/m3 per day'
C      diagCode  = 'SMR     MR     '
C      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C      diagName  = 'SI_lim_lp_field     '
C      diagTitle = 'Large Phytoplankton Silicate Limitation'
C      diagUnits = 'mmol N/m3 per day'
C      diagCode  = 'SMR     MR     '
C      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
C     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C       End New Diagnostics
C       ******************************

      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_GCHEM */

      RETURN
      END
