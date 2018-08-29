C $Header: /u/gcmpack/MITgcm/pkg/gchem/GCHEM_FIELDS.h,v 1.1 2004/11/28 23:48:31 mlosch Exp $
C $Name:  $

#ifdef ALLOW_GCHEM
CBOP
C    !ROUTINE: GCHEM_FIELDS.h
C    !INTERFACE:
 
C    !DESCRIPTION:
C Contains tracer fields specifically for chemical tracers.
C
C  gchemTendency :: 3DxPTRACER_num field that store the tendencies due
C                   to the bio-geochemical model

#ifndef GCHEM_SEPARATE_FORCING
      _RL gchemTendency(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,
     &                  PTRACERS_num)
      COMMON /GCHEM_FIELDS/ 
     &     gchemTendency,
#endif /* GCHEM_SEPARATE_FORCING */

C      _RL sp_npp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL lp_npp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL sz_graze(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL lz_graze(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL pz_graze(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL phy_chl(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL np_frac(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL l_lim_sp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL l_lim_lp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL NO_lim_sp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL NH_lim_sp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL NO_lim_lp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL NH_lim_lp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL SI_lim_lp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      COMMON /GCHEM_FIELDS/
     &     l_lim_sp_field, l_lim_lp_field, phy_chl


CEOP
#endif /* ALLOW_GCHEM */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
