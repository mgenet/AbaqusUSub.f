!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                  !!!
!!! Created by Martin Genet, 2008-2016                               !!!
!!!                                                                  !!!
!!! Laboratoire de Mécanique et de Technologie (LMT), Cachan, France !!!
!!! Lawrence Berkeley National Laboratory, California, USA           !!!
!!! University of California at San Francisco, USA                   !!!
!!! Swiss Federal Institute of Technology (ETH), Zurich, Switzerland !!!
!!! École Polytechnique, Palaiseau, France                           !!!
!!!                                                                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef GUCCIONE_ACTIVE_f
#define GUCCIONE_ACTIVE_f

!                                                           ------------
! ----------------------------------------------------------- INCLUDES -
!                                                           ------------

#include "UTILS/VEC_MAT_TOOLS.f"

!                                                    -------------------
! ---------------------------------------------------- GUCCIONE_ACTIVE -
!                                                    -------------------

SUBROUTINE GUCCIONE_ACTIVE(TA, DTA_DL, T, T0, L, L0, LR, TMAX, B, TRM, TRB)
    IMPLICIT NONE

! -------------------------------------------------------- VARIABLES ---

    DOUBLE PRECISION, INTENT(IN) :: T, T0, L, L0, LR, TMAX, B, TRM, TRB

    DOUBLE PRECISION, INTENT(OUT) :: TA, DTA_DL

    LOGICAL          STRETCH_DEPENDENT_FORCE, STRETCH_DEPENDENT_RELAXATION
    DOUBLE PRECISION M_PI, TA1, DTA1_DL, TR, DTR_DL, W, DW_DTR, DW_DL, TA2, DTA2_DW, DTA2_DL, POS_PART, HEAVISIDE

! --------------------------------------------------- INITIALIZATION ---

!     STRETCH_DEPENDENT_FORCE = .TRUE.
    STRETCH_DEPENDENT_FORCE = .FALSE.

!     STRETCH_DEPENDENT_RELAXATION = .TRUE.
    STRETCH_DEPENDENT_RELAXATION = .FALSE.

    M_PI = 4*ATAN(1.)

! ----------------------------------------------------------------------

    IF (STRETCH_DEPENDENT_FORCE) THEN
        TA1 = TMAX * (1. - EXP(-B * POS_PART(L-L0)))
        DTA1_DL = TMAX * B * HEAVISIDE(L-L0) * EXP(-B * POS_PART(L-L0))
    ELSE
        TA1 = TMAX * (1. - EXP(-B * POS_PART(LR-L0)))
        DTA1_DL = 0.
    ENDIF

    IF (STRETCH_DEPENDENT_RELAXATION) THEN
        TR = TRM * L + TRB
        DTR_DL = TRM
    ELSE
        TR = TRM * LR + TRB
        DTR_DL = 0.
    ENDIF

    IF (T < T0) THEN
        W = M_PI * T / T0
        DW_DTR = 0.
    ELSE IF (T-T0 < TR) THEN
        W = M_PI * (1. + (T-T0) / TR)
        DW_DTR = - M_PI * (T-T0) / TR / TR
    ELSE
        W = M_PI * 2
        DW_DTR = 0.
    ENDIF
    DW_DL = DW_DTR * DTR_DL

    TA2 = (1.-COS(W))/2.
    DTA2_DW = SIN(W)/2.
    DTA2_DL = DTA2_DW * DW_DL

    TA = TA1 * TA2
    DTA_DL = DTA1_DL * TA2 + TA1 * DTA2_DL

    RETURN
END SUBROUTINE GUCCIONE_ACTIVE

#endif
