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

!                                                           ------------
! ----------------------------------------------------------- INCLUDES -
!                                                           ------------

#include "UTILS/VEC_MAT_TOOLS.f"

!                                        -------------------------------
! ---------------------------------------- UANISOHYPERSTR_3D_HOLZAPFEL -
!                                        -------------------------------

SUBROUTINE UANISOHYPERSTR_3D_HOLZAPFEL (EBAR, AJ, UA, DU1, DU2, DU3, TEMP, NOEL, CMNAME, INCMPFLAG, IHYBFLAG, NDI, NSHR, NTENS, NUMSTATEV, STATEV, NUMFIELDV, FIELDV, FIELDVINC, NUMPROPS, PROPS)
    IMPLICIT NONE

!   INPUT VARIABLES
    CHARACTER*80     CMNAME
    INTEGER          NOEL, INCMPFLAG, IHYBFLAG, NDI, NSHR, NTENS, NUMSTATEV, NUMFIELDV, NUMPROPS
    DOUBLE PRECISION EBAR(NTENS), AJ, UA(2), DU1(NTENS+1), DU2((NTENS+1)*(NTENS+2)/2), DU3((NTENS+1)*(NTENS+2)/2), TEMP, STATEV(NUMSTATEV), FIELDV(NUMFIELDV), FIELDVINC(NUMFIELDV), PROPS(NUMPROPS)

!   PASSIVE LAW VARIABLES
    LOGICAL          UNIL
    INTEGER          INDX
    DOUBLE PRECISION A, B, AF, BF, D0, CBAR(NTENS), I1, I4F

!   ACTIVE LAW VARIABLES
    DOUBLE PRECISION M_PI, TA_T, TA_T0, TA_TR, TA_TRM, TA_TRB, TA_W, TA_L, TA_L0, TA_LR, TA_B, TA_TMAX, TA_RSS, TA_RNN, TA, DTA, POS_PART, HEAVISIDE

!     PRINT *, "IHYBFLAG = ", IHYBFLAG
!     PRINT *, "INCMPFLAG = ", INCMPFLAG
!     PRINT *, "EBAR = ", EBAR
!     PRINT *, "AJ = ", AJ

!   UNILATERAL CONDITION
    UNIL = .TRUE.
!     UNIL = .FALSE.

!   INITIALIZE PASSIVE LAW VARIABLES
    A  = PROPS(1)
    B  = PROPS(2)
    AF = PROPS(3)
    BF = PROPS(4)
    D0 = PROPS(5)

    CBAR(1) = 2*EBAR(1)+1
    CBAR(2) = 2*EBAR(2)+1
    CBAR(3) = 2*EBAR(3)+1
    CBAR(4) = 2*EBAR(4)
    CBAR(5) = 2*EBAR(5)
    CBAR(6) = 2*EBAR(6)

    I1  = CBAR(1) + CBAR(2) + CBAR(3)
    I4F = CBAR(1)

!   INITIALIZE ACTIVE LAW VARIABLES
    M_PI     = 4*ATAN(1.)
    TA_TMAX = PROPS( 9)
    TA_RSS  = PROPS(10)
    TA_RNN  = PROPS(11)
    TA_B    = PROPS(12)
    TA_L0   = PROPS(13)
    TA_LR   = PROPS(14)
    TA_T0   = PROPS(17)
    TA_TRM  = PROPS(18)
    TA_TRB  = PROPS(19)

!   ENERGY DENSITY
    UA(1) = (A/2/B) * EXP(B*(I1-3.))
    IF ((.NOT.UNIL).OR.((UNIL).AND.(I4F.GT.1.))) THEN
        UA(1) = UA(1) + (AF/2/BF) * (EXP(BF*(I4F-1.)**2) - 1.)
    ENDIF
    IF (INCMPFLAG .EQ. 0) THEN
        UA(1) = UA(1) + (1./D0) * ((AJ**2-1.)/2 - LOG(AJ))
    ENDIF
!     PRINT *, "UA = ", UA

!   FIRST DERIVATIVE
    DU1 = 0.
    DU1(1) = A * EXP(B*(I1-3.))
    IF ((.NOT.UNIL).OR.((UNIL).AND.(I4F .GT. 1.))) THEN
        DU1(1) = DU1(1) + 2 * AF * (I4F-1.) * EXP(BF*(I4F-1.)**2)
    ENDIF
    DU1(2) = A * EXP(B*(I1-3.))
    DU1(3) = A * EXP(B*(I1-3.))
    IF (INCMPFLAG .EQ. 0) THEN
        DU1(7) = (1./D0) * (AJ - 1./AJ)
    ENDIF
!     PRINT *, "DU1 = ", DU1

!   ACTIVE FORCE
    TA_L  = TA_LR * SQRT(CBAR(1))
    TA_TR = TA_TRM * TA_L + TA_TRB
    TA_T = FIELDV(1)
    TA_T = MODULO(TA_T, TA_T0+TA_TR)
    IF (TA_T <= TA_T0) THEN
        TA_W = M_PI * TA_T/TA_T0
    ELSE
        TA_W = M_PI * (1.+(TA_T-TA_T0)/TA_TR)
    ENDIF
    TA = TA_TMAX * (1.-EXP(-TA_B*POS_PART(TA_L-TA_L0))) * (1.-COS(TA_W)) / 2.
    DU1(1) = DU1(1) + TA
    DU1(2) = DU1(2) + TA * TA_RSS
    DU1(3) = DU1(3) + TA * TA_RNN

!   SECOND DERIVATIVE
    DU2 = 0.
    DU2(INDX(1,1)) = (2*A*B) * EXP(B*(I1-3.))
    IF ((.NOT.UNIL).OR.((UNIL).AND.(I4F .GT. 1.))) THEN
        DU2(INDX(1,1)) = DU2(INDX(1,1)) + (4*AF) * EXP(BF*(I4F-1.)**2) + (8*AF*BF) * (I4F-1.)**2 * EXP(BF*(I4F-1.)**2)
    ENDIF
    DU2(INDX(1,2)) = (2*A*B) * EXP(B*(I1-3.))
    DU2(INDX(1,3)) = (2*A*B) * EXP(B*(I1-3.))
    DU2(INDX(2,2)) = (2*A*B) * EXP(B*(I1-3.))
    DU2(INDX(2,3)) = (2*A*B) * EXP(B*(I1-3.))
    DU2(INDX(3,3)) = (2*A*B) * EXP(B*(I1-3.))
    IF (INCMPFLAG .EQ. 0) THEN
        DU2(28) = (1./D0) * (1. + 1./AJ**2)
    ENDIF
!     PRINT *, "DU2 = ", DU2

!   ACTIVE FORCE
    DTA = TA_TMAX * TA_B*HEAVISIDE(TA_L-TA_L0) * TA_LR/SQRT(CBAR(1)) * EXP(-TA_B*POS_PART(TA_L-TA_L0)) * (1.-COS(TA_W)) / 2.
    IF (TA_T >= TA_T0) THEN
        DTA = DTA - TA_TMAX * (1.-EXP(-TA_B*POS_PART(TA_L-TA_L0))) * (SIN(TA_W)/2.) * (M_PI*(TA_T-TA_T0)/TA_TR/TA_TR) * (TA_TRM*TA_LR/SQRT(CBAR(1))) ! should have a /2 here no?
    ENDIF
    DU2(INDX(1,1)) = DU2(INDX(1,1)) + DTA
    DU2(INDX(2,2)) = DU2(INDX(2,2)) + DTA * TA_RSS
    DU2(INDX(3,3)) = DU2(INDX(3,3)) + DTA * TA_RNN

!   STATE VARIABLES
    STATEV(1) = AJ
    STATEV(2) = TA_L
    STATEV(3) = TA

    RETURN
END SUBROUTINE UANISOHYPERSTR_3D_HOLZAPFEL
