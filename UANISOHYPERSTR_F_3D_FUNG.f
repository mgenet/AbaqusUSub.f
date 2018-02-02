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
#include "CARDIAC_MECHANICS/FUNG_PASSIVE.f"
#include "CARDIAC_MECHANICS/GUCCIONE_ACTIVE.f"
#include "CARDIAC_MECHANICS/ABAQUS_VOLUMETRIC_PENALTY.f"

!                                             --------------------------
! --------------------------------------------- UANISOHYPERSTR_3D_FUNG -
!                                             --------------------------

SUBROUTINE UANISOHYPERSTR_3D_FUNG (EBAR, AJ, UA, DU1, DU2, DU3, TEMP, NOEL, CMNAME, INCMPFLAG, IHYBFLAG, NDI, NSHR, NTENS, NUMSTATEV, STATEV, NUMFIELDV, FIELDV, FIELDVINC, NUMPROPS, PROPS)
    IMPLICIT NONE

!   INPUT VARIABLES
    CHARACTER*80     CMNAME
    INTEGER          NOEL, INCMPFLAG, IHYBFLAG, NDI, NSHR, NTENS, NUMSTATEV, NUMFIELDV, NUMPROPS
    DOUBLE PRECISION EBAR(NTENS), AJ, UA(2), DU1(NTENS+1), DU2((NTENS+1)*(NTENS+2)/2), DU3((NTENS+1)*(NTENS+2)/2), TEMP, STATEV(NUMSTATEV), FIELDV(NUMFIELDV), FIELDVINC(NUMFIELDV), PROPS(NUMPROPS)

!   PASSIVE LAW VARIABLES
    LOGICAL          TOO_LARGE_DEFORMATION
    INTEGER          INDX
    DOUBLE PRECISION C0, B0(NTENS), D0, EXP_TERM, UE, DUE_DE(NTENS), DDUE_DDE(NTENS*(NTENS+1)/2), UJE, DUJE_DJE, DDUJE_DDJE

!   ACTIVE LAW VARIABLES
    DOUBLE PRECISION TA_T, TA_T0, TA_TRM, TA_TRB, TA_TR, TA_W, TA_L, DTAL_DE, TA_L0, TA_LR, TA_B, TA_TMAX, TA_RFF, TA_RSS, TA_RNN, TA, DTA_DTAL

!   GROWTH LAW VARIABLES
    DOUBLE PRECISION LG_DT, LG, DLG, LG_TAU, LG_X, LG_X0, LG_X1, LG_GAMMA, JG, JE, DJE_DAJ, DDJE_DDAJ, POS_PART, HEAVISIDE

!     PRINT *, "IHYBFLAG = ", IHYBFLAG
!     PRINT *, "INCMPFLAG = ", INCMPFLAG
!     PRINT *, "EBAR = ", EBAR
!     PRINT *, "AJ = ", AJ

!   INITIALIZE PASSIVE LAW VARIABLES
    C0    = PROPS(1)
    B0(1) = PROPS(2)
    B0(2) = PROPS(3)
    B0(3) = PROPS(4)
    B0(4) = PROPS(5)
    B0(5) = PROPS(6)
    B0(6) = PROPS(7)
    D0    = PROPS(8)
!     PRINT *, "C0 = ", C0
!     PRINT *, "B0 = ", B0
!     PRINT *, "D0 = ", D0

!   INITIALIZE ACTIVE LAW VARIABLES
    TA_TMAX = PROPS( 9)
    TA_RFF  = PROPS(10)
    TA_RSS  = PROPS(11)
    TA_RNN  = PROPS(12)
    TA_B    = PROPS(13)
    TA_L0   = PROPS(14)
    TA_LR   = PROPS(15)

    TA_T0   = PROPS(17)
    TA_TRM  = PROPS(18)
    TA_TRB  = PROPS(19)

!   INITIALIZE GROWTH LAW VARIABLES
    LG_TAU   = PROPS(25)
    LG_X0    = PROPS(26)
    LG_X1    = PROPS(27)
    LG_GAMMA = PROPS(28)

!   GROWTH
    LG_DT = FIELDVINC(2)
    LG    = STATEV(5)
    LG_X  = SQRT(2.*EBAR(1)+1.) ! Ça devrait être E, là, pas EBAR, n'est-ce pas?
    IF (LG_TAU .NE. 0) THEN
        DLG = (LG_DT/LG_TAU) * (POS_PART(LG_X-LG_X1)**LG_GAMMA - POS_PART(LG_X0-LG_X)**LG_GAMMA)
    ELSE
        DLG = 0.
    ENDIF
    LG = LG + DLG
    JG = (1.+LG)**3
    JE = AJ/JG
    DJE_DAJ = 1./JG
    DDJE_DDAJ = 0.

!   HYPERELASTICITY
    CALL FUNG_PASSIVE(C0, B0, EBAR, UE, DUE_DE, DDUE_DDE, TOO_LARGE_DEFORMATION)

!   ACTIVE CONTRACTION
    TA_T = FIELDV(1)
    TA_L = TA_LR * SQRT(2.*EBAR(1)+1.)
    DTAL_DE = TA_LR / SQRT(2.*EBAR(1)+1.)
    CALL GUCCIONE_ACTIVE(TA, DTA_DTAL, TA_T, TA_T0, TA_L, TA_L0, TA_LR, TA_TMAX, TA_B, TA_TRM, TA_TRB)

!   VOLUMETRIC PENALTY
    IF (INCMPFLAG .EQ. 0) THEN
        CALL ABAQUS_VOLUMETRIC_PENALTY(D0, JE, UJE, DUJE_DJE, DDUJE_DDJE)
    ENDIF

!   ENERGY DENSITY
    UA(1) = UE
    IF (INCMPFLAG .EQ. 0) THEN
        UA(1) = UA(1) + UJE
    ENDIF
!     PRINT *, "UA = ", UA

!   FIRST DERIVATIVE
    DU1 = 0.
    DU1(1) = DUE_DE(1)
    DU1(2) = DUE_DE(2)
    DU1(3) = DUE_DE(3)
    DU1(4) = DUE_DE(4)
    DU1(5) = DUE_DE(5)
    DU1(6) = DUE_DE(6)
    DU1(1) = DU1(1) + TA * TA_RFF
    DU1(2) = DU1(2) + TA * TA_RSS
    DU1(3) = DU1(3) + TA * TA_RNN
    IF (INCMPFLAG .EQ. 0) THEN
        DU1(7) = DU1(7) + DUJE_DJE * DJE_DAJ
    ENDIF
!     PRINT *, "DU1 = ", DU1

!   SECOND DERIVATIVE
    DU2 = 0.
    DU2(INDX(1,1)) = DDUE_DDE(INDX(1,1))
    DU2(INDX(1,2)) = DDUE_DDE(INDX(1,2))
    DU2(INDX(1,3)) = DDUE_DDE(INDX(1,3))
    DU2(INDX(1,4)) = DDUE_DDE(INDX(1,4))
    DU2(INDX(1,5)) = DDUE_DDE(INDX(1,5))
    DU2(INDX(1,6)) = DDUE_DDE(INDX(1,6))
    DU2(INDX(2,2)) = DDUE_DDE(INDX(2,2))
    DU2(INDX(2,3)) = DDUE_DDE(INDX(2,3))
    DU2(INDX(2,4)) = DDUE_DDE(INDX(2,4))
    DU2(INDX(2,5)) = DDUE_DDE(INDX(2,5))
    DU2(INDX(2,6)) = DDUE_DDE(INDX(2,6))
    DU2(INDX(3,3)) = DDUE_DDE(INDX(3,3))
    DU2(INDX(3,4)) = DDUE_DDE(INDX(3,4))
    DU2(INDX(3,5)) = DDUE_DDE(INDX(3,5))
    DU2(INDX(3,6)) = DDUE_DDE(INDX(3,6))
    DU2(INDX(4,4)) = DDUE_DDE(INDX(4,4))
    DU2(INDX(4,5)) = DDUE_DDE(INDX(4,5))
    DU2(INDX(4,6)) = DDUE_DDE(INDX(4,6))
    DU2(INDX(5,5)) = DDUE_DDE(INDX(5,5))
    DU2(INDX(5,6)) = DDUE_DDE(INDX(5,6))
    DU2(INDX(6,6)) = DDUE_DDE(INDX(6,6))
    DU2(INDX(1,1)) = DU2(INDX(1,1)) + DTA_DTAL * DTAL_DE * TA_RFF
    DU2(INDX(1,2)) = DU2(INDX(1,2)) + DTA_DTAL * DTAL_DE * TA_RSS / 2
    DU2(INDX(1,3)) = DU2(INDX(1,3)) + DTA_DTAL * DTAL_DE * TA_RNN / 2
    IF (INCMPFLAG .EQ. 0) THEN
        DU2(28) = DU2(28) + DDUJE_DDJE * DJE_DAJ * DJE_DAJ + DUJE_DJE * DDJE_DDAJ
    ENDIF
!     PRINT *, "DU2 = ", DU2

!   STATE VARIABLES
    STATEV(1) = JE
    STATEV(2) = TA_L
    STATEV(3) = TA
    STATEV(4) = LG_X
    STATEV(5) = LG

!     STATEV( 6) = EBAR(1)
!     STATEV( 7) = EBAR(2)
!     STATEV( 8) = EBAR(3)
!     STATEV( 9) = EBAR(4)
!     STATEV(10) = EBAR(5)
!     STATEV(11) = EBAR(6)

    RETURN
END SUBROUTINE UANISOHYPERSTR_3D_FUNG
