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
#include "CARDIAC_MECHANICS/GUCCIONE_ACTIVE.f"
#include "CARDIAC_MECHANICS/ABAQUS_VOLUMETRIC_PENALTY.f"

!                                        -------------------------------
! ---------------------------------------- UANISOHYPERINV_3D_HOLZAPFEL -
!                                        -------------------------------

SUBROUTINE UANISOHYPERINV_3D_HOLZAPFEL (AINV, UA, ZETA, NFIBERS, NINV, UI1, UI2, UI3, TEMP, NOEL, CMNAME, INCMPFLAG, IHYBFLAG, NUMSTATEV, STATEV, NUMFIELDV, FIELDV, FIELDVINC, NUMPROPS, PROPS)
    IMPLICIT NONE

! -------------------------------------------------- INPUT VARIABLES ---

!   BASIC VARIABLES
    CHARACTER*80     CMNAME
    INTEGER          NFIBERS, NINV, NOEL, INCMPFLAG, IHYBFLAG, NUMSTATEV, NUMFIELDV, NUMPROPS
    DOUBLE PRECISION AINV(NINV), UA(2), ZETA(NFIBERS*(NFIBERS-1)/2), UI1(NINV), UI2(NINV*(NINV+1)/2), UI3(NINV*(NINV+1)/2), TEMP, STATEV(NUMSTATEV), FIELDV(NUMFIELDV), FIELDVINC(NUMFIELDV), PROPS(NUMPROPS)

!   PASSIVE LAW VARIABLES
    LOGICAL          UNIL
    DOUBLE PRECISION A, B, AF, BF, D0, I1BAR, AJ, AJE, DAJE_DAJ, DDAJE_DDAJ, I4FFBAR, UJE, DUJE_DAJE, DDUJE_DDAJE

!   ACTIVE LAW VARIABLES
    DOUBLE PRECISION TA_T, TA_T0, TA_TR, TA_TRM, TA_TRB, TA_W, TA_L, TA_L0, TA_LR, TA_B, TA_TMAX, TA, DTA_DTAL, DTAL_DI4, TA_DT, POS_PART, HEAVISIDE

!   GROWTH LAW VARIABLES
    INTEGER          LAMBDAG_NUM_IT
    DOUBLE PRECISION LAMBDAG_X, LAMBDAG_OLD, LAMBDAG_DT, DLAMBDAG, LAMBDAG_TAU, LAMBDAG_XMIN, LAMBDAG_XMAX, LAMBDAG_GAMMA, LAMBDAG, LAMBDAG_DOT, LAMBDAG_T, LAMBDAG_ERR, LAMBDAG_MIN, LAMBDAG_MAX, LAMBDAG_RELAX, LAMBDAG_RES, LAMBDAG_RES_OLD

! --------------------------------------------------- INITIALIZATION ---

!   FLAGS
!     PRINT *, "IHYBFLAG = ", IHYBFLAG
!     PRINT *, "INCMPFLAG = ", INCMPFLAG

!   UNILATERAL CONDITION
    UNIL = .TRUE.
!     UNIL = .FALSE.

!   PASSIVE LAW VARIABLES
    A  = PROPS(1)
    B  = PROPS(2)
    AF = PROPS(3)
    BF = PROPS(4)
    D0 = PROPS(5)

    I1BAR   = AINV(1)
    I4FFBAR = AINV(4)
    AJ      = AINV(3)
!     PRINT *, "I1BAR = ", I1BAR
!     PRINT *, "I4FFBAR = ", I4FFBAR
!     PRINT *, "AJ = ", AJ

!   INITIALIZE ACTIVE LAW VARIABLES
    TA_TMAX = PROPS( 9)
    TA_B    = PROPS(10)
    TA_L0   = PROPS(11)
    TA_LR   = PROPS(12)
    TA_T0   = PROPS(13)
    TA_TRM  = PROPS(14)
    TA_TRB  = PROPS(15)

    TA_T  = FIELDV(1)
    TA_DT = FIELDVINC(1)
!     PRINT *, "TA_T = ", TA_T
!     PRINT *, "TA_DT = ", TA_DT

!   GROWTH LAW VARIABLES
    LAMBDAG_TAU   = PROPS(17)
    LAMBDAG_XMIN  = PROPS(18)
    LAMBDAG_XMAX  = PROPS(19)
    LAMBDAG_GAMMA = PROPS(20)
    LAMBDAG_MIN   = PROPS(21)
    LAMBDAG_MAX   = PROPS(22)

    LAMBDAG_X = STATEV(4)
    LAMBDAG   = STATEV(5)
!     PRINT *, "LAMBDAG_X = ", LAMBDAG_X
!     PRINT *, "LAMBDAG = ", LAMBDAG

    LAMBDAG_T  = FIELDV(2)
    LAMBDAG_DT = FIELDVINC(2)
!     PRINT *, "LAMBDAG_T = ", LAMBDAG_T
!     PRINT *, "LAMBDAG_DT = ", LAMBDAG_DT

! ------------------------------------------------------ GROWTH LOOP ---

    LAMBDAG_OLD    = LAMBDAG
    LAMBDAG_DOT    = 0.
    LAMBDAG_NUM_IT = 0
    DO
!         PRINT *, "LAMBDAG_NUM_IT = ", LAMBDAG_NUM_IT

!       LAMBDAG
        LAMBDAG = LAMBDAG_OLD + LAMBDAG_DOT*LAMBDAG_DT
!         PRINT *, "LAMBDAG = ", LAMBDAG

!       AJE
        AJE = AJ/(1.+LAMBDAG)**3
        DAJE_DAJ = 1./(1.+LAMBDAG)**3
        DDAJE_DDAJ = 0.
!         PRINT *, "AJE = ", AJE

!       DUJE_DAJE
        IF (INCMPFLAG .EQ. 0) THEN
            CALL ABAQUS_VOLUMETRIC_PENALTY(D0, AJE, UJE, DUJE_DAJE, DDUJE_DDAJE)
        ELSE
            UJE = 0.
            DUJE_DAJE = 0.
            DDUJE_DDAJE = 0.
        ENDIF
!         PRINT *, "DUJE_DAJE = ", DUJE_DAJE
!         PRINT *, "DDUJE_DDAJE = ", DDUJE_DDAJE

!       LAMBDAG_X
        LAMBDAG_X = I1BAR * AJE**(2./3)
!         LAMBDAG_X = SQRT(I4FFBAR)
!         LAMBDAG_X = 3 * AJE / (1+LAMBDAG) * DUJE_DAJE
!         PRINT *, "LAMBDAG_X = ", LAMBDAG_X

!       LAMBDAG_RES
        LAMBDAG_RES_OLD = LAMBDAG_RES
        IF (LAMBDAG_TAU .EQ. 0) THEN
            EXIT
        ENDIF
        LAMBDAG_RES = 1./LAMBDAG_TAU
        IF (LAMBDAG_X > LAMBDAG_XMAX) THEN
            LAMBDAG_RES = LAMBDAG_RES * (LAMBDAG_X-LAMBDAG_XMAX)**LAMBDAG_GAMMA
            IF ((LAMBDAG > 0.) .AND. (LAMBDAG < LAMBDAG_MAX)) THEN
                LAMBDAG_RES = LAMBDAG_RES * (1. - (LAMBDAG/LAMBDAG_MAX)**LAMBDAG_GAMMA)
            ELSEIF (LAMBDAG .GE. LAMBDAG_MAX) THEN
                LAMBDAG_RES = 0.
            ENDIF
        ELSEIF (LAMBDAG_X < LAMBDAG_XMIN) THEN
            LAMBDAG_RES = - LAMBDAG_RES * (LAMBDAG_XMIN-LAMBDAG_X)**LAMBDAG_GAMMA
            IF ((LAMBDAG < 0.) .AND. (LAMBDAG > LAMBDAG_MIN)) THEN
                LAMBDAG_RES = LAMBDAG_RES * (1. - (LAMBDAG/LAMBDAG_MIN)**LAMBDAG_GAMMA)
            ELSEIF (LAMBDAG <= LAMBDAG_MIN) THEN
                LAMBDAG_RES = 0.
            ENDIF
        ELSE
            LAMBDAG_RES = 0.
        ENDIF
        LAMBDAG_RES = LAMBDAG_DOT - LAMBDAG_RES
!         PRINT *, "LAMBDAG_RES = ", LAMBDAG_RES

!       EXIT TEST
        LAMBDAG_ERR = ABS(LAMBDAG_RES)*LAMBDAG_DT
        IF (LAMBDAG_ERR < 1E-6) THEN
            EXIT
        ENDIF

        IF (LAMBDAG_NUM_IT > 100) THEN
            PRINT *, "TOO MANY GROWTH ITERATION…"
            EXIT
        ENDIF

!       RELAXATION
        CALL AITKEN(LAMBDAG_RELAX, LAMBDAG_NUM_IT, LAMBDAG_RES, LAMBDAG_RES_OLD)

!       LAMBDAG_DOT
        LAMBDAG_DOT = LAMBDAG_DOT - LAMBDAG_RELAX * LAMBDAG_RES
!         PRINT *, "LAMBDAG_DOT = ", LAMBDAG_DOT

!       COUNTER
        LAMBDAG_NUM_IT = LAMBDAG_NUM_IT + 1
    ENDDO

! ----------------------------------------------- ACTIVE CONTRACTION ---

    TA_L = TA_LR * SQRT(I4FFBAR)
    DTAL_DI4 = TA_LR / 2 / SQRT(I4FFBAR)
    CALL GUCCIONE_ACTIVE(TA, DTA_DTAL, TA_T, TA_T0, TA_L, TA_L0, TA_LR, TA_TMAX, TA_B, TA_TRM, TA_TRB)

! --------------------------------------------------- ENERGY DENSITY ---

    UA(1) = (A/2/B) * EXP(B*(I1BAR-3.))
    IF ((.NOT. UNIL) .OR. ((UNIL) .AND. (I4FFBAR .GE. 1.))) THEN
        UA(1) = UA(1) + (AF/2/BF) * (EXP(BF*(I4FFBAR-1.)**2) - 1.)
    ENDIF
    IF (INCMPFLAG .EQ. 0) THEN
        UA(1) = UA(1) + UJE
    ENDIF
!     PRINT *, "UA = ", UA

! ------------------------------------------------- FIRST DERIVATIVE ---

    UI1 = 0.
    UI1(1) = (A/2) * EXP(B*(I1BAR-3.))
    IF ((.NOT. UNIL) .OR. ((UNIL) .AND. (I4FFBAR .GE. 1.))) THEN
        UI1(4) = UI1(4) + AF * (I4FFBAR-1.) * EXP(BF*(I4FFBAR-1.)**2)
    ENDIF
    UI1(4) = UI1(4) + TA/2
    IF (INCMPFLAG .EQ. 0) THEN
        UI1(3) = UI1(3) + DUJE_DAJE * DAJE_DAJ
    ENDIF
!     PRINT *, "UI1 = ", UI1

! ------------------------------------------------ SECOND DERIVATIVE ---

    UI2 = 0.
    UI2(1) = (A*B/2) * EXP(B*(I1BAR-3.))
    IF ((.NOT. UNIL) .OR. ((UNIL) .AND. (I4FFBAR .GE. 1.))) THEN
        UI2(10) = UI2(10) + AF * EXP(BF*(I4FFBAR-1.)**2) + 2*AF*BF * (I4FFBAR-1.)**2 * EXP(BF*(I4FFBAR-1.)**2)
    ENDIF
    UI2(10) = UI2(10) + DTA_DTAL/2 * DTAL_DI4
    IF (INCMPFLAG .EQ. 0) THEN
        UI2(6) = UI2(6) + DDUJE_DDAJE * DAJE_DAJ * DAJE_DAJ + DUJE_DAJE * DDAJE_DDAJ
    ENDIF
!     PRINT *, "UI2 = ", UI2

! -------------------------------------------------- STATE VARIABLES ---

    STATEV(1) = AJE
    STATEV(2) = TA_L
    STATEV(3) = TA
    STATEV(4) = LAMBDAG_X
    STATEV(5) = LAMBDAG

! -------------------------------------------------------------- END ---

    RETURN
END SUBROUTINE UANISOHYPERINV_3D_HOLZAPFEL
