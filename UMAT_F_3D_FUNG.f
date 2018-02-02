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

!                                                       ----------------
! --------------------------------------------------------- INCLUDES ---
!                                                       ----------------

#include "UTILS/PRINT_UMAT_DATA.f"
#include "UTILS/VEC_MAT_TOOLS.f"
#include "ALGORITHMS/AITKEN.f"
#include "FINITE TRANSFORMATION/PUSH_FORWARD.f"
#include "FINITE TRANSFORMATION/JAUMANN.f"

!                                                         --------------
! ----------------------------------------------------------- SDVINI ---
!                                                         --------------

! SUBROUTINE SDVINI(STATEV, COORDS, NSTATV, NCRDS, NOEL, NPT, LAYER, KSPT)
!     INCLUDE 'ABA_PARAM.INC'
!
!     INTEGER          NSTATV, NCRDS, NOEL, NPT, LAYER, KSPT
!     DOUBLE PRECISION STATEV(NSTATV), COORDS(NCRDS)
!
! !     STATEV(1) = 1. ! THE_G
! !     STATEV(2) = 1. ! DETFG
! !     STATEV(3) = 1. ! DETFE
! !     STATEV(4) = 1. ! DETF
!
! !     PRINT *, "STATEV = ", STATEV
!
!     RETURN
! END SUBROUTINE SDVINI

!                                                     ------------------
! ------------------------------------------------------- IS_GROWING ---
!                                                     ------------------

LOGICAL FUNCTION IS_GROWING(T, LAMBDAG_TINI, LAMBDAG_TFIN)
    IMPLICIT NONE

    DOUBLE PRECISION LAMBDAG_MAX, T, LAMBDAG_TINI, LAMBDAG_TFIN

!     IF ((INT(MOD(T+3-1e-9, 4.)) == 0)) THEN
!         IS_GROWING = .TRUE.
!     ELSE
!         IS_GROWING = .FALSE.
!     ENDIF
! !     PRINT *, "MOD = ", MOD(T+2-1e-9, 3.)
! !     PRINT *, "IS_GROWING = ", IS_GROWING

!     IF ((INT(MOD(T+2-1e-9, 3.)) == 0)) THEN
!         IS_GROWING = .TRUE.
!     ELSE
!         IS_GROWING = .FALSE.
!     ENDIF
! !     PRINT *, "MOD = ", MOD(T+2-1e-9, 3.)
! !     PRINT *, "IS_GROWING = ", IS_GROWING

    IF ((T .GT. LAMBDAG_TINI) .AND. (T .LE. LAMBDAG_TFIN)) THEN
        IS_GROWING = .TRUE.
    ELSE
        IS_GROWING = .FALSE.
    ENDIF

    RETURN
END FUNCTION IS_GROWING

!                                                       ----------------
! ------------------------------------------------------- UMAT_3D_FUNG -
!                                                       ----------------

SUBROUTINE UMAT_3D_FUNG (STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)
    IMPLICIT NONE

! -------------------------------------------------- INPUT VARIABLES ---

    CHARACTER*80     CMNAME
    INTEGER          NDI, NSHR, NTENS, NSTATV, NPROPS, NOEL, NPT, LAYER, KSPT, KSTEP, KINC
    DOUBLE PRECISION  STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS), SSE, SPD, SCD, RPL, DDSDDT(NTENS), DRPLDE(NTENS), DRPLDT, STRAN(NTENS), DSTRAN(NTENS), TIME(2), DTIME, TEMP, DTEMP, PREDEF, DPRED, PROPS(NPROPS), COORDS(NDI), DROT(NDI,NDI), PNEWDT, CELENT, DFGRD0(NDI,NDI), DFGRD1(NDI,NDI)

!     CALL PRINT_UMAT_DATA(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

! -------------------------------------------------- LOCAL VARIABLES ---

    INTEGER          I, J, K, L, IJ, KL, GET_IJ, GET_I, GET_J
    DOUBLE PRECISION C0, B01111, B02222, B03333, B01212, B01313, B02323, T, DT, IDNDI(NDI,NDI), FT(NDI,NDI), DET33, DETFT, CT(NDI,NDI), ET(NDI,NDI), FE(NDI,NDI), DETFE, CE(NDI,NDI), EE(NDI,NDI), EXP_TERM, PK2_VEC(NTENS), PK2_MAT(NDI,NDI), STRESS_MAT(NDI,NDI), H(NTENS,NTENS)

    LOGICAL          IS_GROWING
    INTEGER          LAMBDAG_NUM_IT
    DOUBLE PRECISION LAMBDAG_TINI, LAMBDAG_TFIN, LAMBDAG_XMIN, LAMBDAG_XMAX, LAMBDAG_GAMMA, LAMBDAG_TAU, LAMBDAG_MIN, LAMBDAG_MAX, LAMBDAG_DOT, LAMBDAG_OLD, LAMBDAG_RES, LAMBDAG_RES_OLD, LAMBDAG, LAMBDAG_X, POS_PART, LAMBDAG_ERR, LAMBDAG_RELAX

! --------------------------------------------------- INITIALIZATION ---

    DT = DTIME
    T  = TIME(2)+DT

    CALL SETID33(IDNDI)
    FT    = DFGRD1
    DETFT = DET33(FT)
    CT    = MATMUL(TRANSPOSE(FT), FT)
    ET    = (CT - IDNDI)/2

    C0     = PROPS(1)
    B01111 = PROPS(2)
    B02222 = PROPS(3)
    B03333 = PROPS(4)
    B01212 = PROPS(5)
    B01313 = PROPS(6)
    B02323 = PROPS(7)

    LAMBDAG_TINI  = PROPS( 9)
    LAMBDAG_TFIN  = PROPS(10)
    LAMBDAG_XMIN  = PROPS(11)
    LAMBDAG_XMAX  = PROPS(12)
    LAMBDAG_GAMMA = PROPS(13)
    LAMBDAG_TAU   = PROPS(14)
    LAMBDAG_MIN   = PROPS(15)
    LAMBDAG_MAX   = PROPS(16)

    LAMBDAG_X = STATEV(1)
    LAMBDAG   = STATEV(2)

! ------------------------------------------------------ GROWTH LOOP ---

    LAMBDAG_OLD    = LAMBDAG
    LAMBDAG_DOT    = 0.
    LAMBDAG_NUM_IT = 0
    DO
!       LAMBDAG
        LAMBDAG = LAMBDAG_OLD + LAMBDAG_DOT*DT
!         PRINT *, "LAMBDAG = ", LAMBDAG

!       FE, CE, EE
        FE = FT / (1.+LAMBDAG)
        CE = MATMUL(TRANSPOSE(FE), FE)
        EE = (CE - IDNDI)/2
!         PRINT *, "FE = ", FE
!         PRINT *, "CE = ", CE
!         PRINT *, "EE = ", EE

!         STRAN = STRAN + DSTRAN
!         PRINT *, "STRAN = ", STRAN

!       EXP_TERM
        EXP_TERM = B01111*EE(1,1)**2 &
                 + B02222*EE(2,2)**2 &
                 + B03333*EE(3,3)**2 &
                 + B01212*EE(1,2)**2 &
                 + B01212*EE(1,2)*EE(2,1) &
                 + B01212*EE(2,1)*EE(1,2) &
                 + B01212*EE(2,1)**2 &
                 + B01313*EE(1,3)**2 &
                 + B01313*EE(1,3)*EE(3,1) &
                 + B01313*EE(3,1)*EE(1,3) &
                 + B01313*EE(3,1)**2 &
                 + B02323*EE(2,3)**2 &
                 + B02323*EE(2,3)*EE(3,2) &
                 + B02323*EE(3,2)*EE(2,3) &
                 + B02323*EE(3,2)**2
!         PRINT *, "EXP_TERM = ", EXP_TERM
        IF (EXP_TERM .GT. 700) THEN
            PRINT *, "STEP ", KSTEP, ", INCREMENT ", KINC, ", ELEMENT ", NOEL, ", INTEGRATION POINT ", NPT, ": TOO LARGE DEFORMATION…"
            EXP_TERM = 700
            PNEWDT = 0.5
            EXIT
        ENDIF
        EXP_TERM = EXP(EXP_TERM)
!         PRINT *, "EXP_TERM = ", EXP_TERM

!       LAMBDAG_X
!         LAMBDAG_X = SQRT(CE(1,1))
        LAMBDAG_X = (1./(1.+LAMBDAG)**3) * C0 * EXP_TERM * (  B01111*EE(1,1)*CT(1,1) &
                                                           +  B02222*EE(2,2)*CT(2,2) &
                                                           +  B03333*EE(3,3)*CT(3,3) &
                                                           +4*B01212*EE(1,2)*CT(1,2) &
                                                           +4*B01313*EE(1,3)*CT(1,3) &
                                                           +4*B02323*EE(2,3)*CT(2,3))

!       LAMBDAG_RES
        IF (IS_GROWING(T, LAMBDAG_TINI, LAMBDAG_TFIN)) THEN
            LAMBDAG_RES_OLD = LAMBDAG_RES
            LAMBDAG_RES = 1./LAMBDAG_TAU
            IF (LAMBDAG_X .GT. LAMBDAG_XMAX) THEN
                LAMBDAG_RES = LAMBDAG_RES * (LAMBDAG_X-LAMBDAG_XMAX)**LAMBDAG_GAMMA
                IF ((LAMBDAG .GT. 0.) .AND. (LAMBDAG .LT. LAMBDAG_MAX)) THEN
                    LAMBDAG_RES = LAMBDAG_RES * (1. - (LAMBDAG/LAMBDAG_MAX)**LAMBDAG_GAMMA)
                ELSEIF (LAMBDAG .GE. LAMBDAG_MAX) THEN
                    LAMBDAG_RES = 0.
                ENDIF
            ELSEIF (LAMBDAG_X .LT. LAMBDAG_XMIN) THEN
                LAMBDAG_RES = - LAMBDAG_RES * (LAMBDAG_XMIN-LAMBDAG_X)**LAMBDAG_GAMMA
                IF ((LAMBDAG .LT. 0.) .AND. (LAMBDAG .GT. LAMBDAG_MIN)) THEN
                    LAMBDAG_RES = LAMBDAG_RES * (1. - (LAMBDAG/LAMBDAG_MIN)**LAMBDAG_GAMMA)
                ELSEIF (LAMBDAG .LE. LAMBDAG_MIN) THEN
                    LAMBDAG_RES = 0.
                ENDIF
            ELSE
                LAMBDAG_RES = 0.
            ENDIF
            LAMBDAG_RES = LAMBDAG_DOT - LAMBDAG_RES
        ELSE
            EXIT
        ENDIF

!       EXIT TEST
        LAMBDAG_ERR = ABS(LAMBDAG_RES)*DT
        IF (LAMBDAG_ERR .LT. 1E-6) THEN
            EXIT
        ENDIF

        IF (LAMBDAG_NUM_IT .GT. 100) THEN
            PRINT *, "STEP ", KSTEP, ", INCREMENT ", KINC, ", ELEMENT ", NOEL, ", INTEGRATION POINT ", NPT, ": TOO MANY GROWTH ITERATION…"
            PNEWDT = 0.5
            EXIT
        ENDIF

!       RELAXATION
        CALL AITKEN(LAMBDAG_RELAX, LAMBDAG_NUM_IT, LAMBDAG_RES, LAMBDAG_RES_OLD)

!       LAMBDAG_DOT
        LAMBDAG_DOT = LAMBDAG_DOT - LAMBDAG_RELAX * LAMBDAG_RES

!       COUNTER
        LAMBDAG_NUM_IT = LAMBDAG_NUM_IT + 1
    ENDDO

! ----------------------------------------------------------- STRESS ---

!   SECOND PIOLA KIRSCHOFF STRESS
    PK2_VEC(1) = C0 * EXP_TERM *   B01111*EE(1,1)
    PK2_VEC(2) = C0 * EXP_TERM *   B02222*EE(2,2)
    PK2_VEC(3) = C0 * EXP_TERM *   B03333*EE(3,3)
    PK2_VEC(4) = C0 * EXP_TERM * 2*B01212*EE(1,2)
    PK2_VEC(5) = C0 * EXP_TERM * 2*B01313*EE(1,3)
    PK2_VEC(6) = C0 * EXP_TERM * 2*B02323*EE(2,3)
!     PRINT *, "PK2_VEC = ", PK2_VEC

!   CAUCHY STRESS
    CALL MATSYM33_FROM_VECCOL6(PK2_MAT, PK2_VEC)
    CALL CAUCHYSTRESS33_FROM_PK2STRESS33(STRESS_MAT, PK2_MAT, FT)
    CALL VECCOL6_FROM_MATSYM33(STRESS, STRESS_MAT)
!     PRINT *, "STRESS = ", STRESS

! --------------------------------------------------------- JACOBIAN ---

!   TANGENT IN TERMS OF SECOND PIOLA KIRSCHOFF STRESS
    H = 0.
    H(1,1) = C0 * EXP_TERM * (B01111+2*B01111*EE(1,1)*B01111*EE(1,1))
    H(1,2) = C0 * EXP_TERM * (        2*B01111*EE(1,1)*B02222*EE(2,2))
    H(1,3) = C0 * EXP_TERM * (        2*B01111*EE(1,1)*B03333*EE(3,3))
    H(1,4) = C0 * EXP_TERM * (        4*B01111*EE(1,1)*B01212*EE(1,2))
    H(1,5) = C0 * EXP_TERM * (        4*B01111*EE(1,1)*B01313*EE(1,3))
    H(1,6) = C0 * EXP_TERM * (        4*B01111*EE(1,1)*B02323*EE(2,3))
    H(2,1) = C0 * EXP_TERM * (        2*B02222*EE(2,2)*B01111*EE(1,1))
    H(2,2) = C0 * EXP_TERM * (B02222+2*B02222*EE(2,2)*B02222*EE(2,2))
    H(2,3) = C0 * EXP_TERM * (        2*B02222*EE(2,2)*B03333*EE(3,3))
    H(2,4) = C0 * EXP_TERM * (        4*B02222*EE(2,2)*B01212*EE(1,2))
    H(2,5) = C0 * EXP_TERM * (        4*B02222*EE(2,2)*B01313*EE(1,3))
    H(2,6) = C0 * EXP_TERM * (        4*B02222*EE(2,2)*B02323*EE(2,3))
    H(3,1) = C0 * EXP_TERM * (        2*B03333*EE(3,3)*B01111*EE(1,1))
    H(3,2) = C0 * EXP_TERM * (        2*B03333*EE(3,3)*B02222*EE(2,2))
    H(3,3) = C0 * EXP_TERM * (B03333+2*B03333*EE(3,3)*B03333*EE(3,3))
    H(3,4) = C0 * EXP_TERM * (        4*B03333*EE(3,3)*B01212*EE(1,2))
    H(3,5) = C0 * EXP_TERM * (        4*B03333*EE(3,3)*B01313*EE(1,3))
    H(3,6) = C0 * EXP_TERM * (        4*B03333*EE(3,3)*B02323*EE(2,3))
    H(4,1) = C0 * EXP_TERM * (        4*B01212*EE(1,2)*B01111*EE(1,1))
    H(4,2) = C0 * EXP_TERM * (        4*B01212*EE(1,2)*B02222*EE(2,2))
    H(4,3) = C0 * EXP_TERM * (        4*B01212*EE(1,2)*B03333*EE(3,3))
    H(4,4) = C0 * EXP_TERM * (B01212+8*B01212*EE(1,2)*B01212*EE(1,2))
    H(4,5) = C0 * EXP_TERM * (        8*B01212*EE(1,2)*B01313*EE(1,3))
    H(4,6) = C0 * EXP_TERM * (        8*B01212*EE(1,2)*B02323*EE(2,3))
    H(5,1) = C0 * EXP_TERM * (        4*B01313*EE(1,3)*B01111*EE(1,1))
    H(5,2) = C0 * EXP_TERM * (        4*B01313*EE(1,3)*B02222*EE(2,2))
    H(5,3) = C0 * EXP_TERM * (        4*B01313*EE(1,3)*B03333*EE(3,3))
    H(5,4) = C0 * EXP_TERM * (        8*B01313*EE(1,3)*B01212*EE(1,2))
    H(5,5) = C0 * EXP_TERM * (B01313+8*B01313*EE(1,3)*B01313*EE(1,3))
    H(5,6) = C0 * EXP_TERM * (        8*B01313*EE(1,3)*B02323*EE(2,3))
    H(6,1) = C0 * EXP_TERM * (        4*B02323*EE(2,3)*B01111*EE(1,1))
    H(6,2) = C0 * EXP_TERM * (        4*B02323*EE(2,3)*B02222*EE(2,2))
    H(6,3) = C0 * EXP_TERM * (        4*B02323*EE(2,3)*B03333*EE(3,3))
    H(6,4) = C0 * EXP_TERM * (        8*B02323*EE(2,3)*B01212*EE(1,2))
    H(6,5) = C0 * EXP_TERM * (        8*B02323*EE(2,3)*B01313*EE(1,3))
    H(6,6) = C0 * EXP_TERM * (B02323+8*B02323*EE(2,3)*B02323*EE(2,3))
!     PRINT *, "H = ", H

!   GROWTH TERMS IN TANGENT
    H = (1./(1.+LAMBDAG)**4) * H

!   TANGENT IN TERMS OF CAUCHY STRESS
    CALL CAUCHYJACOBIAN66_FROM_PK2JACOBIAN66(DDSDDE, H, FT)
!     PRINT *, "DDSDDE = ", DDSDDE

!   JAUMANN TERMS IN TANGENT
    CALL JAUMANN(DDSDDE, STRESS)
!     PRINT *, "DDSDDE = ", DDSDDE

!     DO IJ = 4,6
!     DO KL = 1,6
!         DDSDDE(IJ,KL) = DDSDDE(IJ,KL) / SQRT(2.)
!         DDSDDE(KL,IJ) = DDSDDE(KL,IJ) / SQRT(2.)
!     END DO
!     END DO

! -------------------------------------------------- STATE VARIABLES ---

    STATEV(1) = LAMBDAG_X
    STATEV(2) = LAMBDAG

! -------------------------------------------------------------- END ---

    RETURN
END SUBROUTINE UMAT_3D_FUNG
