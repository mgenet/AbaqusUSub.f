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
#include "FINITE TRANSFORMATION/PUSH_FORWARD.f"
#include "FINITE TRANSFORMATION/JAUMANN.f"

!                                                      -----------------
! ------------------------------------------------------ VUMAT_3D_FUNG -
!                                                      -----------------

SUBROUTINE VUMAT_3D_FUNG (NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL, STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH, PROPS, DENSITY, STRAININC, RELSPININC, TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD, STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD, TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW, STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
    IMPLICIT NONE

! -------------------------------------------------- INPUT VARIABLES ---

    CHARACTER*80 CMNAME
    INTEGER      NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL
    REAL         STEPTIME, TOTALTIME, DT, COORDMP(NBLOCK,*), CHARLENGTH(NBLOCK), PROPS(NPROPS), DENSITY(NBLOCK), STRAININC(NBLOCK,NDIR+NSHR), RELSPININC(NBLOCK,NSHR), TEMPOLD(NBLOCK), STRETCHOLD(NBLOCK,NDIR+NSHR), DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR), FIELDOLD(NBLOCK,NFIELDV), STRESSOLD(NBLOCK,NDIR+NSHR), STATEOLD(NBLOCK,NSTATEV), ENERINTERNOLD(NBLOCK), ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK), STRETCHNEW(NBLOCK,NDIR+NSHR), DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR), FIELDNEW(NBLOCK,NFIELDV), STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV), ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)

! -------------------------------------------------- LOCAL VARIABLES ---

    INTEGER          I, J, K, L, IJ, KL, GET_IJ, GET_I, GET_J
    DOUBLE PRECISION C0, B01111, B02222, B03333, B01212, B01313, B02323, T, DT, IDNDI(NDI,NDI), FT(NDI,NDI), DET33, DETFT, CT(NDI,NDI), ET(NDI,NDI), FE(NDI,NDI), DETFE, CE(NDI,NDI), EE(NDI,NDI), EXP_TERM, PK2_VEC(NTENS), PK2_MAT(NDI,NDI), STRESS_MAT(NDI,NDI), DUMMY, H(NTENS,NTENS)

    INTEGER          IS_GROWING, LAMBDAG_NUM_IT
    DOUBLE PRECISION LAMBDAG_TINI, LAMBDAG_TFIN, LAMBDAG_XMIN, LAMBDAG_XMAX, LAMBDAG_GAMMA, LAMBDAG_TAU, LAMBDAG_MIN, LAMBDAG_MAX, LAMBDAG_DOT, LAMBDAG_OLD, LAMBDAG_RES, LAMBDAG_RES_OLD, LAMBDAG, LAMBDAG_X, POS_PART, LAMBDAG_ERR, LAMBDAG_RELAX

! --------------------------------------------------- INITIALIZATION ---

    C0      = PROPS(1)
    B01111 = PROPS(2)
    B02222 = PROPS(3)
    B03333 = PROPS(4)
    B01212 = PROPS(5)
    B01313 = PROPS(6)
    B02323 = PROPS(7)
    DUMMY   = PROPS(8)

    LAMBDAG_TINI  = PROPS( 9)
    LAMBDAG_TFIN  = PROPS(10)
    LAMBDAG_XMIN    = PROPS(11)
    LAMBDAG_XMAX    = PROPS(12)
    LAMBDAG_GAMMA = PROPS(13)
    LAMBDAG_TAU   = PROPS(14)
    LAMBDAG_MIN   = PROPS(15)
    LAMBDAG_MAX   = PROPS(16)

    CALL SETID33(IDNDI)
!     PRINT *, "IDNDI = ", IDNDI

    DT = DTIME
    T  = TIME(2)+DT

    FT    = DFGRD1
    DETFT = DET33(FT)
    CT    = MATMUL(TRANSPOSE(FT), FT)
    ET    = (CT - IDNDI)/2
!     STATEV( 3) = FT(1,1)
!     STATEV( 4) = FT(1,2)
!     STATEV( 5) = FT(1,3)
!     STATEV( 6) = FT(2,1)
!     STATEV( 7) = FT(2,2)
!     STATEV( 8) = FT(2,3)
!     STATEV( 9) = FT(3,1)
!     STATEV(10) = FT(3,2)
!     STATEV(11) = FT(3,3)

! ----------------------------------------------------------- GROWTH ---

!   LAMBDAG
    IF (LAMBDAG_TAU > 0.) THEN
        LAMBDAG = T/LAMBDAG_TAU
    ELSE
        LAMBDAG = 0.
    ENDIF

!   FE, CE, EE
    FE    = FT / (1.+LAMBDAG)
    DETFE = DET33(FE)
    CE    = MATMUL(TRANSPOSE(FE), FE)
    EE    = (CE - IDNDI)/2
!     PRINT *, "FE = ", FE
!     PRINT *, "CE = ", CE
!     PRINT *, "EE = ", EE

!     STATEV(3) = DETFE
!     STATEV(3) = EE(1,1)

!   EXP_TERM
    EXP_TERM =  B01111*EE(1,1)**2 &
             +  B02222*EE(2,2)**2 &
             +  B03333*EE(3,3)**2 &
             +4*B01212*EE(1,2)**2 &
             +4*B01313*EE(1,3)**2 &
             +4*B02323*EE(2,3)**2
!      PRINT *, "EXP_TERM = ", EXP_TERM
    IF (EXP_TERM > 700) THEN
        PRINT *, "STEP ", KSTEP, ", INCREMENT ", KINC, ", ELEMENT ", NOEL, ", INTEGRATION POINT ", NPT, ": TOO LARGE DEFORMATION…"
        EXP_TERM = 700
        PNEWDT = 0.5
    ENDIF
    EXP_TERM = EXP(EXP_TERM)
!     PRINT *, "EXP_TERM = ", EXP_TERM

!     STATEV(3) = EXP_TERM

! ----------------------------------------------------------- STRESS ---

!   SECOND PIOLA KIRSCHOFF STRESS
    PK2_VEC(1) = C0 * EXP_TERM *   B01111*EE(1,1)
    PK2_VEC(2) = C0 * EXP_TERM *   B02222*EE(2,2)
    PK2_VEC(3) = C0 * EXP_TERM *   B03333*EE(3,3)
    PK2_VEC(4) = C0 * EXP_TERM * 2*B01212*EE(1,2)
    PK2_VEC(5) = C0 * EXP_TERM * 2*B01313*EE(1,3)
    PK2_VEC(6) = C0 * EXP_TERM * 2*B02323*EE(2,3)
!     PRINT *, "PK2_VEC = ", PK2_VEC

!     STATEV(3) = PK2_VEC(4)

!   CAUCHY STRESS
    CALL MATSYM33_FROM_VECCOL6(PK2_MAT, PK2_VEC)
    CALL CAUCHYSTRESS33_FROM_PK2STRESS33(STRESS_MAT, PK2_MAT, FT)
    CALL VECCOL6_FROM_MATSYM33(STRESS, STRESS_MAT)
!     PRINT *, "STRESS = ", STRESS

! -------------------------------------------------- STATE VARIABLES ---

    STATEV(1) = LAMBDAG_X
    STATEV(2) = LAMBDAG

! -------------------------------------------------------------- END ---

    RETURN
END SUBROUTINE UMAT_3D_FUNG
