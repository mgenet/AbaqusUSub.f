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

#include "UTILS/VEC_MAT_TOOLS.f"
#include "FINITE TRANSFORMATION/POLAR_DECOMPOSITION.f"
#include "FINITE TRANSFORMATION/PUSH_FORWARD.f"
#include "FINITE TRANSFORMATION/JAUMANN.f"

!                                              -------------------------
! ---------------------------------------------- UMAT_TEST_OBJECTIVITY -
!                                              -------------------------

SUBROUTINE UMAT_TEST_OBJECTIVITY (STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDETF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)
    IMPLICIT NONE

! -------------------------------------------------- INPUT VARIABLES ---

    CHARACTER*80 CMNAME
    INTEGER NDI, NSHR, NTENS, NSTATV, NPROPS, NOEL, NPT, LAYER, KSPT, KSTEP, KINC
    DOUBLE PRECISION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS), SSE, SPD, SCD, RPL, DDSDDT(NTENS), DRPLDE(NTENS), DRPLDT, STRAN(NTENS), DSTRAN(NTENS), TIME(2), DTIME, TEMP, DTEMP, PREDETF, DPRED, PROPS(NPROPS), COORDS(NDI), DROT(NDI,NDI), PNEWDT, CELENT, DFGRD0(NDI,NDI), DFGRD1(NDI,NDI)

! -------------------------------------------------- LOCAL VARIABLES ---

    DOUBLE PRECISION C0, B0(NTENS), B0_XYZ(NDI,NDI), F_XYZ(NDI,NDI), R_XYZ(NDI,NDI), U_XYZ(NDI,NDI), C_XYZ(NDI,NDI), C_FSN(NDI,NDI), IDNDI(NDI,NDI), E_FSN(NDI,NDI), EXP_TERM, PK2_FSN(NDI,NDI), PK2_XYZ(NDI,NDI), STRESS_XYZ(NDI,NDI), STRESS_FSN(NDI,NDI), H_FSN(NTENS,NTENS), H_XYZ(NTENS,NTENS)

    INTEGER N_ELEMS, I, J, K, L, IJ, KL, GET_I, GET_J, GET_IJ
    PARAMETER (N_ELEMS = 1)
    DOUBLE PRECISION ORIENT(N_ELEMS,9)
    COMMON /ORIENT/ ORIENT

! ------------------------------------------------------------ PRINT ---

!     PRINT *, "KSTEP (step number) =", KINC
!     PRINT *, "KINC (increment number) =", KINC
!     PRINT *, "NOEL (element number) =", NOEL
!     PRINT *, "NPT (integration point number) =", NPT
!     PRINT *, "DROT = ", DROT

! --------------------------------------------------- INITIALISATION ---

    C0    = PROPS(1)
    B0(1) = PROPS(2)
    B0(2) = PROPS(3)
    B0(3) = PROPS(4)
    B0(4) = PROPS(5)
    B0(5) = PROPS(6)
    B0(6) = PROPS(7)

    B0_XYZ(1,1) = ORIENT(NOEL,1)
    B0_XYZ(1,2) = ORIENT(NOEL,2)
    B0_XYZ(1,3) = ORIENT(NOEL,3)
    B0_XYZ(2,1) = ORIENT(NOEL,4)
    B0_XYZ(2,2) = ORIENT(NOEL,5)
    B0_XYZ(2,3) = ORIENT(NOEL,6)
    B0_XYZ(3,1) = ORIENT(NOEL,7)
    B0_XYZ(3,2) = ORIENT(NOEL,8)
    B0_XYZ(3,3) = ORIENT(NOEL,9)

! ----------------------------------------------------------- STRAIN ---

    F_XYZ = DFGRD1
!     PRINT *, "F_XYZ = ", F_XYZ

    CALL POLAR_DECOMPOSITION_RU_33(F_XYZ, R_XYZ, U_XYZ)
!     PRINT *, "R_XYZ = ", R_XYZ
!     PRINT *, "U_XYZ = ", U_XYZ

    C_XYZ = MATMUL(TRANSPOSE(F_XYZ), F_XYZ)
!     PRINT *, "C_XYZ = ", C_XYZ

    C_FSN = MATMUL(MATMUL(B0_XYZ, C_XYZ), TRANSPOSE(B0_XYZ))
!     PRINT *, "C_FSN = ", C_FSN

    CALL SETID33(IDNDI)
    E_FSN = (C_FSN - IDNDI)/2
!     PRINT *, "E_FSN = ", E_FSN

    STATEV(1) = E_FSN(1,1)
    STATEV(2) = E_FSN(2,2)
    STATEV(3) = E_FSN(3,3)
    STATEV(4) = E_FSN(1,2)
    STATEV(5) = E_FSN(1,3)
    STATEV(6) = E_FSN(2,3)

! ----------------------------------------------------------- STRESS ---

!   EXP_TERM
    EXP_TERM = B0(1) * E_FSN(1,1) * E_FSN(1,1) &
             + B0(2) * E_FSN(2,2) * E_FSN(2,2) &
             + B0(3) * E_FSN(3,3) * E_FSN(3,3) &
             + B0(4) * E_FSN(1,2) * E_FSN(1,2) &
             + B0(4) * E_FSN(2,1) * E_FSN(2,1) &
             + B0(4) * E_FSN(1,2) * E_FSN(2,1) &
             + B0(4) * E_FSN(2,1) * E_FSN(1,2) &
             + B0(5) * E_FSN(1,3) * E_FSN(1,3) &
             + B0(5) * E_FSN(3,1) * E_FSN(3,1) &
             + B0(5) * E_FSN(1,3) * E_FSN(3,1) &
             + B0(5) * E_FSN(3,1) * E_FSN(1,3) &
             + B0(6) * E_FSN(2,3) * E_FSN(2,3) &
             + B0(6) * E_FSN(3,2) * E_FSN(3,2) &
             + B0(6) * E_FSN(2,3) * E_FSN(3,2) &
             + B0(6) * E_FSN(3,2) * E_FSN(2,3)
    EXP_TERM = EXP(EXP_TERM)
!     PRINT *, "EXP_TERM = ", EXP_TERM

!   ENERGY
    SSE = (C0/2) * (EXP_TERM - 1.)

!   SECOND PIOLA KIRSCHOFF STRESS
    PK2_FSN(1,1) = C0 * EXP_TERM * B0(1) * E_FSN(1,1)
    PK2_FSN(2,2) = C0 * EXP_TERM * B0(2) * E_FSN(2,2)
    PK2_FSN(3,3) = C0 * EXP_TERM * B0(3) * E_FSN(3,3)
    PK2_FSN(1,2) = C0 * EXP_TERM * B0(4) * E_FSN(1,2) &
                 + C0 * EXP_TERM * B0(4) * E_FSN(2,1)
    PK2_FSN(2,1) = C0 * EXP_TERM * B0(4) * E_FSN(2,1) &
                 + C0 * EXP_TERM * B0(4) * E_FSN(1,2)
    PK2_FSN(1,3) = C0 * EXP_TERM * B0(5) * E_FSN(1,3) &
                 + C0 * EXP_TERM * B0(5) * E_FSN(3,1)
    PK2_FSN(3,1) = C0 * EXP_TERM * B0(5) * E_FSN(3,1) &
                 + C0 * EXP_TERM * B0(5) * E_FSN(1,3)
    PK2_FSN(2,3) = C0 * EXP_TERM * B0(6) * E_FSN(2,3) &
                 + C0 * EXP_TERM * B0(6) * E_FSN(3,2)
    PK2_FSN(3,2) = C0 * EXP_TERM * B0(6) * E_FSN(3,2) &
                 + C0 * EXP_TERM * B0(6) * E_FSN(2,3)
!     PRINT *, "PK2_FSN = ", PK2_FSN

    PK2_XYZ = MATMUL(MATMUL(TRANSPOSE(B0_XYZ), PK2_FSN), B0_XYZ)
!     PRINT *, "PK2_XYZ = ", PK2_XYZ

!   CAUCHY STRESS
    CALL CAUCHYSTRESS33_FROM_PK2STRESS33(STRESS_XYZ, PK2_XYZ, F_XYZ)
    CALL VECCOL6_FROM_MATSYM33(STRESS, STRESS_XYZ)
!     PRINT *, "STRESS = ", STRESS

    STRESS_FSN = MATMUL(MATMUL(MATMUL(MATMUL(B0_XYZ, TRANSPOSE(R_XYZ)), STRESS_XYZ), R_XYZ), TRANSPOSE(B0_XYZ))
!     PRINT *, "STRESS_FSN = ", STRESS_FSN

    STATEV(7)  = STRESS_FSN(1,1)
    STATEV(8)  = STRESS_FSN(2,2)
    STATEV(9)  = STRESS_FSN(3,3)
    STATEV(10) = STRESS_FSN(1,2)
    STATEV(11) = STRESS_FSN(1,3)
    STATEV(12) = STRESS_FSN(2,3)

! ---------------------------------------------------------- TANGENT ---

!   TANGENT IN TERMS OF SECOND PIOLA KIRSCHOFF STRESS
    H_FSN = 0.
    H_FSN(1,1) = C0 * EXP_TERM * (B0(1) + B0(1) * E_FSN(1,1) * (2 * B0(1) * E_FSN(1,1)                         ))
    H_FSN(1,2) = C0 * EXP_TERM * (        B0(1) * E_FSN(1,1) * (2 * B0(2) * E_FSN(2,2)                         ))
    H_FSN(1,3) = C0 * EXP_TERM * (        B0(1) * E_FSN(1,1) * (2 * B0(3) * E_FSN(3,3)                         ))
    H_FSN(1,4) = C0 * EXP_TERM * (        B0(1) * E_FSN(1,1) * (2 * B0(4) * E_FSN(1,2) + 2 * B0(4) * E_FSN(2,1)))
    H_FSN(1,5) = C0 * EXP_TERM * (        B0(1) * E_FSN(1,1) * (2 * B0(5) * E_FSN(1,3) + 2 * B0(5) * E_FSN(3,1)))
    H_FSN(1,6) = C0 * EXP_TERM * (        B0(1) * E_FSN(1,1) * (2 * B0(6) * E_FSN(2,3) + 2 * B0(6) * E_FSN(3,2)))

    H_FSN(2,1) = H_FSN(1,2)
    H_FSN(2,2) = C0 * EXP_TERM * (B0(2) + B0(2) * E_FSN(2,2) * (2 * B0(2) * E_FSN(2,2)                         ))
    H_FSN(2,3) = C0 * EXP_TERM * (        B0(2) * E_FSN(2,2) * (2 * B0(3) * E_FSN(3,3)                         ))
    H_FSN(2,4) = C0 * EXP_TERM * (        B0(2) * E_FSN(2,2) * (2 * B0(4) * E_FSN(1,2) + 2 * B0(4) * E_FSN(2,1)))
    H_FSN(2,5) = C0 * EXP_TERM * (        B0(2) * E_FSN(2,2) * (2 * B0(5) * E_FSN(1,3) + 2 * B0(5) * E_FSN(3,1)))
    H_FSN(2,6) = C0 * EXP_TERM * (        B0(2) * E_FSN(2,2) * (2 * B0(6) * E_FSN(2,3) + 2 * B0(6) * E_FSN(3,2)))

    H_FSN(3,1) = H_FSN(1,3)
    H_FSN(3,2) = H_FSN(2,3)
    H_FSN(3,3) = C0 * EXP_TERM * (B0(3) + B0(3) * E_FSN(3,3) * (2 * B0(3) * E_FSN(3,3)                         ))
    H_FSN(3,4) = C0 * EXP_TERM * (        B0(3) * E_FSN(3,3) * (2 * B0(4) * E_FSN(1,2) + 2 * B0(4) * E_FSN(2,1)))
    H_FSN(3,5) = C0 * EXP_TERM * (        B0(3) * E_FSN(3,3) * (2 * B0(5) * E_FSN(1,3) + 2 * B0(5) * E_FSN(3,1)))
    H_FSN(3,6) = C0 * EXP_TERM * (        B0(3) * E_FSN(3,3) * (2 * B0(6) * E_FSN(2,3) + 2 * B0(6) * E_FSN(3,2)))

    H_FSN(4,1) = H_FSN(1,4)
    H_FSN(4,2) = H_FSN(2,4)
    H_FSN(4,3) = H_FSN(3,4)
    H_FSN(4,4) = C0 * EXP_TERM * (B0(4) + B0(4) * E_FSN(1,2) * (2 * B0(4) * E_FSN(1,2) + 2 * B0(4) * E_FSN(2,1)))
    H_FSN(4,5) = C0 * EXP_TERM * (        B0(4) * E_FSN(1,2) * (2 * B0(5) * E_FSN(1,3) + 2 * B0(5) * E_FSN(3,1)))
    H_FSN(4,6) = C0 * EXP_TERM * (        B0(4) * E_FSN(1,2) * (2 * B0(6) * E_FSN(2,3) + 2 * B0(6) * E_FSN(3,2)))

    H_FSN(5,1) = H_FSN(1,5)
    H_FSN(5,2) = H_FSN(2,5)
    H_FSN(5,3) = H_FSN(3,5)
    H_FSN(5,4) = H_FSN(4,5)
    H_FSN(5,5) = C0 * EXP_TERM * (B0(5) + B0(5) * E_FSN(1,3) * (2 * B0(5) * E_FSN(1,3) + 2 * B0(5) * E_FSN(3,1)))
    H_FSN(5,6) = C0 * EXP_TERM * (        B0(5) * E_FSN(1,3) * (2 * B0(6) * E_FSN(2,3) + 2 * B0(6) * E_FSN(3,2)))

    H_FSN(6,1) = H_FSN(1,6)
    H_FSN(6,2) = H_FSN(2,6)
    H_FSN(6,3) = H_FSN(3,6)
    H_FSN(6,4) = H_FSN(4,6)
    H_FSN(6,5) = H_FSN(5,6)
    H_FSN(6,6) = C0 * EXP_TERM * (B0(6) + B0(6) * E_FSN(2,3) * (2 * B0(6) * E_FSN(2,3) + 2 * B0(6) * E_FSN(3,2)))
!     PRINT *, "H_FSN = ", H_FSN

    H_XYZ = 0.
    DO IJ = 1,6
    DO KL = 1,6
        DO I = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
            H_XYZ(IJ,KL) = H_XYZ(IJ,KL) &
                         + H_FSN(GET_IJ(3,I,J), GET_IJ(3,K,L)) * B0_XYZ(GET_I(3,IJ),I) &
                                                               * B0_XYZ(GET_J(3,IJ),J) &
                                                               * B0_XYZ(GET_I(3,KL),K) &
                                                               * B0_XYZ(GET_J(3,KL),L)
        ENDDO
        ENDDO
        ENDDO
        ENDDO
    ENDDO
    ENDDO

!   TANGENT IN TERMS OF CAUCHY STRESS
    CALL CAUCHYJACOBIAN66_FROM_PK2JACOBIAN66(DDSDDE, H_XYZ, F_XYZ)
!     PRINT *, "DDSDDE = ", DDSDDE

!   JAUMANN TERMS IN TANGENT
    CALL JAUMANN(DDSDDE, STRESS)
!     PRINT *, "DDSDDE = ", DDSDDE

! -------------------------------------------------------------- END ---

    RETURN
END SUBROUTINE UMAT_TEST_OBJECTIVITY
