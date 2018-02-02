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

    INTEGER I, J, K_STATEV
    DOUBLE PRECISION C0, B0_FFFF, B0_SSSS, B0_NNNN, B0_FSFS, B0_FNFN, B0_SNSN, IDNDI(NDI,NDI), F(NDI,NDI), C(NDI,NDI), E(NDI,NDI), EXP_TERM, PK2_VEC(NTENS), PK2_MAT(NDI,NDI), STRESS_MAT(NDI,NDI), H(NTENS,NTENS), R(NDI, NDI), U(NDI, NDI), DFGRD0_INV(NDI,NDI), DF(NDI,NDI), DR(NDI,NDI), DU(NDI,NDI), ROT(NDI,NDI)

! --------------------------------------------------- INITIALISATION ---

    C0      = PROPS(1)
    B0_FFFF = PROPS(2)
    B0_SSSS = PROPS(3)
    B0_NNNN = PROPS(4)
    B0_FSFS = PROPS(5)
    B0_FNFN = PROPS(6)
    B0_SNSN = PROPS(7)

! ------------------------------------------------------------ PRINT ---

!     PRINT *, "KSTEP (step number) =", KINC
!     PRINT *, "KINC (increment number) =", KINC
!     PRINT *, "NOEL (element number) =", NOEL
!     PRINT *, "NPT (integration point number) =", NPT

!     PRINT *, "DROT = ", DROT

!     IF (KINC .EQ. 1) THEN
!         CALL SETID33(ROT)
!     ELSE
!         K_STATEV = 1
!         DO I = 1,3
!             DO J = 1,3
!                 ROT(I,J) = STATEV(K_STATEV)
!                 K_STATEV = K_STATEV + 1
!             ENDDO
!         ENDDO
!         ROT = MATMUL(DROT, ROT)
!     ENDIF
!     K_STATEV = 1
!     DO I = 1,3
!         DO J = 1,3
!             STATEV(K_STATEV) = ROT(I,J)
!             K_STATEV = K_STATEV + 1
!         ENDDO
!     ENDDO
!     PRINT *, "ROT = ", ROT

!     CALL POLAR_DECOMPOSITION_RU_33(DFGRD1, R, U)
!     PRINT *, "F = ", DFGRD1
!     PRINT *, "R = ", R
!     PRINT *, "U = ", U
!     CALL INV33(DFGRD0, DFGRD0_INV)
!     DF = MATMUL(DFGRD1, DFGRD0_INV)
!     CALL POLAR_DECOMPOSITION_RU_33(DF, DR, DU)
!     PRINT *, "DF = ", DF
!     PRINT *, "DR = ", DR
!     PRINT *, "DU = ", DU

! ----------------------------------------------------------- STRAIN ---

    F = DFGRD1
    C = MATMUL(TRANSPOSE(F), F)

    CALL SETID33(IDNDI)
    E = (C - IDNDI)/2

! ----------------------------------------------------------- STRESS ---

!   EXP_TERM
    EXP_TERM =  B0_FFFF*E(1,1)**2 &
             +  B0_SSSS*E(2,2)**2 &
             +  B0_NNNN*E(3,3)**2 &
             +4*B0_FSFS*E(1,2)**2 &
             +4*B0_FNFN*E(1,3)**2 &
             +4*B0_SNSN*E(2,3)**2
!     PRINT *, "EXP_TERM = ", EXP_TERM
    EXP_TERM = EXP(EXP_TERM)
!     PRINT *, "EXP_TERM = ", EXP_TERM

!   SECOND PIOLA KIRSCHOFF STRESS
    PK2_VEC(1) = C0 * EXP_TERM *   B0_FFFF*E(1,1)
    PK2_VEC(2) = C0 * EXP_TERM *   B0_SSSS*E(2,2)
    PK2_VEC(3) = C0 * EXP_TERM *   B0_NNNN*E(3,3)
    PK2_VEC(4) = C0 * EXP_TERM * 2*B0_FSFS*E(1,2)
    PK2_VEC(5) = C0 * EXP_TERM * 2*B0_FNFN*E(1,3)
    PK2_VEC(6) = C0 * EXP_TERM * 2*B0_SNSN*E(2,3)
!     PRINT *, "PK2_VEC = ", PK2_VEC

!   CAUCHY STRESS
    CALL MATSYM33_FROM_VECCOL6(PK2_MAT, PK2_VEC)
    CALL CAUCHYSTRESS33_FROM_PK2STRESS33(STRESS_MAT, PK2_MAT, F)
    CALL VECCOL6_FROM_MATSYM33(STRESS, STRESS_MAT)
!     PRINT *, "STRESS = ", STRESS

! ---------------------------------------------------------- TANGENT ---

!   TANGENT IN TERMS OF SECOND PIOLA KIRSCHOFF STRESS
    H = 0.
    H(1,1) = C0 * EXP_TERM * (B0_FFFF+2*B0_FFFF*E(1,1)*B0_FFFF*E(1,1))
    H(1,2) = C0 * EXP_TERM * (        2*B0_FFFF*E(1,1)*B0_SSSS*E(2,2))
    H(1,3) = C0 * EXP_TERM * (        2*B0_FFFF*E(1,1)*B0_NNNN*E(3,3))
    H(1,4) = C0 * EXP_TERM * (        4*B0_FFFF*E(1,1)*B0_FSFS*E(1,2))
    H(1,5) = C0 * EXP_TERM * (        4*B0_FFFF*E(1,1)*B0_FNFN*E(1,3))
    H(1,6) = C0 * EXP_TERM * (        4*B0_FFFF*E(1,1)*B0_SNSN*E(2,3))
    H(2,1) = C0 * EXP_TERM * (        2*B0_SSSS*E(2,2)*B0_FFFF*E(1,1))
    H(2,2) = C0 * EXP_TERM * (B0_SSSS+2*B0_SSSS*E(2,2)*B0_SSSS*E(2,2))
    H(2,3) = C0 * EXP_TERM * (        2*B0_SSSS*E(2,2)*B0_NNNN*E(3,3))
    H(2,4) = C0 * EXP_TERM * (        4*B0_SSSS*E(2,2)*B0_FSFS*E(1,2))
    H(2,5) = C0 * EXP_TERM * (        4*B0_SSSS*E(2,2)*B0_FNFN*E(1,3))
    H(2,6) = C0 * EXP_TERM * (        4*B0_SSSS*E(2,2)*B0_SNSN*E(2,3))
    H(3,1) = C0 * EXP_TERM * (        2*B0_NNNN*E(3,3)*B0_FFFF*E(1,1))
    H(3,2) = C0 * EXP_TERM * (        2*B0_NNNN*E(3,3)*B0_SSSS*E(2,2))
    H(3,3) = C0 * EXP_TERM * (B0_NNNN+2*B0_NNNN*E(3,3)*B0_NNNN*E(3,3))
    H(3,4) = C0 * EXP_TERM * (        4*B0_NNNN*E(3,3)*B0_FSFS*E(1,2))
    H(3,5) = C0 * EXP_TERM * (        4*B0_NNNN*E(3,3)*B0_FNFN*E(1,3))
    H(3,6) = C0 * EXP_TERM * (        4*B0_NNNN*E(3,3)*B0_SNSN*E(2,3))
    H(4,1) = C0 * EXP_TERM * (        4*B0_FSFS*E(1,2)*B0_FFFF*E(1,1))
    H(4,2) = C0 * EXP_TERM * (        4*B0_FSFS*E(1,2)*B0_SSSS*E(2,2))
    H(4,3) = C0 * EXP_TERM * (        4*B0_FSFS*E(1,2)*B0_NNNN*E(3,3))
    H(4,4) = C0 * EXP_TERM * (B0_FSFS+8*B0_FSFS*E(1,2)*B0_FSFS*E(1,2))
    H(4,5) = C0 * EXP_TERM * (        8*B0_FSFS*E(1,2)*B0_FNFN*E(1,3))
    H(4,6) = C0 * EXP_TERM * (        8*B0_FSFS*E(1,2)*B0_SNSN*E(2,3))
    H(5,1) = C0 * EXP_TERM * (        4*B0_FNFN*E(1,3)*B0_FFFF*E(1,1))
    H(5,2) = C0 * EXP_TERM * (        4*B0_FNFN*E(1,3)*B0_SSSS*E(2,2))
    H(5,3) = C0 * EXP_TERM * (        4*B0_FNFN*E(1,3)*B0_NNNN*E(3,3))
    H(5,4) = C0 * EXP_TERM * (        8*B0_FNFN*E(1,3)*B0_FSFS*E(1,2))
    H(5,5) = C0 * EXP_TERM * (B0_FNFN+8*B0_FNFN*E(1,3)*B0_FNFN*E(1,3))
    H(5,6) = C0 * EXP_TERM * (        8*B0_FNFN*E(1,3)*B0_SNSN*E(2,3))
    H(6,1) = C0 * EXP_TERM * (        4*B0_SNSN*E(2,3)*B0_FFFF*E(1,1))
    H(6,2) = C0 * EXP_TERM * (        4*B0_SNSN*E(2,3)*B0_SSSS*E(2,2))
    H(6,3) = C0 * EXP_TERM * (        4*B0_SNSN*E(2,3)*B0_NNNN*E(3,3))
    H(6,4) = C0 * EXP_TERM * (        8*B0_SNSN*E(2,3)*B0_FSFS*E(1,2))
    H(6,5) = C0 * EXP_TERM * (        8*B0_SNSN*E(2,3)*B0_FNFN*E(1,3))
    H(6,6) = C0 * EXP_TERM * (B0_SNSN+8*B0_SNSN*E(2,3)*B0_SNSN*E(2,3))
!     PRINT *, "H = ", H

!   TANGENT IN TERMS OF CAUCHY STRESS
    CALL CAUCHYJACOBIAN66_FROM_PK2JACOBIAN66(DDSDDE, H, F)
!     PRINT *, "DDSDDE = ", DDSDDE

!   JAUMANN TERMS IN TANGENT
    CALL JAUMANN(DDSDDE, STRESS)
!     PRINT *, "DDSDDE = ", DDSDDE

! -------------------------------------------------------------- END ---

    RETURN
END SUBROUTINE UMAT_TEST_OBJECTIVITY
