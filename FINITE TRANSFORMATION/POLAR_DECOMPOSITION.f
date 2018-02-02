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

#ifndef POLAR_DECOMPOSITION_F
#define POLAR_DECOMPOSITION_F

#include "UTILS/VEC_MAT_TOOLS.f"

!                                                -----------------------
! ------------------------------------------------ POLAR_DECOMPOSITION -
!                                                -----------------------

SUBROUTINE POLAR_DECOMPOSITION_RU_33(F, R, U)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: F(3,3)
    DOUBLE PRECISION, INTENT(OUT) :: R(3,3), U(3,3)

    INTEGER          I, INFO
    DOUBLE PRECISION C(3,3), WW(3,3), W(3), WORK(8), C_DIAG(3,3), U_DIAG(3,3), UINV(3,3)

!     PRINT *, "F = ", F

    C = MATMUL(TRANSPOSE(F), F)
!     PRINT *, "C = ", C

    WW = C
    CALL DSYEV('V', 'U', 3, WW, 3, W, WORK, 8, INFO)

!     PRINT *, "INFO = ", INFO
!     PRINT *, "WW = ", WW
!     PRINT *, "W = ", W

    C_DIAG = MATMUL(MATMUL(TRANSPOSE(WW), C), WW)
!     PRINT *, "C_DIAG = ", C_DIAG

!     PRINT *, "U_DIAG = ", U_DIAG
    DO I = 1,3
        U_DIAG(I,I) = SQRT(C_DIAG(I,I))
    ENDDO
!     PRINT *, "U_DIAG = ", U_DIAG

    U = MATMUL(MATMUL(WW, U_DIAG), TRANSPOSE(WW))
!     PRINT *, "U = ", U

    CALL INV33(U, UINV)
    R = MATMUL(F, UINV)
!     PRINT *, "R = ", R

    RETURN
END SUBROUTINE POLAR_DECOMPOSITION_RU_33

#endif