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

#ifndef PUSH_FORWARD_f
#define PUSH_FORWARD_f

!                                                -----------------------
! ------------------------------------------------ PUSH STRESS FORWARD -
!                                                -----------------------

SUBROUTINE CAUCHYSTRESS33_FROM_PK2STRESS33(CAUCHY_STRESS, PK2_STRESS, F)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: PK2_STRESS(3,3), F(3,3)
    DOUBLE PRECISION, INTENT(OUT) :: CAUCHY_STRESS(3,3)

    DOUBLE PRECISION DET33

    CAUCHY_STRESS = MATMUL(MATMUL(F, PK2_STRESS), TRANSPOSE(F)) / DET33(F)

    RETURN
END SUBROUTINE CAUCHYSTRESS33_FROM_PK2STRESS33

!                                              -------------------------
! ---------------------------------------------- PUSH JACOBIAN FORWARD -
!                                              -------------------------

SUBROUTINE CAUCHYJACOBIAN66_FROM_PK2JACOBIAN66(CAUCHY_JACOBIAN, PK2_JACOBIAN, F)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: PK2_JACOBIAN(6,6), F(3,3)
    DOUBLE PRECISION, INTENT(OUT) :: CAUCHY_JACOBIAN(6,6)

    INTEGER          IJ, KL, I, J, K, L, GET_IJ, GET_I, GET_J
    DOUBLE PRECISION DET33

    CAUCHY_JACOBIAN = 0.
    DO IJ = 1,6
    DO KL = 1,6
        DO I = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
            CAUCHY_JACOBIAN(IJ,KL) = CAUCHY_JACOBIAN(IJ,KL) &
                                   + PK2_JACOBIAN(GET_IJ(3,I,J), GET_IJ(3,K,L)) * F(GET_I(3,IJ),I) &
                                                                                * F(GET_J(3,IJ),J) &
                                                                                * F(GET_I(3,KL),K) &
                                                                                * F(GET_J(3,KL),L)
        END DO
        END DO
        END DO
        END DO
    END DO
    END DO
    CAUCHY_JACOBIAN = CAUCHY_JACOBIAN / DET33(F)

    RETURN
END SUBROUTINE CAUCHYJACOBIAN66_FROM_PK2JACOBIAN66

#endif