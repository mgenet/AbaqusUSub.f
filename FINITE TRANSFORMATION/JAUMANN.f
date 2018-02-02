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

#ifndef JAUMANN_F
#define JAUMANN_F

!                                                            -----------
! ------------------------------------------------------------ JAUMANN -
!                                                            -----------

SUBROUTINE JAUMANN(DDSDDE, STRESS)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: STRESS(6)
    DOUBLE PRECISION, INTENT(OUT) :: DDSDDE(6,6)

    INTEGER          IJ, KL, GET_I, GET_J
    DOUBLE PRECISION STRESS33(3,3), ID33(3,3)

!     DDSDDE(1,1) = DDSDDE(1,1) + 2*STRESS(1)
!     DDSDDE(2,2) = DDSDDE(2,2) + 2*STRESS(2)
!     DDSDDE(3,3) = DDSDDE(3,3) + 2*STRESS(3)
!     DDSDDE(4,4) = DDSDDE(4,4) + (STRESS(1)+STRESS(2))/2
!     DDSDDE(5,5) = DDSDDE(5,5) + (STRESS(1)+STRESS(3))/2
!     DDSDDE(6,6) = DDSDDE(6,6) + (STRESS(2)+STRESS(3))/2
!     DDSDDE(1,4) = DDSDDE(1,4) +   STRESS(4)
!     DDSDDE(1,5) = DDSDDE(1,5) +   STRESS(5)
!     DDSDDE(2,4) = DDSDDE(2,4) +   STRESS(4)
!     DDSDDE(2,6) = DDSDDE(2,6) +   STRESS(6)
!     DDSDDE(3,5) = DDSDDE(3,5) +   STRESS(5)
!     DDSDDE(3,6) = DDSDDE(3,6) +   STRESS(6)
!     DDSDDE(4,5) = DDSDDE(4,5) +   STRESS(6)/2
!     DDSDDE(4,6) = DDSDDE(4,6) +   STRESS(5)/2
!     DDSDDE(5,6) = DDSDDE(5,6) +   STRESS(4)/2
!     DDSDDE(4,1) = DDSDDE(1,4)
!     DDSDDE(5,1) = DDSDDE(1,5)
!     DDSDDE(4,2) = DDSDDE(2,4)
!     DDSDDE(6,2) = DDSDDE(2,6)
!     DDSDDE(5,3) = DDSDDE(3,5)
!     DDSDDE(6,3) = DDSDDE(3,6)
!     DDSDDE(5,4) = DDSDDE(4,5)
!     DDSDDE(6,4) = DDSDDE(4,6)
!     DDSDDE(6,5) = DDSDDE(5,6)

    CALL MATSYM33_FROM_VECCOL6(STRESS33, STRESS)
    CALL SETID33(ID33)

    DO IJ = 1,6
    DO KL = 1,6
        DDSDDE(IJ,KL) = DDSDDE(IJ,KL) + (ID33(GET_I(3,IJ),GET_I(3,KL)) * STRESS33(GET_J(3,IJ),GET_J(3,KL)) &
                                       + ID33(GET_J(3,IJ),GET_J(3,KL)) * STRESS33(GET_I(3,IJ),GET_I(3,KL)) &
                                       + ID33(GET_J(3,IJ),GET_I(3,KL)) * STRESS33(GET_I(3,IJ),GET_J(3,KL)) &
                                       + ID33(GET_I(3,IJ),GET_J(3,KL)) * STRESS33(GET_J(3,IJ),GET_I(3,KL))) / 2
    END DO
    END DO

    RETURN
END SUBROUTINE JAUMANN

#endif