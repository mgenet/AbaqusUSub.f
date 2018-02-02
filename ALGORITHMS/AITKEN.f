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

!                                                             ----------
! ------------------------------------------------------------- AITKEN -
!                                                             ----------

SUBROUTINE AITKEN(RELAX, NUM_IT, RES, RES_OLD)
    IMPLICIT NONE

    INTEGER, INTENT(IN) ::          NUM_IT
    DOUBLE PRECISION, INTENT(IN) :: RES, RES_OLD

    DOUBLE PRECISION, INTENT(OUT) :: RELAX

    IF ((NUM_IT .EQ. 0) .OR. (RES == RES_OLD)) THEN
        RELAX = 1.
    ELSE
        RELAX = - RELAX * RES_OLD / (RES - RES_OLD)
    ENDIF

    RETURN
END SUBROUTINE AITKEN

