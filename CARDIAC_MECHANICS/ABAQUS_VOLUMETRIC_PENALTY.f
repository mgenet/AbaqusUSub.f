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

#ifndef ABAQUS_VOLUMETRIC_PENALTY_f
#define ABAQUS_VOLUMETRIC_PENALTY_f

!                                                           ------------
! ----------------------------------------------------------- INCLUDES -
!                                                           ------------

!                                          -----------------------------
! ------------------------------------------ ABAQUS_VOLUMETRIC_PENALTY -
!                                          -----------------------------

SUBROUTINE ABAQUS_VOLUMETRIC_PENALTY(D, J, UJ, DUJ_DJ, DDUJ_DDJ)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: J, D
    DOUBLE PRECISION, INTENT(OUT) :: UJ, DUJ_DJ, DDUJ_DDJ

    UJ = (1./D) * ((J*J-1.)/2 - LOG(J))
    DUJ_DJ = (1./D) * (J - 1./J)
    DDUJ_DDJ = (1./D) * (1. + 1./J/J)

    RETURN
END SUBROUTINE ABAQUS_VOLUMETRIC_PENALTY

#endif
