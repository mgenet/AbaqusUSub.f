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

#include "UANISOHYPERINV_F_3D_HOLZAPFEL.f"

!                                                    -------------------
! ---------------------------------------------------- UANISOHYPER_INV -
!                                                    -------------------

SUBROUTINE UANISOHYPER_INV (AINV, UA, ZETA, NFIBERS, NINV, UI1, UI2, UI3, TEMP, NOEL, CMNAME, INCMPFLAG, IHYBFLAG, NUMSTATEV, STATEV, NUMFIELDV, FIELDV, FIELDVINC, NUMPROPS, PROPS)
!     INCLUDE 'ABA_PARAM.INC'

    CHARACTER*80     CMNAME
    INTEGER          NFIBERS, NINV, NOEL, INCMPFLAG, IHYBFLAG, NUMSTATEV, NUMFIELDV, NUMPROPS
    DOUBLE PRECISION AINV(NINV), UA(2), ZETA(NFIBERS*(NFIBERS-1)/2), UI1(NINV), UI2(NINV*(NINV+1)/2), UI3(NINV*(NINV+1)/2), TEMP, STATEV(NUMSTATEV), FIELDV(NUMFIELDV), FIELDVINC(NUMFIELDV), PROPS(NUMPROPS)

    IF (CMNAME .EQ. 'UANISOHYPERINV_F_3D_HOLZAPFEL') THEN
        CALL UANISOHYPERINV_3D_HOLZAPFEL(AINV, UA, ZETA, NFIBERS, NINV, UI1, UI2, UI3, TEMP, NOEL, CMNAME, INCMPFLAG, IHYBFLAG, NUMSTATEV, STATEV, NUMFIELDV, FIELDV, FIELDVINC, NUMPROPS, PROPS)
    END IF

    RETURN
END SUBROUTINE UANISOHYPER_INV
