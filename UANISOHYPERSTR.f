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

#include "UANISOHYPERSTR_F_3D_FUNG.f"
#include "UANISOHYPERSTR_F_3D_HOLZAPFEL.f"

!                                                 ----------------------
! ------------------------------------------------- UANISOHYPER_STRAIN -
!                                                 ----------------------

SUBROUTINE UANISOHYPER_STRAIN (EBAR, AJ, UA, DU1, DU2, DU3, TEMP, NOEL, CMNAME, INCMPFLAG, IHYBFLAG, NDI, NSHR, NTENS, NUMSTATEV, STATEV, NUMFIELDV, FIELDV, FIELDVINC, NUMPROPS, PROPS)
!     INCLUDE 'ABA_PARAM.INC'

    CHARACTER*80     CMNAME
    INTEGER          NOEL, INCMPFLAG, IHYBFLAG, NDI, NSHR, NTENS, NUMSTATEV, NUMFIELDV, NUMPROPS
    DOUBLE PRECISION EBAR(NTENS), AJ, UA(2), DU1(NTENS+1), DU2((NTENS+1)*(NTENS+2)/2), DU3((NTENS+1)*(NTENS+2)/2), TEMP, STATEV(NUMSTATEV), FIELDV(NUMFIELDV), FIELDVINC(NUMFIELDV), PROPS(NUMPROPS)

    IF (CMNAME .EQ. 'UANISOHYPERSTR_F_3D_FUNG') THEN
        CALL UANISOHYPERSTR_3D_FUNG (EBAR, AJ, UA, DU1, DU2, DU3, TEMP, NOEL, CMNAME, INCMPFLAG, IHYBFLAG, NDI, NSHR, NTENS, NUMSTATEV, STATEV, NUMFIELDV, FIELDV, FIELDVINC, NUMPROPS, PROPS)
    ELSE IF (CMNAME .EQ. 'UANISOHYPERSTR_F_3D_HOLZAPFEL') THEN
        CALL UANISOHYPERSTR_3D_HOLZAPFEL (EBAR, AJ, UA, DU1, DU2, DU3, TEMP, NOEL, CMNAME, INCMPFLAG, IHYBFLAG, NDI, NSHR, NTENS, NUMSTATEV, STATEV, NUMFIELDV, FIELDV, FIELDVINC, NUMPROPS, PROPS)
    ENDIF

    RETURN
END SUBROUTINE UANISOHYPER_STRAIN

