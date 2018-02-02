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

#ifndef PRINT_UMAT_DATA_f
#define PRINT_UMAT_DATA_f

!                                                           ------------
! ----------------------------------------------------------- INCLUDES -
!                                                           ------------

!                                                    -------------------
! ---------------------------------------------------- PRINT_UMAT_DATA -
!                                                    -------------------

SUBROUTINE PRINT_UMAT_DATA(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

    CHARACTER*80     CMNAME
    INTEGER          NDI, NSHR, NTENS, NSTATV, NPROPS, NOEL, NPT, LAYER, KSPT, KSTEP, KINC
    DOUBLE PRECISION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS), SSE, SPD, SCD, RPL, DDSDDT(NTENS), DRPLDE(NTENS), DRPLDT, STRAN(NTENS), DSTRAN(NTENS), TIME(2), DTIME, TEMP, DTEMP, PREDEF, DPRED, PROPS(NPROPS), COORDS(3), DROT(3,3), PNEWDT, CELENT, DFGRD0(3,3), DFGRD1(3,3)

    PRINT *, "CMNAME (MATERIAL NAME) = ", CMNAME
    PRINT *, "KSTEP (STEP NUMBER) = ", KSTEP
    PRINT *, "KINC (INCREMENT NUMBER) = ", KINC
    PRINT *, "NOEL (ELEMENT NUMBER) = ", NOEL
    PRINT *, "NPT (INTEGRATION POINT NUMBER) = ", NPT
    PRINT *, "NPROPS (NUMBER OF PROPERTY VARIABLES) = ", NPROPS
    PRINT *, "PROPS (PROPERTY VARIABLES) = ", PROPS
    PRINT *, "NSTATV (NUMBER OF STATE VARIABLES) = ", NSTATV
    PRINT *, "STATEV (STATE VARIABLES) = ", STATEV
    PRINT *, "STRESS (CAUCHY STRESS AT THE BEGINING OF THE INCREMENT) = ", STRESS
    PRINT *, "STRAN (TOTAL STRAIN AT THE BEGINING OF THE INCREMENT) = ", STRAN
    PRINT *, "DSTRAN (TOTAL STRAIN INCREMENT FOR THE ITERATION) = ", DSTRAN

END SUBROUTINE PRINT_UMAT_DATA

#endif