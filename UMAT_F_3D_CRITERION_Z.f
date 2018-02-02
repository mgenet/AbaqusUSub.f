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

C
C REQUESTED SUBROUTINES
C
      #include '../ABAQUS_TOOLS/TOOLS.f'
C
      #include 'DAMAGE/DAMAGE_3D_Z_HOOKE.f'
      #include 'DAMAGE/DAMAGE_3D_Z_FORCES.f'
C
C ABAQUS UMAT_3D_CRITERION_Z SUBROUTINE
C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     3 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      #include 'ABA_PARAM.INC'
C
      CHARACTER CMNAME(80)
      INTEGER NTENS, NSTATV, NPROPS, NOEL, NPT, KSTEP, KINC
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),
     2 PREDEF,DPRED,PROPS(NPROPS),COORDS(3),DROT(3,3),PNEWDT,
     3 DFGRD0(3,3),DFGRD1(3,3)
C
      DOUBLE PRECISION STRANI(NTENS), E1, E2, E3, nu12, nu13, nu23,
     1 G12, G13, G23, sigmac3, d3, dmax
C
      E1 = PROPS(1)
      E2 = PROPS(2)
      E3 = PROPS(3)
      nu12 = PROPS(4)
      nu13 = PROPS(5)
      nu23 = PROPS(6)
      G12 = PROPS(7)
      G13 = PROPS(8)
      G23 = PROPS(9)
C
      sigmac3 = PROPS(10)
C
      d3 = STATEV(1)
C
      CALL DAMAGE_Z_HOOKE(DDSDDE, E1, E2, E3, nu12, nu13, nu23,
     1 G12, G13, G23, d3)
      CALL ADD_VEC(STRANI, STRAN, DSTRAN, NTENS)
      CALL PROD_MAT_VEC(STRESS, DDSDDE, STRANI, NTENS)
C
      dmax = 1. - 1e-6
      IF ((d3 < dmax) .AND. (STRESS(3) > sigmac3)) THEN
        d3 = dmax
        STATEV(1) = d3
        CALL DAMAGE_Z_HOOKE(DDSDDE, E1, E2, E3, nu12, nu13, nu23,
       1 G12, G13, G23, d3)
        CALL PROD_MAT_VEC(STRESS, DDSDDE, STRANI, NTENS)
      ENDIF
C
      RETURN
      END












