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
      #include 'DAMAGE_LAWS/DAMAGE_LAW1.f'
C
C ABAQUS UMAT_3D_DAMAGE_Z SUBROUTINE
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
      DOUBLE PRECISION TSTRAN(NTENS), E1, E2, E3, nu12, nu13, nu23,
     1 G12, G13, G23, Y30, Y31, d3, Y3, dmax
C
!       PRINT *, "NOEL (element number) =", NOEL
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
      Y30 = PROPS(10)
      Y31 = PROPS(11)
C
      d3 = STATEV(1)
C
      CALL ADD_VEC(TSTRAN, STRAN, DSTRAN, NTENS)
      CALL DAMAGE_3D_Z_FORCES(Y3, TSTRAN, E1, E2, E3, nu12, nu13,
     1 nu23, G12, G13, G23, d3)
C
!       PRINT *, "Y3 =", Y3
C
      CALL DAMAGE_LAW1(d3, Y3, Y30, Y31)
      dmax = 1. - 1e-6
      IF (d3 > dmax) THEN
        d3 = dmax
      ENDIF
C
!       PRINT *, "d3 =", d3
C
      IF (d3 > STATEV(1)) THEN
        STATEV(1) = d3
      ELSE
        d3 = STATEV(1)
      ENDIF
C
!       PRINT *, "d3 =", d3
C
      CALL DAMAGE_3D_Z_HOOKE(DDSDDE, E1, E2, E3, nu12, nu13, nu23,
     1 G12, G13, G23, d3)
      CALL PROD_MAT_VEC(STRESS, DDSDDE, TSTRAN, NTENS)
C
      RETURN
      END












