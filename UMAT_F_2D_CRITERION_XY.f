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
      #include 'DAMAGE/DAMAGE_2D_XY_HOOKE.f'
C
C ABAQUS UMAT_2D_CRITERION_XY SUBROUTINE
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
     1 G12, G13, G23, sigmac1, sigmac2, d1, d2, dmax
C
!       PRINT *, "NOEL (element number) =", NOEL
!       PRINT *, "NPT (integration point number) =", NPT
!       PRINT *, "KSTEP (step number) =", KSTEP
!       PRINT *, "KINC (increment number) =", KINC
!       PRINT *, "NSTATV (number of state variables) =", NSTATV
!       PRINT *, "STATEV (state variables)"
!       DO I = 1, NSTATV
!           PRINT *, STATEV(I)
!       END DO
!       PRINT *, "STRESS (total strain at the begining of the increment)"
!       DO I = 1, NTENS
!           PRINT *, STRESS(I)
!       END DO
!       PRINT *, "STRAN (total strain at the begining of the increment)"
!       DO I = 1, NTENS
!           PRINT *, STRAN(I)
!       END DO
!       PRINT *, "DSTRAN (total strain increment for the iteration)"
!       DO I = 1, NTENS
!           PRINT *, DSTRAN(I)
!       END DO
!       PRINT *, "STRANI (total strain at the begining of the iteration)"
!       CALL ADD_VEC(STRANI, STRAN, DSTRAN, NTENS)
!       DO I = 1, NTENS
!           PRINT *, STRANI(I)
!       END DO
C
      E1 = PROPS(1)
      E2 = PROPS(2)
      nu12 = PROPS(3)
      G12 = PROPS(4)
C
      sigmac1 = PROPS(5)
      sigmac2 = PROPS(6)
C
      d1 = STATEV(1)
      d2 = STATEV(2)
C
!       PRINT *, "d1 =", d1
!       PRINT *, "d2 =", d2
C
      CALL DAMAGE_2D_XY_HOOKE(DDSDDE, E1, E2, nu12, G12, d1, d2)
      CALL ADD_VEC(STRANI, STRAN, DSTRAN, NTENS)
      CALL PROD_MAT_VEC(STRESS, DDSDDE, STRANI, NTENS)
C
      dmax = 1. - 1e-6
C
      IF ((d1 < dmax) .AND. (STRESS(1) > sigmac1)) THEN
        d1 = dmax
        STATEV(1) = d1
        CALL DAMAGE_2D_XY_HOOKE(DDSDDE, E1, E2, nu12, G12, d1, d2)
        CALL PROD_MAT_VEC(STRESS, DDSDDE, STRANI, NTENS)
      ENDIF
C
      IF ((d2 < dmax) .AND. (STRESS(2) > sigmac2)) THEN
        d2 = dmax
        STATEV(2) = d2
        CALL DAMAGE_2D_XY_HOOKE(DDSDDE, E1, E2, nu12, G12, d1, d2)
        CALL PROD_MAT_VEC(STRESS, DDSDDE, STRANI, NTENS)
      ENDIF
C
!       PRINT *, "d1 =", d1
!       PRINT *, "d2 =", d2
!       PRINT *, "STRESS"
!       DO I = 1, NTENS
!           PRINT *, STRESS(I)
!       END DO
C
      RETURN
      END












