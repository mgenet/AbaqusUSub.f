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
      #include 'DAMAGE/DAMAGE_3D_XY1_HOOKE.f'
!       #include 'DAMAGE/DAMAGE_3D_XY1_HOOKE_DIRECT.f'
      #include 'DAMAGE/DAMAGE_3D_XY1_FORCES.f'
C
      #include 'DAMAGE_LAWS/DAMAGE_LAW2.f'
C
C ABAQUS UMAT_3D_DAMAGE_XY1 SUBROUTINE
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
      DOUBLE PRECISION TSTRAN(NTENS), E1, E2, E3, nu12,
     1 nu13, nu23, G12, G13, G23, Y10, Y20, d1, d2, Y1, Y2
C
!       PRINT *, "KINC (increment number) =", KINC
!       PRINT *, "NOEL (element number) =", NOEL
!       PRINT *, "NPT (integration point number) =", NPT
!       PRINT *, "KSTEP (step number) =", KSTEP
!       PRINT *, "NSTATV (number of state variables) =", NSTATV
!       PRINT *, "STATEV (state variables)"
!       DO I = 1, NSTATV
!           PRINT *, STATEV(I)
!       END DO
!       PRINT *, "STRESS (stress at the begining of the increment)"
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
!       PRINT *, "TSTRAN (total strain at the begining of the iteration)"
!       CALL ADD_VEC(TSTRAN, STRAN, DSTRAN, NTENS)
!       DO I = 1, NTENS
!           PRINT *, TSTRAN(I)
!       END DO
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
      Y10 = PROPS(10)
      Y20 = PROPS(11)
C
      d1 = STATEV(1)
      d2 = STATEV(2)
C
!       PRINT *, "d1 =", d1
!       PRINT *, "d2 =", d2
C
      CALL ADD_VEC(TSTRAN, STRAN, DSTRAN, NTENS)
      CALL DAMAGE_3D_XY1_FORCES(Y1, Y2, TSTRAN, E1, E2, E3, nu12, nu13,
     1 nu23, G12, G13, G23, d1, d2)
      CALL DAMAGE_LAW2(d1, Y1, Y10)
      CALL DAMAGE_LAW2(d2, Y2, Y20)
C
      IF (d1 > STATEV(1)) THEN
        STATEV(1) = d1
      ELSE
        d1 = STATEV(1)
      ENDIF
      IF (d2 > STATEV(2)) THEN
        STATEV(2) = d2
      ELSE
        d2 = STATEV(2)
      ENDIF
C
!       PRINT *, "d1 =", d1
!       PRINT *, "d2 =", d2
C
      CALL DAMAGE_3D_XY1_HOOKE(DDSDDE, E1, E2, E3, nu12, nu13, nu23,
     1 G12, G13, G23, d1, d2)
      CALL PROD_MAT_VEC(STRESS, DDSDDE, TSTRAN, NTENS)
C
!       PRINT *, "STRESS (stress at the end of the iteration)"
!       DO I = 1, NTENS
!           PRINT *, STRESS(I)
!       END DO
C
!       PRINT *, "DDSDDE"
!       DO I = 1, 6
!         DO J = 1, I
!           PRINT *, 'K(', I, ' ,', J, ' ) =', DDSDDE(I, J)
!         ENDDO
!       ENDDO
C
      RETURN
      END












