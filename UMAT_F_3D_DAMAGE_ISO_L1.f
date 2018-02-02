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
      #include 'DAMAGE/DAMAGE_3D_ISO_HOOKE.f'
      #include 'DAMAGE/DAMAGE_3D_ISO_FORCES.f'
C
      #include 'DAMAGE_LAWS/DAMAGE_LAW1.f'
C
C ABAQUS UMAT_3D_DAMAGE_ISO SUBROUTINE
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
      DOUBLE PRECISION TSTRAN(NTENS), E, nu, Y0, Y1, d, Y, dmax
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
      E = PROPS(1)
      nu = PROPS(2)
C
      Y0 = PROPS(3)
      Y1 = PROPS(4)
C
      d = STATEV(1)
C
!       PRINT *, "d =", d
C
      CALL ADD_VEC(TSTRAN, STRAN, DSTRAN, NTENS)
      CALL DAMAGE_3D_ISO_FORCES(Y, TSTRAN, E, nu, d)
      CALL DAMAGE_LAW1(d, Y, Y0, Y1)
C
      dmax = 1. - 1e-6
      IF (d > dmax) THEN
        d = dmax
      ENDIF
C
      IF (d > STATEV(1)) THEN
        STATEV(1) = d
      ELSE
        d = STATEV(1)
      ENDIF
C
!       PRINT *, "d =", d
C
      CALL DAMAGE_3D_ISO_HOOKE(DDSDDE, E, nu, d)
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












