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
      #include '../TOOLS_FORTRAN/TOOLS.f'
C
      #include 'ELASTICITY/ELASTICITY_3D_ISO_HOOKE.f'
!       #include 'ELASTICITY/ELASTICITY_3D_ISO_HOOKE_DIRECT.f'
C
C ABAQUS UMAT_3D_ELASTICITY_ISO SUBROUTINE
C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1 DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     2 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     3 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
        parameter (nprecd=1) !if gfortran compilation
!         #include 'ABA_PARAM.INC' !if abaqus compilation
C
      CHARACTER CMNAME(80)
      INTEGER NTENS, NSTATV, NPROPS, NOEL, NPT, KSTEP, KINC
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),
     2 PREDEF,DPRED,PROPS(NPROPS),COORDS(3),DROT(3,3),PNEWDT,
     3 DFGRD0(3,3),DFGRD1(3,3)
C
      INTEGER I, J
      DOUBLE PRECISION STRANI(NTENS), E, nu
C
!       PRINT *, "INPUT VALUES"
!       PRINT *, "CMNAME (material's name) = ", CMNAME
!       PRINT *, "NOEL (element number) =", NOEL
!       PRINT *, "NPT (integration point number) =", NPT
!       PRINT *, "COORDS (coords of the integration point)"
!       DO I = 1, 3
!           PRINT *, COORDS(I)
!       END DO
!       PRINT *, "NPROPS (number of material constants) =", NPROPS
!       PRINT *, "PROPS (material constants)"
!       DO I = 1, NPROPS
!           PRINT *, PROPS(I)
!       END DO
!       PRINT *, "KSTEP (step number) =", KSTEP
!       PRINT *, "KINC (increment number) =", KINC
!       PRINT *, "TIME (step) =", TIME(1)
!       PRINT *, "TIME (total) =", TIME(2)
!       PRINT *, "DTIME (time increment) =", DTIME
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
!       PRINT *, "STRANI (total strain at the begining of the iteration)"
!       CALL ADD_VEC(STRANI, STRAN, DSTRAN, NTENS)
!       DO I = 1, NTENS
!           PRINT *, STRAN(I)
!       END DO
!       PRINT *, "LAYER (layer number) =", LAYER
!       PRINT *, "KSPT (section point within the current layer) =", KSPT
!       PRINT *, "SSE (specific elastic strain) =", SSE
!       PRINT *, "SPD (specific plastic dissipation) =", SPD
!       PRINT *, "SCD (specific creep dissipation) =", SCD
!       PRINT *, RPL
!       PRINT *, DDSDDT
!       PRINT *, DRPLDEDRPLDE
!       PRINT *, DRPLDTDRPLDT
!       PRINT *, TEMP
!       PRINT *, DTEMP
!       PRINT *, "PREDEF =", PREDEF(1)
!       PRINT *, "DPRED =", DPRED(1)
!       PRINT *, "NDI (number of direct stress components) =", NDI
!       PRINT *, "NSHR (number of engineering shear stress) =", NSHR
!       PRINT *, "NTENS (size of the stress or strain component) =", NTENS
!       PRINT *, "DROT (rotational increment matrix)"
!       DO I = 1, 3
!           DO J = 1, 3
!               PRINT *, DROT(I, J)
!           END DO
!       END DO
!       PRINT *, "CELENT (characteristic element length) =", CELENT
!       PRINT *, "DFGRD0 (deformation gradient at the beginning of the
!      1 increment)"
!       DO I = 1, 3
!           DO J = 1, 3
!               PRINT *, DFGRD0(I, J)
!           END DO
!       END DO
!       PRINT *, "DFGRD1 (deformation gradient at the end of the
!      1 increment)"
!       DO I = 1, 3
!           DO J = 1, 3
!               PRINT *, DFGRD1(I, J)
!           END DO
!       END DO
C
!       PRINT *, "CALCULATIONS"
      E = PROPS(1)
!       PRINT *, "E =", E
      nu = PROPS(2)
!       PRINT *, "nu =", nu
C
      CALL ADD_VEC(STRANI, STRAN, DSTRAN, NTENS)
      CALL ELASTICITY_3D_ISO_HOOKE(DDSDDE, E, nu)
      CALL PROD_MAT_VEC(STRESS, DDSDDE, STRANI, NTENS)
C
      PRINT *, "OUTPUT VALUES"
      PRINT *, "STRESS (stress at the end of the increment)"
      DO I = 1, NTENS
          PRINT *, STRESS(I)
      END DO
      PRINT *, "DDSDDE (jacobian matrix)"
      DO I = 1, NTENS
          DO J = 1, NTENS
              PRINT *, DDSDDE(I, J)
          END DO
      END DO
      PRINT *, "PNEWDT (suggested new time increment) =", PNEWDT
C
      RETURN
      END
C
C
C
