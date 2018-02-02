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
      #include 'ELASTICITY/ELASTICITY_3D.f'
C
C ABAQUS UMAT_EXT_ELASTICITY SUBROUTINE
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
      INTEGER I, J
      DOUBLE PRECISION STRANI(NTENS), E, nu
C
!       PRINT *, "KINC (increment number) =", KINC
!       PRINT *, "NOEL (element number) =", NOEL
!       PRINT *, "NPT (integration point number) =", NPT
!       PRINT *, "KSTEP (step number) =", KSTEP
C
      CALL ADD_VEC(STRANI, STRAN, DSTRAN, NTENS)
      CALL ELASTICITY_3D(DDSDDE, NOEL, NTENS)
      CALL PROD_MAT_VEC(STRESS, DDSDDE, STRANI, NTENS)
C
      RETURN
      END
C
C
C
