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
      #include 'ELASTICITY/ELASTICITY_3D_ORTHO_HOOKE.f'
!       #include 'ELASTICITY/ELASTICITY_3D_ORTHO_HOOKE_DIRECT.f'
C
C ABAQUS UMAT_3D_ELASTICITY_ORTHO SUBROUTINE
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
     1 G12, G13, G23
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
      CALL ADD_VEC(STRANI, STRAN, DSTRAN, NTENS)
      CALL ELASTICITY_3D_ORTHO_HOOKE(DDSDDE, E1, E2, E3, nu12, nu13, nu23,
     1 G12, G13, G23)
      CALL PROD_MAT_VEC(STRESS, DDSDDE, STRANI, NTENS)
C
      RETURN
      END
C
C
C
