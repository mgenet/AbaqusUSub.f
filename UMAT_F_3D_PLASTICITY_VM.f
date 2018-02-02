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

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      #include 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),PROPS(NPROPS)
C
      DIMENSION STRESSOLD(NTENS),DSTATEV(NTENS)
C
      real E, nu, G, lambda, mu, dlambda
C
      PRINT *, "\n INPUT VALUES"
      PRINT *, "NOEL (element number) = ", NOEL
      PRINT *, "NPT (integration point number) = ", NPT
      PRINT *, "KSTEP (step number) = ", KSTEP
      PRINT *, "KINC (increment number) = ", KINC
      PRINT *, "NSTATV (number of state variables) = ", NSTATV
      PRINT *, "STATEV (state variables)"
      DO I = 1, NSTATV
          PRINT *, STATEV(I)
      END DO
      PRINT *, "STRAN (total strain)"
      DO I = 1, NTENS
          PRINT *, STRAN(I)
      END DO
      PRINT *, "DSTRAN (total strain increment)"
      DO I = 1, NTENS
          PRINT *, DSTRAN(I)
      END DO
C
!       E1 = PROPS(1)
!       E2 = PROPS(2)
!       E3 = PROPS(3)
!       nu12 = PROPS(4)
!       nu13 = PROPS(5)
!       nu23 = PROPS(6)
!       G12 = PROPS(7)
!       G13 = PROPS(8)
!       G23 = PROPS(9)
!       sigmac = PROPS(10)
C
!       epsilonp11 = STATEV(1)
!       epsilonp22 = STATEV(2)
!       epsilonp33 = STATEV(3)
!       epsilonp12 = STATEV(4)
!       epsilonp13 = STATEV(5)
!       epsilonp23 = STATEV(6)
C
      E = PROPS(1)
      PRINT *, "E = ", E
      nu = PROPS(2)
      PRINT *, "nu = ", nu
      G = PROPS(1) / 2 / (1 + PROPS(2))
      PRINT *, "G = ", G
      sigmac = PROPS(3)
C
      lambda = (E * nu) / ((1 + nu) * (1 - (2*nu)))
      PRINT *, "lambda = ", lambda
      mu = E / (2 * (1 + nu))
      PRINT *, "mu = ", mu
C
      DDSDDE(1,1) = lambda + 2 * mu
      DDSDDE(1,2) = lambda
      DDSDDE(1,3) = lambda
      DDSDDE(1,4) = 0.
      DDSDDE(1,5) = 0.
      DDSDDE(1,6) = 0.
      DDSDDE(2,1) = lambda
      DDSDDE(2,2) = lambda + 2 * mu
      DDSDDE(2,3) = lambda
      DDSDDE(2,4) = 0.
      DDSDDE(2,5) = 0.
      DDSDDE(2,6) = 0.
      DDSDDE(3,1) = lambda
      DDSDDE(3,2) = lambda
      DDSDDE(3,3) = lambda + 2 * mu
      DDSDDE(3,4) = 0.
      DDSDDE(3,5) = 0.
      DDSDDE(3,6) = 0.
      DDSDDE(4,1) = 0.
      DDSDDE(4,2) = 0.
      DDSDDE(4,3) = 0.
      DDSDDE(4,4) = 2 * mu
      DDSDDE(4,5) = 0.
      DDSDDE(4,6) = 0.
      DDSDDE(5,1) = 0.
      DDSDDE(5,2) = 0.
      DDSDDE(5,3) = 0.
      DDSDDE(5,4) = 0.
      DDSDDE(5,5) = 2 * mu
      DDSDDE(5,6) = 0.
      DDSDDE(6,1) = 0.
      DDSDDE(6,2) = 0.
      DDSDDE(6,3) = 0.
      DDSDDE(6,4) = 0.
      DDSDDE(6,5) = 0.
      DDSDDE(6,6) = 2 * mu
      sigmac = PROPS(3)
C
!       DDSDDE(1,1) = (nu23**2*E1**2*E3**2-E1**2*E2*E3)/(nu23**2*E1*
!      1 E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*
!      2 E3+nu13**2*E1**2*E2)
!       DDSDDE(1,2) = -((nu12*E1*E2**2+nu13*nu23*E1**2*E2)*E3)/(nu23**2*
!      1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*
!      2 E3+nu13**2*E1**2*E2)
!       DDSDDE(1,3) = -(nu12*nu23*E1*E2*E3**2+nu13*E1**2*E2*E3)/(nu23**2*
!      1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
!      2 E1**2*E2)
!       DDSDDE(1,4) = 0.
!       DDSDDE(1,5) = 0.
!       DDSDDE(1,6) = 0.
!
!       DDSDDE(2,1) = -((nu12*E1*E2**2+nu13*nu23*E1**2*E2)*E3)/(nu23**2*
!      1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
!      2 E1**2*E2)
!       DDSDDE(2,2) = -(E1*E2**2*E3-nu13**2*E1**2*E2**2)/(nu23**2*E1*
!      1 E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
!      2 E1**2*E2)
!       DDSDDE(2,3) = -(nu23*E1*E2*E3**2+nu12*nu13*E1*E2**2*E3)/(nu23**2*
!      1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
!      2 E1**2*E2)
!       DDSDDE(2,4) = 0.
!       DDSDDE(2,5) = 0.
!       DDSDDE(2,6) = 0.
!
!       DDSDDE(3,1) = -(nu12*nu23*E1*E2*E3**2+nu13*E1**2*E2*E3)/(nu23**2*
!      1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
!      2 E1**2*E2)
!       DDSDDE(3,2) = -(nu23*E1*E2*E3**2+nu12*nu13*E1*E2**2*E3)/(nu23**2*
!      1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
!      2 E1**2*E2)
!       DDSDDE(3,3) = ((nu12**2*E2**2-E1*E2)*E3**2)/(nu23**2*E1*
!      1 E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
!      2 E1**2*E2)
!       DDSDDE(3,4) = 0.
!       DDSDDE(3,5) = 0.
!       DDSDDE(3,6) = 0.
!
!       DDSDDE(4,1) = 0.
!       DDSDDE(4,2) = 0.
!       DDSDDE(4,3) = 0.
!       DDSDDE(4,4) = 2 * G12
!       DDSDDE(4,5) = 0.
!       DDSDDE(4,6) = 0.
!
!       DDSDDE(5,1) = 0.
!       DDSDDE(5,2) = 0.
!       DDSDDE(5,3) = 0.
!       DDSDDE(5,4) = 0.
!       DDSDDE(5,5) = 2 * G13
!       DDSDDE(5,6) = 0.
!
!       DDSDDE(6,1) = 0.
!       DDSDDE(6,2) = 0.
!       DDSDDE(6,3) = 0.
!       DDSDDE(6,4) = 0.
!       DDSDDE(6,5) = 0.
!       DDSDDE(6,6) = 2 * G23
C
      DO I = 1, NTENS
        STRESSOLD(I) = STRESS(I)
        DO J = 1, NTENS
          STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
        ENDDO
      ENDDO
C
!       q = ((STRESS(3)+DSTRESS(3))-((STRESS(3)+DSTRESS(3))+
!      1 (STRESS(2)+DSTRESS(2))+(STRESS(1)+DSTRESS(1)))/3)**2
!       q = q + ((STRESS(2)+DSTRESS(2))-((STRESS(3)+DSTRESS(3))+
!      1 (STRESS(2)+DSTRESS(2))+(STRESS(1)+DSTRESS(1)))/3)**2
!       q = q + ((STRESS(1)+DSTRESS(1))-((STRESS(3)+DSTRESS(3))+
!      1 (STRESS(2)+DSTRESS(2))+(STRESS(1)+DSTRESS(1)))/3)**2
!       q = q + (STRESS(6)+DSTRESS(6))**2+(STRESS(5)+DSTRESS(5))**2+
!      1 (STRESS(6)+DSTRESS(6))**2+(STRESS(4)+DSTRESS(4))**2+
!      2 (STRESS(5)+DSTRESS(5))**2+(STRESS(4)+DSTRESS(4))**2
!       PRINT *, q
C
      q = (STRESS(3)-(STRESS(3)+STRESS(2)+STRESS(1))/3)**2+
     1 (STRESS(2)-(STRESS(3)+STRESS(2)+STRESS(1))/3)**2+
     2 (STRESS(1)-(STRESS(3)+STRESS(2)+STRESS(1))/3)**2+
     3 2*STRESS(6)**2+2*STRESS(5)**2+2*STRESS(4)**2
!       PRINT *, q
C
      IF (sqrt(3*q/2) - sigmac > 0) THEN
10      IF (abs(sqrt(3*q/2) - sigmac) / sigmac > 1e-6) THEN
          PRINT *, "q = ", sqrt(3*q/2)
          PRINT *, "sigmac = ", sigmac
          dlambda = (sqrt(3*q/2) - sigmac) / 2 / G / sqrt(3*q/2)
          PRINT *, "dlambda = ", dlambda
C
          DSTATEV(1)=dlambda*(STRESSOLD(1)-(STRESSOLD(3)+
         1 STRESSOLD(2)+STRESSOLD(1))/3)
          DSTATEV(2)=dlambda*(STRESSOLD(2)-(STRESSOLD(3)+
         2 STRESSOLD(2)+STRESSOLD(1))/3)
          DSTATEV(3)=dlambda*(STRESSOLD(3)-(STRESSOLD(3)+
         3 STRESSOLD(2)+STRESSOLD(1))/3)
          DSTATEV(4)=dlambda*STRESSOLD(4)
          DSTATEV(5)=dlambda*STRESSOLD(5)
          DSTATEV(6)=dlambda*STRESSOLD(6)
C
!         PRINT *, DSTATEV(1)
!         PRINT *, DSTATEV(2)
!         PRINT *, DSTATEV(3)
!         PRINT *, DSTATEV(4)
!         PRINT *, DSTATEV(5)
!         PRINT *, DSTATEV(6)
C
!         STATEV(1)=STATEV(1)+dlambda*(STRESSOLD(1)-(STRESSOLD(3)+
!        1 STRESSOLD(2)+STRESSOLD(1))/3)
!         STATEV(2)=STATEV(2)+dlambda*(STRESSOLD(2)-(STRESSOLD(3)+
!        2 STRESSOLD(2)+STRESSOLD(1))/3)
!         STATEV(3)=STATEV(3)+dlambda*(STRESSOLD(3)-(STRESSOLD(3)+
!        3 STRESSOLD(2)+STRESSOLD(1))/3)
!         STATEV(4)=STATEV(4)+dlambda*STRESSOLD(4)
!         STATEV(5)=STATEV(5)+dlambda*STRESSOLD(5)
!         STATEV(6)=STATEV(6)+dlambda*STRESSOLD(6)
C
          STATEV(1) = STATEV(1) + DSTATEV(1)
          STATEV(2) = STATEV(2) + DSTATEV(2)
          STATEV(3) = STATEV(3) + DSTATEV(3)
          STATEV(4) = STATEV(4) + DSTATEV(4)
          STATEV(5) = STATEV(5) + DSTATEV(5)
          STATEV(6) = STATEV(6) + DSTATEV(6)
C
!         PRINT *, STATEV(1)
!         PRINT *, STATEV(2)
!         PRINT *, STATEV(3)
!         PRINT *, STATEV(4)
!         PRINT *, STATEV(5)
!         PRINT *, STATEV(6)
C
          DO I = 1, NTENS
            STRESSOLD(I) = STRESS(I)
            DO J = 1, NTENS
              STRESS(I) = STRESS(I)-DDSDDE(I,J)*DSTATEV(J)
            ENDDO
          ENDDO
C
          q = (STRESS(3)-(STRESS(3)+STRESS(2)+STRESS(1))/3)**2+
         1 (STRESS(2)-(STRESS(3)+STRESS(2)+STRESS(1))/3)**2+
         2 (STRESS(1)-(STRESS(3)+STRESS(2)+STRESS(1))/3)**2+
         3 2*STRESS(6)**2+2*STRESS(5)**2+2*STRESS(4)**2
C
        GOTO 10
        ENDIF
      ENDIF
C
!       PRINT *, "OUTPUT VALUES"
!       PRINT *, "STRESS (stress at the end of the increment)"
!       DO I = 1, NTENS
!           PRINT *, STRESS(I)
!       END DO
C
!       PRINT *, "DDSDDE (jacobian matrix)"
!       DO I = 1, NTENS
!           DO J = 1, NTENS
!               PRINT *, DDSDDE(I, J)
!           END DO
!       END DO
C
      RETURN
      END












