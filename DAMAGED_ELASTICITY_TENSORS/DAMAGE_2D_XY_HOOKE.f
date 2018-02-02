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
C DAMAGE_2D_XY_HOOKE SUBROUTINE
C
      SUBROUTINE DAMAGE_2D_XY_HOOKE(K, E1, E2, nu12, G12, d1, d2)
C
      INTEGER I, J, INFO
      DOUBLE PRECISION K(3,3), E1, E2, nu12, G12, d1, d2
C
!       PRINT *, "E1 =", E1
!       PRINT *, "E2 =", E2
!       PRINT *, "nu12 =", nu12
!       PRINT *, "G12 =", G12
!       PRINT *, "d1 =", d1
!       PRINT *, "d2 =", d2
C
      K(1,1) = 1./E1/(1-d1)
      K(1,2) = -nu12/E1
      K(1,3) = 0.

      K(2,1) = -nu12/E1
      K(2,2) = 1./E2/(1-d2)
      K(2,3) = 0.

      K(3,1) = 0.
      K(3,2) = 0.
      K(3,3) = 1./2./G12
C
      CALL DPOTRF('U', 3, K, 3, INFO)
      CALL DPOTRI('U', 3, K, 3, INFO)
      DO I = 1, 3
        DO J = 1, I
          K(I, J) = K(J, I)
!           PRINT *, 'K(', I, ' ,', J, ' ) =', K(J, I)
        ENDDO
      ENDDO
C
      END
C
C
C
