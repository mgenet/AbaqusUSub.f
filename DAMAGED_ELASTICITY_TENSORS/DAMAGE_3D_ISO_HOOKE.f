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
C DAMAGE_3D_ISO_HOOKE SUBROUTINE
C
      SUBROUTINE DAMAGE_3D_ISO_HOOKE(K, E, nu, d)
C
      INTEGER I, J, INFO
      DOUBLE PRECISION K(6,6), E, nu, d
C
      K(1,1) = 1./E/(1-d)
      K(1,2) = -nu/E/(1-d)
      K(1,3) = -nu/E/(1-d)
      K(1,4) = 0.
      K(1,5) = 0.
      K(1,6) = 0.

      K(2,1) = -nu/E/(1-d)
      K(2,2) = 1./E/(1-d)
      K(2,3) = -nu/E/(1-d)
      K(2,4) = 0.
      K(2,5) = 0.
      K(2,6) = 0.

      K(3,1) = -nu/E/(1-d)
      K(3,2) = -nu/E/(1-d)
      K(3,3) = 1./E/(1-d)
      K(3,4) = 0.
      K(3,5) = 0.
      K(3,6) = 0.

      K(4,1) = 0.
      K(4,2) = 0.
      K(4,3) = 0.
      K(4,4) = (1.+nu)/E/(1-d)
      K(4,5) = 0.
      K(4,6) = 0.

      K(5,1) = 0.
      K(5,2) = 0.
      K(5,3) = 0.
      K(5,4) = 0.
      K(5,5) = (1.+nu)/E/(1-d)
      K(5,6) = 0.

      K(6,1) = 0.
      K(6,2) = 0.
      K(6,3) = 0.
      K(6,4) = 0.
      K(6,5) = 0.
      K(6,6) = (1.+nu)/E/(1-d)
C
      CALL DPOTRF('U', 6, K, 6, INFO)
      CALL DPOTRI('U', 6, K, 6, INFO)
      DO I = 1, 6
        DO J = 1, I
          K(I, J) = K(J, I)
        ENDDO
      ENDDO
C
      END
C
C
C
