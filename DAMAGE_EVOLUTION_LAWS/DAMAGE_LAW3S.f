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
C DAMAGE_LAW3S SUBROUTINE
C
      SUBROUTINE DAMAGE_LAW3S(d, Y, Y0, Y1, dseuil)
C
      DOUBLE PRECISION d, Y, Y0, Y1, dseuil
C
      IF (sqrt(Y) <= sqrt(Y0)) THEN
        d = 0.
      ELSEIF (sqrt(Y) > sqrt(Y0)) THEN
        d1 = dseuil * (sqrt(Y1) - sqrt(Y0)) / (2.*sqrt(Y1) - sqrt(Y0))
        IF (sqrt(Y) <= sqrt(Y1)) THEN
          d = d1 * (sqrt(Y) - sqrt(Y0)) / (sqrt(Y1) - sqrt(Y0))
        ELSEIF (sqrt(Y) > sqrt(Y1)) THEN
          d = dseuil - ((dseuil - d1) * sqrt(Y1) / sqrt(Y));
        ENDIF
      ENDIF
C
      END
C
C
C
