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
C DAMAGE_LAW1 SUBROUTINE
C
      SUBROUTINE DAMAGE_LAW1(d, Y, Y0, Y1)
C
      DOUBLE PRECISION d, Y, Y0, Y1
C
      IF (sqrt(Y) < sqrt(Y0)) THEN
        d = 0.
      ELSEIF (sqrt(Y) > sqrt(Y1)) THEN
        d = 1.
      ELSE
        d = (sqrt(Y) - sqrt(Y0)) / (sqrt(Y1) - sqrt(Y0))
      ENDIF
C
      END
C
C
C
