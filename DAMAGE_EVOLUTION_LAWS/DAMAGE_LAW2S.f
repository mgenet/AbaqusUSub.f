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
C DAMAGE_LAW2S SUBROUTINE
C
      SUBROUTINE DAMAGE_LAW2S(d, Y, Y0, dseuil)
C
      DOUBLE PRECISION d, Y, Y0, dseuil
C
      IF (sqrt(Y) < sqrt(Y0)) THEN
        d = 0.
      ELSE
        d = dseuil * (sqrt(Y) - sqrt(Y0)) / sqrt(Y)
      ENDIF
C
      END
C
C
C
