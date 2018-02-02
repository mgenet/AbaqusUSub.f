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
      #include 'DAMAGE_LAWS/DAMAGE_LAW3S.f'
C
C DAMAGE_LAW3 SUBROUTINE
C
      SUBROUTINE DAMAGE_LAW3(d, Y, Y0, Y1)
C
      DOUBLE PRECISION d, Y, Y0, Y1, dseuil
C
      dseuil = 1.
      CALL DAMAGE_LAW3S(d, Y, Y0, Y1, dseuil)
C
      END
C
C
C
