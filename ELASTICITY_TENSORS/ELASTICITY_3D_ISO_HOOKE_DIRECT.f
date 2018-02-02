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
C ELASTICITY_3D_ISO_HOOKE_DIRECT SUBROUTINE
C
      SUBROUTINE ELASTICITY_3D_ISO_HOOKE(K, E, nu)
C
      DOUBLE PRECISION K(6,6), E, nu, lambda, mu
C
      lambda = (E * nu) / ((1 + nu) * (1 - (2*nu)))
      mu = E / (2 * (1 + nu))
C
      K(1,1) = lambda + 2 * mu
      K(1,2) = lambda
      K(1,3) = lambda
      K(1,4) = 0.
      K(1,5) = 0.
      K(1,6) = 0.
      K(2,1) = lambda
      K(2,2) = lambda + 2 * mu
      K(2,3) = lambda
      K(2,4) = 0.
      K(2,5) = 0.
      K(2,6) = 0.
      K(3,1) = lambda
      K(3,2) = lambda
      K(3,3) = lambda + 2 * mu
      K(3,4) = 0.
      K(3,5) = 0.
      K(3,6) = 0.
      K(4,1) = 0.
      K(4,2) = 0.
      K(4,3) = 0.
      K(4,4) = 2 * mu
      K(4,5) = 0.
      K(4,6) = 0.
      K(5,1) = 0.
      K(5,2) = 0.
      K(5,3) = 0.
      K(5,4) = 0.
      K(5,5) = 2 * mu
      K(5,6) = 0.
      K(6,1) = 0.
      K(6,2) = 0.
      K(6,3) = 0.
      K(6,4) = 0.
      K(6,5) = 0.
      K(6,6) = 2 * mu
C
      END
C
C
C
