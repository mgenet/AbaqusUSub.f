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
      #include 'ELASTICITY/ELASTICITY_3D_ORTHO_HOOKE.f'
C
C ELASTICITY_3D_ISO_HOOKE SUBROUTINE
C
      SUBROUTINE ELASTICITY_3D_ISO_HOOKE(K, E, nu)
C
      DOUBLE PRECISION K(6,6), E, nu, G
C
      G = E / 2. / (1+nu)
C
      CALL ELASTICITY_3D_ORTHO_HOOKE(K, E, E, E, nu, nu, nu, G, G, G)
C
      END
C
C
C
