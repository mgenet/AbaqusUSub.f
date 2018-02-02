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
C ELASTICITY_3D_ORTHO_HOOKE_DIRECT SUBROUTINE
C
      SUBROUTINE ELASTICITY_3D_ORTHO_HOOKE(K, E1, E2, E3, nu12, nu13,
     1 nu23, G12, G13, G23)
C
      DOUBLE PRECISION K(6,6), E1, E2, E3, nu12, nu13, nu23, G12, G13,
     1 G23
C
      K(1,1) = (nu23**2*E1**2*E3**2-E1**2*E2*E3)/(nu23**2*E1*
     1 E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*
     2 E3+nu13**2*E1**2*E2)
      K(1,2) = -((nu12*E1*E2**2+nu13*nu23*E1**2*E2)*E3)/(nu23**2*
     1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*
     2 E3+nu13**2*E1**2*E2)
      K(1,3) = -(nu12*nu23*E1*E2*E3**2+nu13*E1**2*E2*E3)/(nu23**2*
     1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
     2 E1**2*E2)
      K(1,4) = 0.
      K(1,5) = 0.
      K(1,6) = 0.

      K(2,1) = -((nu12*E1*E2**2+nu13*nu23*E1**2*E2)*E3)/(nu23**2*
     1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
     2 E1**2*E2)
      K(2,2) = -(E1*E2**2*E3-nu13**2*E1**2*E2**2)/(nu23**2*E1*
     1 E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
     2 E1**2*E2)
      K(2,3) = -(nu23*E1*E2*E3**2+nu12*nu13*E1*E2**2*E3)/(nu23**2*
     1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
     2 E1**2*E2)
      K(2,4) = 0.
      K(2,5) = 0.
      K(2,6) = 0.

      K(3,1) = -(nu12*nu23*E1*E2*E3**2+nu13*E1**2*E2*E3)/(nu23**2*
     1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
     2 E1**2*E2)
      K(3,2) = -(nu23*E1*E2*E3**2+nu12*nu13*E1*E2**2*E3)/(nu23**2*
     1 E1*E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
     2 E1**2*E2)
      K(3,3) = ((nu12**2*E2**2-E1*E2)*E3**2)/(nu23**2*E1*
     1 E3**2+(nu12**2*E2**2+(2*nu12*nu13*nu23-1)*E1*E2)*E3+nu13**2*
     2 E1**2*E2)
      K(3,4) = 0.
      K(3,5) = 0.
      K(3,6) = 0.

      K(4,1) = 0.
      K(4,2) = 0.
      K(4,3) = 0.
      K(4,4) = 2 * G12
      K(4,5) = 0.
      K(4,6) = 0.

      K(5,1) = 0.
      K(5,2) = 0.
      K(5,3) = 0.
      K(5,4) = 0.
      K(5,5) = 2 * G13
      K(5,6) = 0.

      K(6,1) = 0.
      K(6,2) = 0.
      K(6,3) = 0.
      K(6,4) = 0.
      K(6,5) = 0.
      K(6,6) = 2 * G23
C
      END
C
C
C
