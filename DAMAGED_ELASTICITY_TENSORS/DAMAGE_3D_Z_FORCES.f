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
C DAMAGE_3D_Z_FORCES SUBROUTINE
C
      SUBROUTINE DAMAGE_3D_Z_FORCES(Y3, EPS, E1, E2, E3, nu12, nu13,
     1 nu23, G12, G13, G23, d3)
C
      DOUBLE PRECISION Y3, EPS(6), E1, E2, E3, nu12, nu13, nu23,
     1 G12, G13, G23, d3
C
      Y3 = -2*EPS(2)*EPS(3)*nu13*nu12**3*E2**2
      Y3 = Y3 +2*EPS(1)*E1**2*EPS(3)*nu13
      Y3 = Y3 +EPS(1)**2*E1**2*nu23**2*nu12**2
      Y3 = Y3 +EPS(2)**2*nu13**2*E2**2*nu12**2
      Y3 = Y3 +2*EPS(2)*EPS(3)*nu23*E1**2
      Y3 = Y3 -2*EPS(3)**2*E1*nu12**2*E2
      Y3 = Y3 +EPS(1)**2*E1**2*nu13**2
      Y3 = Y3 +EPS(2)**2*E1**2*nu23**2
      Y3 = Y3 +EPS(3)**2*nu12**4*E2**2
      Y3 = Y3 +2*EPS(1)**2*E1**2*nu12*nu23*nu13
      Y3 = Y3 +2*EPS(1)*E1**2*EPS(2)*nu12*nu23**2
      Y3 = Y3 +2*EPS(1)*E1*EPS(2)*E2*nu12*nu13**2
      Y3 = Y3 -2*EPS(1)*E1*EPS(3)*E2*nu12**3*nu23
      Y3 = Y3 -2*EPS(1)*E1*EPS(3)*E2*nu13*nu12**2
      Y3 = Y3 +2*EPS(1)*E1**2*EPS(2)*nu23*nu13
      Y3 = Y3 +2*EPS(1)*E1*EPS(2)*nu23*nu13*E2*nu12**2
      Y3 = Y3 +2*EPS(1)*E1**2*EPS(3)*nu12*nu23
      Y3 = Y3 +2*EPS(2)**2*E2*E1*nu12*nu23*nu13
      Y3 = Y3 +EPS(3)**2*E1**2
      Y3 = Y3 -2*EPS(2)*EPS(3)*nu23*E1*nu12**2*E2
      Y3 = Y3 +2*EPS(2)*EPS(3)*nu13*nu12*E2*E1
      Y3 = Y3 * (1./2.) * E3 * E2**2
      Y3 = Y3 /(E1*E2-E1*nu23**2*E3+E1*nu23**2*E3*d3-nu12**2*E2**2-
     2 2*nu12*E2*nu23*nu13*E3+2*nu12*E2*nu23*nu13*E3*d3-
     3 nu13**2*E2*E3+nu13**2*E2*E3*d3)**2
C
      END
C
C
C
