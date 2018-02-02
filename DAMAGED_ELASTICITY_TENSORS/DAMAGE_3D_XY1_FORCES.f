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
C DAMAGE_3D_XY1_FORCES SUBROUTINE
C
      SUBROUTINE DAMAGE_3D_XY1_FORCES(Y1, Y2, EPS, E1, E2, E3, nu12,
     1 nu13, nu23, G12, G13, G23, d1, d2)
C
      DOUBLE PRECISION Y1, Y2, EPS(6), E1, E2, E3, nu12, nu13, nu23,
     1 G12, G13, G23, d1, d2
C
      Y1 = EPS(2)**2*E2**4*nu12**2
      Y1 = Y1 +2*EPS(1)**2*E1**2*E2*nu23**2*E3*d2
      Y1 = Y1 +2*EPS(1)*E1*EPS(2)*E2**3*nu12
      Y1 = Y1 +EPS(1)**2*E1**2*E2**2
      Y1 = Y1 +4*EPS(1)*E1*EPS(3)*E2*E3**2*nu12*nu23**3*d2
      Y1 = Y1 -2*EPS(1)*E1*EPS(2)*E2**2*nu12*nu23**2*E3
      Y1 = Y1 +4*EPS(1)*E1*EPS(2)*E2**2*nu12*nu23**2*E3*d2
      Y1 = Y1 +2*EPS(1)*E1*EPS(3)*E2*E3**2*nu13*nu23**2*d2
      Y1 = Y1 +2*EPS(1)*E1*EPS(3)*E2**2*E3*nu12*nu23
      Y1 = Y1 -2*EPS(1)*E1*EPS(3)*E2*E3**2*nu12*nu23**3
      Y1 = Y1 -2*EPS(1)*E1*EPS(3)*E2*E3**2*nu13*nu23**2
      Y1 = Y1 -2*EPS(1)*E1*EPS(3)*E2**2*E3*nu12*nu23*d2
      Y1 = Y1 -2*EPS(1)*E1*EPS(3)*E2*E3**2*nu12*nu23**3*d2**2
      Y1 = Y1 +2*EPS(1)*E1*EPS(3)*E2**2*E3*nu13
      Y1 = Y1 -2*EPS(1)**2*E1**2*nu23**4*E3**2*d2
      Y1 = Y1 +EPS(1)**2*E1**2*nu23**4*E3**2*d2**2
      Y1 = Y1 -2*EPS(1)**2*E1**2*E2*nu23**2*E3
      Y1 = Y1 -2*EPS(1)*E1*EPS(2)*E2**2*nu12*d2**2*nu23**2*E3
      Y1 = Y1 +4*EPS(1)*E1*EPS(2)*E2*nu23**3*nu13*E3**2*d2
      Y1 = Y1 -2*EPS(1)*E1*EPS(2)*E2*nu23**3*nu13*E3**2*d2**2
      Y1 = Y1 -2*EPS(1)*E1*EPS(2)*E2**3*nu12*d2
      Y1 = Y1 +2*EPS(1)*E1*EPS(2)*E2**2*nu23*nu13*E3
      Y1 = Y1 -2*EPS(1)*E1*EPS(2)*E2**2*nu23*nu13*E3*d2
      Y1 = Y1 -2*EPS(1)*E1*EPS(2)*E2*nu23**3*nu13*E3**2
      Y1 = Y1 +2*EPS(2)*E2**2*EPS(3)*E3**2*nu23*nu13**2
      Y1 = Y1 -4*EPS(2)*E2**3*EPS(3)*E3*nu23*nu12**2*d2
      Y1 = Y1 -2*EPS(2)**2*E2**4*nu12**2*d2
      Y1 = Y1 +EPS(3)**2*E2**2*E3**2*nu13**2
      Y1 = Y1 +EPS(2)**2*E2**4*nu12**2*d2**2
      Y1 = Y1 +2*EPS(2)**2*E2**3*nu12*nu23*nu13*E3
      Y1 = Y1 -4*EPS(2)**2*E2**3*nu12*nu23*nu13*E3*d2
      Y1 = Y1 +2*EPS(2)*E2**3*EPS(3)*E3*nu13*nu12
      Y1 = Y1 +2*EPS(2)*E2**2*EPS(3)*E3**2*nu13*nu12*nu23**2
      Y1 = Y1 -4*EPS(2)*E2**2*EPS(3)*E3**2*nu13*nu12*nu23**2*d2
      Y1 = Y1 +2*EPS(2)*E2**3*EPS(3)*E3*nu23*nu12**2
      Y1 = Y1 +EPS(2)**2*E2**2*nu13**2*E3**2*nu23**2
      Y1 = Y1 -2*EPS(2)**2*E2**2*nu13**2*E3**2*nu23**2*d2
      Y1 = Y1 +2*EPS(2)*E2**3*EPS(3)*E3*nu23*nu12**2*d2**2
      Y1 = Y1 +2*EPS(2)**2*E2**3*nu12*nu23*nu13*E3*d2**2
      Y1 = Y1 +2*EPS(2)*E2**2*EPS(3)*E3**2*nu13*nu12*nu23**2*d2**2
      Y1 = Y1 +EPS(2)**2*E2**2*nu13**2*E3**2*nu23**2*d2**2
      Y1 = Y1 -2*EPS(2)*E2**2*d2*EPS(3)*E3**2*nu23*nu13**2
      Y1 = Y1 -2*EPS(2)*E2**3*d2*EPS(3)*E3*nu13*nu12
      Y1 = Y1 -2*EPS(3)**2*E2**2*E3**2*nu12*nu23*nu13*d2
      Y1 = Y1 -2*EPS(3)**2*E2**2*E3**2*nu12**2*nu23**2*d2
      Y1 = Y1 +EPS(3)**2*E2**2*E3**2*nu12**2*nu23**2
      Y1 = Y1 +EPS(3)**2*E2**2*E3**2*nu12**2*nu23**2*d2**2
      Y1 = Y1 +2*EPS(3)**2*E2**2*E3**2*nu12*nu23*nu13
      Y1 = Y1 +EPS(1)**2*E1**2*nu23**4*E3**2
      Y1 = Y1 * (1./2.) * E1
      Y1 = Y1 / (E1*E2-E1*nu23**2*E3+E1*nu23**2*E3*d2-nu12**2*E2**2+
     2 nu12**2*E2**2*d2+nu12**2*E2**2*d1-nu12**2*E2**2*d1*d2-
     3 2*nu12*E2*nu23*nu13*E3+2*nu12*E2*nu23*nu13*E3*d2+
     4 2*nu12*E2*nu23*nu13*E3*d1-2*nu12*E2*nu23*nu13*E3*d1*d2-
     5 nu13**2*E2*E3+nu13**2*E2*E3*d1)**2
C
      Y2 = 2*EPS(2)*E2*EPS(3)*E3*nu23*E1**2
      Y2 = Y2 +EPS(2)**2*E2**2*nu13**4*E3**2*d1**2
      Y2 = Y2 -2*EPS(2)**2*E2**2*E1*nu13**2*E3
      Y2 = Y2 +2*EPS(2)**2*E2**2*E1*nu13**2*E3*d1
      Y2 = Y2 -2*EPS(2)*E2**2*EPS(3)*E3**2*nu13**3*nu12
      Y2 = Y2 +EPS(2)**2*E2**2*E1**2
      Y2 = Y2 +2*EPS(1)*E1**2*EPS(3)*E3**2*nu13*nu23**2
      Y2 = Y2 +EPS(1)**2*E1**2*nu12**2*E2**2
      Y2 = Y2 +EPS(3)**2*E3**2*E1**2*nu23**2
      Y2 = Y2 +2*EPS(3)**2*E3**2*E1*nu12*E2*nu23*nu13
      Y2 = Y2 -2*EPS(3)**2*E3**2*E1*nu12*E2*nu23*nu13*d1
      Y2 = Y2 +2*EPS(1)*E1**2*EPS(2)*E2**2*nu12
      Y2 = Y2 -2*EPS(3)**2*E3**2*nu12**2*E2**2*nu13**2*d1
      Y2 = Y2 +EPS(3)**2*E3**2*nu12**2*E2**2*nu13**2
      Y2 = Y2 +EPS(3)**2*E3**2*nu12**2*E2**2*nu13**2*d1**2
      Y2 = Y2 -4*EPS(1)*E1*EPS(3)*nu12*nu23*E2*E3**2*nu13**2*d1
      Y2 = Y2 +2*EPS(1)*E1*EPS(3)*E2**2*E3*nu13*nu12**2
      Y2 = Y2 -2*EPS(1)**2*E1**2*nu23**2*E3**2*nu13**2*d1
      Y2 = Y2 +EPS(1)**2*E1**2*nu23**2*E3**2*nu13**2
      Y2 = Y2 -2*EPS(1)**2*E1**2*nu12**2*E2**2*d1
      Y2 = Y2 +2*EPS(1)**2*E1**2*nu12*E2*nu23*nu13*E3
      Y2 = Y2 -4*EPS(1)**2*E1**2*nu12*E2*nu23*nu13*E3*d1
      Y2 = Y2 +2*EPS(1)*E1**2*EPS(3)*nu12*nu23*E2*E3
      Y2 = Y2 +2*EPS(1)*E1*EPS(3)*nu12*nu23*E2*E3**2*nu13**2
      Y2 = Y2 -4*EPS(1)*E1*EPS(3)*E2**2*E3*nu13*nu12**2*d1
      Y2 = Y2 +2*EPS(1)*E1**2*EPS(2)*E2*nu23*nu13*E3
      Y2 = Y2 -2*EPS(1)*E1*EPS(2)*E2*nu23*nu13**3*E3**2
      Y2 = Y2 +4*EPS(1)*E1*EPS(2)*E2*nu23*nu13**3*E3**2*d1
      Y2 = Y2 -2*EPS(1)*E1*EPS(2)*E2**2*nu12*nu13**2*E3
      Y2 = Y2 +4*EPS(1)*E1*EPS(2)*E2**2*nu12*nu13**2*E3*d1
      Y2 = Y2 +2*EPS(1)**2*E1**2*nu12*E2*nu23*nu13*E3*d1**2
      Y2 = Y2 +EPS(2)**2*E2**2*nu13**4*E3**2
      Y2 = Y2 -2*EPS(1)*E1*EPS(2)*E2*nu23*nu13**3*E3**2*d1**2
      Y2 = Y2 -2*EPS(1)*E1*EPS(2)*E2**2*nu12*nu13**2*E3*d1**2
      Y2 = Y2 -2*EPS(1)*E1**2*d1*EPS(3)*nu12*nu23*E2*E3
      Y2 = Y2 -2*EPS(1)*E1**2*d1*EPS(3)*E3**2*nu13*nu23**2
      Y2 = Y2 +2*EPS(1)*E1*EPS(3)*nu12*nu23*E2*E3**2*nu13**2*d1**2
      Y2 = Y2 +2*EPS(1)*E1*EPS(3)*E2**2*E3*nu13*nu12**2*d1**2
      Y2 = Y2 +4*EPS(2)*E2**2*EPS(3)*E3**2*nu13**3*nu12*d1
      Y2 = Y2 -2*EPS(2)*E2*EPS(3)*E3**2*nu23*E1*nu13**2
      Y2 = Y2 -2*EPS(2)*E2**2*EPS(3)*E3**2*nu13**3*nu12*d1**2
      Y2 = Y2 +2*EPS(2)*E2*EPS(3)*E3**2*nu23*E1*nu13**2*d1
      Y2 = Y2 -2*EPS(2)*E2**2*EPS(3)*E3*nu13*nu12*d1*E1
      Y2 = Y2 +2*EPS(2)*E2**2*EPS(3)*E3*nu13*nu12*E1
      Y2 = Y2 -2*EPS(2)**2*E2**2*nu13**4*E3**2*d1
      Y2 = Y2 +EPS(1)**2*E1**2*nu23**2*E3**2*nu13**2*d1**2
      Y2 = Y2 +EPS(1)**2*E1**2*nu12**2*E2**2*d1**2
      Y2 = Y2 -2*EPS(1)*E1**2*d1*EPS(2)*E2**2*nu12
      Y2 = Y2 -2*EPS(1)*E1**2*d1*EPS(2)*E2*nu23*nu13*E3
      Y2 = Y2 * (1./2.) * E2
      Y2 = Y2 / (E1*E2-E1*nu23**2*E3+E1*nu23**2*E3*d2-nu12**2*E2**2+
     2 nu12**2*E2**2*d2+nu12**2*E2**2*d1-nu12**2*E2**2*d1*d2-
     3 2*nu12*E2*nu23*nu13*E3+2*nu12*E2*nu23*nu13*E3*d2+
     4 2*nu12*E2*nu23*nu13*E3*d1-2*nu12*E2*nu23*nu13*E3*d1*d2-
     5 nu13**2*E2*E3+nu13**2*E2*E3*d1)**2
C
      END
C
C
C
