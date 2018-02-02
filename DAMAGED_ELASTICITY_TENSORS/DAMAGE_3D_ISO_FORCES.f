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
C DAMAGE_3D_ISO_FORCES SUBROUTINE
C
      SUBROUTINE DAMAGE_3D_ISO_FORCES(Y, EPS, E, nu, d)
C
      DOUBLE PRECISION Y, EPS(6), E, nu, d
C
      Y = E*(nu*EPS(1)**2-2*EPS(2)*EPS(3)*nu+2*EPS(6)**2*nu+nu*EPS(2)**2+
     1 2*EPS(5)**2*nu+nu*EPS(3)**2-2*EPS(1)*EPS(2)*nu+2*EPS(4)**2*nu-
     2 2*EPS(1)*EPS(3)*nu-EPS(6)**2-EPS(1)**2-EPS(5)**2-EPS(2)**2-
     3 EPS(3)**2-EPS(4)**2)/((1+nu)*(2*nu-1))/2
C
      END
C
C
C
