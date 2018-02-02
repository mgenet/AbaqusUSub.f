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

#ifndef VEC_MAT_TOOLS_f
#define VEC_MAT_TOOLS_f

!                                                          -------------
! ---------------------------------------------------------- HEAVISIDE -
!                                                          -------------

DOUBLE PRECISION FUNCTION HEAVISIDE(X)
    IMPLICIT NONE
    DOUBLE PRECISION X
    IF (X .GE. 0.) THEN
        HEAVISIDE = 1.
    ELSE
        HEAVISIDE = 0.
    ENDIF
    RETURN
END FUNCTION HEAVISIDE

!                                                           ------------
! ----------------------------------------------------------- POS_PART -
!                                                           ------------

DOUBLE PRECISION FUNCTION POS_PART(X)
    IMPLICIT NONE
    DOUBLE PRECISION X
    IF (X .GE. 0) THEN
        POS_PART = X
    ELSE
        POS_PART = 0.
    ENDIF
    RETURN
END FUNCTION POS_PART

!                                                           ------------
! ----------------------------------------------------------- NEG_PART -
!                                                           ------------

DOUBLE PRECISION FUNCTION NEG_PART(X)
    IMPLICIT NONE
    DOUBLE PRECISION X
    IF (X .LE. 0) THEN
        NEG_PART = X
    ELSE
        NEG_PART = 0.
    ENDIF
    RETURN
END FUNCTION NEG_PART

!                                                              ---------
! -------------------------------------------------------------- DET22 -
!                                                              ---------

DOUBLE PRECISION FUNCTION DET22(M)
    IMPLICIT NONE
    DOUBLE PRECISION M(2,2)
    DET22 = M(1,1)*M(2,2)-M(2,1)*M(1,2)
    RETURN
END FUNCTION DET22

!                                                              ---------
! -------------------------------------------------------------- DET33 -
!                                                              ---------

DOUBLE PRECISION FUNCTION DET33(M)
    IMPLICIT NONE
    DOUBLE PRECISION M(3,3)
    DET33 = M(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2)) &
          - M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1)) &
          + M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
    RETURN
END FUNCTION DET33

!                                                              ---------
! -------------------------------------------------------------- COF33 -
!                                                              ---------

SUBROUTINE COF33(M, COFM)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: M(3,3)
    DOUBLE PRECISION, INTENT(OUT) :: COFM(3,3)
    COFM(1,1) =        M(2,2)*M(3,3) - M(2,3)*M(3,2)
    COFM(1,2) = (-1.)*(M(2,1)*M(3,3) - M(3,1)*M(2,3))
    COFM(1,3) =        M(2,1)*M(3,2) - M(2,2)*M(3,1)

    COFM(2,1) = (-1.)*(M(1,2)*M(3,3) - M(1,3)*M(3,2))
    COFM(2,2) =        M(1,1)*M(3,3) - M(1,3)*M(3,1)
    COFM(2,3) = (-1.)*(M(1,1)*M(3,2) - M(1,2)*M(3,1))

    COFM(3,1) =        M(1,2)*M(2,3) - M(2,2)*M(1,3)
    COFM(3,2) = (-1.)*(M(1,1)*M(2,3) - M(2,1)*M(1,3))
    COFM(3,3) =        M(1,1)*M(2,2) - M(1,2)*M(2,1)
    RETURN
END SUBROUTINE COF33

!                                                              ---------
! -------------------------------------------------------------- ADJ33 -
!                                                              ---------

SUBROUTINE ADJ33(M, ADJM)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: M(3,3)
    DOUBLE PRECISION, INTENT(OUT) :: ADJM(3,3)
    CALL COF33(M, ADJM)
    ADJM = TRANSPOSE(ADJM)
    RETURN
END SUBROUTINE ADJ33

!                                                              ---------
! -------------------------------------------------------------- INV33 -
!                                                              ---------

SUBROUTINE INV33(M, INVM)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: M(3,3)
    DOUBLE PRECISION, INTENT(OUT) :: INVM(3,3)
    DOUBLE PRECISION DET33
    CALL ADJ33(M, INVM)
    INVM = INVM / DET33(M)
    RETURN
END SUBROUTINE INV33

!                                                            -----------
! ------------------------------------------------------------ SETID22 -
!                                                            -----------

SUBROUTINE SETID22(M)
    IMPLICIT NONE
    DOUBLE PRECISION M(2,2)
    M      = 0.
    M(1,1) = 1.
    M(2,2) = 1.
    RETURN
END SUBROUTINE SETID22

!                                                            -----------
! ------------------------------------------------------------ SETID33 -
!                                                            -----------

SUBROUTINE SETID33(M)
    IMPLICIT NONE
    DOUBLE PRECISION M(3,3)
    M      = 0.
    M(1,1) = 1.
    M(2,2) = 1.
    M(3,3) = 1.
    RETURN
END SUBROUTINE SETID33

!                                              -------------------------
! ---------------------------------------------- MATSYM33_FROM_VECCOL6 -
!                                              -------------------------

SUBROUTINE MATSYM33_FROM_VECCOL6(MATSYM33, VECCOL6)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: VECCOL6(6)
    DOUBLE PRECISION, INTENT(OUT) :: MATSYM33(3,3)
    MATSYM33(1,1) = VECCOL6(1)
    MATSYM33(1,2) = VECCOL6(4)
    MATSYM33(1,3) = VECCOL6(5)
    MATSYM33(2,1) = VECCOL6(4)
    MATSYM33(2,2) = VECCOL6(2)
    MATSYM33(2,3) = VECCOL6(6)
    MATSYM33(3,1) = VECCOL6(5)
    MATSYM33(3,2) = VECCOL6(6)
    MATSYM33(3,3) = VECCOL6(3)
!     PRINT *, VECCOL6
!     PRINT *, MATSYM33
    RETURN
END SUBROUTINE MATSYM33_FROM_VECCOL6

SUBROUTINE VECCOL6_TO_MATSYM33(VECCOL6, MATSYM33)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: VECCOL6(6)
    DOUBLE PRECISION, INTENT(OUT) :: MATSYM33(3,3)
    CALL MATSYM33_FROM_VECCOL6(MATSYM33, VECCOL6)
END SUBROUTINE VECCOL6_TO_MATSYM33

!                                              -------------------------
! ---------------------------------------------- VECCOL6_FROM_MATSYM33 -
!                                              -------------------------

SUBROUTINE VECCOL6_FROM_MATSYM33(VECCOL6, MATSYM33)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: MATSYM33(3,3)
    DOUBLE PRECISION, INTENT(OUT) :: VECCOL6(6)
    VECCOL6(1) = MATSYM33(1,1)
    VECCOL6(2) = MATSYM33(2,2)
    VECCOL6(3) = MATSYM33(3,3)
    VECCOL6(4) = MATSYM33(1,2)
    VECCOL6(5) = MATSYM33(1,3)
    VECCOL6(6) = MATSYM33(2,3)
    RETURN
END SUBROUTINE VECCOL6_FROM_MATSYM33

SUBROUTINE MATSYM33_TO_VECCOL6(MATSYM33, VECCOL6)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: MATSYM33(3,3)
    DOUBLE PRECISION, INTENT(OUT) :: VECCOL6(6)
    CALL VECCOL6_FROM_MATSYM33(VECCOL6, MATSYM33)
END SUBROUTINE MATSYM33_TO_VECCOL6

!                                                             ----------
! ------------------------------------------------------------- GET_IJ -
!                                                             ----------

INTEGER FUNCTION GET_IJ(NDIM, I, J)
    IMPLICIT NONE
    INTEGER NDIM, I, J
!     PRINT *, "NDIM = ", NDIM
!     PRINT *, "I = ", I
!     PRINT *, "J = ", J
    IF (NDIM.EQ.2) THEN
        IF ((I.EQ.1).AND.(J.EQ.1)) THEN
            GET_IJ = 1
        ELSE IF ((I.EQ.2).AND.(J.EQ.2)) THEN
            GET_IJ = 2
        ELSE IF ((I.EQ.1).AND.(J.EQ.2)) THEN
            GET_IJ =3
        ELSE IF ((I.EQ.2).AND.(J.EQ.1)) THEN
            GET_IJ = 3
        ELSE
            PRINT *, "I and J must be between 1 and 2. Aborting."
            STOP
        END IF
    ELSE IF (NDIM.EQ.3) THEN
        IF ((I.EQ.1).AND.(J.EQ.1)) THEN
            GET_IJ = 1
        ELSE IF ((I.EQ.2).AND.(J.EQ.2)) THEN
            GET_IJ = 2
        ELSE IF ((I.EQ.3).AND.(J.EQ.3)) THEN
            GET_IJ = 3
        ELSE IF ((I.EQ.1).AND.(J.EQ.2)) THEN
            GET_IJ = 4
        ELSE IF ((I.EQ.2).AND.(J.EQ.1)) THEN
            GET_IJ = 4
        ELSE IF ((I.EQ.1).AND.(J.EQ.3)) THEN
            GET_IJ = 5
        ELSE IF ((I.EQ.3).AND.(J.EQ.1)) THEN
            GET_IJ = 5
        ELSE IF ((I.EQ.2).AND.(J.EQ.3)) THEN
            GET_IJ = 6
        ELSE IF ((I.EQ.3).AND.(J.EQ.2)) THEN
            GET_IJ = 6
        ELSE
            PRINT *, "I and J must be between 1 and 3. Aborting."
            STOP
        END IF
    ELSE
        PRINT *, "NDIM must be 2 or 3. Aborting."
        STOP
    END IF
!     PRINT *, "GET_IJ = ", GET_IJ
    RETURN
END FUNCTION GET_IJ

!                                                              ---------
! -------------------------------------------------------------- GET_I -
!                                                              ---------

INTEGER FUNCTION GET_I(NDIM, IJ)
    IMPLICIT NONE
    INTEGER NDIM, IJ
    IF (NDIM.EQ.2) THEN
        IF (IJ.EQ.1) THEN
            GET_I = 1
        ELSE IF (IJ.EQ.2) THEN
            GET_I = 2
        ELSE IF (IJ.EQ.3) THEN
            GET_I = 1
        ELSE
            PRINT *, "IJ must be between 1 and 3. Aborting."
            STOP
        END IF
    ELSE IF (NDIM.EQ.3) THEN
        IF (IJ.EQ.1) THEN
            GET_I = 1
        ELSE IF (IJ.EQ.2) THEN
            GET_I = 2
        ELSE IF (IJ.EQ.3) THEN
            GET_I = 3
        ELSE IF (IJ.EQ.4) THEN
            GET_I = 1
        ELSE IF (IJ.EQ.5) THEN
            GET_I = 1
        ELSE IF (IJ.EQ.6) THEN
            GET_I = 2
        ELSE
            PRINT *, "IJ must be between 1 and 6. Aborting."
            STOP
        END IF
    ELSE
        PRINT *, "NDIM must be 2 or 3. Aborting."
        STOP
    END IF
    RETURN
END FUNCTION GET_I

!                                                              ---------
! -------------------------------------------------------------- GET_J -
!                                                              ---------

INTEGER FUNCTION GET_J(NDIM, IJ)
    IMPLICIT NONE
    INTEGER NDIM, IJ
    IF (NDIM.EQ.2) THEN
        IF (IJ.EQ.1) THEN
            GET_J = 1
        ELSE IF (IJ.EQ.2) THEN
            GET_J = 2
        ELSE IF (IJ.EQ.3) THEN
            GET_J = 2
        ELSE
            PRINT *, "IJ must be between 1 and 3. Aborting."
            STOP
        END IF
    ELSE IF (NDIM.EQ.3) THEN
        IF (IJ.EQ.1) THEN
            GET_J = 1
        ELSE IF (IJ.EQ.2) THEN
            GET_J = 2
        ELSE IF (IJ.EQ.3) THEN
            GET_J = 3
        ELSE IF (IJ.EQ.4) THEN
            GET_J = 2
        ELSE IF (IJ.EQ.5) THEN
            GET_J = 3
        ELSE IF (IJ.EQ.6) THEN
            GET_J = 3
        ELSE
            PRINT *, "IJ must be between 1 and 6. Aborting."
            STOP
        END IF
    ELSE
        PRINT *, "NDIM must be 2 or 3. Aborting."
        STOP
    END IF
    RETURN
END FUNCTION GET_J

!                                                               --------
! --------------------------------------------------------------- INDX -
!                                                               --------

INTEGER FUNCTION INDX(I, J)
    IMPLICIT NONE
    INTEGER I, J, II, JJ
    II = MIN(I,J)
    JJ = MAX(I,J)
    INDX = II + JJ*(JJ-1)/2
    RETURN
END FUNCTION INDX

#endif