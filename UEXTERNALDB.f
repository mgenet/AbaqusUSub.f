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

!                                                           ------------
! ----------------------------------------------------------- INCLUDES -
!                                                           ------------

!                                                               --------
! --------------------------------------------------------------- UMAT -
!                                                               --------

SUBROUTINE UEXTERNALDB (LOP, LRESTART, TIME, DTIME, KSTEP, KINC)

! -------------------------------------------------------- VARIABLES ---

    INTEGER, INTENT(IN) :: LOP, LRESTART, KSTEP, KINC
    DOUBLE PRECISION, INTENT(IN) :: TIME(2), DTIME

    INTEGER N_ELEMS
    CHARACTER*80 ORIENT_FILE

! --------------------------------------------------- SET PARAMETERS ---

!     PARAMETER (N_ELEMS = 1)
!     PARAMETER (ORIENT_FILE="/home/genet/BOULOT/RESEARCH/Abaqus/ObjectivityTesting/orient-XYZ.inp")

    PARAMETER (N_ELEMS = 1)
    PARAMETER (ORIENT_FILE="/home/genet/BOULOT/RESEARCH/Abaqus/ObjectivityTesting/orient-YZX.inp")

!     PARAMETER (N_ELEMS = 1)
!     PARAMETER (ORIENT_FILE="/home/genet/BOULOT/RESEARCH/Abaqus/ObjectivityTesting/orient-ZXY.inp")

! ----------------------------------------------- VARIABLES (CONT'D) ---

    DOUBLE PRECISION ORIENT(N_ELEMS,9)
    COMMON /ORIENT/ ORIENT

! ------------------------------------------------- READ ORIENT FILE ---

    IF (LOP .EQ. 0) THEN
        OPEN (UNIT = 1, FILE = ORIENT_FILE)
        DO K_ELEM = 1,N_ELEMS
            READ (1,*) ORIENT(K_ELEM,1), ORIENT(K_ELEM,2), ORIENT(K_ELEM,3), ORIENT(K_ELEM,4), ORIENT(K_ELEM,5), ORIENT(K_ELEM,6), ORIENT(K_ELEM,7), ORIENT(K_ELEM,8), ORIENT(K_ELEM,9)
        ENDDO
!         PRINT *, "ORIENT = ", ORIENT
    ENDIF

! -------------------------------------------------------------- END ---

    RETURN
END SUBROUTINE UEXTERNALDB
