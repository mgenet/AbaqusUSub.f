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
C ELASTICITY_3D SUBROUTINE
C
      SUBROUTINE ELASTICITY_3D(K, NOEL, NTENS)
C
      CHARACTER*100 STEL, FSTEL, FILENAME
      INTEGER NOEL, NTENS, I, J
      DOUBLE PRECISION K(NTENS, NTENS)
C
!       PRINT *, NOEL
      WRITE(STEL, *) NOEL
      IF (NOEL < 10) THEN
        WRITE(FSTEL, *) '00000', STEL(12:12)
      ELSE IF (NOEL < 100) THEN
        WRITE(FSTEL, *) '0000', STEL(11:12)
      ELSE IF (NOEL < 1000) THEN
        WRITE(FSTEL, *) '000', STEL(10:12)
      ELSE IF (NOEL < 10000) THEN
        WRITE(FSTEL, *) '00', STEL(9:12)
      ELSE IF (NOEL < 100000) THEN
        WRITE(FSTEL, *) '0', STEL(8:12)
      ELSE IF (NOEL < 1000000) THEN
        WRITE(FSTEL, *) '', STEL(7:12)
      ENDIF
!       PRINT *, FSTEL
!       WRITE(FILENAME, *) '/home/genet/SCIENCE/UCSB/CALCULS/',
!       1 'bm_heterogeneous/ELEMENTS_EFF_MED_PROPERTIES/', FSTEL(2:7)
      WRITE(FILENAME, *) '/u/binary_model/CALCULS/',
     1 'binary_model/heterogeneous/ELEMENTS_EFF_MED_PROPERTIES/',
     2 FSTEL(2:7)
      FILENAME = FILENAME(2:)
!       PRINT *, FILENAME
C
      OPEN(UNIT=1, FILE=FILENAME)
C
      DO I = 1, NTENS
        DO J = 1, I
          READ(1, *) K(I, J)
          K(J, I) = K(I, J)
!           PRINT *, 'K(', I, ' ,', J, ' ) =', K(I, J)
        ENDDO
      ENDDO
C
      CLOSE(UNIT=1)
C
!       CALL DPOTRF('U', 6, K, 6, INFO)
!       CALL DPOTRI('U', 6, K, 6, INFO)
C
!       DO I = 1, 6
!         DO J = 1, I
!           K(I, J) = K(J, I)
! !           PRINT *, 'K(', I, ' ,', J, ' ) =', K(J, I)
!         ENDDO
!       ENDDO
C
      END
C
C
C
