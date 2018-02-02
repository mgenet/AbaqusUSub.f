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

#include "VUMAT_3D_FUNG.f"

!                                                               --------
! --------------------------------------------------------------- UMAT -
!                                                               --------

SUBROUTINE VUMAT (NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL, STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH, PROPS, DENSITY, STRAININC, RELSPININC, TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD, STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD, TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW, STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)

!     INCLUDE 'VABA_PARAM.INC'

    CHARACTER*80 CMNAME
    INTEGER      NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL
    REAL         STEPTIME, TOTALTIME, DT, COORDMP(NBLOCK,*), CHARLENGTH(NBLOCK), PROPS(NPROPS), DENSITY(NBLOCK), STRAININC(NBLOCK,NDIR+NSHR), RELSPININC(NBLOCK,NSHR), TEMPOLD(NBLOCK), STRETCHOLD(NBLOCK,NDIR+NSHR), DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR), FIELDOLD(NBLOCK,NFIELDV), STRESSOLD(NBLOCK,NDIR+NSHR), STATEOLD(NBLOCK,NSTATEV), ENERINTERNOLD(NBLOCK), ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK), STRETCHNEW(NBLOCK,NDIR+NSHR), DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR), FIELDNEW(NBLOCK,NFIELDV), STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV), ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)

    IF (CMNAME .EQ. 'VUMAT_F_3D_FUNG') THEN
        CALL VUMAT_3D_FUNG (NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL, STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH, PROPS, DENSITY, STRAININC, RELSPININC, TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD, STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD, TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW, STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
    ENDIF

    RETURN
END SUBROUTINE VUMAT



