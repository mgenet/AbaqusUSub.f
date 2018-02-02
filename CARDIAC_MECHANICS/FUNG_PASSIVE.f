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

#ifndef FUNG_PASSIVE_f
#define FUNG_PASSIVE_f

!                                                           ------------
! ----------------------------------------------------------- INCLUDES -
!                                                           ------------

#include "UTILS/VEC_MAT_TOOLS.f"

!                                                       ----------------
! ------------------------------------------------------- FUNG_PASSIVE -
!                                                       ----------------

SUBROUTINE FUNG_PASSIVE(C0, B0, E, U, DU_DE, DDU_DDE, TOO_LARGE_DEFORMATION)
    IMPLICIT NONE

! -------------------------------------------------------- VARIABLES ---

    DOUBLE PRECISION, INTENT(IN) :: C0, B0(6), E(6)

    LOGICAL, INTENT(OUT) :: TOO_LARGE_DEFORMATION
    DOUBLE PRECISION, INTENT(OUT) :: U, DU_DE(6), DDU_DDE(21)

    INTEGER          INDX
    DOUBLE PRECISION EXP_TERM

! ----------------------------------------------------------------------

!   EXP TERM
    EXP_TERM =  B0(1)*E(1)**2 &
             +  B0(2)*E(2)**2 &
             +  B0(3)*E(3)**2 &
             +4*B0(4)*E(4)**2 & ! 20150122: OK, Il n'y a pas les deux dans EBAR
             +4*B0(5)*E(5)**2 & ! 20150122: OK, Il n'y a pas les deux dans EBAR
             +4*B0(6)*E(6)**2   ! 20150122: OK, Il n'y a pas les deux dans EBAR
!     PRINT *, "EXP_TERM = ", EXP_TERM
    IF (EXP_TERM > 700.) THEN
        PRINT *, "TOO LARGE DEFORMATION…"
        PRINT *, "E = ", E
        PRINT *, "EXP_TERM = ", EXP_TERM
        EXP_TERM = 700.
        TOO_LARGE_DEFORMATION = .TRUE.
    ELSE
        TOO_LARGE_DEFORMATION = .FALSE.
    ENDIF
    EXP_TERM = EXP(EXP_TERM)
!     PRINT *, "EXP_TERM = ", EXP_TERM

!   ENERGY DENSITY
    U = (C0/2) * (EXP_TERM - 1.)
!     PRINT *, "U = ", U

!   FIRST DERIVATIVE
    DU_DE = 0.
    DU_DE(1) = C0 * EXP_TERM *   B0(1)*E(1)
    DU_DE(2) = C0 * EXP_TERM *   B0(2)*E(2)
    DU_DE(3) = C0 * EXP_TERM *   B0(3)*E(3)
    DU_DE(4) = C0 * EXP_TERM * 2*B0(4)*E(4) ! 20150122: Comprends pas, ça devrait être 4 là non? Mais si je mets 4 alors je n'obtiens pas le même résultat qu'avec la loi de Fung d'Abaqus...
    DU_DE(5) = C0 * EXP_TERM * 2*B0(5)*E(5) ! 20150612: Non, en fait c'est bien un 2 car en vrai l'énergie contient les termes E12**2, E21**2, E12*E21 & E21*E12
    DU_DE(6) = C0 * EXP_TERM * 2*B0(6)*E(6)

!   SECOND DERIVATIVE
    DDU_DDE = 0.
    DDU_DDE(INDX(1,1)) = C0 * EXP_TERM * (B0(1)+2*B0(1)*E(1)*B0(1)*E(1))
    DDU_DDE(INDX(1,2)) = C0 * EXP_TERM * (      2*B0(1)*E(1)*B0(2)*E(2))
    DDU_DDE(INDX(1,3)) = C0 * EXP_TERM * (      2*B0(1)*E(1)*B0(3)*E(3))
    DDU_DDE(INDX(1,4)) = C0 * EXP_TERM * (      4*B0(1)*E(1)*B0(4)*E(4))
    DDU_DDE(INDX(1,5)) = C0 * EXP_TERM * (      4*B0(1)*E(1)*B0(5)*E(5))
    DDU_DDE(INDX(1,6)) = C0 * EXP_TERM * (      4*B0(1)*E(1)*B0(6)*E(6))
    DDU_DDE(INDX(2,2)) = C0 * EXP_TERM * (B0(2)+2*B0(2)*E(2)*B0(2)*E(2))
    DDU_DDE(INDX(2,3)) = C0 * EXP_TERM * (      2*B0(2)*E(2)*B0(3)*E(3))
    DDU_DDE(INDX(2,4)) = C0 * EXP_TERM * (      4*B0(2)*E(2)*B0(4)*E(4))
    DDU_DDE(INDX(2,5)) = C0 * EXP_TERM * (      4*B0(2)*E(2)*B0(5)*E(5))
    DDU_DDE(INDX(2,6)) = C0 * EXP_TERM * (      4*B0(2)*E(2)*B0(6)*E(6))
    DDU_DDE(INDX(3,3)) = C0 * EXP_TERM * (B0(3)+2*B0(3)*E(3)*B0(3)*E(3))
    DDU_DDE(INDX(3,4)) = C0 * EXP_TERM * (      4*B0(3)*E(3)*B0(4)*E(4))
    DDU_DDE(INDX(3,5)) = C0 * EXP_TERM * (      4*B0(3)*E(3)*B0(5)*E(5))
    DDU_DDE(INDX(3,6)) = C0 * EXP_TERM * (      4*B0(3)*E(3)*B0(6)*E(6))
    DDU_DDE(INDX(4,4)) = C0 * EXP_TERM * (B0(4)+8*B0(4)*E(4)*B0(4)*E(4))
    DDU_DDE(INDX(4,5)) = C0 * EXP_TERM * (      8*B0(4)*E(4)*B0(5)*E(5))
    DDU_DDE(INDX(4,6)) = C0 * EXP_TERM * (      8*B0(4)*E(4)*B0(6)*E(6))
    DDU_DDE(INDX(5,5)) = C0 * EXP_TERM * (B0(5)+8*B0(5)*E(5)*B0(5)*E(5))
    DDU_DDE(INDX(5,6)) = C0 * EXP_TERM * (      8*B0(5)*E(5)*B0(6)*E(6))
    DDU_DDE(INDX(6,6)) = C0 * EXP_TERM * (B0(6)+8*B0(6)*E(6)*B0(6)*E(6))

    RETURN
END SUBROUTINE FUNG_PASSIVE

#endif
