!! Klausmeier model Fortran File!! (adapted from cusp demo)

! Subroutine FUNC sets up the equations!

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

  DOUBLE PRECISION Z, X, W, RA1, RA2, RA3, IA1, IA2, IA3, L, nu, c, a, m, gamma, Rlam, Ilam, Rlamg, Ilamg, Rlamgg, Ilamgg, RAG1, RAG2, RAG3, IAG1, IAG2, IAG3, RAGG1, RAGG2, RAGG3, IAGG1, IAGG2, IAGG3, mu

  ! Define the variables we take derivates of !
  X = U(1)
  Z = U(2)
  W = U(3)

  RA1 = U(4)
  RA2 = U(5)
  RA3 = U(6)
  IA1 = U(7)
  IA2 = U(8)
  IA3 = U(9)

  !RAG1 = U(10)
  !RAG2 = U(11)
  !RAG3 = U(12)
  IAG1 = U(10)
  IAG2 = U(11)
  IAG3 = U(12)

  RAGG1 = U(13)
  RAGG2 = U(14)
  RAGG3 = U(15)
  !IAGG1 = U(19)
  !IAGG2 = U(20)
  !IAGG3 = U(21)


  ! Define the parameters in the model !
  nu = PAR(1)
  c = PAR(2)
  a = PAR(3)
  m = PAR(4)
  gamma = PAR(5)
  Rlam = PAR(6)
  Ilam = PAR(7)
  Rlamg = PAR(8)
  Ilamg = PAR(9)
  Rlamgg = PAR(10)
  L = PAR(11)
  Ilamgg = PAR(12)
  mu = PAR(13)



  ! Define the RHS of the equations we are modelling !
  F(1) = L*Z
  F(2) = L*(m*X - W*X*X - c*Z)
  F(3) = L*(W*X*X + W - a)/(nu + c)

  !alpha
  F(4) = IA1*gamma + L*RA2
  F(5) = -IA1*Ilam*L + IA2*gamma + L*RA1*Rlam - 2*L*RA1*W*X + L*RA1*m - L*RA2*c - L*RA3*X**2
  F(6) = -IA3*Ilam*L/(c + nu) + IA3*gamma + 2*L*RA1*W*X/(c + nu) + L*RA3*Rlam/(c + nu) + L*RA3*X**2/(c + nu) + L*RA3/(c + nu)
  F(7) = IA2*L - RA1*gamma
  F(8) = IA1*L*Rlam - 2*IA1*L*W*X + IA1*L*m - IA2*L*c - IA3*L*X**2 + Ilam*L*RA1 - RA2*gamma
  F(9) = 2*IA1*L*W*X/(c + nu) + IA3*L*Rlam/(c + nu) + IA3*L*X**2/(c + nu) + IA3*L/(c + nu) + Ilam*L*RA3/(c + nu) - RA3*gamma

  !alpha_gamma
  !F(10) = IA1 + IAG1*gamma + L*RAG2
  !F(11) = -IA1*Ilamg*L + IA2 - IAG1*Ilam*L + IAG2*gamma + L*RA1*Rlamg + L*RAG1*Rlam - 2*L*RAG1*W*X + L*RAG1*m - L*RAG2*c - L*RAG3*X**2
  !F(12) = -IA3*Ilamg*L/(c + nu) + IA3 - IAG3*Ilam*L/(c + nu) + IAG3*gamma + L*RA3*Rlamg/(c + nu) + 2*L*RAG1*W*X/(c + nu) + L*RAG3*Rlam/(c + nu) + L*RAG3*X**2/(c + nu) + L*RAG3/(c + nu)
  F(10) = IAG2*L - mu*RA1
  F(11) = -2*IAG1*L*W*X + IAG1*L*m - IAG2*L*c - IAG3*L*X**2 + Ilamg*L*RA1 - mu*RA2
  F(12) = 2*IAG1*L*W*X/(c + nu) + IAG3*L*X**2/(c + nu) + IAG3*L/(c + nu) + Ilamg*L*RA3/(c + nu) - mu*RA3

  !alpha_gammagamma
  F(13) = 2*mu*IAG1 + L*RAGG2
  F(14) = -2*mu*IAG1*Ilamg*L + 2*mu*IAG2 + L*RA1*Rlamgg + 2*L*RAG1*Rlamg - 2*L*RAGG1*W*X + L*RAGG1*m - L*RAGG2*c - L*RAGG3*X**2
  F(15) = -2*mu*IAG3*Ilamg*L/(c + nu) + 2*mu*IAG3 + L*RA3*Rlamgg/(c + nu) + 2*L*RAG3*Rlamg/(c + nu) + 2*L*RAGG1*W*X/(c + nu) + L*RAGG3*X**2/(c + nu) + L*RAGG3/(c + nu)
  !F(19) = IAGG2*L - 2*RAG1 - RAGG1*gamma
  !F(20) = IA1*L*Rlamgg + 2*IAG1*L*Rlamg + IAGG1*L*Rlam - 2*IAGG1*L*W*X + IAGG1*L*m - IAGG2*L*c - IAGG3*L*X**2 + Ilam*L*RAGG1 + 2*Ilamg*L*RAG1 + Ilamgg*L*RA1 - 2*RAG2 - RAGG2*gamma
  !F(21) = IA3*L*Rlamgg/(c + nu) + 2*IAG3*L*Rlamg/(c + nu) + 2*IAGG1*L*W*X/(c + nu) + IAGG3*L*Rlam/(c + nu) + IAGG3*L*X**2/(c + nu) + IAGG3*L/(c + nu) + Ilam*L*RAGG3/(c + nu) + 2*Ilamg*L*RAG3/(c + nu) + Ilamgg*L*RA3/(c + nu) - 2*RAG3 - RAGG3*gamma
  

END SUBROUTINE FUNC


! Subroutine STPNT sets up the starting point of the continuation algorithm !

SUBROUTINE STPNT(NDIM,U,PAR,T)
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION L
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T

  PAR(1) = 182.5d0   ! nu
  PAR(2) = 0.8d0       ! c
  PAR(3) = 2.5d0     ! a  
  PAR(4) = 0.45d0    ! m
  PAR(5) = 0.0d0     ! gamma
  PAR(6) = 0.0d0     ! real part of lambda
  PAR(7) = 0.0d0     ! imaginary part of lambda
  PAR(11)=201.58d0   ! wavelength / period

  U(1) = 5.369d0     !X initial steady state   
  U(2) = 0.0d0       !Z initial steady state
  U(3) = 0.0838d0    !W initial steady state

  !we write over these when intialising the continuation
  U(4) = 0.0d0       !RA1 inital value
  U(5) = 0.0d0       !RA2 inital value
  U(6) = 0.0d0       !RA3 inital value
  U(7) = 0.0d0       !IM1 inital value
  U(8) = 0.0d0       !IM2 intial value
  U(9) = 0.0d0       !IM3 intial value

  U(10) = 0.0d0       !IAG1 inital value
  U(11) = 0.0d0       !IAG2 inital value
  U(12) = 0.0d0       !IAG3 inital value
  U(13) = 0.0d0       !RAGG1 inital value
  U(14) = 0.0d0       !RAGG2 inital value
  U(15) = 0.0d0       !RAGG3 inital value




END SUBROUTINE STPNT


! Subroutine PLVS is what we want to track over time. !

SUBROUTINE PVLS(NDIM,U,PAR)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

  DOUBLE PRECISION GETP


END SUBROUTINE PVLS




SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
  DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
  INTEGER J


    do j=1,15
      fb(j) = U0(j) - U1(j)
    end do


END SUBROUTINE BCND

SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*), U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
  DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)
  INTEGER J

  fi(1) = -1.0
  fi(2) = 0.0
  fi(3) = 0.0
  fi(4) = 0.0
  fi(5) = 0.0

  !orthoganality and phase fixing for alpha
  do j=4,6
    fi(1) = fi(1) + U(j)*U(j) + U(j+3)*U(j+3)
    fi(2) = fi(2) + UOLD(j)*U(j+3) - UOLD(j+3)*U(j)
  end do

  !phase fixing for V bar bar
  do j=1,3
    fi(3) = fi(3) + UPOLD(j)*(UOLD(j) - U(j))
  end do

  !orthogonality for alpha_gamma
  do j=4,6
    fi(4) = fi(4) + U(j)*U(j+6)
  end do

  !orthogonality for alpha_gammagamma
  do j=4,6
    fi(5) = fi(5) + U(j)*U(j+9)
  end do



END SUBROUTINE ICND

SUBROUTINE FOPT 
END SUBROUTINE FOPT

