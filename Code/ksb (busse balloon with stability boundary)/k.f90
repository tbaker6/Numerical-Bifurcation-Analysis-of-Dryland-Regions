!! Klausmeier model Fortran File!! (adapted from cusp demo)

! Subroutine FUNC sets up the equations!

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

  DOUBLE PRECISION Z, X, W, nu, c, a, m

  ! Define the variables we take derivates of !
  X = U(1)
  Z = U(2)
  W = U(3)

  ! Define the parameters in the model !
  nu = PAR(1)
  c = PAR(2)
  a = PAR(3)
  m = PAR(4)

  ! Define the RHS of the equations we are modelling !
  F(1)= Z
  F(2) = m*X - W*X*X - c*Z
  F(3) = (W*X*X + W - a)/(nu + c)
  

END SUBROUTINE FUNC


! Subroutine STPNT sets up the starting point of the continuation algorithm !

SUBROUTINE STPNT(NDIM,U,PAR,T)
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T

  PAR(1) = 182.5d0   ! nu
  PAR(2) = 0.8d0     ! c
  PAR(3) = 2.5d0     ! a  -- start near, but below, the reported Hopf (~2.78)
  PAR(4) = 0.45d0    ! m

  U(1) = 5.369d0    
  U(2) = 0.0d0       
  U(3) = 0.0838d0   

END SUBROUTINE STPNT


! Subroutine PLVS is what we want to track over time. !

SUBROUTINE PVLS(NDIM,U,PAR)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

  DOUBLE PRECISION GETP

   !PAR(5) = PAR(2)

END SUBROUTINE PVLS


SUBROUTINE BCND 
END SUBROUTINE BCND

SUBROUTINE ICND 
END SUBROUTINE ICND

SUBROUTINE FOPT 
END SUBROUTINE FOPT

