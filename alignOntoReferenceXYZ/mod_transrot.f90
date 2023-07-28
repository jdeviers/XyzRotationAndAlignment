MODULE mod_transrot
  USE mod_precision
  implicit none

  CONTAINS

  SUBROUTINE TRANSLATE_N5_TO_ORIGIN(NBATOMS,XYZ,XYZ_translated)
    implicit none

    INTEGER              :: NBATOMS,ALLOCSTAT
    REAL(dp),ALLOCATABLE :: XYZ(:,:),XYZ_translated(:,:)

    IF (.NOT. ALLOCATED(XYZ_translated)) ALLOCATE(XYZ_translated(NBATOMS,3),stat=ALLOCSTAT)
    IF (ALLOCSTAT .NE. 0) STOP '**** ALLOC. FAILED: NOT ENOUGH MEMORY ****'

    WRITE(*,'(A,I3,A)') '****  Using N5 index = ',N5_index,' ****'

    IF (XYZ(N5_index,1) .NE. 0.D0) XYZ_translated(:,1) = XYZ(:,1) - XYZ(N5_index,1)
    IF (XYZ(N5_index,2) .NE. 0.D0) XYZ_translated(:,2) = XYZ(:,2) - XYZ(N5_index,2)
    IF (XYZ(N5_index,3) .NE. 0.D0) XYZ_translated(:,3) = XYZ(:,3) - XYZ(N5_index,3)

  END SUBROUTINE TRANSLATE_N5_TO_ORIGIN

  FUNCTION Rx(alpha)
    implicit none

    REAL(dp) :: alpha
    REAL(dp) :: Rx(3,3)

    Rx = reshape((/1.0D0,     0.0D0,      0.0D0, &
                   0.0D0,COS(alpha),-SIN(alpha), &
                   0.0D0,SIN(alpha), COS(alpha)/), (/3,3/),ORDER=(/2,1/))
  END FUNCTION Rx

  FUNCTION Ry(beta)
    implicit none

    REAL(dp) :: beta
    REAL(dp) :: Ry(3,3)

    Ry = reshape((/ COS(beta),0.0D0,SIN(beta), &
                        0.0D0,1.0D0,    0.0D0, &
                   -SIN(beta),0.0D0,COS(beta)/), (/3,3/),ORDER=(/2,1/))
  END FUNCTION Ry

  FUNCTION Rz(gamma)
    implicit none

    REAL(dp) :: gamma
    REAL(dp) :: Rz(3,3)

    Rz = reshape((/COS(gamma),-SIN(gamma),0.0D0, &
                   SIN(gamma), COS(gamma),0.0D0, &
                        0.0D0,      0.0D0,1.0D0/), (/3,3/),ORDER=(/2,1/))
  END FUNCTION Rz

END MODULE mod_transrot