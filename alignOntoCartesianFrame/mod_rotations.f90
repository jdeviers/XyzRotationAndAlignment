MODULE mod_rotations
  USE mod_precision
  implicit none


  CONTAINS

  SUBROUTINE CC_ROTATION_X(alpha,NBATOMS,XYZ,XYZ_rotated)
! Right-hand rotation about x axis
    implicit none

    INTEGER(i4)          :: i,NBATOMS
    REAL(dp)             :: Rx(3,3),slice(3,1)
    REAL(dp)             :: alpha
    REAL(dp),ALLOCATABLE :: XYZ(:,:),XYZ_rotated(:,:)

    Rx = reshape((/1.0D0,     0.0D0,      0.0D0, &
                   0.0D0,COS(alpha),-SIN(alpha), &
                   0.0D0,SIN(alpha), COS(alpha)/), (/3,3/),ORDER=(/2,1/))

    IF (.NOT. ALLOCATED(XYZ_rotated)) ALLOCATE(XYZ_rotated(NBATOMS,3))

    DO i=1,NBATOMS
      slice = RESHAPE(XYZ(i,:),(/3,1/))
      XYZ_rotated(i,:) = RESHAPE(MATMUL(Rx,slice),(/3/))
    END DO

  END SUBROUTINE CC_ROTATION_X

! ------------

  SUBROUTINE CC_ROTATION_Y(beta,NBATOMS,XYZ,XYZ_rotated)
! Right-hand rotation about y axis
    implicit none

    INTEGER(i4)          :: i,NBATOMS
    REAL(dp)             :: Ry(3,3),slice(3,1)
    REAL(dp)             :: beta
    REAL(dp),ALLOCATABLE :: XYZ(:,:),XYZ_rotated(:,:)

    Ry = reshape((/ COS(beta),0.0D0,SIN(beta), &
                        0.0D0,1.0D0,    0.0D0, &
                   -SIN(beta),0.0D0,COS(beta)/), (/3,3/),ORDER=(/2,1/))

    IF (.NOT. ALLOCATED(XYZ_rotated)) ALLOCATE(XYZ_rotated(NBATOMS,3))

    DO i=1,NBATOMS
      slice = RESHAPE(XYZ(i,:),(/3,1/))
      XYZ_rotated(i,:) = RESHAPE(MATMUL(Ry,slice),(/3/))
    END DO

  END SUBROUTINE CC_ROTATION_Y

! ----------

  SUBROUTINE CC_ROTATION_Z(gamma,NBATOMS,XYZ,XYZ_rotated)
! Right-hand rotation about z axis
    implicit none

    INTEGER(i4) :: i,NBATOMS
    REAL(dp) :: Rz(3,3),slice(3,1)
    REAL(dp) :: gamma
    REAL(dp),ALLOCATABLE :: XYZ(:,:),XYZ_rotated(:,:)

    Rz = reshape((/COS(gamma),-SIN(gamma),0.0D0, &
                   SIN(gamma), COS(gamma),0.0D0, &
                        0.0D0,      0.0D0,1.0D0/), (/3,3/),ORDER=(/2,1/))

    IF (.NOT. ALLOCATED(XYZ_rotated)) ALLOCATE(XYZ_rotated(NBATOMS,3))

    DO i=1,NBATOMS
      slice = RESHAPE(XYZ(i,:),(/3,1/))
      XYZ_rotated(i,:) = RESHAPE(MATMUL(Rz,slice),(/3/))
    END DO

  END SUBROUTINE CC_ROTATION_Z

END MODULE mod_rotations