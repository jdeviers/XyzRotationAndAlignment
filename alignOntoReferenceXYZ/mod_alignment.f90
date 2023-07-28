MODULE mod_alignment
  USE mod_precision
  USE mod_procedures
  USE mod_transrot
  implicit none

  CONTAINS

! ----------

  SUBROUTINE ALIGN(NBATOMS,XYZ_REF,XYZ_TOROTATE,XYZ_ROTATED,ROTATIONS,ROTMAT)
!
! -- Loops over the 3 axes
!
    implicit none

    INTEGER               :: i,NBATOMS,pass
    INTEGER               :: ALLOCSTAT
    REAL(dp)              :: RMSD
    REAL(dp)              :: ROTMAT(3,3)      
    REAL(dp),ALLOCATABLE  :: ROTATIONS(:,:)
    REAL(dp),ALLOCATABLE  :: XYZ_REF(:,:),XYZ_TOROTATE(:,:),XYZ_ROTATED(:,:)

! -- Allocations
    ALLOCATE(ROTATIONS(10,3),STAT=ALLOCSTAT)
    IF (ALLOCSTAT .NE. 0) STOP '**** ALLOC. FAILED: NOT ENOUGH MEMORY ****'
    ALLOCATE(XYZ_ROTATED(NBATOMS,3),STAT=ALLOCSTAT)
    IF (ALLOCSTAT .NE. 0) STOP '**** ALLOC. FAILED: NOT ENOUGH MEMORY ****'

! -- Initial assignments
    ROTATIONS = 0.D0
    XYZ_ROTATED = XYZ_TOROTATE
    pass=0

! -- Perform the 3 rotations
    DO WHILE (RMSD .GT. RMSD_thresh)
      DO i=1,3
        CALL OPT_ROT(NBATOMS,i,XYZ_REF,XYZ_ROTATED,ROTATIONS,RMSD,pass)
      END DO
      pass = pass+1
      IF (pass .EQ. 10) EXIT
    END DO

! -- Build the ROTMAT from the 3 rotations
    CALL BUILD_ROTMAT(ROTATIONS,ROTMAT)

  END SUBROUTINE ALIGN

! ----------

  SUBROUTINE OPT_ROT(NBATOMS,axis,XYZ_REF,XYZ_TOROTATE,ROTATIONS,RMSD,pass)
!
! -- Performs the rotation, locates the best angle
!
    implicit none

    INTEGER               :: NBATOMS,axis,i,pass
    INTEGER               :: ALLOCSTAT
    REAL(dp)              :: RMSD
    REAL(dp)              :: PW_RMSD(1000)
    REAL(dp)              :: ROT(3,3)
    REAL(dp),PARAMETER    :: stepwise_angle = (2.D0*pi)/1000.D0
    REAL(dp),ALLOCATABLE  :: ROTATIONS(:,:)
    REAL(dp),ALLOCATABLE  :: XYZ_REF(:,:),XYZ_TOROTATE(:,:),WORKING_XYZ(:,:)

! -- Allocations
    ALLOCATE(WORKING_XYZ(NBATOMS,3),STAT=ALLOCSTAT)
    IF (ALLOCSTAT .NE. 0) STOP '**** ALLOC. FAILED: NOT ENOUGH MEMORY ****'

! -- Initial assignments
    WORKING_XYZ = XYZ_TOROTATE
    PW_RMSD(:) = 1000.D0
    
    IF (axis .EQ. 1) ROT = Rx(stepwise_angle)
    IF (axis .EQ. 2) ROT = Ry(stepwise_angle)
    IF (axis .EQ. 3) ROT = Rz(stepwise_angle)

    DO i=1,1000
      CALL DO_ROTATION(NBATOMS,ROT,WORKING_XYZ,WORKING_XYZ)
      CALL PAIRWISE_RMSD(NBATOMS,XYZ_REF,WORKING_XYZ,PW_RMSD(i))
    END DO

    ROTATIONS(pass,axis) = MINLOC(PW_RMSD,1) * stepwise_angle
    WRITE(*,'(A,F6.3,A,I1,A,F6.3,A)') ' * Rotating by ',ROTATIONS(pass,axis),' rad. along axis ',axis, &
                                      ' for P_RMSD = ',MINVAL(PW_RMSD),' Angstrom.'

    RMSD = MINVAL(PW_RMSD)

    IF (axis .EQ. 1) ROT = Rx(ROTATIONS(pass,axis))
    IF (axis .EQ. 2) ROT = Ry(ROTATIONS(pass,axis))
    IF (axis .EQ. 3) ROT = Rz(ROTATIONS(pass,axis))

    CALL DO_ROTATION(NBATOMS,ROT,XYZ_TOROTATE,XYZ_TOROTATE)

  END SUBROUTINE OPT_ROT

! ----------

  SUBROUTINE DO_ROTATION(NBATOMS,MAT,XYZ_IN,XYZ_OUT)
    implicit none

    INTEGER               :: NBATOMS,i
    INTEGER               :: ALLOCSTAT
    REAL(dp)              :: slice(3)
    REAL(dp)              :: MAT(3,3)
    REAL(dp),ALLOCATABLE  :: XYZ_IN(:,:),XYZ_OUT(:,:)

! -- Allocations
    ALLOCSTAT = 0
    IF (.NOT. ALLOCATED(XYZ_OUT)) ALLOCATE(XYZ_OUT(NBATOMS,3),STAT=ALLOCSTAT)
    IF (ALLOCSTAT .NE. 0) STOP '**** ALLOC. FAILED: NOT ENOUGH MEMORY ****'

    DO i=1,NBATOMS
      slice = RESHAPE(XYZ_IN(i,:),(/3/))                ! prepare slice in column vector for matmul
      XYZ_OUT(i,:) = RESHAPE(MATMUL(MAT,slice),(/3/))
    END DO

  END SUBROUTINE DO_ROTATION

  SUBROUTINE PAIRWISE_RMSD(NBATOMS,XYZ_1,XYZ_2,PW_RMSD)
!
! -- Contains scoring function for rotation angle
!
    implicit  none

    INTEGER               :: NBATOMS,i
    REAL(dp)              :: PW_RMSD
    REAL(dp),ALLOCATABLE  :: XYZ_1(:,:),XYZ_2(:,:)

! -- Initial assignments
    PW_RMSD = 0.D0

    DO i = 1,NBATOMS
      PW_RMSD = PW_RMSD + NORM2(XYZ_1(i,:) - XYZ_2(i,:))
    END DO

    PW_RMSD = PW_RMSD/NBATOMS

  END SUBROUTINE PAIRWISE_RMSD

! ----------

  SUBROUTINE BUILD_ROTMAT(ROTATIONS,ROTMAT)
! -- Builds the rotation matrix, applies it onto the geom to test it is correct.
    implicit none

    INTEGER              :: i
    REAL(dp)             :: alpha,beta,gamma
    REAL(dp)             :: Rx(3,3),Ry(3,3),Rz(3,3),ID(3,3)
    REAL(dp)             :: ROTMAT(3,3)
    REAL(dp),ALLOCATABLE :: ROTATIONS(:,:)

    alpha=0.d0 ; beta=0.d0 ; gamma = 0.d0
    DO i=1,UBOUND(ROTATIONS,1)
      alpha = alpha + ROTATIONS(i,1)
      beta  = beta  + ROTATIONS(i,2)
      gamma = gamma + ROTATIONS(i,3)
    END DO

    ID = 0.D0
    DO i=1,3
      ID(i,i) = 1.D0
    END DO

    Rx = reshape((/1.0D0,     0.0D0,      0.0D0, &
                   0.0D0,COS(alpha),-SIN(alpha), &
                   0.0D0,SIN(alpha), COS(alpha)/), (/3,3/),ORDER=(/2,1/))

    Ry = reshape((/ COS(beta),0.0D0,SIN(beta), &
                        0.0D0,1.0D0,    0.0D0, &
                   -SIN(beta),0.0D0,COS(beta)/), (/3,3/),ORDER=(/2,1/))

    Rz = reshape((/COS(gamma),-SIN(gamma),0.0D0, &
                   SIN(gamma), COS(gamma),0.0D0, &
                        0.0D0,      0.0D0,1.0D0/), (/3,3/),ORDER=(/2,1/))


    ROTMAT = MATMUL(Rx,ID)
    ROTMAT = MATMUL(Ry,ROTMAT)
    ROTMAT = MATMUL(Rz,ROTMAT)

  END SUBROUTINE

END MODULE mod_alignment
