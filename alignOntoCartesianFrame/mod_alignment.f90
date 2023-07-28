MODULE mod_alignment
  USE mod_precision
  USE mod_procedures
  USE mod_rotations
  implicit none

  CONTAINS

  SUBROUTINE ALIGN(NBATOMS,initial_XYZ,final_XYZ,rotations,ROTMAT)
! -- N10 -> N5 must be along x (towards positive)
! -- C1 -> C2 must be along y (towards positive)
! -- If these conditions are satisfied, then the z-alignment is also sorted out.
    implicit none

    INTEGER(i4)           :: i,NBATOMS
    REAL(dp)              :: alignment = 1.D0, range = 360.D0 
    REAL(dp)              :: ROTMAT(3,3)      
    REAL(dp),ALLOCATABLE  :: rotations(:)
    REAL(dp),ALLOCATABLE  :: initial_XYZ(:,:),final_XYZ(:,:)

! -- Initialising the incrementally rotated matrix
    ALLOCATE(rotations(3))
    rotations = (/0.D0,0.D0,0.D0/)
    ALLOCATE(final_XYZ(NBATOMS,3))
    final_XYZ = initial_XYZ

    DO i=1,3
      CALL OPT_ROT(range,i,NBATOMS,final_XYZ,rotations,alignment)
    END DO
    WRITE(*,'(A,F12.6,A,F12.6,A,F12.6)') '  --->  Range: ',range,', angle increment: ',range/1000.D0,' deg., alignment: ',alignment

    CALL BUILD_ROTMAT(rotations,ROTMAT)

  END SUBROUTINE ALIGN

! ----------

  SUBROUTINE OPT_ROT(range,axis,NBATOMS,XYZ,rotations,alignment)
! If range is set to 30, then XYZ will be rotated from -15 deg (2.*pi*(-15.)/360.) to +15 deg (2.*pi*15./360.)
! over 100 steps, i.e in increments of (2.*pi*range/360.)/100.
! An initial rotation of (-(2.*pi*(range/2.)/360.) - increment) will be performed, 
! so the system explores angles around the previous aligned geom, and not just above.
! Alignment is done here one-shot because to build a rotation matrix, the order of transformations matter,
! so the rotations() vector should become a tensor, and the final rotmat could be built in a stepwise manner.
! This is unnecessary: giving a small enough angle pace works fine (at the cost of memory and runtime, but fine here).
    implicit none

    INTEGER(i4)           :: i,NBATOMS,axis
    INTEGER(i4),PARAMETER :: N10 = 255, N5 = 246, C1 = 253, C2 = 254
    REAL(dp)              :: range,initial_rot,stepwise_rot,do_rotation,alignment
    REAL(dp)              :: N10_N5(3),C1_C2(3),unit_N10_N5(3),unit_C1_C2(3)
    REAL(dp)              :: to_maximise(1000)
    REAL(dp),ALLOCATABLE  :: rotations(:)
    REAL(dp),ALLOCATABLE  :: XYZ(:,:),working_XYZ(:,:)

    ALLOCATE(working_XYZ(NBATOMS,3))
    working_XYZ = XYZ

    initial_rot = (-2.D0*pi*(range/2.D0) / 360.D0) - (2.D0*range / 360.D0)/1000.D0
    stepwise_rot = (2.D0*pi*(range/360.D0))/1000.D0

    IF (axis .EQ. 1) CALL CC_ROTATION_X(initial_rot,NBATOMS,working_XYZ,working_XYZ)
    IF (axis .EQ. 2) CALL CC_ROTATION_Y(initial_rot,NBATOMS,working_XYZ,working_XYZ)
    IF (axis .EQ. 3) CALL CC_ROTATION_Z(initial_rot,NBATOMS,working_XYZ,working_XYZ)

    DO i=1,1000
      IF (axis .EQ. 1) CALL CC_ROTATION_X(stepwise_rot,NBATOMS,working_XYZ,working_XYZ)
      IF (axis .EQ. 2) CALL CC_ROTATION_Y(stepwise_rot,NBATOMS,working_XYZ,working_XYZ)
      IF (axis .EQ. 3) CALL CC_ROTATION_Z(stepwise_rot,NBATOMS,working_XYZ,working_XYZ)

! -- tmp vars to make the DOT_PROD call more legible
      N10_N5 = ( RESHAPE(working_XYZ(N5,:),(/3/)) - RESHAPE(working_XYZ(N10,:),(/3/)) )
      C1_C2  = ( RESHAPE(working_XYZ(C2,:),(/3/)) - RESHAPE(working_XYZ(C1,:),(/3/)) )
      unit_N10_N5 = N10_N5 / NORM2(N10_N5)
      unit_C1_C2  = C1_C2 / NORM2(C1_C2)

      to_maximise(i) = DOT_PRODUCT(unit_N10_N5,(/1.D0,0.D0,0.D0/))+DOT_PRODUCT(unit_C1_C2,(/0.D0,1.D0,0.D0/))
!      WRITE(*,'(3F8.4)') (initial_rot+i*stepwise_rot),to_maximise(i)
    END DO

! -- Identify and perform the rotation that yields the best alignment
    do_rotation = initial_rot + MAXLOC(to_maximise,1)*stepwise_rot
    alignment = MAXVAL(to_maximise)

    IF (axis .EQ. 1) THEN
      WRITE(*,'(2(A,F8.4))') " * Must rotate about x by : ",do_rotation," radians for an alignment of ", MAXVAL(to_maximise)
      CALL CC_ROTATION_X(do_rotation,NBATOMS,XYZ,XYZ)
      rotations(1) = rotations(1) + do_rotation
    ELSE IF (axis .EQ. 2) THEN
      WRITE(*,'(2(A,F8.4))') " * Must rotate about y by : ",do_rotation," radians for an alignment of ", MAXVAL(to_maximise)
      CALL CC_ROTATION_Y(do_rotation,NBATOMS,XYZ,XYZ)
      rotations(2) = rotations(2) + do_rotation
    ELSE IF (axis .EQ. 3) THEN
      WRITE(*,'(2(A,F8.4))') " * Must rotate about z by : ",do_rotation," radians for an alignment of ", MAXVAL(to_maximise)
      CALL CC_ROTATION_Z(do_rotation,NBATOMS,XYZ,XYZ)
      rotations(3) = rotations(3) + do_rotation
    END IF

    DEALLOCATE(working_XYZ)
  END SUBROUTINE OPT_ROT

! ----------

  SUBROUTINE BUILD_ROTMAT(rotations,ROTMAT)
! -- Builds the rotation matrix, applies it onto the geom to test it is correct.
    implicit none

    INTEGER(i4)          :: i
    REAL(dp)             :: alpha,beta,gamma
    REAL(dp)             :: Rx(3,3),Ry(3,3),Rz(3,3),ID(3,3)
    REAL(dp)             :: ROTMAT(3,3)
    REAL(dp),ALLOCATABLE :: rotations(:)

    alpha = rotations(1)
    beta  = rotations(2)
    gamma = rotations(3)

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