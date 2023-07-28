PROGRAM superimpose
  USE mod_precision
  USE mod_procedures
  USE mod_transrot
  USE mod_alignment
  implicit none
!
! The vector rotations(alpha,beta,gamma) keeps track of cumulated rotations. 
! Can finally reconstruct a rotation matrix one from it once the molecule is aligned.
! This method is safer than multiplying successive rotation matrices. 
!

! -- Ref geom var
  REAL(dp),ALLOCATABLE         :: xyz_ref(:,:)
  CHARACTER(LEN=50)            :: xyzfile_ref
! -- To be rotated geom var
  REAL(dp),ALLOCATABLE         :: xyz_torotate(:,:)
  CHARACTER(LEN=50)            :: xyzfile_torotate
! -- Common variables
  INTEGER                      :: nb_atoms
  CHARACTER(LEN=4),ALLOCATABLE :: atom_names(:)
! -- Output variables
  REAL(dp)                     :: rot_matrix(3,3)
  REAL(dp),ALLOCATABLE         :: rotations(:,:)
  REAL(dp),ALLOCATABLE         :: xyz_rotated(:,:)
  CHARACTER(LEN=50)            :: outfile


! -- Read atomic positions from xyz file, translate N5 to origin, update positions
  CALL GETARG(1,xyzfile_ref)
  CALL READ_XYZ(xyzfile_ref,nb_atoms,atom_names,xyz_ref)

  CALL TRANSLATE_N5_TO_ORIGIN(nb_atoms,xyz_ref,xyz_ref)
  CALL MODIFY_FILENAME(xyzfile_ref,outfile,'_trans.xyz')
  CALL WRITE_XYZ(outfile,nb_atoms,atom_names,xyz_ref)

  CALL GETARG(2,xyzfile_torotate)
  CALL READ_XYZ(xyzfile_torotate,nb_atoms,atom_names,xyz_torotate)

  CALL TRANSLATE_N5_TO_ORIGIN(nb_atoms,xyz_torotate,xyz_torotate)
  CALL MODIFY_FILENAME(xyzfile_torotate,outfile,'_trans.xyz')
  CALL WRITE_XYZ(outfile,nb_atoms,atom_names,xyz_torotate)

! -- Rotate atomic positions to align 
  CALL ALIGN(nb_atoms,xyz_ref,xyz_torotate,xyz_rotated,rotations,rot_matrix)
  CALL MODIFY_FILENAME(xyzfile_torotate,outfile,'_rot.xyz')
  CALL WRITE_XYZ(outfile,nb_atoms,atom_names,xyz_rotated)

  WRITE(*,'(/,A,3F8.2,/)') 'Rotations: ',rotations(:,:)

! -- Write rotation matrix to file
  CALL MODIFY_FILENAME(xyzfile_torotate,outfile,'_rotmat.dat')
  CALL WRITE_MAT_TO_FILE(outfile,rot_matrix)

  DEALLOCATE(rotations)
  DEALLOCATE(atom_names,xyz_ref,xyz_torotate,xyz_rotated)
END PROGRAM superimpose
