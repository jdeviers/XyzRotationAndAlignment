PROGRAM reorient
  USE mod_precision
  USE mod_procedures
  USE mod_translation
  USE mod_rotations
  USE mod_alignment
  implicit none
!
! The vector rotations(alpha,beta,gamma) keeps track of cumulated rotations. 
! Can finally reconstruct a rotation matrix from it once the molecule is aligned.
! This method is safer than multiplying successive rotation matrices. 
!
  INTEGER(i4)                  :: nb_atoms
  INTEGER(i4),PARAMETER        :: N10_index = 255
  REAL(dp)                     :: rot_matrix(3,3)
  REAL(dp),ALLOCATABLE         :: rotations(:)
  REAL(dp),ALLOCATABLE         :: xyz_ini(:,:),xyz_trans(:,:),xyz_rot(:,:)
  CHARACTER(LEN=4),ALLOCATABLE :: atom_names(:) 
  CHARACTER(LEN=50)            :: xyzfile,outfile

! -- Read atomic positions from xyz file
  CALL GETARG(1,xyzfile)
  CALL READ_XYZ(xyzfile,nb_atoms,atom_names,xyz_ini)

! -- Translate atomic positions so N10 (index 255 for FADH) is at the origin
  CALL TRANSLATE_TO_ORIGIN(xyz_ini,nb_atoms,N10_index,xyz_trans) 

! -- Rotate atomic positions to align 
  CALL ALIGN(nb_atoms,xyz_trans,xyz_rot,rotations,rot_matrix)
  CALL MODIFY_FILENAME(xyzfile,outfile,'_rot.xyz')
  CALL WRITE_XYZ(outfile,nb_atoms,atom_names,xyz_rot)

  WRITE(*,'(/,A,3F8.2,/)') 'Rotations: ',rotations(:)

! -- Write rotation matrix to file
  CALL MODIFY_FILENAME(xyzfile,outfile,'_rotmat.dat')
  CALL WRITE_MAT_TO_FILE(outfile,rot_matrix)

  DEALLOCATE(rotations)
  DEALLOCATE(atom_names,xyz_ini,xyz_trans,xyz_rot)
END PROGRAM reorient
