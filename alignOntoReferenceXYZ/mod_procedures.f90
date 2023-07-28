MODULE mod_procedures
  USE mod_precision
  implicit none

  CONTAINS

  SUBROUTINE READ_XYZ(INFILE,NBATOMS,ATNAMES,XYZ)
    implicit none

    INTEGER                      :: i,io,ALLOCSTAT,NBATOMS
    REAL(dp),ALLOCATABLE         :: XYZ(:,:)
    CHARACTER(LEN=*)             :: INFILE
    CHARACTER(LEN=4),ALLOCATABLE :: ATNAMES(:) 

    OPEN(10,file=INFILE,status='OLD',action='READ',iostat=io)
    IF (io .NE. 0) STOP "Failed to open coordinate file."
    READ(10,*) NBATOMS
    READ(10,*)
    
    IF (.NOT. ALLOCATED(XYZ)) ALLOCATE(XYZ(NBATOMS,3),stat=ALLOCSTAT)
    IF (ALLOCSTAT .NE. 0) STOP '**** ALLOC. FAILED: NOT ENOUGH MEMORY ****'

    IF (.NOT. ALLOCATED(ATNAMES)) ALLOCATE(ATNAMES(NBATOMS),stat=ALLOCSTAT)
    IF (ALLOCSTAT .NE. 0) STOP '**** ALLOC. FAILED: NOT ENOUGH MEMORY ****'
    
    DO i=1,NBATOMS
      READ(10,*) ATNAMES(i), XYZ(i,:)
    END DO
    CLOSE(10)
  END SUBROUTINE READ_XYZ

! ----------

  SUBROUTINE MODIFY_FILENAME(INFILE,OUTFILE,TO_APPEND)
    implicit none

    INTEGER          :: cut
    CHARACTER(LEN=*) :: INFILE,OUTFILE,TO_APPEND

    WRITE(*,'(A)') INFILE
    cut = (INDEX(INFILE,'.xyz') - 1)

    OUTFILE = INFILE(1:cut)//TO_APPEND
    WRITE(*,'(A)') OUTFILE

  END SUBROUTINE MODIFY_FILENAME

! ----------

  SUBROUTINE WRITE_XYZ(OUTFILE,NBATOMS,ATNAMES,XYZ)
    implicit none

    INTEGER                      :: i,io,NBATOMS
    REAL(dp),ALLOCATABLE         :: XYZ(:,:)
    CHARACTER(LEN=4),ALLOCATABLE :: ATNAMES(:)
    CHARACTER(LEN=*)             :: OUTFILE

    OPEN(10,file=OUTFILE,status='UNKNOWN',action='WRITE',iostat=io)
    IF (io .NE. 0) STOP "**** Failed to open outfile for writing. ****"
    WRITE(10,'(I4)') NBATOMS
    WRITE(10,'(A)') 'Translated/rotated coordinates'
    DO i=1,NBATOMS
      WRITE(10,'(A,3F12.6)') ATNAMES(i),XYZ(i,:)
    END DO
    WRITE(10,*)
    CLOSE(10)

  END SUBROUTINE WRITE_XYZ

! ----------

  SUBROUTINE WRITE_TRAJ(OUTFILE,NBATOMS,ATNAMES,XYZ)
    implicit none

    INTEGER                      :: i,io,NBATOMS
    REAL(dp),ALLOCATABLE         :: XYZ(:,:)
    CHARACTER(LEN=4),ALLOCATABLE :: ATNAMES(:)
    CHARACTER(LEN=*)             :: OUTFILE

    OPEN(10,file=OUTFILE,status='UNKNOWN',action='WRITE',position='APPEND',iostat=io)
    IF (io .NE. 0) STOP "**** Failed to open trajfile for writing. ****"
    WRITE(10,'(I4)') NBATOMS
    WRITE(10,'(A)') 'Translated/rotated coordinates'
    DO i=1,NBATOMS
      WRITE(10,'(A,3F12.6)') ATNAMES(i),XYZ(i,:)
    END DO
    CLOSE(10)

  END SUBROUTINE WRITE_TRAJ

! ----------

  SUBROUTINE WRITE_MAT_TO_FILE(OUTFILE,MAT)
    implicit none

    INTEGER          :: i,io
    REAL(dp)         :: MAT(3,3)
    CHARACTER(LEN=*) :: OUTFILE

    OPEN(10,file=OUTFILE,status='UNKNOWN',action='WRITE',iostat=io)
    IF (io .NE. 0) STOP "**** Failed to open outfile for writing. ****"
    DO i=1,3
      WRITE(10,'(3F12.6)') MAT(i,:)
    END DO
    CLOSE(10)

  END SUBROUTINE WRITE_MAT_TO_FILE

END MODULE mod_procedures