MODULE mod_translation
  USE mod_precision
  implicit none

  CONTAINS

  SUBROUTINE TRANSLATE_TO_ORIGIN(XYZ,NBATOMS,AT_INDEX,XYZ_translated)
    implicit none

    INTEGER(i4)          :: NBATOMS,AT_INDEX
    REAL(dp),ALLOCATABLE :: XYZ(:,:),XYZ_translated(:,:)

    ALLOCATE(XYZ_translated(NBATOMS,3))
! -- Just in case N10 is already at origin:
    XYZ_translated = XYZ

! -- And this bit does the translation of needed:
    IF (XYZ(AT_INDEX,1) .NE. 0.0D0) XYZ_translated(:,1) = XYZ(:,1) - XYZ(AT_INDEX,1)  
    IF (XYZ(AT_INDEX,2) .NE. 0.0D0) XYZ_translated(:,2) = XYZ(:,2) - XYZ(AT_INDEX,2) 
    IF (XYZ(AT_INDEX,3) .NE. 0.0D0) XYZ_translated(:,3) = XYZ(:,3) - XYZ(AT_INDEX,3) 

  END SUBROUTINE TRANSLATE_TO_ORIGIN

END MODULE mod_translation
