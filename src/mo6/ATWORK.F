      MODULE ATHEAP 
C 
C  This module sets up the large DBAT and ATHBER arrays 
C  to allocatable status. 
C 
      REAL, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: DBAT 
      REAL, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: ATHBER 
C 
      END MODULE 
C 
C---------------------------------------------------------------------- 
C 
      SUBROUTINE ATALLOC 
C 
C  This routine allocates the large DBAT and ATHBER arrays 
C 
      USE ATHEAP 
C 
      ALLOCATE (DBAT(25,24,5,6,6,28)) 
      ALLOCATE (ATHBER(2,25,2,4,28,50)) 
C 
      RETURN 
      END 
C 
C---------------------------------------------------------------------- 
C 
      SUBROUTINE ATDEALLOC 
C 
C  This routine allocates the large DBAT and ATHBER arrays 
C 
      USE ATHEAP 
C 
      DEALLOCATE (DBAT) 
      DEALLOCATE (ATHBER) 
C 
      RETURN 
      END 
