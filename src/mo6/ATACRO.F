      SUBROUTINE ATACRO(CAT,IV) 
C 
C  Ths subroutine computes the Air Toxic conversion factors 
C  for Acrolein. 
C 
C  Changes:  (Last change first) 
C 
C  18 Sep 01 AIR Task 15-22: Comments 
C  12 Jul 01 AIR Task 19: New Subroutine. 
C 
C  Called by:  ATCALX 
C 
C  Calls:  None 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX2.I'    ! ATACROLINE 
      INCLUDE 'ATOX4.I'    ! ATCONV 
C 
C Input arguments 
C 
      INTEGER IV 
      REAL    CAT(3) 
C 
C  Compute the Acrolein fractions, weighted by catalyst. 
C  (Note: At this time, the fractions are independent of 
C  low or high emitter class). Acrolein is TOG based, 
C  so IVT is always 2. 
C 
      ATCONV(1,2,6) = ATACROLEIN(IV,1) *  CAT(1)       !Normals 
     *              + ATACROLEIN(IV,2) * (CAT(2) + CAT(3)) 
C 
      ATCONV(2,2,6) = ATCONV(1,2,6)                    !Highs 
C       
      RETURN 
      END 
       
