      SUBROUTINE ATACRO(CAT,IV) 
C 
C  Ths subroutine computes the Air Toxic conversion factors 
C  for Acrolein. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX2.I'    ! ATACROLINE 
      INCLUDE 'ATOX4.I'    ! ATCONV 
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
       
      RETURN 
      END 
       
