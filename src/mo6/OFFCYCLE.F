      REAL FUNCTION OFFCYCLE(INH,IP,BER) 
C 
C    Computes the off-cycle emission rate for light duty gasoline 
C    vehicles and trucks.  The off-cycle emission rate is added to 
C    the base emissions prior to applying the speed correction. 
C 
C    Called by: SPEED_CF 
C 
C    Calls: 
C 
C 
C    Changes: (Last change first) 
C 
C    1 May 00 @EPA-ddj     1-020, removed DATA statements 
C    7 Jan 00 @DynTel-ddj 1-020, New Function 
C 
C 
C  Input on call: 
C    argument list: INH,IP,BER 
C 
C    common blocks: 
C 
C 
C  Output on return: 
C    common blocks: 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C   P1      R    Array of linear terms off-cycle coefficients 
C   P2      R    Array of square terms off-cycle coefficients 
C 
      IMPLICIT NONE 
C 
      INTEGER, INTENT(IN) :: INH 
      INTEGER, INTENT(IN) :: IP 
C 
      INTEGER, PARAMETER  :: CO   = 2 
      INTEGER, PARAMETER  :: NOX  = 3 
      INTEGER, PARAMETER  :: HIGH = 2 
C 
      REAL, INTENT(IN)    :: BER 
      REAL, DIMENSION(3)  :: P1 = (/ 0.305,    0.984,    0.332/) 
      REAL, DIMENSION(3)  :: P2 = (/-0.02492, -0.07638, -0.04745/) 
C 
      IF(IP .EQ. CO .AND. INH .EQ. HIGH) THEN         ! CO emitter, no 
         OFFCYCLE = 0.0                               ! correction. 
      ELSE  IF(IP .EQ. NOX .AND. BER .GE. 3.50) THEN 
            OFFCYCLE = 0.581                          ! Cap NOx. 
          ELSE 
             OFFCYCLE = BER * ( P1(IP) + BER*P2(IP) ) ! All other. 
      END IF 
C 
      OFFCYCLE = MAX(OFFCYCLE,0.0) 
C 
      RETURN 
      END FUNCTION 
