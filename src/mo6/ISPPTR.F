      FUNCTION ISPPTR(MY,IP,IV,ICASE) 
C 
C  ISPPTR selects the pointer to the bag speed normalization coefficients in 
C  arrays SK1/SK2C/SKG/SK2N & the scf's in SCUNA1 for LDGV, LDGT, & MC. 
C 
C  Called by BIGCFX and BIGSAL. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  23 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: MY,IP,IV,ICASE 
C    common blocks: 
C    /REGION/ IREJN 
C    /SPEED2/ MYGSP1,MYGSMC,MYGSP2 
C    /SPEED4/ MAXSP1,LBLAST 
C    /SPEED5/ MYHSPC,MAXHSP 
C 
C  Output on return: 
C 
C    function: ISPPTR 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  IUB1ST   I    Upper Bound my of 1ST myg of new LDGV/T equation case 
C 
C  Notes: 
C 
C  Error code 61 causes QUITER to print an error message referencing the 
C  calling pointer function and the failure to find a myg pointer.  QUITER 
C  then stops execution.  Hence, control never returns to the calling function. 
C  No explicit branch after the call out of the function is needed. 
C 
C 23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C           certification naming conventions. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'REGION.I' 
      INCLUDE 'SPEED2.I' 
      INCLUDE 'SPEED4.I' 
      INCLUDE 'SPEED5.I' 
C 
      INTEGER ICASE, MY, IV, IG, IUB1ST, IP, ISPPTR 
C 
C  Branch to myg index array lookup for the scf coefficients case requested. 
C 
      GOTO(10,20,30,40),ICASE 
C 
C  Case 1: "Old" (= Not changed since MOBILE2) scf coefficients for the fifth 
C          degree polynomial equation for LDGV/T.  These are dependent on the 
C          my, altitude and vehicle class. 
C 
   10 DO 15 IG=1,MAXSP1 
      IF(MY.LE.MYGSP1(IG,IREJN,IV)) GOTO 90 
   15 CONTINUE 
C 
      IG=1 
      CALL QUITER(0.,12,61,MY) 
C 
C  Case 2: "Old" (= Not changed since MOBILE2) scf coefficients for MC. 
C          These are dependent on the my, altitude and pollutant type. 
C 
   20 DO 25 IG=1,MAXSP1 
      IF(MY.LE.MYGSMC(IG,IREJN,IP)) GOTO 90 
   25 CONTINUE 
C 
      IG=1 
      CALL QUITER(0.,13,61,MY) 
C 
C  Case 3: coefficients for new for MOBILE4 A/X+B & A+B*X+C*X**2 scf equations. 
C          for LDGV/T.  These are dependent on my, pollutant type and vehicle 
C          class. 
C 
   30 IUB1ST=MYGSP2(IP,IV) 
      IF(MY.LE.IUB1ST) IG=1 
      IF(MY.GT.IUB1ST.AND.MY.LE.LBLAST) IG=MY-IUB1ST+1 
      IF(MY.GT.LBLAST) IG=LBLAST-IUB1ST+1 
      GOTO 90 
C 
C  Case 4: Subcase of 3 - a different set of A/X + B coefficients is used 
C          when the inputted speed exceeds FTP speed, the pollutant is HC or CO, 
C          the vehicle type is LDGV/T and the my is at or beyond a specified 
C          year.  This second set of new coefficients vary with my and iv. 
C 
   40 IG=MY-MYHSPC(1,IV)+1 
      IF(MY.GT.MYHSPC(2,IV)) IG=MAXHSP(IV) 
C 
C  Assign the myg cell pointer. 
C 
   90 ISPPTR=IG 
C 
      RETURN 
      END 
