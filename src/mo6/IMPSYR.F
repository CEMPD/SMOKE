      INTEGER FUNCTION IMPSYR(ICY,MY,IV) 
C 
C  IMPSYR finds the year that I/M coverage begins for a given model year 
C  and vehicle type. This may involve a change in I/M programs as long as 
C  the programs operate continuously (whenever any program ends, another 
C  program starts the following year). 
C 
C  Called by PCLEFT 
C 
C  Calls FINDIM 
C 
C  Changes: (Last change first) 
C 
C  20 Sep 00 @EPA- BG Added constant IEXHAUST=1 for call to FINDIM 
C 15 Sep 00 AIR Task 08: Updated FINDIM calls to select exhaust I/M. 
C  2 Jun 99 @DynTel-MLA 2-000  Account for the I/M program's grace period. 
C 10 Feb 99 @DynTel-MLA 2-663  This is a new routine for Mobile6. 
C 
C 
C  Input on call: 
C 
C    common blocks: 
C    /IMPAR1/ ICYIM,MODYR 
C    /VVNAMS/ VVGAS     (IVTYPE.I) 
C    /VVTYPS/ VTGASBUS  (IVTYPE.I) 
C 
C  Output on return: 
C    argument list: 
C    function FINDIM 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------ 
C  ICOV     I    A calendar year covered by some I/M program 
C  IMPGM    I    I/M program number 
C  IVTAM    I    Veh. type (mapped into range of vehs. subject to I/M) 
C 
C  Notes: 
C 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IVTYPE.I' 
C 
      INTEGER FINDIM 
C 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: MY 
C 
      INTEGER IYCOV 
      INTEGER IMPGM 
      INTEGER MYCOV 
      INTEGER, SAVE :: IEXHAUST=1 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      IMPSYR = 0 
C 
      IYCOV = ICY 
      MYCOV = MY 
C 
   15 IMPGM = FINDIM(IYCOV,MYCOV,IV,IEXHAUST) 
C 
C     If there was a program operating in the specified calendar year, 
C     look for an earlier program covering the year just prior to the 
C     start year of the current program. 
C 
      IF (IMPGM.NE.0) THEN 
C 
        IMPSYR = ICYIM(1,IMPGM) 
        IYCOV = ICYIM(1,IMPGM) - 1 
        IF (MYCOV.GT.IYCOV) MYCOV = IYCOV 
        GOTO 15 
C 
      END IF  
C 
   99 CONTINUE 
      RETURN 
      END 
