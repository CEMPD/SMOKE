      SUBROUTINE GETDSF(NSCEN,INERR,RC) 
C 
C  GETDSF reads Diesel sales fractions. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  05 Jul 00 @EPA-ddj bug#169 Added bounds checks for diesel fractions 
C  08 Jun 00 AIR Task 03: Removed nonblock DO constructs. 
C  25 May 2000 @EPA-elg Inputs are now Free Format rather than Fixed Format. 
C  25 Jan 2000 @EPA-djb changes in arrays. 
C  08 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    common blocks: 
C    /IOUCOM/ IOUGEN 
C    /MAXIMA/ MAXYRS 
C 
C  Output on return: 
C 
C    common block: 
C    /FLAGS1/ DSFLAG,GLDSFLAG 
C    /REGISF/ USRGSF 
C 
C  Local variable / array dictionary: 
C 
C  Name    Type              Description 
C  ------  ----  ---------------------------------------------------- 
C   JDX     I     Year index for diesel sales fractions 
C   IVTEMP  I     vehicle class index for class groupings 
C 
C  Notes: 
C 
C  GETDSF is new for Mobile6. This code used to reside in RDCITY. 
C     
      IMPLICIT NONE  
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'REGISF.I' 
C 
      INTEGER RC 
C 
      INTEGER JDX 
      INTEGER IVTEMP 
      INTEGER NSCEN, INERR 
C 
      RC = 0 
C 
C  Optional user-input of gas/diesel sales fractions (as diesel!) 
C  In Mobile5, this subroutine would be called because DSFLAG was 2. In 
C  Mobile6, a labeled input line instructed us to call this routine and 
C  set DSFLAG to 2. 
C 
      DSFLAG = 2 
C 
C  Read diesel sales fraction of model year sales. 
C  Read rows in JDX order (youngest age to oldest age). 
C  Vehicle class groupings are: 
C 
C     1 : LDV 
C     2 : LDT1 
C     3 : LDT2 
C     4 : LDT3 
C     5 : LDT4 
C     6 : HDV2B 
C     7 : HDV3 
C     8 : HDV4 
C     9 : HDV5 
C    10 : HDV6 
C    11 : HDV7 
C    12 : HDV8A 
C    13 : HDV8B 
C    14 : HDBS (school bus) 
C 
      DO IVTEMP=1,14 
        READ(IOUGEN,*,ERR=95,END=96) 
     *      (USRGSF(JDX,IVTEMP),JDX=1,25) 
      END DO 
C 
      CALL QUITER(0.0,0,614,INERR) 
C 
      DO JDX=1,MAXYRS 
        DO IVTEMP=1,14 
          IF(USRGSF(JDX,IVTEMP).LT.0.0  .OR. 
     &       USRGSF(JDX,IVTEMP).GT.1.0) THEN 
          CALL QUITER(USRGSF(JDX,IVTEMP),0,88,INERR) 
          USRGSF(JDX,IVTEMP)=0.0 
          END IF 
        END DO 
      END DO 
C 
C     Copy the data to the global data area, if we have not 
C     already started reading input for some scenario. 
C 
      IF (NSCEN.EQ.0) THEN 
        GLDSFLAG = 2 
        GLUSGF = USRGSF 
      ENDIF 
C 
      GOTO 99 
C 
C  Error or EOF on any attempted read. 
C 
   95 RC = 2 
   96 RC = RC - 1 
C 
   98 CALL QUITER(0.,0,138,INERR) 
C 
   99 RETURN 
      END 
