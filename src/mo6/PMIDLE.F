      SUBROUTINE PMIDLE(ICY) 
C 
C  This subroutine returns the PM idle emission factor (in g/hr) 
C  for calendar year ICY. 
C 
C Changes: 
C 
C 07 May 01 AIR Task 11: New Subroutine for MOBILE6.1. Updated 
C           code to handle all MY/JDX. Removed code which 
C           incorrectly added g/mi values to g/hr emission rates. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Called by EFCALXP5. 
C 
C  Calls DIESLF 
C 
C  Input on call: 
C    common blocks: /CARBON/ SOFPNT 
C 
C 
C  Output on return: CIDLE 
C 
C     Local Variable Dictionary: 
C 
C  Variable    Type      Description 
C ----------------------------------------------------------------------- 
C   CIDLE      R         Idle emission factor, g/hr 
C   IDLPNT     I         MY pointer 
C   IDX        I         Model year indicator 
C   IPTR       I         MY pointer 
C   IV         I         Vehicle class 
C   IVPART     i         PART5 vehicle class 
C   MY         I         Model year 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IDLE.I'    ! NIDLYR, IDLYRS, BASIDL 
      INCLUDE 'IVTYPE.I'  ! MAXVEH, VDM6PRT, VTHDDV2B, VTHDDV8B, VTURBAN, VTSCHOOL 
      INCLUDE 'MAXIMA.I'  ! MAXIH, MAXYRS 
      INCLUDE 'MYRCAL.I'  ! TF 
      INCLUDE 'PART3.I'   ! VPM, EFPM 
      INCLUDE 'PSIZE.I'   ! PSDIE 
      INCLUDE 'VMXCOM.I'  ! VMTMIX 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
C 
C  Declare local variables/arrays. 
C 
      INTEGER IDLPNT 
      INTEGER IDX 
      INTEGER IPTR 
      INTEGER IV 
      INTEGER IVPART 
      INTEGER JDX 
      INTEGER MY 
      REAL    CIDLE 
C 
C  Loop for each vehicle class 
C 
      DO IV=1,MAXVEH 
C 
C  PM idle emissions are only available for diesel 2B-8B, and 
C  urban/school busses 
C 
      IF(IV.GE.VTHDDV2B .AND. IV.LE.VTHDDV8B .OR. 
     *   IV.EQ.VTURBAN  .OR.  IV.EQ.VTSCHOOL) THEN 
C 
C  Skip if VMTMIX is zero 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
C 
C  Map the MOBILE6 vehicle classes into PART5 diesel space. 
C 
      IVPART = VDM6PRT(IV) - 7 
C 
C  Loop for all ages 
C 
        DO JDX=1,MAXYRS 
C 
C  Set IDX 
C 
          IDX=(MAXYRS+1)-JDX 
C 
C  If travel fraction for this age/vehicle is 0, skip. 
C 
          IF(TF(IDX,IV) .EQ. 0.) CYCLE 
C 
C  Set model year 
C 
          MY=ICY-JDX+1 
C 
C  Determine the model year pointer 
C 
      DO 10 IPTR=1,NIDLYR 
      IF(MY .LE. IDLYRS(IPTR)) GOTO 11 
   10 CONTINUE 
 
   11 IDLPNT=IPTR 
C 
C     Multiply by the particle size cutoff, PSDIE: 
C 
      CIDLE = BASIDL(IDLPNT,IVPART) * PSDIE 
C 
C  Store in arrays, weighted by travel fractions 
C  (NOTE: values are NOT stored in database arrays) 
C 
          EFPM(11,IV) = EFPM(11,IV) + CIDLE * TF(IDX,IV) 
C 
        END DO !JDX 
C 
        VPM(11) = VPM(11) + EFPM(11,IV) * VMTMIX(IV) 
C 
      ENDIF 
C 
      END DO !IV 
C 
  99  RETURN 
      END 
