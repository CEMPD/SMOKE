      SUBROUTINE STOERQ(INERR) 
C 
C  STOERQ attempts to store the user-supplied specifications for new FTP basic 
C  emission rate equations and when successful (i.e., when room found in the 
C  user supplied rate arrays). 
C 
C  Called by GETERQ. 
C 
C  Calls QUITER. 
C   
C  Changes: (Last change first) 
C 
C  02 Aug 00 @EPA-bag Bug62  Added conditional to skip any post-1980 LDG 
C           parameters entered with label BASIC EFS. 
C   3 Jul 00 @EPA-ddj Bug62  Added conditional to skip any 1981-1993 LDG 
C           parameters entered with label BASIC EFS. 
C  23 Aug 99 @DynTel-BG 2-660 replaced constants for ISR=1(=RUN) 
C            and ISR=2(=START). 
C  17 Feb 1999 @DynTel-rwk 2-668 Pre-81 code, removed MOBILE5 50K 
C                                deterioration code. 
C  08 Jan 1999 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: INERR 
C    common blocks: 
C    /BASEQ4/ NEWOML,NEWPAR,BERNEW,NEWCT 
C    /BASEQ6/ MAXERU 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /BASEQ4/ NEWFIT 
C    /BASEQ5/ ERUZML,ERUDR,ESUZML,ESUDR 
C    /BASEQ6/ MYGERU,NUMERU 
C 
C  Notes: 
C 
C 13 Apr 98  @DynTel-RJD  Added arrays ESUZML and ESUDR to accommodate 
C             starts emissions. 
C 
C  None. 
C 
      IMPLICIT NONE 
      INCLUDE 'BASEQ4.I' 
      INCLUDE 'BASEQ5.I' 
      INCLUDE 'BASEQ6.I' 
      INCLUDE 'IVTYPE.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(INOUT) :: INERR 
C 
C     Declare local variables 
C 
      INTEGER IR 
      INTEGER IV 
      INTEGER ISR 
      INTEGER IFAC 
      INTEGER IERNEW 
      INTEGER IP 
      INTEGER, PARAMETER :: RUN=1,START=2 
C 
      DO IFAC = 1, NEWCT 
        NEWFIT(IFAC)=0 
        IR=NEWPAR(1,IFAC) 
        IV=NEWPAR(2,IFAC) 
        IP=NEWPAR(3,IFAC) 
C 
C     skip if inputs are for post-1980 light duty gas vehicles 
C     get_8193, get94e, and get94p deal with those vehicles 
C 
        IF((VVLDGAS(IV)==1) .AND. 
     &     (NEWPAR(4,IFAC)>1980 .OR. NEWPAR(5,IFAC)>1980) ) CYCLE 
C 
          IF (NEWOML(IFAC)=="R".OR.NEWOML(IFAC)=="r") THEN 
            ISR=RUN 
          ELSE IF (NEWOML(IFAC)=="S".OR.NEWOML(IFAC)=="s") THEN 
            ISR=START 
          END IF 
C 
C  Check if there's room to store equation IFAC's parameters. 
C 
          IERNEW=NUMERU(IP,IV,IR,ISR)+1 
          IF(IERNEW > MAXERU) THEN 
            CALL QUITER(0.,IFAC,51,INERR) 
          ELSE 
C 
C  There is room. 
C 
            NUMERU(IP,IV,IR,ISR)=IERNEW 
            MYGERU(IERNEW,1,IP,IV,IR,ISR)=NEWPAR(4,IFAC) 
            MYGERU(IERNEW,2,IP,IV,IR,ISR)=NEWPAR(5,IFAC) 
            IF (ISR == RUN) THEN 
              ERUZML(IERNEW,IP,IV,IR)=BERNEW(1,IFAC) 
              ERUDR(IERNEW,IP,IV,IR)=BERNEW(2,IFAC) 
            ELSE IF (ISR == START) THEN 
              ESUZML(IERNEW,IP,IV,IR)=BERNEW(1,IFAC) 
              ESUDR(IERNEW,IP,IV,IR) =BERNEW(2,IFAC) 
            END IF 
            NEWFIT(IFAC)=1 
C 
          END IF    ! IERNEW greater than MAXERU 
      END DO 
C 
      RETURN 
      END 
