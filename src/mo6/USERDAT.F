      MODULE USERDAT 
C 
C  This module contains subroutines and functions involved in the reading and subsequent 
C  processing of User input Basic Emission Rates. 
C 
C 
C  Changes: (Last change first) 
C 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  02 Oct 00 @EPA-djb Bug230 Changed write to * to IOUOUT. 
C  07 Jul 00 @AIR- DK    changed case of subroutine statements to all caps 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  03 Aug 99 @Dyntel nh 2-660. Changed organization of linked lists 
C            obviating use of Count_user. This changed GET_USER, 
C            PR_SR_P, PR_AVG_P, CHK_USER, CHK_AVG,DEALLOC_8193 and 
C            deleted COUNT_USER 
C 
C   1 Jul 99 @Dyntel nh 2-660. Modified CHKLAB calls in USER_FILE to conform 
C            to user defined type changes in LABELS.FOR 
C  28 Jun 99 @Dyntel-NH 2-000  
C 
C subroutines/functions list: 
C 
C   Name     Type              Description 
C  ------    ----  ----------------------------------------------------- 
C GET_USER   SBR   reads and stores the user input records in 
C                       appropriate data structures. 
C PR_SR_P    SBR   prints the stored user-input running/start 
C                       records. This is a diagnostic function and its 
C                       invocation is commented out. 
C PR_AVG_P   SBR   same as PR_SR_P but for average records. 
C CHK_USER   SBR   checks for the existence of a relevant user-input 
C                      running/start record and if it exists calculates the 
C                      emission factor. 
c CHK_AVG    SBR   checkes for the existence of a relevant stored 
C                     user-input average record and if it exists returns a 
C                     pointer to these data. 
C HNUSER     RFCN  calculates the emissions when invoked by CHK_USER. 
C USER_file   SBR   tests for the existence of the user input file 
C                     referenced in the input file and then tests the label 
C                     integrity. 
C DEALLOC_8193 SBR Deallocates user input arrays 
C 
C These routines use allocatable arrays, user defined types and pointers. 
C SR_REC is a user defined type which holds the variables included in a USER 
C input record which if applicable will be used in place of the hard coded 
C Block Data. Additionally it includes a pointer to the next similar record 
C 
C USER DEFINED TYPES 
C 
C   TYPE Name  MEMBERS  MEMBER ATTRIBUTE            Description 
C    ------    -------  ----------------    ------------------------------------- 
C    SR_REC 
C               MYF      I                   First model Year. 
C               MYL      I                   Final Model Year 
C               ZML      R                   Zero mile 
C               DR       R                   Deterioration Rate. 
C               next     P                   Pointer to next record. 
C               SR_RECP  A                   SR_RECP( 
C                                            array of all normal/high user 
C                                            input records. If user input 
C                                            apply to a run the pointer in 
C                                            SR_P will point to the first 
C                                            record if any satisfying the 
C                                            subscript conditions. (otherwise 
C                                            it will not be associated) 
C 
C   AVG_REC                                  user defined type which holds the data from a user input record 
C                                            and a pointer to the next similar record. 
C               MYF      I                   First Model Year 
C               MYL      I                   Final Model Year 
C               ZML      R                   Zero Mile 
C               ADJ      R                   Adjustment 
C               DR1      R                   First Deterioration Rate 
C               DR2      R                   Second Deterioration Rate 
C               CPT1     R                   First Critical Point 
C               CPT2     R                   Second Critical Point 
C               NEXT     P                   Next record pointer 
C 
C 
C   SR1                                       a user defined type consisting 
C                                             only of pointers to type SR_REC. 
C                                             This ploy is necessary because 
C                                             a quirk of FORTRAN 90 disallows a 
C                                             straight array of pointers but 
C                                             does allow an array of User 
c                                             Defined Type which can be pointers! 
C 
C               HEAD     P 
C               TAIL P 
C 
C               SR_P     A                     SR_P(IP,IV,IFDS,ISR,INH) 
C                                              IP=3 Pollutant 
C                                              IV=5 Vehicle 
C                                              IFDS=3 Fuel Delivery System 
C                                              ISR=2 Start/Running 
C                                              INH=2 Normal/High 
C 
C   AVG1                                      the corresponding structure for 
C                                             average records. 
C               HEAD     P 
C               TAIL P 
C 
C               AVG_P    A                     AVG_P(IP,IV,IFDS) 
C                                              IP3 Pollutant 
C                                              IV=5 Vehicle 
C                                              IFDS=3 Fuel delivery System 
C 
      IMPLICIT NONE 
      TYPE SR_REC 
        INTEGER ::MYF,MYL 
        REAL :: ZML,DR 
         TYPE (SR_REC), POINTER::NEXT 
      END TYPE SR_REC 
 
      TYPE AVG_REC 
        INTEGER ::MYF,MYL 
        REAL ::ZML,ADJ,DR1,DR2,CPT1,CPT2 
         TYPE (AVG_REC),POINTER::NEXT 
      END TYPE AVG_REC 
 
      TYPE SR1 
        TYPE (SR_REC),POINTER :: HEAD,TAIL 
      END TYPE SR1 
 
      TYPE AVG1 
        TYPE (AVG_REC),POINTER :: HEAD,TAIL 
      END TYPE AVG1 
 
      TYPE (SR1),ALLOCATABLE :: SR_P(:,:,:,:,:) 
      TYPE (SR_REC),POINTER :: SRP,PNULL 
      TYPE (SR_REC) :: SR 
      TYPE (AVG1),ALLOCATABLE :: AVG_P(:,:,:) 
      TYPE (AVG_REC),POINTER :: AVGP,PNULLAVG 
      TYPE(AVG_REC) :: AVG 
      CHARACTER*80 FILENAME 
 
      CONTAINS 
      SUBROUTINE GET_USER() 
C 
C GETUSER reads and stores the user input records in appropriate 
C data structures. 
C 
C Changes 
C 
C 5-8-99 @Dyntel NH Request 2-660. Improved organization of record storage. 
C 6-7-99 @Dyntel NH Request 2-660 
C 
C Called by GET_8193 
C 
C       Inputs on call 
C                       NONE 
C       Outputs on return 
C                        sr_p and avg_p arrays and associated linked lists 
C                        populated. 
C       Local array subscripts: 
C                       None 
C 
C OUTPUT DATA STRUCTURES. 
C 
C TYPE (SR1),ALLOCATABLE :: SR_P(:,:,:,:,:) 
C TYPE (AVG1),ALLOCATABLE :: AVG_P(:,:,:,:,:) 
C 
C       Local variable/array dictionary 
C 
C       Name    Type            Description 
C       ----------------------------------- 
C       aline   Character*89    A line of user file is read into aline 
C       i       Integer         A general loop index 
C       k       Integer         Index for the sr_recp array 
C       k2      Integer         Index for the avg_recp array 
C       ip      Integer         Pollutant read from user supplied record 
C       iv      Integer         Vehicle read from user supplied record 
C       ifds    Integer         Fuel delivery system from user supplied record 
C       myf     Integer         First model year read from user supplied record 
c       myl     Integer         Final model year read from user supplied record 
C       isr     Integer         1 for running record 2 for start record 
C       inh     Integer         1 for Normal record, 2 for high record 
C       source  Integer         A logical IO unit. 
C 
C       Notes 
C 
C Get_user reads the user supplied input file and stores the records. The head pointer in 
C sr_p is set to point to the first record with matching indices for the 
C start/running records and the head pointer in avg_p to the first average record 
C Each of the corresponding tail pointers point to the last record and are only 
C used during population. 
 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
       
      INTEGER, EXTERNAL :: JUNIT 
       
      CHARACTER*89 ALINE 
      CHARACTER*1 ISRCH 
      INTEGER K,K2,IP,IV,IFDS,MYF,MYL,ISR,INH,SOURCE 
      REAL ZML,DR,ADJ,DR1,DR2,CPT1,CPT2 
       
ccs   Add declarations of M6DIR and M6FILE 
      CHARACTER(LEN=200) M6DIR 
      CHARACTER(LEN=280) M6FILE 
 
      SOURCE=JUNIT() 
 
ccs   All M6 input files live in SMK_M6PATH 
      CALL GETENV( "SMK_M6PATH", M6DIR ) 
      M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME 
         
      OPEN(FILE=M6FILE,UNIT=SOURCE,ACTION='READ',STATUS='OLD') 
      K=0 
      K2=0 
      NULLIFY(PNULL) 
      NULLIFY(PNULLAVG) 
      DO 
         READ(SOURCE,100,END=1000)ALINE 
         SELECT CASE (ALINE(1:1) ) 
             CASE ('N','n','H','h') 
               SELECT CASE (ALINE(1:1)) 
                   CASE ('N','n') 
                     INH=1 
                   CASE ('H','h') 
                     INH=2 
                   CASE DEFAULT 
               END SELECT 
               READ(ALINE,110)IP,IV,IFDS,MYF,MYL,ZML,DR,ISRCH 
               SELECT CASE (ISRCH) 
                   CASE ('S','s') 
                     ISR=2 
                   CASE ('R','r') 
                     ISR=1 
                   CASE DEFAULT 
                     PRINT*,"Bad record final must be S or R" 
                     STOP! add error handling 
               END SELECT 
               K=K+1 
 
              IF(.NOT.ASSOCIATED(SR_P(IP,IV,IFDS,ISR,INH)%HEAD))THEN 
                ALLOCATE(SRP) 
                SR=SR_REC(MYF,MYL,ZML,DR,PNULL) 
                SRP=SR 
                SR_P(IP,IV,IFDS,ISR,INH)%HEAD=>SRP 
                SR_P(IP,IV,IFDS,ISR,INH)%TAIL=>SRP 
              ELSE 
                SRP=>SR_P(IP,IV,IFDS,ISR,INH)%TAIL 
                ALLOCATE(SRP%NEXT) 
                SRP=>SRP%NEXT 
                SR=SR_REC(MYF,MYL,ZML,DR,PNULL) 
                SRP=SR 
                SR_P(IP,IV,IFDS,ISR,INH)%TAIL=>SRP 
              ENDIF 
             CASE ('A','a') 
               READ(ALINE,120)IP,IV,IFDS,MYF,MYL,ZML,ADJ,DR1,DR2,CPT1, 
     C            CPT2 
               K2=K2+1 
               IF (IV>1) THEN 
                IV=2 
               END IF 
              IF(.NOT.ASSOCIATED(AVG_P(IP,IV,IFDS)%HEAD))THEN 
                ALLOCATE(AVGP) 
                AVG=AVG_REC(MYF,MYL,ZML,ADJ,DR1,DR2,CPT1,CPT2,PNULLAVG) 
                AVGP=AVG 
                AVG_P(IP,IV,IFDS)%HEAD=>AVGP 
                AVG_P(IP,IV,IFDS)%TAIL=>AVGP 
              ELSE 
                AVGP=>AVG_P(IP,IV,IFDS)%TAIL 
                ALLOCATE(AVGP%NEXT) 
                AVGP=>AVGP%NEXT 
                AVG=AVG_REC(MYF,MYL,ZML,ADJ,DR1,DR2,CPT1,CPT2,PNULLAVG) 
                AVGP=AVG 
                AVG_P(IP,IV,IFDS)%TAIL=>AVGP 
              ENDIF 
             CASE DEFAULT 
         END SELECT 
      END DO 
1000  CONTINUE 
      CLOSE(SOURCE) 
100   FORMAT(A80) 
110   FORMAT(1X,3I2,2I5,1X,F9.5,1X,F10.6,1X,A1) 
120   FORMAT(1X,3I2,2I5,4(1X,F8.6),2(1X,F8.4)) 
      END SUBROUTINE 
 
 
      SUBROUTINE PR_SR_P() 
C 
C pr_sr_p prints the stored user supplied starts/running records. 
C 
C Changes 
C 
C 02 Oct 00 @EPA-djb Bug230 Changed write to * to IOUOUT. 
C 5 Aug 99  @Dyntel NH Request 2-660. Changes related to modified storage. 
C 17 Jun 99 @Dyntel NH Request 2-660 
C 
C Called by GET_8193 
C 
C       Local variable/array dictionary 
C 
C       Name    Type            Description 
C       ----------------------------------- 
C Local Array subscripts 
C 
C       chsr    Character*1     2 maps index into representative character 
C       chnh    Character*1     2 maps index into representative character 
C 
C Local Variables 
C       ip      Integer         Pollutant index 
C       iv      Integer         Vehicle index 
C       ifds    Integer         Fuel Delivery index 
C       isr     Integer         Start/Running index 
C       inh     Integer         Normal/Running index 
C 
C       Module variables (See declarations above) 
C 
C INPUT ARGUMENTS 
C                NONE 
C OUTPUT ARGUMENTS 
C                None. 
C OUTPUT DATA STRUCTURES. 
C                None 
C 
C       Local variable/array dictionary 
C 
C       Name    Type            Description 
C       ----------------------------------- 
C       ip      Integer         Pollutant read from user supplied record 
C       iv      Integer         Vehicke read from user supplied record 
C       ifds    Integer         Fuel delivery system from user supplied record 
C       myf     Integer         First model year read from user supplied record 
c       myl     Integer         Final model year read from user supplied record 
C       isr     Integer         1 for running record 2 for start record 
C       inh     Integer         1 for Normal record, 2 for high record 
C 
C 
C If invoked this subroutine will print the stored start/running records to 
C the default output 
C 
      INCLUDE 'IOUCOM.I' 
      INTEGER IP,IV,IFDS,ISR,INH 
      CHARACTER(LEN=1), DIMENSION(2) :: CHSR=(/'R','S'/) 
      CHARACTER(LEN=1), DIMENSION(2) :: CHNH=(/'N','H'/) 
C 
C For each entry in the start/running record pointer array test to for 
C association. If associated print the record and follow the chain until 
C not associated. 
C 
      DO INH=1,2 
         DO ISR=1,2 
            DO IFDS=1,3 
               DO IV=1,5 
                  DO IP=1,3 
                     SRP=>SR_P(IP,IV,IFDS,ISR,INH)%HEAD 
C 
C SRP is set to point to the record pointed to in the sr_p array for the corresponding 
C INH,ISR,IFDS,IV,IP. If this pointer is ASSOCIATED (with a record) that record is printed 
C The pointer chain is then followed (and printed) until a non associated pointer is found. 
C 
                        DO WHILE(ASSOCIATED(SRP)) 
                                WRITE(IOUOUT,100)CHNH(INH), 
     *                          IP,IV,IFDS,SRP%MYF,SRP%MYL,SRP%ZML, 
     *                                  SRP%DR,CHSR(ISR) 
                                SRP=>SRP%NEXT 
                        END DO 
                  END DO 
               END DO 
            END DO 
         END DO 
      END DO 
100   FORMAT(1X,A1,3I2,2I6,F8.4,F10.6,A2) 
      END SUBROUTINE 
 
      SUBROUTINE PR_AVG_P() 
C 
C pr_avg_p prints the stored user supplied average running records. 
C 
C Changes 
C 
C 02 Oct 00 @EPA-djb Bug230 Changed write to * to IOUOUT. 
C 5 Aug 99  @Dyntel NH Request 2-660. Minor change related to storage. 
C 17 Jun 99 @Dyntel NH Request 2-660 
C 
C Called by GET_8193 
C 
C Module variables (See declarations above) 
C 
C       Name 
C 
C       avg_p 
C       avgp 
C 
C       Local variable/array dictionary 
C 
C       Name    Type            Description 
C       ----------------------------------- 
C 
C       ip      Integer         Pollutant index 
C       iv      Integer         Vehicle index 
C       ifds    Integer         Fuel Delivery index 
C 
      INCLUDE 'IOUCOM.I' 
      INTEGER IP,IV,IFDS 
 
      DO IFDS=1,3 
        DO IV=1,2 
                DO IP=1,3 
                        AVGP=>AVG_P(IP,IV,IFDS)%HEAD 
C 
C AVGP is set to point to the record pointed to in the avg_p array for the 
C corresponding IFDS,IV,IP. If this pointer is ASSOCIATED (with a record) that 
C record is printed.The pointer chain is then followed (and printed) until a 
C non associated pointer is found. 
C 
                DO WHILE(ASSOCIATED(AVGP)) 
                        WRITE(IOUOUT,100) 
     C                  IP,IV,IFDS,AVGP%MYF,AVGP%MYL,AVGP%ZML,AVGP%ADJ, 
     C                 AVGP%DR1,AVGP%DR2,AVGP%CPT1,AVGP%CPT2 
                                AVGP=>AVGP%NEXT 
                        END DO 
                END DO 
        END DO 
      END DO 
100   FORMAT(1X,3I2,2I6,4(1X,F8.6),2(1X,F8.4)) 
      END SUBROUTINE 
 
       SUBROUTINE CHK_USER(MY,IP,IV,IFDS,ISR,INH,KMILES,HNBEFUSR,EXISTS) 
C 
C chk_user determines if a user supplied starts/running record is stored which 
C matches the model year, pollutant, vehicle, fuels delivery system, 
C starts/running and normal/high. If a record is available it is used to 
C calculate the emissions. 
C 
C CALLED BY EF81_93 
C 
C Calls hnbefusr 
C Changes 
C 
C 5 Aug 99  @Dyntel NH Request 2-660. Minor change related to storage. 
C 17 Jun 99 @Dyntel NH Request 2-660 
C 
C INPUT ARGUMENTS 
C                MY Model Year 
C                IP Pollutant 
C                IV Vehicle 
c                IFDS Fuel Delivery System 
C                ISR Start/Running 
C                INH Normal/High 
C                KMILES Mileage 
C Output Arguments 
C                hnbefusr Calculated emissions. 
C                exists True if relevant record found and emisdsions calculated 
C 
C Module Variables (see declaration above) 
C 
C       Name 
C       ---- 
C       sr_p 
c       srp 
C       Local variable/array dictionary 
C 
C       Name    Type            Description 
C       ----------------------------------- 
CC      my      Integer         Model year 
C       ip      Integer         Pollutant index 
C       iv      Integer         Vehicle index 
C       ifds    Integer         Fuel Delivery index 
C       isr     Integer         Start/running 
C       inh     Integer         Normal/High 
C       key     Integer         Used to distiguish case as function 
C                               of isr,inh,ip 
C       kmiles  Real            Mileage 
C       hnbefusr Real           Function calculates emissions for user input 
C       exists  Logical         Relevant user record is present. 
C       Inrange Logical         The user myf and myl satisfy model year 
 
       IMPLICIT NONE 
       INTEGER MY,IP,IV,IFDS,ISR,INH 
       REAL KMILES,HNBEFUSR 
       LOGICAL EXISTS,INRANGE 
       EXISTS=.FALSE. 
C 
C If sr_p is not allocated there are no user input records so return 
C with exists= false 
C 
        IF(.NOT.ALLOCATED(SR_P))RETURN 
        SRP=>SR_P(IP,IV,IFDS,ISR,INH)%HEAD 
                DO WHILE(ASSOCIATED(SRP)) 
C 
C There is a record so test the model year for applicability. 
C 
                    INRANGE=MY.GE.SRP%MYF .AND. MY.LE.SRP%MYL 
                    IF (INRANGE) THEN 
C 
C The record is for current model year so calculate emissions and return 
C 
                        HNBEFUSR=HNUSER(KMILES,SRP%ZML,SRP%DR, 
     *                          INH,ISR,IP) 
                        EXISTS=.TRUE. 
                        RETURN 
                    ELSE 
                        EXISTS=.FALSE. 
                    END IF 
C 
C Get the pointer to the next record and if it is associated continue search. 
C 
                    SRP=>SRP%NEXT 
                END DO 
100    FORMAT(I5,5I2,2I5,2(1X,F8.4)) 
       END SUBROUTINE 
 
       SUBROUTINE CHK_AVG(IFDS,IV,IP,MY,USER,EXISTS) 
C 
C chk_avg determines if a user supplied record matching the model year, fuel 
C delivery system, vehicle, and pollutant. If such a record is available the 
C subroutine returns a pointer to the record to the calling program 
C 
C Changes 
C    05 Aug 99 @Dyntel 2-660. Cosmetic changes. Comments etc 
C 
C Called by BEREF 
C 
C INPUT ARGUMENTS 
C                MY Model Year 
C                IP Pollutant 
C                IV Vehicle 
c                IFDS Fuel Delivery System 
C Output Arguments 
C                exists True if relevant record found and emisdsions calculated/ 
C                user; Pointer to record to be used for calculation. 
C 
C Module variables 
C 
C  AVG_P 
C 
C       Local variable/array dictionary 
C 
C       Name    Type            Description 
C       ----------------------------------- 
C       my      Integer         Model year 
C       ip      Integer         Pollutant index 
C       iv      Integer         Vehicle index 
C       ifds    Integer         Fuel Delivery index 
C       user    Pointer         Points to user supplied record. 
C       exists  Logical         Set true if user supplied record exists. 
C       inrange Logical         Set true if model year satisfied by user 
c                                supplied record 
      IMPLICIT NONE 
      INTEGER IFDS,IV,IP,MY 
      LOGICAL EXISTS,INRANGE 
      TYPE (AVG_REC),POINTER :: USER 
 
      EXISTS=.FALSE. 
C 
C If avg_p is not allocated there are no user records so return 
C 
      IF (.NOT.ALLOCATED (AVG_P)) THEN 
       RETURN 
      END IF 
      USER=>AVG_P(IP,IV,IFDS)%HEAD 
C 
C There may be a suitable user record so test model year 
C 
       DO WHILE (ASSOCIATED(USER)) 
           INRANGE=MY.GE.USER%MYF .AND. MY.LE.USER%MYL 
           IF (INRANGE) THEN 
C 
C Found a user record return to caller and use it 
C 
               EXISTS=.TRUE. 
               RETURN 
           ELSE 
              EXISTS=.FALSE. 
           END IF 
C 
C Follow link to next record 
C 
           USER=>USER%NEXT 
       END DO 
 
      END  SUBROUTINE 
 
 
      REAL FUNCTION HNUSER(KMILES,ZML,DR,INH,ISR,IP) 
C 
C HNUSER calculates the emissions using a user supplied record 
C 
C Changes 
C  03 AUG 99 @Dyntel NH. Changed case construct to if then else. 
C Called by chk_user 
C 
C INPUT ARGUMENTS 
C               KMILES mileage 
C               ZML Zero Mile 
C               DR Deterioration Factor 
C               inh =1 if NORMAL =2 if HIGH 
C               isr =1 if RUN =2 if START 
C 
C OUTPUT 
C       FUNCTION value 
C 
      IMPLICIT NONE 
      REAL KMILES,ZML,DR 
      INTEGER INH,ISR,IP 
      INTEGER NOX,RUN,NORMAL 
      PARAMETER (NOX=3,RUN=1,NORMAL=1) 
 
      SELECT CASE (IP) 
          CASE (NOX) 
            IF (INH.EQ.NORMAL .AND. ISR.EQ.RUN) THEN 
                HNUSER=ZML+DR*KMILES 
            ELSE IF (ISR.EQ.RUN) THEN 
                HNUSER=ZML 
            ELSE 
                HNUSER=ZML+DR*KMILES 
            END IF 
          CASE DEFAULT 
            IF (INH.EQ.NORMAL .AND. ISR.EQ.RUN) THEN 
                HNUSER=ZML+DR*KMILES 
            ELSE IF (INH.EQ.NORMAL) THEN 
                HNUSER=ZML+DR*KMILES 
            ELSE IF (ISR.EQ.RUN) THEN 
                HNUSER=ZML 
            ELSE 
                HNUSER=ZML 
            END IF 
      END SELECT 
      END FUNCTION 
 
 
       SUBROUTINE DEALLOC_8193 
C 
C deallocates storage for user supplied 81 93 emission records. 
C 
C Called by REINIT 
C 
C CALLS none 
C 
C Changes 
C 03 Aug 99 @Dyntel NH 2-660. Changed to deallocate linked lists. 
C 
C 27 Jul 99 @Dyntel NH 2-660. Added to reinitialize after using user supplied 
c           81 83 emission records 
C 
C       Local variable/array dictionary 
C 
C       Name    Type            Description 
C       ----------------------------------- 
C       ip      Integer         Pollutant index 
C       iv      Integer         Vehicle index 
C       ifds    Integer         Fuel Delivery index 
C       isr     Integer         Start/running 
C       inh     Integer         Normal/High 
C       SRPX    type(sr_rec),pointer  POINTER to a START/RUNNING record 
C       AVGX    type(avg_rec),pointer POINTER to an AVERAGE RUNNING record 
 
       TYPE (SR_REC),POINTER :: SRPX 
       TYPE (AVG_REC),POINTER :: AVGX 
       INTEGER DEALLOC 
       INTEGER INH 
       INTEGER ISR 
       INTEGER IFDS 
       INTEGER IP 
       INTEGER IV 
 
       IF (ALLOCATED(SR_P)) THEN 
 
       DO INH=1,2 
         DO ISR=1,2 
            DO IFDS=1,3 
               DO IV=1,5 
                  DO IP=1,3 
                      SRP =SR_P(IP,IV,IFDS,ISR,INH)%HEAD 
                      DO WHILE(ASSOCIATED(SRP)) 
                        SRPX=>SRP 
                        SRP=>SRP%NEXT 
                        DEALLOCATE(SRPX) 
                      END DO 
                  END DO 
               END DO 
            END DO 
         END DO 
       END DO 
       DEALLOCATE (SR_P,STAT=DEALLOC) 
       IP=DEALLOC!SIMPLY TO AVOID WARNING. 
       END IF 
 
       IF (ALLOCATED(AVG_P)) THEN 
 
       DO IFDS=1,3 
        DO IV=1,2 
          DO IP=1,3 
               AVGP =AVG_P(IP,IV,IFDS)%HEAD 
               DO WHILE(ASSOCIATED(AVGP)) 
                  AVGX=>AVGP 
                  AVGP=>AVGP%NEXT 
                  DEALLOCATE(AVGX) 
               END DO 
          END DO 
        END DO 
       END DO 
       DEALLOCATE (AVG_P,STAT=DEALLOC) 
       END IF 
       RETURN 
       END SUBROUTINE 
 
       SUBROUTINE INIT_SR_AVG 
       IMPLICIT NONE 
C 
C Initializes data structures for storage of user supplied 81 93 emission records 
C 
       INTEGER INH,ISR,IFDS,IV,IP 
 
       ALLOCATE(SR_P(3,5,3,2,2)) !ALLOCATE STORAGE FOR NORMAL START POINTERS. 
       ALLOCATE(AVG_P(3,2,3))  !ALLOCATE STORAGE FOR AVERAGE RUNNING RECORDS. 
       DO INH=1,2 
         DO ISR=1,2 
            DO IFDS=1,3 
               DO IV=1,5 
                  DO IP=1,3 
                      NULLIFY (SR_P(IP,IV,IFDS,ISR,INH)%HEAD) 
                      NULLIFY (SR_P(IP,IV,IFDS,ISR,INH)%TAIL) 
                  END DO 
               END DO 
            END DO 
         END DO 
       END DO 
 
       DO IFDS=1,3 
        DO IV=1,2 
          DO IP=1,3 
            NULLIFY(AVG_P(IP,IV,IFDS)%HEAD) 
            NULLIFY(AVG_P(IP,IV,IFDS)%TAIL) 
          END DO 
        END DO 
       END DO 
       END SUBROUTINE 
 
      SUBROUTINE USER_FILE(INPREC,INERR,RC) 
C 
C user_files checks for the existence of the user file named in the run 
C input file 
C 
      USE LABELS 
      IMPLICIT NONE 
 
      INCLUDE 'IOUCOM.I' 
C  
      INTEGER, EXTERNAL :: JUNIT    
 
      CHARACTER   INPREC*(*) 
      INTEGER     INERR 
      INTEGER     RC 
      INTEGER     BER 
      CHARACTER   LABEL*19 
      INTEGER     LABTYPE 
      INTEGER     SECTION 
      CHARACTER   M6REC*150 
      CHARACTER   STDLAB*40 
      INTEGER     SOURCE 
       
ccs   Add declarations of M6DIR and M6FILE 
      CHARACTER(LEN=200) M6DIR 
      CHARACTER(LEN=280) M6FILE 
C 
C     FIND OUT WHETHER WE ARE READING FROM AN EXTERNAL DATA FILE. 
C                  
      IF (INPREC(20:20).EQ.'@') THEN 
        SOURCE = IOUGEN 
        M6REC = INPREC 
      ELSE 
C 
C       READ THE DATA FROM AN EXTERNAL DATA FILE. 
C 
        SOURCE = JUNIT() 
C 
C       FIND THE TOKEN NUMBER OF THE RECORD LABEL THAT MUST APPEAR 
C       ON A RECORD SOMEWHERE IN THE HEADER OF THE EXTERNAL DATA FILE. 
C 
        CALL CHKLAB('81-93 LDG EFS',STDLAB,BER,SECTION) 
C 
C       OPEN THE EXTERNAL DATA FILE. 
C 
ccs        CALL UCCOMP(INPREC(21:), FILENAME) 
        CALL LJCOMP(INPREC(21:), FILENAME) 
 
ccs     All M6 input files live in SMK_M6PATH 
        CALL GETENV( "SMK_M6PATH", M6DIR ) 
        M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME 
         
        OPEN(SOURCE,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     *       IOSTAT=RC,ERR=80) 
        WRITE(IOUREP,120) FILENAME 
  120   FORMAT(/'* Reading non-default BASIC EMISSION RATES from ', 
     &          'the following external',/'* data file: ',A) 
C 
C       READ THE REQUIRED HEADER LINE FROM THE FILE. 
C 
        CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
        IF (RC.NE.0) GOTO 85 
C 
C       THE NEXT TWO ASSIGNMENT STATEMENTS ARE PRESENT SOLELY TO SUPPRESS 
C       COMPILER GENERATED WARNINGS. THEY HAVE NO FUNCTIONAL IMPORTANCE. 
C 
        LABEL = STDLAB(1:1) 
        LABTYPE = SECTION 
C 
C       CHECK FOR A VALID LABELED INPUT RECORD. 
C 
        CALL UCCOMP(M6REC(1:19),LABEL) 
        CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
        IF (LABTYPE.NE.BER) GOTO 90 
      ENDIF 
      CLOSE(SOURCE) 
      RETURN 
 
C 
C     ERROR OPENING THE EXTERNAL DATA FILE. 
C 
   80 WRITE(IOUERR,180) FILENAME 
  180 FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
      INERR = INERR + 1 
      GOTO 99 
C 
C     ERROR READING THE HEADER OF THE EXTERNAL DATA FILE. 
C 
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,185) FILENAME 
  185     FORMAT(/'*** ERROR: no header record was found in ', 
     &            'the following external data file:',/'***',A) 
          RC = 2 
          INERR = INERR + 1 
        CASE (1:) 
          WRITE(IOUERR,186) FILENAME 
  186     FORMAT(/'*** ERROR reading the header of the ', 
     &            'following external data file:',/'*** ',A) 
        CASE DEFAULT 
      END SELECT 
      GOTO 99 
C 
C     INVALID LABEL ON THE HEADER RECORD. 
C 
   90 WRITE(IOUERR,190) M6REC 
  190 FORMAT(/'*** ERROR: the following record is not allowed in ', 
     &        'the header of the external data file:', 
     &      /'*** ',A) 
      RC = 3 
      INERR = INERR + 1 
C 
   99 RETURN 
C 
       END SUBROUTINE 
 
      END MODULE 
