C 
C  IMPAR5 contains several I/M program variables. 
C 
C  Changes: (Last change first) 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  ATPUEF    R    User-supplied or default value for the effectiveness 
C                 of test-and-repair ATP programs. 
C  CRHDGV    R    HDGV I/M credit adj. for annual, bienn. HC, CO, NOx 
C  DCNTSV    R    Holds the program default values for DISCNT. 
C  DISCNT    R    The discounts to be applied to the I/M program credits 
C                 if the I/M program is test-and-repair (not test-only). 
C  LBIM4P    I    lower bound for Tech4+ model year group. 
C  MYGIM2    I    model year groups for Tech 1 and 2 (Tech3=2, for I/M) 
C  PPGUEF    R    User-supplied or default value for the effectiveness 
C                 of test-and-repair purge/pressure programs. 
C 
C  Notes: 
C 
C 1 Feb 96 @DynTel-MLA 2-611  Added IP dimension to DISCNT. Added new array 
C    DCNTSV with the same dimensions as DISCNT. This holds the default 
C    values for the I/M program discounts. It is used to reinitialize  
C    DISCNT at the beginning of each new run. 
C 23 Feb 96 @DynTel-MLA 2-611  Added new common variables ATPUEF and PPGUEF 
C    which hold the user-supplied values for the effectiveness of the ATP 
C    and the Purge/Pressure programs. These values are read from the 
C    Program Effectiveness Record, which appears in the one-time data 
C    section if a flag is set on the I/M Program Descriptive Record. 
C  3 Dec 96 DynTel@MLA Expanded array MYGIM2 to accomodate the new LDGT 
C           vehicle categories. 
C 
C  The data is initialized in BD_IM2(MYGIM2,LBIM4P), GETEFR.FOR(ATPUEF,PPGUEF), 
C 
      INTEGER MYGIM2,LBIM4P 
      REAL CRHDGV,DISCNT,DCNTSV,ATPUEF,PPGUEF 
      COMMON /IMPAR5/ MYGIM2(3,5),LBIM4P,CRHDGV(3,2), 
     &                DISCNT(3,3),DCNTSV(3,3),ATPUEF,PPGUEF 
