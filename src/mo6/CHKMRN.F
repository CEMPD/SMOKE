      SUBROUTINE CHKMRN(INERR,RC) 
C 
C  CHKMRN checks for consistency between user-supplied mileage accrumulation 
C  rates and registration distributions. Whether or not the distributions are 
C  replaced, they are normalized. 
C 
C  Called by CHKINP 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  12 Apr 00 @EPA-djb Add check for HDDV mileage inputs and adapt to new 
C            array size for JULMYR. 
C  24 Jan 00 @DynTel-MLA  1-011  Changed VTLDDT to VTLDDT12 because 
C            LDDT was split into 2 categories and names have changed. 
C  31 Aug 98 @DynTel-ZK 2-000 Explicit typing  
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: INERR 
C    common blocks: 
C    /FLAGS2/ MYMRFG 
C    /MAXIMA/ MAXYRS 
C    /VVNAMS/ VVHDDV  (IVTYPE.I) 
C    /VVTYPS/ VTLDDT12, VTLDDV, VTLDGT1, VTLDGV  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /MYRSAV/ AMAR,JULMYR 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  SUMMYR   R    sum of myr for IV; used to normalize user entered JULMYR 
C 
C  Notes: 
C 
C  This subroutine is new for Mobile6. 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRSAV.I' 
C 
      INTEGER  RC, IV, INERR, JDX 
      REAL  SUMMYR 
C 
      RC = 0 
C 
C  Check for user supplied my mileage accumulation rates or registration 
C  distributions. 
C 
      IF(MYMRFG.EQ.1) GOTO 65 
C      
   45 DO 55 IV=1,MAXVEH 
        DO 50 JDX=1,MAXYRS 
          IF(AMAR(JDX,IV).LE.0.0) 
     *      CALL QUITER(AMAR(JDX,IV),IV,21,INERR) 
   50   CONTINUE 
   55 CONTINUE 
C 
C  Normalize JULMYR. 
C 
   65 DO 80 IV=1,16 
        SUMMYR=0.0 
        DO 70 JDX=1,MAXYRS 
          SUMMYR=JULMYR(JDX,IV)+SUMMYR 
   70   CONTINUE 
        IF(SUMMYR.LE.0.00001) THEN 
          CALL QUITER(SUMMYR,IV,22,INERR) 
          SUMMYR=1. 
        ENDIF 
        DO 75 JDX=1,MAXYRS 
          JULMYR(JDX,IV)=JULMYR(JDX,IV)/SUMMYR 
   75   CONTINUE 
   80 CONTINUE 
C 
   99 RETURN 
      END 
