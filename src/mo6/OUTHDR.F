      SUBROUTINE OUTHDR(IOUOUT,MIN,NFILE,NRUN) 
C 
C  OUTHDR writes the output header on unit IOUOUT.  
C 
C  Called by MOBILE. 
C 
C  Changes: (Last change first) 
C 
C  09-Dec-1998 @DynTel-ZK Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: IOUOUT,NFILE,NRUN 
C    common block: /VDATA / VERSION 
C 
C  Notes: 
C  12-Dec-1994 @ CSC-rwk Req. 2-000 Include file VDATA was added to make it 
C  easier to change the name of the version of the model.  This was done for quality 
C  control. 
C  15-Oct-1997 @ DynTel-gjr 2-000 Changed VERSON to VERSION. 
C   3 Jun 98 @DynTel-MLA 2-659 In Mobile6, the header does not contain a 
C            PROJID. Users may supply their own titles by echoing comment 
C            lines that are placed at the beginning of their input files.  
C            Additional comments may appear anywhere in the input file. 
C 
      IMPLICIT NONE 
      INCLUDE 'VDATA.I' 
      INTEGER IOUOUT  
      INTEGER   NFILE 
      INTEGER   NRUN 
C 
      CHARACTER BUFFER1*90 
      CHARACTER BUFFER2*70 
      CHARACTER FNAME*40 
      CHARACTER MIN*(*) 
C 
CC    WRITE(BUFFER1,120) NFILE,NRUN 
CC120 FORMAT('Processing input for file ',I9,', run ',I9,'.') 
CC    CALL LJCOMP(BUFFER1,BUFFER2) 
C      

ccs      CALL UCCOMP(MIN,FNAME) 
      CALL LJCOMP(MIN,FNAME)

      WRITE(BUFFER1,120) FNAME,NFILE,NRUN 
  120 FORMAT('Input file: ',A40,' (file ',I9,', run ',I9,').') 
      CALL LJCOMP(BUFFER1,BUFFER2) 
C 
      WRITE(IOUOUT,200) VERSION, BUFFER2 
  200 FORMAT('**************************************************', 
     &       '*************************', 
     &    2(/'* ',A70,T75,'*'), 
     &      /'**************************************************', 
     &       '*************************') 
CC200 FORMAT('** ',T4,A80) 
C 
      RETURN 
      END 
