C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETTTC(INERR,*)
C
C  GETTTC Form TTCDAT array that contains Technical Training and
C         Certification Benefits using default values or values from TTC.D
C         file. 
C
C  Called by ONESEC.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /TTCINF/ TTCDEF,ITTCFG
C    /IOUCOM/  IOUTTC              
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C     /TTCINF/ TTCDAT,TTCSML
C
C  Notes:
C
C  April-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Subroutine was developed in accordance with this task request.  
C 
C 
      INCLUDE 'TTCINF.I'
      INCLUDE 'IOUCOM.I'
c
cmh Add definition of FILNAME, ENVL
      CHARACTER*132 FILNAME, ENVL
C
      TTCSML=' '
      IF ( ITTCFG .LT. 5) THEN
C
C Restore the default values       
C
        DO 10 I=1,18
        DO 10 J=1,3
        DO 10 K=1,24
	   TTCDAT(K,J,I)=TTCDEF(K,J,I)
  10    CONTINUE
      ELSE
C
C  Read data file TTC.D and load values into real array TTCDAT
C
cmh Add path to input file name
c
         CALL GETENV("EMS_LOC",ENVL)
         FILNAME = ENVL(1:LNBLNK(ENVL))//'/TTC.D'
         OPEN (IOUPFF, FILE=FILNAME, ERR=98)
cmh      OPEN (IOUTTC, FILE='TTC.D',ERR=98)
c
         READ (IOUTTC,101,ERR=98) TTCSML
 101     FORMAT (A80)
         DO 20 I = 1,18
         DO 20 J = 1,3
           READ (IOUTTC,100,ERR=98) (TTCDAT(K,J,I), K = 1,24)
 100       FORMAT (24F5.3)
  20     CONTINUE
         CLOSE(IOUTTC)
      ENDIF
      RETURN	 
C
  98  CALL QUITER(0.,0,163,INERR)
  99  RETURN1
      END
