      SUBROUTINE CEVMPD 
C 
C  CEVMPD calculates miles per day by the type and relative (to ICY) age of 
C  the vehicle.  The formula uses the REGMOD computed TFMYM (January 1 
C  mileage accumulation rates).  As with TYMYM and BMYTPD, one set 
C  of rates is calculated for use on all scenarios of the run. 
C 
C  Called by REGMOD. 
C 
C  Changes : 
C 
C   7 Sep 99 @DynTel-ddj 2-710 Added option that sets the miles per day 
C            array to the FTP value of 7.5 when FTP MILES PER DAY flag 
C            is not equal to 0. 
C  18 May 99 @DynTel-HQ 2-650 removed TPDCOM.I 
C  26 Aug 98 @DynTel-ZK Explicit Typing 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    common blocks: 
C    /MAXIMA/ MAXYRS 
C    /MYRCAL/ TFMYM 
C 
C  Output on return: 
C 
C    common blocks: 
C    /CEVBMY/ BMYMPD 
C 
C  Notes: 
C 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6. 
      IMPLICIT NONE  
      INCLUDE 'CEVBMY.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRCAL.I' 
       
      INTEGER IV, JDX 
C 
      DO 10 IV = 1, MAXVEH 
        IF(FTP_INPUT_FLG.EQ.0) THEN 
          DO JDX = 1, MAXYRS 
            BMYMPD(JDX,IV) = (TFMYM(JDX,IV)*100000)/365. 
          END DO 
        ELSE 
C 
C     FTP input option, miles/day set equal to 7.5 miles. 
C 
          DO JDX = 1, MAXYRS 
            BMYMPD(JDX,IV) = 7.5 
          END DO 
        END IF 
   10 CONTINUE 
C 
      RETURN 
      END 
