C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTDTG
C
C  OUTDTG writes a message that indicates that the effects of the
C  detergent gasoline additive rulemaking have been altered by the
C  user. The values provided by the user are printed as well.
C
C  Called by OUTHD3 and OUTHD4
C
C  Calls None
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C
C    /DTGGAS/ DGPISY,DGPIEY,DGPHSF,FIDFRC,FIDEFF,IVDFRC,IVDEFF
C    /IOUCOM/ IOUREP
C
C
C  Notes:
C
C  September 15 1994, @CSC-pme Request 2-464. Subroutine created for
C  Detergent Gas User Output in Mobile5b.
C  07-June-1995 @ CSC-wcs Request # 0-448 Mobile5b testing (fix echoed
C    parameters so that their order matches that of the labels above)
C  13-Sep-1996 @DynTel-yc Request #2-630 The detergent gas phase-in effect
C               was updated by refining phase-in period
C
      INCLUDE 'DTGGAS.I'
      INCLUDE 'IOUCOM.I'
C
C     Write message.
C
      WRITE(IOUREP,200)
 200  FORMAT(/,' ','Effects of the Detergent Gasoline Additive ',
     *           'Rulemaking have been altered by user:',/)
C
C     Write parameters' header.
C
      WRITE(IOUREP,210)
C
 210  FORMAT(' ','PIEY  Pct    FI %   FIHC   FICO   FINO  ',
     *                       ' IV %   IVHC   IVCO   IVNO',/,
     *       ' ','---- ------ ------ ------ ------ ------ ',
     *                       '------ ------ ------ ------')
C
C         Write parameter values.
C
C  07-June-1995 @ CSC-wcs Request # 0-448 Mobile5b testing (fix echoed
C    parameters so that their order matches that of the labels above)
C
      WRITE(IOUREP,220) DGPIEY,DGPHSF*100,FIDFRC*100,
     &                 FIDEFF(1)*100,FIDEFF(2)*100,FIDEFF(3)*100,
     &                 IVDFRC*100,
     &                 IVDEFF(1)*100,IVDEFF(2)*100,IVDEFF(3)*100
C
 220  FORMAT(' ',I4,9(1X,F6.2),/)
C
      RETURN
      END
