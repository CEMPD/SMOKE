C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTDT5(ICY)
C
C  OUTDT5 writes scenario emission factor results by model year for selected
C  vehicle types on unit IOUREP using output format (OUTFMT) 5.
C
C  Called by OUTPUT.
C
C  Calls OUTBY1, OUTBY2, OUTBY3, AND OUTBY4.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /EVAL/   MEVAL,MONTH
C    /FLAGS4/ PRTFLG
C    /IOUCOM/ IOUREP
C
C  Notes:
C
C  (1) OUTDT5 was added for MOBILE4 In-House Version 02 (By Model Year Tables).
C
C  (2) OUTDT5 prints by model year results only for those pollutants picked
C      by the user via PRTFLG.
C
C  (3) OUTDT5 is called after the OUTFMT=3 output (scenario header + output
C      in 112 column report format) is printed.  Starting on a new page, two
C      By Model Year tables will now be printed on each page thereafter, until
C      all picked vehicle types' results are written.  The maximum space
C      required by OUTFMT = 5 for a scenario is thus 5 pages (OUTFMT = 3
C      output plus eight By Model Year tables).
C
C  (4) Diesel classes (LDDV, LDDT, HDDV) do not have tampering offsets,
C      evaporative HC, running loss HC or refueling loss HC.  Motorcycles do
C      not have any of these cases either, except for evaporative HC.  The
C      fields for these pollutants for these vehicle classes are eliminated
C      in the By Model Year Table(s).
C
C
      INCLUDE 'EVAL.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'IOUCOM.I'
C
C
      IF(MEVAL.NE.1) THEN
      WRITE(IOUREP,5) MONTH(MEVAL)
   5  FORMAT('1','By-Model-Year output is not available for ',A4,
     *' 1st runs.')
      RETURN
      ENDIF
C
      GOTO(10,20,30,40),PRTFLG
C
C  HC requested.
C
   10 CALL OUTBY1(ICY)
      RETURN
C
C  CO requested.
C
   20 CALL OUTBY2(ICY)
      RETURN
C
C  NOX requested.
C
   30 CALL OUTBY3(ICY)
      RETURN
C
C  HC, CO & NOX requested.
C
   40 CALL OUTBY4(ICY)
C
      RETURN
      END
