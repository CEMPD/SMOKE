C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTHD1
C
C  OUTHD1 writes the run header on unit IOUREP using output format (OUTFMT) 1.
C  The format is intended for postprocessing of the data.  Hence, the
C  header consists of the field descriptors of the data lines.  Each scenario
C  generates 4 identically formatted lines, 1 each for evap HC, total HC,
C  CO and NOx.  The evap HC line omits the last field (All Veh C.E.F.).
C
C  Called by OUTPUT.
C
C  Input on call:
C
C    common blocks:
C    /CITRV1/ RVPBAS,RVPIUS,IUSESY
C    /FLAGS2/ IMFLAG
C    /FLAGS3/ ATPFLG,OUTFMT
C    /IOUCOM/ IOUREP
C    /MAXIMA/ MAXVEH
C    /OXY1/ SHRMKT,OXYCNT,IGASHW
C    /RFORM1/ RFGFLG
C    /STRING/ NOYES
C    /TEMPS/ TEMMIN,TEMMAX
C
C  Output on return:
C
C    common block: /REGION/ INITPR
C
C  Local array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IMF      I    indicator for I/M program no=1 yes=2
C  NAMVEH   C*8  names of vehicle classes for field headings
C  NAMSPD  C*14  names of vehicle classes for field headings (speeds)
C  NAMCEF  C*11  names of vehicle classes for field headings (EF's)
C  NAMMIX  C*12  names of vehicle classes for field headings (Mix)
C
C
C  OUTFMT 1 contains the same set of fields as MOBILE3's OUTFMT 2, except that
C
C  (1) The fields are all on 1 record, instead of split across 3 records,
C      for a given scenario and pollutant.  Field print order is the same.
C
C  (2) The emission factors are printed with 3, instead of 2, decimal places.
C
C  (3) The scenario id label (new for MOBILE4) is appended as the last field.
C
C  OUTFMT 1 is the "long" numeric format, with a record length of 212.
C  OUTFMT 6 was added for MOBILE5
C  Sep-13-1994 @ CSC-bsg Request 2-000, 2-392 added include files
C
C
      INCLUDE 'CITRV1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'OXY1.I'
      INCLUDE 'REGION.I'
      INCLUDE 'RFORM1.I'
      INCLUDE 'STRING.I'
      INCLUDE 'TEMPS.I'
C
      CHARACTER*8  NAMVEH(9)
      CHARACTER*14 NAMSPD(9)
      CHARACTER*11 NAMCEF(9)
      CHARACTER*12 NAMMIX(9)
C
C
      DATA NAMVEH/
     *   'LDGV    ', 'LDGT1   ', 'LDGT2   ', 'LDGT    ', 'HDGV    ',
     *   'LDDV    ', 'LDDT    ', 'HDDV    ', ' MC     '/
C
      DATA NAMSPD/
     *   '"LDGV  speed",', '"LDGT1 speed",', '"LDGT2 speed",',
     *   '"LDGT  speed",', '"HDGV  speed",', '"LDDV  speed",',
     *   '"LDDT  speed",', '"HDDV  speed",', '"MC    speed",'/
C
      DATA NAMCEF/
     *   '"LDGV  ef",', '"LDGT1 ef",', '"LDGT2 ef",',
     *   '"LDGT  ef",', '"HDGV  ef",', '"LDDV  ef",',
     *   '"LDDT  ef",', '"HDDV  ef",', '"MC    ef",'/
C
      DATA NAMMIX/
     *   '"LDGV  mix",', '"LDGT1 mix",', '"LDGT2 mix",',
     *   '"LDGT  mix",', '"HDGV  mix",', '"LDDV  mix",',
     *   '"LDDT  mix",', '"HDDV  mix",', '"MC    mix",'/
C
      IMF=1
      IF(IMFLAG.NE.1) IMF=2
      INITPR=0
C
      MAXV=MAXVEH+1
      IF(OUTFMT.EQ.1)
     *   WRITE(IOUREP,200)
     *   NAMVEH,(NAMVEH(IV),IV=1,3), (NAMVEH(IV),IV=5,MAXV)
C
  200 FORMAT(
     *' ','R',73X,'P'/
     *' ','e',44X,'Amb.',20X,'Alt. o',28X,'Composite Emission Factors',
     *        45X,'Vehicle Mix'/
     *' ','g CY',13X,'Vehicle Speeds',14X,'Tmp   Cold/Hot Start  in Ft.'
     *    ,' l',9(3X,A5),' All Veh',8(1X,A5),' Scenario Title'/
     *' ','- -- ',39('-'),' ---  ',16('-'),' ------ - ',79('-'),1X,
     *    47('-'),1X,16('-'))
C
      IF(OUTFMT.EQ.6.AND.RFGFLG.EQ.1)
     *   WRITE(IOUREP,210) NOYES(IMF),NOYES(ATPFLG),NOYES(RFGFLG),
     *   TEMMIN,TEMMAX,RVPBAS,RVPIUS,IUSESY,SHRMKT,OXYCNT,NOYES(IGASHW),
     *   NAMCEF,(NAMSPD(IV),IV=1,3),(NAMSPD(IV),IV=5,MAXV),
     *  (NAMMIX(IV),IV=1,3),(NAMMIX(IV),IV=5,MAXV)
C
  210 FORMAT(/'"I/M Program:","",',A3,
     * ',"","Anti-tam Program:","",',A3,
     * ',"","Reformulated Gas:","",',A3/
     * '"Minimum Temp:","",',F5.1,',"","Maximum Temp:","",',F5.1,
     * ',"","Period 1 RVP:","",',F5.1,',"","Period 2 RVP:","",',F5.1,
     * ',"","Period 2 start year:","","",',I4/
     * '"Gasoline Market Share:","","",',F5.3,
     * ',"","Ether Blend Market Share:","","",',F5.3,
     * ',"","Alcohol Blend Market Share:","","",',F5.3/
     * 5('"",'),'"Ether Blend Oxygen Content:","","",',F5.3,
     * ',"","Alcohol Blend Oxygen Content:","","",',F5.3,
     * ',"","Alcohol Blend RVP Waiver:","","",',A3/
     * '" ",', '" ",', '"Composite Emission Factors",', 9('"",'),
     * '"Operating Percentages",', 3('"",'),
     * '"Vehicle Speed",', 9('"",'), '"Vehicle Mix"'/
     * '"Pollutant",', '"Cal. Year",', 9(A11), '"All Veh",',
     * '"Amb. Temp.",', '"Cold",', '"Hot",', '"Stabilized",',
     * '"Alt. (Ft)",', '"Region",', 8(A14),
     * 8(A12), '"Scenario Title"'/)
C
      IF(OUTFMT.EQ.6.AND.RFGFLG.EQ.2)
     *   WRITE(IOUREP,220) NOYES(IMF),NOYES(ATPFLG),NOYES(RFGFLG),
     *   TEMMIN,TEMMAX,RVPBAS,RVPIUS,IUSESY,
     *   NAMCEF,(NAMSPD(IV),IV=1,3),(NAMSPD(IV),IV=5,MAXV),
     *  (NAMMIX(IV),IV=1,3),(NAMMIX(IV),IV=5,MAXV)
C
  220 FORMAT(/'"I/M Program:","",',A3,
     * ',"","Anti-tam Program:","",',A3,
     * ',"","Reformulated Gas:","",',A3/
     * '"Minimum Temp:","",',F5.1,',"","Maximum Temp:","",',F5.1,
     * ',"","Period 1 RVP:","",',F5.1,',"","Period 2 RVP:","",',F5.1,
     * ',"","Period 2 start year:","","",',I4/
     * '" ",', '" ",', '"Composite Emission Factors",', 9('"",'),
     * '"Operating Percentages",', 3('"",'),
     * '"Vehicle Speed",', 9('"",'), '"Vehicle Mix"'/
     * '"Pollutant",', '"Cal. Year",', 9(A11), '"All Veh",',
     * '"Amb. Temp.",', '"Cold",', '"Hot",', '"Stabilized",',
     * '"Alt. (Ft)",', '"Region",', 8(A14),
     * 8(A12), '"Scenario Title"'/)
      RETURN
      END
