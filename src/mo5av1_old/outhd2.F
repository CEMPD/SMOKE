C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTHD2
C
C  OUTHD2 writes the run header on unit IOUREP using output format (OUTFMT) 2.
C  The format is intended for postprocessing of the data.  Hence, the
C  header consists of the field descriptors of the data lines.  Each scenario
C  generates 4 identically formatted lines, 1 each for evap HC, total HC,
C  CO and NOx.  The evap HC line prints blanks for the AllVeh C.E.F.
C
C  Called by OUTPUT.
C
C  Input on call:
C
C    common blocks:
C    /IOUCOM/ IOUREP
C    /MAXIMA/ MAXVEH
C
C  Output on return:
C
C    common block: /REGION/ INITPR
C
C  Local array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  NAMVEH   C*8  names of vehicle classes for field headings
C
C  Notes:
C
C  OUTFMT 2 contains the same set of fields as MOBILE3's OUTFMT 1, except that
C
C  (1) The speed vector is omitted.
C
C  (2) The emission factors and vmt mix fractions are printed for all vehicle
C      classes are printed.
C
C  OUTFMT 2 is the "short" numeric format.  With a record length of 130
C  columns, each record will fit entirely onto an output page without
C  truncation.
C
C  OUTHD2 was added for MOBILE4.
C
C
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
C
      CHARACTER*8 NAMVEH(9)
C
      DATA NAMVEH/
     *   'LDGV    ','LDGT1   ','LDGT2   ','LDGT    ','HDGV    ',
     *   'LDDV    ','LDDT    ','HDDV    ',' MC     '/
C
      INITPR=0
C
      MAXV=MAXVEH+1
      WRITE(IOUREP,200)
     *   NAMVEH,(NAMVEH(IV),IV=1,3), (NAMVEH(IV),IV=5,MAXV)
  200 FORMAT(
     *'0','R',27X,'P'/
     *' ','e    Amb.',19X,'o',22X,'Composite Emission Factors',38X,
     *                                               'Vehicle Mix'/
     *' ','g CY Tmp   Cold/Hot Start   l',9(2X,A5),' AllVeh',1X,8A5/
     *' ','- -- ---- ----------------- - ',69('-'),1X,39('-'))
C
      RETURN
      END
