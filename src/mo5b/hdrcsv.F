C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE HDRCSV
C
C  HDRCSV makes the header for ouptut into a CSV format and in
C  one single record.
C
C  Called by GETHOR, and GETVER.
C
C  Local array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  CYLBL   C*12  Calendar year header
C  DPRLBL  C*31  Right half of program header
C  EVMLBL  C*13  Evaluation month header
C  POLBL   C*12  Pollution header
C  PRJLBL  C*25  Project identification header
C  VHLBL   C*12  Vehicle header
C
C
C  OUTFMT 6 causes the same header to be written once.  This was
C  created Feb 24, 1993.
C  18-June-1993 @ ARC-bsg Subtask 222 CSV version of In-House MOBILE5a
C
      CHARACTER*25 PRJLBL
      CHARACTER*13 EVMLBL
      CHARACTER*12 POLBL,CYLBL,VHLBL
      CHARACTER*31 DPRLBL
C
      INCLUDE 'IOUCOM.I'
C
C      1234567890123456789012345678901234567890
      DATA PRJLBL/
     *'"Project Identification",'/
      DATA EVMLBL/
     *'"Eval Month",'/
      DATA POLBL/
     *'"Pollutant",'/
      DATA CYLBL/
     *'"Cal. Year",'/
      DATA VHLBL/
     *'"Vehicle  ",'/
      DATA DPRLBL/
     *' Descriptive Parameter Record",'/
      WRITE(IOUREP,200) PRJLBL,EVMLBL,CYLBL,POLBL,VHLBL,DPRLBL,DPRLBL,
     *DPRLBL,DPRLBL,DPRLBL
  200 FORMAT(A25,A13,3A12,'"Emission Factor",',
     *'"VMT","Speed","Temperatures","RVP info","Op Mode",',
     *'"I/M Program #1',A30,'"I/M Program #2',A30,
     *'"Anti-Tampering Program',A30,
     *'"Pressure Functionality Check Program',A30,
     *'"Purge Functionality Check Program',A30,
     *'"RFG","Oxy Record","Region","Tier 1(I/M)","LEV (I/M)",',
     *'"ZEV%","Last Field","Evap Phase-in Schedule","Stage II VRS",',
     *'"Onboard VRS","Trip Distribution","Scenario"')
      RETURN
      END
