C Version @(#)$Id$ $Source$ $Date$ 

c  MOBILE5a.01 released March 29, 1993
c
C  MOBILE5a is the fifth version of the FORTRAN program implementation of the
C  Mobile Source Emission Model. It updates and replaces MOBILE4.1. The model
C  estimates HC, CO and NOx exhaust and HC evaporative emission factors by
C  motor vehicles.  For information on using the program and on the
C  differences between MOBILE4.1 and MOBILE5a, consult the User's Guide to
C  MOBILE4.1, EPA-AA-TEB-91-01, and the User's Guide to MOBILE5,
C  EPA-AA-TEB-92-01.  Otherwise questions should be directed to:
C
C         Office of Mobile Sources
C         Emission Planning and Strategies Division
C         2565 Plymouth Road
C         Ann Arbor, Mich. 48105
c
c---------------------------------------------------------------------------------------
C PLEASE NOTE THAT THIS VERSION OF MOBILE5  HAS BEEN SPECIFICALLY CUSTOMIZED
C TO OPERATE WITH THE LAKE MICHIGAN OZONE STUDY (LMOS) EMISSIONS MODELING
C SYSTEM (EMS).
C
c  This copy of MOBILE5 modified by Radian Corporation for LMOS MOVEM.
c  The program has been modified to calculate and report diurnal
c  evaporative emissions separately in grams/mile.  All modifications can
c  be located by searching for the string "LMOS Modification".  Key new
c  output variables are EFDNL (diurnal evaporative emission factor) and
c  EFEVND (total non-diurnal evaporative emission factor which includes
c  hot soak, crankcase, running losses, and resting losses).  The values for
c  these new variables are stored in a new common block named LMOS.
c
c  A new flag is added to the end of the scenario record set by the controlling
c  SAS program whether the scenario record is to generate diurnal emission factors
c  or non-diurnal emission factors.  This flag, ndi_di, is stored in a new common block
c  SCENE3.
c
c
c CHANGES ARE IDENTIFIED BY THE WORDS 'LMOS Modification',
c CHANGES WERE MADE TO THIS PROGRAM DURING THE MONTHS OF FEBRUARY, 1993 THROUGH
c APRIL, 1993.  THIS VERSION OF THE PROGRAM IS BASED ON THE MARCH 26,1993 RELEASE.
c NOTE THAT EARLIER VERSIONS OF MOBILE5A HAD MORE COMPLETE DOCUMENTATION OF VARIABLES
c IN THE HEADER.
c----------------------------------------------------------------------------------------
C
      PROGRAM DRIVER
C
C  DRIVER (as of 3/24/92):
C
C  DRIVER demonstrates how to call the subroutine version of
C  MOBILE5
C
C  Calls MOBILE.
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  BLANK    C    blank filename detector
C  INTER    C    flag for interactive mode ('Y') vs. file mode ('N')
C  M5CALL   I    loop index for MOBILE5 call loop
C
C
      COMMON /IOUCOM/ IOUIMD, IOUGEN, IOUREP, IOUERR, IOUASK
c
c  LMOS Modification 04/08/93
c
c  Add definition of ENVL and M5STAT variables
      CHARACTER*132 ENVL, M5STAT
c
c  LMOS Modification 04/08/93
c  Extend length of M5IN and M5OUT filenames
c
c      CHARACTER*40 M5IN, M5OUT, BLANK, BATIN
      CHARACTER*132 M5IN, M5OUT, BLANK, BATIN
      CHARACTER*1 INTER, BATCH
      DATA BLANK /'                                        '/
C
c LMOS Modification 04/08/93
c
c Comment out options for interactive or batch mode input
c Assume non-interactive input and one input file (i.e., non-batch)
c
c
cC Change for different computers
cC on the MAC IOUGEN = 9 and on the IBM IOUGEN = 5 and on MTS IOUGEN = 1
cC on the MAC IOUREP = 6 and on the IBM IOUREP = 6 and on MTS IOUREP = 7
cC on the MAC IOUASK = 9 and on the IBM IOUASK = 0 and on MTS IOUASK = 6
cC on the MAC IOUERR = 6 and on the IBM IOUERR = 6 and on MTS IOUERR = 7
cC
c      IOUGEN = 5
c      IOUREP = 6
c      IOUASK = 0
c      IOUERR = 6
c      IOBAT = 19
cC
c 310  CONTINUE
cC
c   10 WRITE (IOUASK,20)
c   20 FORMAT (///' ','Welcome to MOBILE5a (26-Mar-93) ',///,
c     * '&','Do you want to use the interactive input mode (Y/N)?')
c      READ (IOUASK,30) INTER
c   30 FORMAT (A1)
cC
c      IF (INTER .EQ. 'y') INTER = 'Y'
c      IF (INTER .EQ. 'n') INTER = 'N'
c      IF (INTER .NE. 'Y' .AND. INTER .NE. 'N') GO TO 310
cC  IBM-PC Mode
c      IF (INTER .EQ. 'Y'.AND.IOUGEN.EQ.5) IOUGEN=0
cC  MTS-Mainframe Mode
c      IF (INTER .EQ. 'Y'.AND.IOUGEN.EQ.1) IOUGEN=5
c      BATCH = 'N'
c      IF (INTER .EQ. 'N') THEN
cC
c 320    CONTINUE
c        WRITE (IOUASK,40)
c   40   FORMAT (/'&', 'Do you want to use the batch file input mode (Y/N
c     1)?')
c        READ (IOUASK,30) BATCH
cC
c        IF (BATCH .EQ. 'y') BATCH = 'Y'
c        IF (BATCH .EQ. 'n') BATCH = 'N'
c        IF (BATCH .NE. 'Y' .AND. BATCH .NE. 'N') GO TO 320
c      END IF
cC
c      IF (INTER .EQ. 'Y') GO TO 90
c      IF (BATCH .EQ. 'Y') GO TO 110
cC
c   50 WRITE (IOUASK,60)
c   60 FORMAT (/'&', 'Please enter the MOBILE5 input filename:')
c      READ (IOUASK,70,END=270) M5IN
c
c  LMOS Modification 04/08/93
c
c  Note that this next format statement is left uncommented since
c  it is accessed by an active portion of the subroutine
c
   70 FORMAT (A40)
cC
c      IF (M5IN .NE. BLANK) GO TO 90
c      WRITE (IOUASK,80)
c   80 FORMAT (' Invalid filename, please try again.'/' ')
c      GO TO 50
cC
c   90 IF(INTER.EQ.'N'.AND.IOUGEN.EQ.9) IOUGEN=5
c      IF(INTER.EQ.'N') OPEN(IOUGEN,FILE=M5IN,STATUS='OLD')
c      WRITE (IOUASK,100)
c  100 FORMAT (/'&', 'Please enter the MOBILE5 output filename:')
c      READ (IOUASK,70,END=270) M5OUT
cC
c      IF (M5OUT .NE. BLANK) GO TO 130
c      WRITE (IOUASK,80)
c      GO TO 90
cC
c  110 WRITE (IOUASK,120)
c  120 FORMAT (/'&', 'Please enter the MOBILE5 batch filename:')
c      READ (IOUASK,70,END=270) BATIN
cC
c      IF (BATIN .NE. BLANK) GO TO 150
c      WRITE (IOUASK,80)
c      GO TO 110
cC
c  130 WRITE (IOUASK,140)
c  140 FORMAT (/)
C
c
c  LMOS Modification 04/08/93
c
c  Assign INTER = 'N' (non-interactive input) and
c  BATCH = 'N' (non-batch input)
c
      INTER = 'N'
      BATCH = 'N'
c
c  Assign output unit numbers.
c
      IOUASK = 99
      IOUGEN = 7
      IOUREP = 8
      IOUERR = 9
      IOUIMD = 10
c
c  Get the EMS environment variables  where the MOBILE5 input, output
c  I/M credit, and alternate LEV input files will be stored.
c  Input and output files will be stored in EMS_LOC.
c
      CALL GETENV("EMS_LOC",ENVL)
      M5IN = ENVL(1:LNBLNK(ENVL))//'/M5IN.DAT'
      M5OUT = ENVL(1:LNBLNK(ENVL))//'/M5OUT.DAT'
      M5STAT = ENVL(1:LNBLNK(ENVL))//'/M5_STATUS_REPORT'
c
      OPEN(IOUGEN,FILE=M5IN,STATUS='OLD')
      OPEN(IOUREP,FILE=M5OUT,STATUS='UNKNOWN')
      OPEN(IOUASK,FILE=M5STAT,STATUS='UNKNOWN')
C
c    Comment out opening of file M5OUT
c
c      OPEN (IOUREP,FILE=M5OUT)
c
c  Remainder of DRIVER should not be used since we have assigned INTER = 'N'
c  and BATCH = 'N'
C
  150 IF (BATCH .EQ. 'Y') OPEN (IOBAT,FILE=BATIN,STATUS='OLD')
C
      IF (INTER .EQ. 'Y') WRITE (IOUASK,160)
  160 FORMAT (//' ', 'Please enter the PROMPT flag...'/'&',
     1       '(Enter 2 for Vertical or 4 for Horizontal Prompting)')
C
      NFILE = 0
C
  170 IF (BATCH .EQ. 'Y') THEN
        IF(INTER.EQ.'N'.AND.IOUGEN.EQ.9) IOUGEN=5
        IF (NFILE .GT. 0) CLOSE (IOUGEN)
        IF (NFILE .GT. 0) CLOSE (IOUREP)
        READ (IOBAT,70,END=230) M5IN, M5OUT
        OPEN (IOUGEN,FILE=M5IN)
        OPEN (IOUREP,FILE=M5OUT)
        NFILE = NFILE + 1
        WRITE (IOUASK,180) NFILE, M5IN
  180   FORMAT (/' ', 'Processing batch file #', I4, ': ', A40)
      END IF
C
      M5CALL = 0
      INERR = 0
C
  190 M5CALL = M5CALL + 1
      CALL MOBILE(INERR,*210)
      WRITE (IOUASK,200) M5CALL, INERR
  200 FORMAT (' ', 'Run #', I3, '  INERR =', I2/)
      GO TO 190
C
  210 WRITE (IOUASK,220)
  220 FORMAT (' ', 'DRIVER calls completed.'/)
C
      IF (BATCH .EQ. 'Y') GO TO 170
C
  230 IF (BATCH .EQ. 'Y') WRITE (IOUASK,240) NFILE
  240 FORMAT (//' ', I4, ' Batch file(s) processed')
C
      IF (IOUASK .EQ. 9) THEN
        WRITE (IOUASK,250)
  250   FORMAT (' Hit Return... then click cancel')
        READ (IOUASK,260)
  260   FORMAT (1X)
      END IF
C
  270 STOP
      END
