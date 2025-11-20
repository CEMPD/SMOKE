
        PROGRAM SMKREPORT

C***********************************************************************
C  subroutine body starts at line 129
C
C  DESCRIPTION:
C    The SMKREPORT routine create emissions and activity reports for one
C    major source category at a time (area, mobile, or point). It permits
C    the user to control the columns and rows in the report through a series
C    of instructions.  These reports allow users to quality assure emissions
C    inventories by comparison and analysis of reports from different stages
C    of SMOKE processing. The outputs can be read into the Java Analysis and
C    Report Tool (JART).
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Revised 7/2003 by A. Holland
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$
C
C***********************************************************************

C...........   MODULES for public variables
C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: QAFMTL3, NMATX, SSFLAG, SLFLAG, NREPORT,
     &                      RPT_, AFLAG, ALLRPT

C.........  This module contains report arrays for each output bin
        USE MODREPBN, ONLY: NSPCIN

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID

C...........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, INVPIDX, NCHARS, JSCC

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER(2)      CRLF
        LOGICAL           ENVYN
        INTEGER           PROMPTFFILE
        INTEGER           SECSDIFF

        EXTERNAL  CRLF, ENVYN, PROMPTFFILE, SECSDIFF

C...........   LOCAL PARAMETERS
        CHARACTER(50), PARAMETER ::
     &  CVSW = '$Name SMOKEv5.2.1_Sep2025$' ! CVS release tag

C...........   Gridding Matrix
        INTEGER, ALLOCATABLE :: GMAT( : ) ! Contiguous gridding matrix

C...........   Speciation matrices
        INTEGER, ALLOCATABLE :: SLMAT( :,: ) ! mole-based
        INTEGER, ALLOCATABLE :: SSMAT( :,: ) ! mass-based

C...........   File units and logical/physical names
        INTEGER :: ADEV = 0   !  ASCII elevated file
        INTEGER :: CDEV = 0   !  reports configuration file
        INTEGER :: EDEV = 0   !  elevated source ID file
        INTEGER :: GDEV = 0   !  gridding supplemental file
        INTEGER :: LDEV = 0   !  log-device
        INTEGER :: NDEV = 0   !  SCC descriptions
        INTEGER :: NIDEV = 0  !  SIC descriptions
        INTEGER :: NMDEV = 0  !  MACT descriptions
        INTEGER :: NNDEV = 0  !  NAICS descriptions
        INTEGER :: NODEV = 0  !  ORIS descriptions
        INTEGER :: MODEV = 0  !  src mapping file
        INTEGER :: PDEV = 0   !  speciation supplemental file
        INTEGER :: NPDEV = 0  !  GSPRO descriptions
        INTEGER :: RDEV(3) = ( / 0,0,0 / ) !  ASCII reports from Cntlmat program
        INTEGER :: SDEV = 0   !  ASCII inven input file
        INTEGER :: TDEV = 0   !  temporal supplemental files
        INTEGER :: YDEV = 0   !  country/state/county names file

        INTEGER, ALLOCATABLE :: ODEV( : )   !  output file unit numbers

        CHARACTER(16)  :: ANAME  = ' '   !  logical name for ASCII inven input
        CHARACTER(16)  :: ENAME  = ' '   !  logical name for I/O API inven input
        CHARACTER(16)  :: CUNAME = ' '   !  multiplicative control matrix input
        CHARACTER(16)  :: GNAME  = ' '   !  gridding matrix input
        CHARACTER(16)  :: LNAME  = ' '   !  layer fractions input file
        CHARACTER(16)  :: PRNAME = ' '   !  projection matrix input
        CHARACTER(16)  :: SLNAME = ' '   !  speciation matrix input
        CHARACTER(16)  :: SSNAME = ' '   !  speciation matrix input
        CHARACTER(16)  :: TNAME  = ' '   !  hourly emissions input file

        CHARACTER(300) :: FNAME = ' '    !  output physical/logical file name
        CHARACTER(300) :: PNAME = ' '    !  previous output file name

C...........   Other local variables
        INTEGER      I, J,  K, L, N       ! indices and counters

        INTEGER      HWID                 ! header width
        INTEGER      IOS                  ! i/o status
        INTEGER      EDIDX                ! ending index of loop
        INTEGER   :: GDIM    = 0          ! dimension of contiguous gridding mat
        INTEGER   :: NSLIN   = 1          ! no. mole input speciation variables
        INTEGER   :: NSSIN   = 1          ! no. mass input speciation variables

        REAL         RNFILES              ! real number of files per report
        REAL         RNSECT               ! real number of sections per report

        LOGICAL       :: EFLAG    = .FALSE.     ! true: error found
        LOGICAL       :: ZEROFLAG = .FALSE.     ! true: report zero values

        CHARACTER(300)     MESG             !  message buffer
        CHARACTER(QAFMTL3) OUTFMT           !  data output format string

        CHARACTER(16) :: PROGNAME = 'SMKREPORT' ! program name

C***********************************************************************
C   begin body of program SMKREPORT

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Prompt for and open REPCONFIG file
        CDEV = PROMPTFFILE(
     &           'Enter logical name for the REPORT CONFIGURATION file',
     &           .TRUE., .TRUE., 'REPCONFIG', PROGNAME )

C.........  Scan report configuration file to determine input file types,
C           get global file flags and settings, and determine maximum
C           values for use in memory allocation.
        CALL SCANREPC( CDEV )

C.........  Get environment variable settings
        ZEROFLAG = ENVYN( 'REPORT_ZERO_VALUES', 'Leave entries ' //
     &                    'with values equal to zero in reports',
     &                    .FALSE., IOS )

C.........  Prompt for and open all other input files
        CALL OPENREPIN( ENAME, ANAME, CUNAME, GNAME, LNAME,
     &                  PRNAME, SLNAME, SSNAME, TNAME, RDEV,
     &                  SDEV, GDEV, PDEV, TDEV, EDEV, YDEV, NDEV,
     &                  NIDEV, NPDEV, ADEV, NMDEV, NNDEV, NODEV )

C.........  Read and store all report instructions
        CALL RDRPRTS( CDEV )

C.........  Build pollutant, activity, emis-type, and species indices to output
C           data columns based on the selected output data from the reports
C.........  Index arrays are stored in the bin module
        CALL BLDREPIDX( SLNAME, SSNAME )

C.........  Allocate memory for gridding matrix (even if not used so that it
C           can be passed through subroutines)
        GDIM = NGRID + 2*NMATX
        ALLOCATE( GMAT( GDIM ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GMAT', PROGNAME )

C.........  Allocate memory for speciation matrices (even if no speciation
C           so that arrays can be passed through subroutines).
        N = 1
        IF( SLFLAG .OR. SSFLAG ) N = NSRC

        IF( SLFLAG ) NSLIN = NSPCIN
        ALLOCATE( SLMAT( N, NSLIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SLMAT', PROGNAME )

        IF( SSFLAG ) NSSIN = NSPCIN
        ALLOCATE( SSMAT( N, NSSIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SSMAT', PROGNAME )

C.........  Read one-time input file data
        CALL RDREPIN( NSLIN, NSSIN, RDEV, SDEV, GDEV, PDEV, TDEV,
     &                EDEV, YDEV, NDEV, NIDEV, NPDEV, NMDEV, NNDEV,
     &                NODEV, ADEV, ENAME, CUNAME, GNAME, LNAME,
     &                PRNAME, SLNAME, SSNAME, GMAT( 1 ),
     &                GMAT( NGRID+1 ), GMAT( NGRID+NMATX+1 ),
     &                SSMAT, SLMAT )

C.........  Preprocess the country/state/county data
c note: Could add routine to reduce list of co/st/cy data to just records
c    n: selected.

C.........  Preprocess the inventory data
c note: Could add routine to reduce source to just records
c    n: selected across all groups.

C.........  Read and store all group definitions
        CALL RDGRPS( CDEV )

C.........  Loop through reports
        DO N = 1, NREPORT

            RPT_ = ALLRPT( N )

C............  Determine number of output files/sections per report

            IF( RPT_%RPTNVAR .GT. RPT_%NUMDATA ) THEN
                RPT_%RPTNVAR = RPT_%NUMDATA
            END IF

            IF( RPT_%RPTMODE .EQ. 1 ) THEN

                RPT_%NUMSECT = 1

                RNFILES = REAL( RPT_%NUMDATA ) / REAL( RPT_%RPTNVAR )

                IF( RNFILES .LT. 1.0 ) THEN

                    RPT_%NUMFILES = 1

                ELSE

                    RPT_%NUMFILES = INT( RNFILES )

                    IF( RNFILES .GT. RPT_%NUMFILES ) THEN
                        RPT_%NUMFILES = RPT_%NUMFILES + 1
                    END IF

                END IF

            ELSE IF( RPT_%RPTMODE .EQ. 2 ) THEN

                RPT_%NUMFILES = 1

                RNSECT = REAL( RPT_%NUMDATA ) / REAL( RPT_%RPTNVAR )

                IF( RNSECT .LT. 1.0 ) THEN

                    RPT_%NUMSECT = 1

                ELSE

                    RPT_%NUMSECT = INT( RNSECT )

                    IF( RNSECT .GT. RPT_%NUMSECT ) THEN
                        RPT_%NUMSECT = RPT_%NUMSECT + 1
                    END IF

                END IF

            ELSE

                RPT_%NUMFILES = 1
                RPT_%NUMSECT = 1

            END IF

            ALLRPT( N )%NUMFILES = RPT_%NUMFILES
            ALLRPT( N )%NUMSECT  = RPT_%NUMSECT
            ALLRPT( N )%RPTNVAR  = RPT_%RPTNVAR


            WRITE( MESG,94010 )
     &             '***** CHECKING INPUTS FOR REPORT', N, ' *****'
            CALL M3MSG2( MESG )

C.............  QA reports configuration file settings
            CALL QAREPIN( N, IOS )

C.............  Skip report if errors are found
            IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) '***** SKIPPING REPORT', N, ' *****'
                CALL M3MSG2( MESG )
                CYCLE
            END IF

C.............  Write message to log and standard output for report that is
C               being processed

C.............  Get file name
            FNAME = RPT_%OFILENAM

C.............  If current file is different than previous
            IF( FNAME .NE. PNAME ) THEN

                IF( ALLOCATED( ODEV ) ) DEALLOCATE( ODEV )

C................  Allocate output file number array
                ALLOCATE( ODEV( RPT_%NUMFILES ), STAT=IOS )
                CALL CHECKMEM( IOS, 'ODEV', PROGNAME )
                ODEV = 0

C.................  When not first report...
                IF( N .GT. 1 ) THEN

C.....................  Add Metadata to previous file if current file number is
C                       different from previous file number
c                    CALL WRMETADAT( FDEV )
c                    note: Need to write this

C.....................  Close output file(s)
                    DO I = 1, RPT_%NUMFILES
                        CLOSE( ODEV( I ) )
                    END DO

                END IF

C.................  Open new output file if current file number is different
C                   previous file number.
                CALL OPENREPOUT( FNAME, ODEV, MODEV )

            END IF


            MESG = BLANK10 // 'Selecting records...'
            CALL M3MSG2( MESG )

C.............  Select inventory records
            CALL SELECTSRC( N )

C.............  Apply gridding information
            IF( RPT_%USEGMAT ) THEN
                CALL REPMRGGRD( N, GMAT( 1 ), GMAT( NGRID+1 ),
     &                          GMAT( NGRID+NMATX+1 ), EFLAG   )
            END IF

C.............  Skip remainder of report if error found so far
            IF( EFLAG ) CYCLE

            MESG = BLANK10 // 'Aggregating output records...'
            CALL M3MSG2( MESG )

C.............  Assign bin numbers to selected records
            CALL ASGNBINS( N )

            MESG = BLANK10 //
     &             'Reading emissions data and writing report...'
            CALL M3MSG2( MESG )

C.............  Update inventory input names and units, depending on status of
C               average day emissions.
            INVPIDX = 0
            IF ( RPT_%AVEDAY ) INVPIDX = 1
            IF( .NOT. AFLAG ) CALL GETSINFO( ENAME )
            IF( JSCC .GT. 0 ) NCHARS = NCHARS - 1  ! duplicate of rdrepin.f

C.............  Determine input units and create conversion factors
            CALL REPUNITS( N )

C.............  If report is multifile or database
            IF( RPT_%RPTMODE .EQ. 1 .OR. RPT_%RPTMODE .EQ. 3
     &          .OR. RPT_%RPTMODE .EQ. 0 ) THEN
                EDIDX = RPT_%NUMFILES

C.............  If report is multisection
            ELSE
                EDIDX = RPT_%NUMSECT

            END IF

            DO I = 1, EDIDX

                J = I
                IF( RPT_%RPTMODE .EQ. 2 ) J = 1

C.............  Write report header
                CALL WRREPHDR( ODEV( J ), N, I, HWID, OUTFMT )

C.............  Loop through time steps (if any) and sum emissions into bins
C               for the appropriate time resolution...

C.............  For mole-based speciation...
                IF( RPT_%USESLMAT ) THEN
                    CALL GENRPRT( ODEV( J ), N, ADEV, MODEV,  ENAME,
     &                     TNAME, LNAME, OUTFMT, SLMAT, ZEROFLAG,
     &                     EFLAG )

C.............  For mass-based and no speciation
                ELSE
                    CALL GENRPRT( ODEV( J ), N, ADEV, MODEV, ENAME,
     &                     TNAME, LNAME, OUTFMT, SSMAT, ZEROFLAG,
     &                     EFLAG )
                END IF

C.............  Save file number to use in next iteration
                PNAME  = FNAME

            END DO    ! end loop over files/sections

        END DO   ! end loop over reports

C.........  Completion with errors
        IF( EFLAG ) THEN
            MESG = 'Problem creating reports'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Normal completion
        ELSE
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )
        END IF

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END PROGRAM SMKREPORT


