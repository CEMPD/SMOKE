
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
C
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C  
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C  
C See file COPYRIGHT for conditions of use.
C  
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C  
C env_progs@mcnc.org
C  
C Pathname: $Source$
C Last updated: $Date$ 
C  
C***********************************************************************

C...........   MODULES for public variables

C.........  This module contains Smkreport-specific settings
        USE MODREPRT

C.........  This module contains report arrays for each output bin
        USE MODREPBN

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

C...........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2       CRLF
        INTEGER           PROMPTFFILE
        INTEGER           SECSDIFF

        EXTERNAL  CRLF, PROMPTFFILE, SECSDIFF

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   Gridding Matrix
        INTEGER, ALLOCATABLE :: GMAT( : ) ! Contiguous gridding matrix

C...........   Speciation matrices
        INTEGER, ALLOCATABLE :: SLMAT( :,: ) ! mole-based
        INTEGER, ALLOCATABLE :: SSMAT( :,: ) ! mass-based

C...........   File units and logical/physical names
        INTEGER         CDEV    !  reports configuration file
        INTEGER         EDEV    !  elevated source ID file
        INTEGER         GDEV    !  gridding supplemental file
        INTEGER         LDEV    !  log-device
        INTEGER         NDEV    !  SCC descriptions
        INTEGER         ODEV    !  output file unit number
        INTEGER         PDEV    !  speciation supplemental file
        INTEGER         SDEV    !  ASCII inven input file
        INTEGER         TDEV    !  temporal supplemental files
        INTEGER         YDEV    !  country/state/county names file

        CHARACTER*16  :: ANAME  = ' '   !  logical name for ASCII inven input 
        CHARACTER*16  :: ENAME  = ' '   !  logical name for I/O API inven input
        CHARACTER*16  :: GNAME  = ' '   !  gridding matrix input
        CHARACTER*16  :: LNAME  = ' '   !  layer fractions input file
        CHARACTER*16  :: SLNAME = ' '   !  speciation matrix input
        CHARACTER*16  :: SSNAME = ' '   !  speciation matrix input
        CHARACTER*16  :: TNAME  = ' '   !  hourly emissions input file

        CHARACTER*300 :: FNAME = ' '    !  output physical/logical file name
        CHARACTER*300 :: PNAME = ' '    !  previous output file name

C...........   Other local variables
        INTEGER      I, K, L, N           ! indices and counters

        INTEGER      HWID                 ! header width
        INTEGER      IOS                  ! i/o status
        INTEGER   :: GDIM    = 0          ! dimension of contiguous gridding mat
        INTEGER   :: NSLIN   = 1          ! no. mole input speciation variables
        INTEGER   :: NSSIN   = 1          ! no. mass input speciation variables

        LOGICAL       :: EFLAG = .FALSE.        ! true: error found

        CHARACTER*300          MESG             !  message buffer
        CHARACTER(LEN=QAFMTL3) OUTFMT           !  data output format string

        CHARACTER*16 :: PROGNAME = 'SMKREPORT' ! program name

C***********************************************************************
C   begin body of program SMKREPORT
        
        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
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

C.........  Prompt for and open all other input files
        CALL OPENREPIN( ENAME, ANAME, GNAME, LNAME, SLNAME, SSNAME, 
     &                  TNAME, SDEV, GDEV, PDEV, TDEV, EDEV, YDEV, NDEV)

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
        CALL RDREPIN( GDIM, NSLIN, NSSIN, SDEV, GDEV, PDEV, TDEV, EDEV, 
     &                YDEV, NDEV, ENAME, GNAME, LNAME, SLNAME, SSNAME, 
     &                GMAT( 1 ), GMAT( NGRID+1 ),
     &                GMAT( NGRID+NMATX+1 ), SSMAT, SLMAT )

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

C.................  When not first report...
                IF( N .GT. 1 ) THEN

C.....................  Add Metadata to previous file if current file number is
C                       different from previous file number
c                    CALL WRMETADAT( FDEV )
c                    note: Need to write this

C.....................  Close output file
                    CLOSE( ODEV )

                END IF

C.................  Open new output file if current file number is different 
C                   previous file number.
                CALL OPENREPOUT( FNAME, ODEV )

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
C               ozone-season emissions.
            IF( .NOT. DESC3( ENAME ) ) THEN

                L = LEN_TRIM( ENAME )
                MESG = 'Could not get description of file "' //
     &                 ENAME( 1:L ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE
                INVPIDX = 0
                IF ( RPT_%O3SEASON ) INVPIDX = 1
                CALL GETSINFO

            END IF

C.............  Determine input units and create conversion factors
            CALL REPUNITS( N )

C.............  Write report header
            CALL WRREPHDR( ODEV, N, HWID, OUTFMT )

C.............  Loop through time steps (if any) and sum emissions into bins
C               for the appropriate time resolution...

C.............  For mole-based speciation...
            IF( RPT_%USESLMAT ) THEN
                CALL GENRPRT( ODEV, N, HWID, ENAME, TNAME, LNAME, 
     &                        OUTFMT, SLMAT, EFLAG )

C.............  For mass-based and no speciation
            ELSE
                CALL GENRPRT( ODEV, N, HWID, ENAME, TNAME, LNAME, 
     &                        OUTFMT, SSMAT, EFLAG )
            END IF

C.............  Save file number to use in next iteration
            PNAME  = FNAME

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
        
 
