	
        PROGRAM CNTLMAT

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Copied from CTLPMAT.F version 4.4 by M. Houyoux 3/99
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C  
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C**********************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

C.........  This module contains the speciation profiles
        USE MODSPRO

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters        
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         ENVYN
        INTEGER         GETIFDSC
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        
        EXTERNAL        ENVYN, GETIFDSC, PROMPTFFILE, PROMPTMFILE

C...........  LOCAL PARAMETERS and their descriptions:

        INTEGER     , PARAMETER :: MXCHRS  = 7
        INTEGER     , PARAMETER :: NPACKET = 7

        CHARACTER*50, PARAMETER :: SCCSW = '@(#)$Id$'
        CHARACTER*20, PARAMETER :: PKTLIST( NPACKET ) = 
     &                          (  / 'CTG                 ',
     &                               'CONTROL             ',
     &                               'ALLOWABLE           ',
     &                               'ADD                 ',
     &                               'REACTIVITY          ',
     &                               'PROJECT AMS         ',
     &                               'PROJECT PTS         '  / )

C...........   LOCAL VARIABLES and their descriptions:

C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C...........   Arrays for each of the packets
        INTEGER      PKTCNT( NPACKET )  ! count of records per packet
        INTEGER      PKTBEG( NPACKET )  ! first line of packet in input file
        INTEGER      XRFCNT( NPACKET )  ! cross-reference info memory needs

C...........   Allocatable local arrays
        
        LOGICAL, ALLOCATABLE :: TFLAG( : )  !  flags:  track these sources

C...........   Logical names and unit numbers
        INTEGER         CDEV         !  control file unit no.
        INTEGER         LDEV         !  log file unit no.
        INTEGER         RDEV         !  report file unit no.
        INTEGER         SDEV         !  ASCII part of inventory unit no.
        INTEGER         TDEV         !  tracking file unit no.

        CHARACTER*16    ANAME   ! logical name for additive control matrix
        CHARACTER*16    ENAME   ! logical name for i/o api inventory input file
        CHARACTER*16    INAME   ! logical name for ASCII inventory input file
        CHARACTER*16    MNAME   ! logical name for multiplicative control matrix
        CHARACTER*16    PNAME   ! logical name for projection matrix
        CHARACTER*16    SNAME   ! logical name for ascii inventory input file

C...........   Other local variables
        INTEGER         CPYEAR       !  control packet year to project to
        INTEGER         IOS          !  I/O status
        INTEGER         ENLEN        !  length of the emissions inven name
        INTEGER         NSRC         !  number of sources
        INTEGER         PYEAR        !  projected year of inventory
        INTEGER         SYEAR        !  year for projecting from

        LOGICAL      :: CFLAG   = .FALSE.  ! true: control cntls in use
        LOGICAL      :: DFLAG   = .FALSE.  ! true: additive cntls in use
        LOGICAL      :: EFLAG   = .FALSE.  ! true: error has occurred
        LOGICAL      :: GFLAG   = .FALSE.  ! true: ctg cntls in use
        LOGICAL      :: JFLAG   = .FALSE.  ! true: projections in use
        LOGICAL      :: KFLAG   = .FALSE.  ! true: tracking file in use
        LOGICAL      :: LFLAG   = .FALSE.  ! true: allowable cntls in use
        LOGICAL      :: RFLAG   = .FALSE.  ! true: reactivty cntls in use
        LOGICAL      :: OFLAG   = .FALSE.  ! true: create report
        LOGICAL      :: YFLAG   = .FALSE.  ! true: projection entries have years

        CHARACTER*7     ACTION             ! buffer for PKTLOOP action
        CHARACTER*300   MESG               ! message buffer
        
        CHARACTER*16  :: PROGNAME = 'CNTLMAT'   !  program name
        
C***********************************************************************
C   begin body of program CNTLMAT
        
        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Get environment variable values...
C.........  Get type of projection entries: with year or without it (EPS)
        YFLAG = ENVYN( 'PROJECTION_YR_SPEC', 
     &                 'Projection entries in year-specific format',
     &                 .TRUE., IOS )

        TFLAG = ENVYN( 'CONTROL_REPORT', 
     &                 'Output a controls report file', .FALSE., IOS )

        KFLAG = ENVYN( 'CONTROL_TRACKING', 
     &                 'Use a special file to track specific sources',
     &                 .FALSE., IOS )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........   Get file names and open files
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for the I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )
        ENLEN = LEN_TRIM( ENAME )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for the ASCII INVENTORY file',
     &           .TRUE., .TRUE., ANAME, PROGNAME )

        CDEV = PROMPTFFILE( 
     &           'Enter logical name for ASCII CONTROL PACKETS file',
     &           .TRUE., .TRUE., 'GCNTL', PROGNAME )


C.........  Get header description of inventory file, error if problem
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             ENAME( 1:ENLEN ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
        ELSE

            CALL GETSINFO

C.............  Store non-category-specific header information
            NSRC = NROWS3D

C.............  Determine if file is a base or future-year inventory file
            BYEAR = GETIFDSC( FDESC3D, '/BASE YEAR/'     , .TRUE.  )
            PYEAR = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .FALSE. )

C.............  Set starting year on which to base possible projections created
C               in this program
            SYEAR = BYEAR
            IF( PYEAR .GT. 0 ) SYEAR = PYEAR

        ENDIF

C.........  Read inventory source characteristics
        CALL M3MSG2( 'Reading in inventory file...' )

C.........  Set inventory variables to read for all source categories
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'CSCC'
        IVARNAMS( 3 ) = 'INVYR'

C.........  Allocate memory for and read required inventory characteristics
        IF( CATEGORY .EQ. 'AREA' ) THEN

            CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, 3, IVARNAMS )

        ELSE IF( CATEGORY .EQ. 'MOBILE' ) THEN

            IVARNAMS( 4 ) = 'CVTYPE'
            CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, 4, IVARNAMS )

        ELSE IF( CATEGORY .EQ. 'POINT' ) THEN

            IVARNAMS( 4 ) = 'ISIC'
            IVARNAMS( 5 ) = 'CSOURC'

            CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, 5, IVARNAMS )

        END IF

C.........  Build unique lists of SCCs per SIC from the inventory arrays
        CALL GENUSLST

C.........  Allocate memory for control packet information in input file.
        CALL ALOCPKTS( CDEV, SYEAR, NPACKET, PKTLIST, CPYEAR, PKTCNT, 
     &                 PKTBEG, XRFCNT )

C.........  Set the flags that indicate which packets are valid
        GFLAG = ( PKTCNT( 1 ) .GT. 0 )
        CFLAG = ( PKTCNT( 2 ) .GT. 0 )
        LFLAG = ( PKTCNT( 3 ) .GT. 0 )
        DFLAG = ( PKTCNT( 4 ) .GT. 0 )
        RFLAG = ( PKTCNT( 5 ) .GT. 0 )
        JFLAG = ( PKTCNT( 6 ) .GT. 0 .OR. PKTCNT( 7 ) .GT. 0 )

C.........  Process packets: this means read packet, sort it, group it into 
C           grouped, x-ref structures in MODXREF, and assign to sources, and
C           write outputs for non-pollutant specific packets.
C.........  For control matrices that depend on pollutants, temporary files
C           will be written if there is more than one pollutant group, and
C           these will be used to store the control data index information for
C           each packet type while determining the pollutants to use in opening
C           the final output files.

        ACTION = 'PROCESS'
        CALL PKTLOOP( CDEV, NSRC, CPYEAR, NPACKET, ACTION, ENAME, 
     &                PKTCNT, PKTBEG, PKTLIST, XRFCNT )

C..........  Process control matrices that depend on pollutants...

C.........  Multiplicative matrix
        IF( CFLAG .OR. GFLAG .OR. LFLAG ) THEN

C..............  Open control matrix
C            CALL OPENCMAT( ENAME, NPOLMULT, 'MULTIPLICATIVE', PNAMMULT )

C            CALL GENMULTC( )  ! Post-processes temporary packet index by source
C                                files

C STOPPED HERE: Need to write opencmat, genmultc, genaddc, report post-processor
        END IF

C.........  Additive matrix
        IF( DFLAG ) THEN

C            CALL OPENCMAT( ENAME, NPOLADD, 'ADDITIVE', PNAMADD )

C            CALL GENADDC( )

        END IF

C.........  Open final report file

C.........  Post-process temporary report file to create final report file

C.........  Successful completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C NOTE: Prune format statements, if possible

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A, ':', I7, ',', I8, ',', I4.4, ',', I3.3, ' -- ', 
     &          3( F8.5, '~>', F8.5, : '; ' )  )

93020   FORMAT( A, ':', I4, ' to ', I4, 1X, I5.5, 1X, I4,
     &          3( F8.5, '~>', F8.5, : '; ' )  )

93030   FORMAT( A, ':', 1X, I5.5, 1X, I4,
     &          3( F8.5, '~>', F8.5, : '; ' )  )


C...........   Internal buffering formats............ 94xxx

94910   FORMAT( A, :, ' at line', I7, :,  '; IOSTAT=', I7 )

94920   FORMAT( A, ' "', A, :,  '" at line', I7, :, '; IOSTAT=', I7 )

94930   FORMAT( A )

        END PROGRAM CNTLMAT

