	
        PROGRAM CNTLMAT

C***********************************************************************
C  program body starts at line 136
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
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C  
C See file COPYRIGHT for conditions of use.
C  
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C  
C smoke@emc.mcnc.org
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

C.........This module is required by the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters        
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         GETIFDSC
        INTEGER         PROMPTFFILE
        
        EXTERNAL        CRLF, ENVYN, GETIFDSC, PROMPTFFILE

C...........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

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
        INTEGER         ATMPDEV      !  file unit no. for tmp ADD file
        INTEGER         CDEV         !  control file unit no.
        INTEGER         CTMPDEV      !  file unit no. for tmp CTL file
        INTEGER         GTMPDEV      !  file unit no. for tmp CTG file
        INTEGER         IDEV         !  tmp unit number if inven is map-formatted
        INTEGER         LDEV         !  log file unit no.
        INTEGER         LTMPDEV      !  file unit no. for tmp ALW file
        INTEGER         SDEV         !  ASCII part of inventory unit no.
        INTEGER         TDEV         !  tracking file unit no.

        CHARACTER*16    ANAME   ! logical name for ASCII inventory input file
        CHARACTER*16    ENAME   ! logical name for i/o api inventory input file
        CHARACTER*16    INAME   ! tmp logical name for inven file of unknown fmt
        CHARACTER*16    MNAME   ! logical name for multiplicative control matrix
        CHARACTER*16    PNAME   ! logical name for projection matrix
        CHARACTER*16    SNAME   ! logical name for ascii inventory input file

C...........   Other local variables
        INTEGER         CPYEAR       !  control packet year to project to
        INTEGER         IOS          !  I/O status
        INTEGER         ENLEN        !  length of the emissions inven name
        INTEGER         NCPE         !  no control packet entries
        INTEGER         NINVARR      !  number inventory variables to input
        INTEGER      :: PYEAR   = 0  !  projected year of inventory
        INTEGER         SYEAR        !  year for projecting from

        LOGICAL      :: CFLAG   = .FALSE.  ! true: control cntls in use
        LOGICAL      :: DFLAG   = .FALSE.  ! true: additive cntls in use
        LOGICAL      :: EFLAG   = .FALSE.  ! true: error has occurred
        LOGICAL      :: GFLAG   = .FALSE.  ! true: ctg cntls in use
        LOGICAL      :: JFLAG   = .FALSE.  ! true: projections in use
        LOGICAL      :: KFLAG   = .FALSE.  ! true: tracking file in use
        LOGICAL      :: LFLAG   = .FALSE.  ! true: allowable cntls in use
        LOGICAL      :: RFLAG   = .FALSE.  ! true: reactivty cntls in use
        LOGICAL      :: SFLAG   = .FALSE.  ! true: EMS-95 fmt controls
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
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Get environment variable values...
C.........  Get type of projection entries: with year or without it (EPS)
        YFLAG = ENVYN( 'PROJECTION_YR_SPEC', 
     &                 'Projection entries in year-specific format',
     &                 .TRUE., IOS )

        KFLAG = ENVYN( 'CONTROL_TRACKING', 
     &                 'Use a special file to track specific sources',
     &                 .FALSE., IOS )

C.........  Warning if tracking file use is attempted
        IF( KFLAG ) THEN
            MESG = 'WARNING: Specific source tracking has not been ' //
     &             'implemented. ' // CRLF() // BLANK10 //
     &             'CONTROL_TRACKING variable will have no effect.'
            CALL M3MSG2( MESG )
        END IF

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open input I/O API and ASCII files
        MESG= 'Enter logical name for the I/O API or MAP INVENTORY file'
        CALL PROMPTWHAT( MESG, FSREAD3, .TRUE., .TRUE., ENAME,
     &                   PROGNAME, INAME, IDEV )

C.........  If input file is ASCII format, then open and read map 
C           file to check files, sets environment for ENAME, opens 
C           files, stores the list of physical file names for the 
C           pollutant files in the MODINFO module, and stores the map
C           file switch in MODINFO as well.
        IF( IDEV .GT. 0 ) THEN

            CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

C.........  Otherwise, open separate I/O API and ASCII files that
C           do not store the pollutants as separate 
        ELSE
            ENAME = INAME
            SDEV = PROMPTFFILE( 
     &             'Enter logical name for the ASCII INVENTORY file',
     &             .TRUE., .TRUE., ANAME, PROGNAME )
        END IF

        CDEV = PROMPTFFILE( 
     &           'Enter logical name for ASCII CONTROL PACKETS file',
     &           .TRUE., .TRUE., 'GCNTL', PROGNAME )

C.........  Store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
        CALL GETSINFO( ENAME )

C.........  Check for future-year inventory file
        PYEAR = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .FALSE. )

C.........  Set starting year on which to base possible projections created
C           in this program
        SYEAR = BYEAR
        IF( PYEAR .GT. 0 ) SYEAR = PYEAR

C.........  Set inventory variables to read for all source categories
        IVARNAMS( 1 ) = 'INVYR'
        IVARNAMS( 2 ) = 'CSCC'
        IVARNAMS( 3 ) = 'CSOURC'

C.........  Set inventory variables to read for specific source categories
        IF( CATEGORY .EQ. 'AREA' ) THEN
            NINVARR = 3

        ELSE IF( CATEGORY .EQ. 'MOBILE' ) THEN
            NINVARR = 4
            IVARNAMS( 4 ) = 'CVTYPE'

        ELSE IF( CATEGORY .EQ. 'POINT' ) THEN
            NINVARR = 4
            IVARNAMS( 4 ) = 'ISIC'

        END IF

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Build unique lists of SCCs per SIC from the inventory arrays
        CALL GENUSLST

C.........  Initialize arrays
        PKTCNT = 0  ! array
        PKTBEG = 0  ! array
        XRFCNT = 0  ! array

C.........  Allocate memory for control packet information in input file.
        CALL ALOCPKTS( CDEV, SYEAR, CPYEAR, PKTCNT, 
     &                 PKTBEG, XRFCNT )

C.........  Set the flags that indicate which packets are valid
        GFLAG = ( PKTCNT( 1 ) .GT. 0 )
        CFLAG = ( PKTCNT( 2 ) .GT. 0 )
        LFLAG = ( PKTCNT( 3 ) .GT. 0 )
        DFLAG = ( PKTCNT( 4 ) .GT. 0 )
        RFLAG = ( PKTCNT( 5 ) .GT. 0 )
        JFLAG = ( PKTCNT( 6 ) .GT. 0 )
        SFLAG = ( PKTCNT( 7 ) .GT. 0 )

C.........  Cannot have CONTROL and EMS_CONTROL packet in same inputs
        IF( CFLAG .AND. SFLAG ) THEN
           MESG = 'CONTROL and EMS_CONTROL packets cannot be ' //
     &            'in the same input file'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Process packets: this means read packet, sort it, group it into 
C           grouped, x-ref structures in MODXREF, and assign to sources, and
C           write outputs for non-pollutant specific packets.
C.........  For control matrices that depend on pollutants, temporary files
C           will be written if there is more than one pollutant group, and
C           these will be used to store the control data index information for
C           each packet type while determining the pollutants to use in opening
C           the final output files.

        ACTION = 'PROCESS'
        CALL PKTLOOP( CDEV, ATMPDEV, CTMPDEV, GTMPDEV, LTMPDEV,  
     &                CPYEAR, ACTION, ENAME, PKTCNT, PKTBEG, XRFCNT )

C.........  Process control matrices that depend on pollutants...

C.........  Multiplicative matrix
        IF( CFLAG .OR. GFLAG .OR. LFLAG .OR. SFLAG ) THEN

C.............  Write-out control matrix
            NCPE = MAX( PKTCNT( 2 ), PKTCNT( 7 ) )
            CALL GENMULTC( ATMPDEV, CTMPDEV, GTMPDEV, LTMPDEV,
     &                     NCPE, PYEAR, ENAME, MNAME, CFLAG, GFLAG,
     &                     LFLAG, SFLAG )
        END IF

C.........  Post-process temporary files to create final report file
        CALL WCNTLREP( ATMPDEV, CTMPDEV, GTMPDEV, LTMPDEV )

C.........  Successful completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END PROGRAM CNTLMAT

