
        SUBROUTINE RDGRPS( )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     The RDGRPS routine reads the groups from the REPCONFIG file, interprets
C     the entries, and stores the fully detailed group information in a table
C     for each type of group.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 7/2000 by M Houyoux
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

C.........  MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains Smkreport-specific settings
        USE MODREPRT

C.........  This module contains report arrays for each output bin
        USE MODREPBN

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  EXTERNAL FUNCTIONS and their descriptions:
c        INTEGER    STR2INT

c        EXTERNAL   STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV       ! output file unit number
        INTEGER     , INTENT (IN) :: RCNT       ! report count

C...........   Local allocatable arrays...

C...........   Region group input arrays
        INTEGER               , ALLOCATABLE :: REGRAW  ( : )
        LOGICAL               , ALLOCATABLE :: REGSTAT ( : )
        CHARACTER(LEN=LENLAB3), ALLOCATABLE :: REGNLABL( : )

C...........   Subgrid input arrays
        INTEGER               , ALLOCATABLE :: SBGRAW  ( : )
        LOGICAL               , ALLOCATABLE :: SBGSTAT ( : )
        CHARACTER(LEN=LENLAB3), ALLOCATABLE :: SUBGLABL( : )

C...........   Other local arrays
        INTEGER   NCNT( NALLPCKT )   ! Count of groups defined by SELECT 
        
C...........   Other local variables

        LOGICAL      :: EFLAG = .FALSE.     ! true: error found

        CHARACTER*16 :: PROGNAME = 'RDGRPS' ! program name

C***********************************************************************
C   begin body of subroutine RDGRPS

C.........  Allocate memory for reading and storing raw information from 
C           the input file.
        ALLOCATE( REGNLABL( NREGRAW ), STAT=IOS )
        CALL CHECKMEM( IOS, 'REGNLABL', PROGNAME )
        ALLOCATE( SUBGLABL( NSBGRAW ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUBGLABL', PROGNAME )
        
        REGNLABL = ' '  ! array
        SUBGLABL = ' '  ! array

C.........  Read in defined group labels
c        CALL READ_GROUPS( FDEV, 'DEFINED LABELS' )

C.........  Initialize number of groups of each type defined with SELECT
        NCNT = 0   ! array
        
C.........  Read in inline group labels and compare to defined groups.  If not
C           defined, try to match as country, state, or county and store
C           additional names.
        CALL READ_GROUPS( FDEV, 'SELECT LABELS' )

C.........  If error was found so far, abort
        IF( EFLAG ) THEN
            MESG = 'Problem with groups in input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set the number for each type of group
        NREGNGRP = NREGRAW + NCNT( REG_IDX )
        NSUBGRID = NSBGRAW + NCNT( SBG_IDX )

C.........  If no groups are defined, leave the subroutine
        IF( NREGNGRP .EQ. 0 .AND.
     &      NSUBGRID .EQ. 0       ) RETURN

C.........  The maximum number of records per group has been set in SCANREPC 
C           and needs to be at least 1
        MXGRPREC = MAX( MXGRPREC, 1 )

C.........  Allocate memory for raw group information. The total number of
C           groups is the sum of the defined groups and the unmatched and valid
C           in-line groups.
C.........  Reallocate memory for group labels so that labels can be reset with
C           the inline labels as well.
        IF( NREGNGRP .GT. 0 ) THEN
            DEALLOCATE( REGNLABL )
            ALLOCATE( REGNLABL( NREGRAW ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REGNLABL', PROGNAME )
            ALLOCATE( REGRAW( MXGRPREC,NREGNGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REGRAW', PROGNAME )
            ALLOCATE( REGSTAT( MXGRPREC,NREGNGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REGSTAT', PROGNAME )

            REGNLABL = ' '     ! array
            REGRAW   = 0       ! array
            REGSTAT  = .TRUE.  ! array (default is include)

        END IF

C.........  Same notes as above, but do for subgrids.
        IF( NSUBGRID .GT. 0 ) THEN
            DEALLOCATE( SUBGLABL )
            ALLOCATE( SUBGLABL( NSUBGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SUBGLABL', PROGNAME )
            ALLOCATE( SBGRAW( MXGRPREC,NSUBGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SBGRAW', PROGNAME )
            ALLOCATE( SBGSTAT( MXGRPREC,NSUBGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SBGSTAT', PROGNAME )

            SBGLABL  = ' '     ! array
            SBGRAW   = ' '     ! array
            SBGSTAT  = .TRUE.  ! array

        END IF

c Note:        STOPPED HERE 

C.........  Read raw group information.  "Raw" means the entries as they
C           appear in the input file, but not converted to the data structures
C           needed for further processing.
        CALL READ_GROUPS( FDEV, 'STORE' )

C.........  Convert entries to data structures needed for further processing.

c            Need tmp array to store flags for unused group assignments
            
C.........  Deallocate raw group information
        IF( ALLOCATED( REGRAW ) ) DEALLOCATE( REGRAW, REGSTAT )
        IF( ALLOCATED( SBGRAW ) ) DEALLOCATE( SBGRAW, SBGSTAT )

C.........  
c Loop through file to 
c determine the SELECT Region statements that do not use defined regions 
c     for memory allocation of unsorted
c Ignore group assignments that do not match defined groups or 
c     characteristics of the inventory
c Flag those records that do not match so they can be skipped later
c Count how many SELECT regions are there



c Assign memory for unsorted, full group list (with SELECT-specific 
c entries as well as standard groups)

C Loop through file to store raw group information, including:
c   include/exclude status

c      
c	it will require a count of countys by state, and by country
c	the group counts will be of records that will be EXCLUDED.  So, the
c       initial value will include all counties in the inventory, and the
c          number of records excluded will be subtracted

c NOTE: all SELECT REGION <Region code> entries will be converted to a group
c    n: by this routine, and given a special name that uses the Region code.

c        NREGRAW is the number of raw input region groups from SCANRPT
c        NSBGRAW is the number of raw input subgrids from SCANRPT

c Allocate memory for the unsorted group arrays

c Read the groups and populate the unsorted group arrays and group labels
c	group names are case insensitive	
c	call routine to scan out comments
c	store the line numbers of the input file for each group
c 	the group arrays are treated differently for different types of groups.
c 	store include and exclude records as separate entries
c
c 	GROUP REGION:
c	  Intialize list of FIPS codes with flags indicating included or
c 	    excluded, and previous operation.
c	  Based on group, go through whole list of inventory region codes and 
C           set up group.
c	  After all Region codes selected, create sorted list of only those that
c	    have NOT been selected by the group and store this with the group
c           label
c	  also store the actual number in the group (the dimensioned value will
c	    be a max based on the number of counties per state, etc.)
c	  Warning if group record is not in the inventory
c	  Warning if included and then excluded or vice-versa

c	SUBGRID:
c	  Create list of x-y grid cells to INCLUDE for each subgrid
c	  Warning if cells do not exist in domain
c 	  cells output to list must be sorted
c	  

c allocate memory for the sorted group arrays
c	use one 2-d array for each type of group.  These will be indexed
c	   by their group number.
cVALIDGRD for grid cells
cEXLCDRGN for region

c store the sorted groups
c Subgrid will store cell numbers that are included
c Region will store Regions that are excluded.


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS
 
C.............  This subprogram provides reads the group definitions
C               in various ways, as indicated by the main program
C               argument.
            SUBROUTINE READ_GROUPS( FDEV, READTYPE )

C.............  External functions
            LOGICAL     CHKINT
            INTEGER     GETNLIST
            INTEGER     INDEX1
            INTEGER     STR2INT

            EXTERNAL    CHKINT, GETNLIST, INDEX1, STR2INT

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: FDEV     ! input file unit
            CHARACTER(*), INTENT (IN) :: READTYPE   ! Reading type

C.............   Local parameters   
            INTEGER, PARAMETER :: MXSEG = 10

C.............  Subprogram local arrays
            CHARACTER*256 SEGMENT( MXSEG )

C.............  Local variables
            INTEGER       I, J, L, N       ! counters and indices

            INTEGER       FIP      ! tmp region code
            INTEGER       IOS      ! i/o status
            INTEGER       IREC     ! line number
            INTEGER       NS       ! no. segments in line
            INTEGER       RCNT     ! record count

            CHARACTER*300 BUFFER   ! tmp line buffer as uppercase
            CHARACTER*300 LINE     ! tmp line buffer
            CHARACTER*300 MESG     ! mesg buffer

            CHARACTER(LEN=LENLAB3) :: PREGNNAM
            CHARACTER(LEN=LENLAB3) :: PSUBGNAM

C----------------------------------------------------------------------

C.............  Rewind input file
            REWIND( FDEV )

C.............  Loop though file to store local array of labeled group names
            DO I = 1, NLINE_RC
            
                READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
                IREC = IREC + 1

                IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                'I/O error', IOS, 
     &                'reading report configuration file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Skip blank lines and comment lines
                IF( BLKORCMT( LINE ) ) CYCLE

C.................  Screen for appended comments and remove them
                CALL RMCOMMNT( '##', LINE )

C.................  Left-justify and convert line to upper case
                BUFFER = ADJUSTL( LINE )
                CALL UPCASE( BUFFER )

C.................  Initialize segment from previous iteration
                SEGMENT( 1:NS ) = ' '

C.................  Parse line into segments
                L = LEN_TRIM( BUFFER )
                NS = GETNLIST( L, BUFFER )
                IF( NS .GT. MXSEG ) NS = MXSEG
                CALL PARSLINE( BUFFER, NS, SEGMENT )

C.................  Interpret line of code.  Set global variables in MODREPRT.
                CALL PRCLINRC( IREC, NS, BUFFER, SEGMENT )

C.................  If a new group...
                SELECT CASE( READTYPE )

C.................  Read the defined groups and store the labels
                CASE( 'DEFINED LABELS' )

C.....................  Get current count of current packet 
                    RCNT = PKT_COUNT( PKT_IDX )
                
                    SELECT CASE( PKT_IDX )

C.....................  Store region group label
                    CASE( REG_IDX )
                        IF( LIN_DEFGRP ) REGNLABL( RCNT ) = GRP_LABEL
                
C.....................  Store subgrid label in local array
                    CASE( SBG_IDX )
                        IF( LIN_DEFGRP ) SUBGLABL( RCNT ) = GRP_LABEL
                
                    END SELECT

C.................  Count SELECT statements that do not use defined 
C                   regions.  Determine those that coorespond to valid 
C                   entries and those that should be ignored.
                CASE( 'SELECT LABELS' )
            
C.....................  Skip if report section not started yet.
                    IF( .NOT. INREPORT ) CYCLE

C.....................  A region is being selected
                    IF( LREGION .AND. RPT_%REGNNAM .NE. PREGNNAM ) THEN

C.........................  Search for region name in defined regions
                        J = INDEX1( RPT_%REGNNAM, NREGRAW, REGNLABL )
                    
C.........................  If region name is not found, then try to convert
C                           to a region code and compare with the inventory.
                        IF( J .LE. 0 ) THEN

C.............................  Check if label is an integer
                            IF( CHKINT( RPT_%REGNNAM ) ) THEN

C.................................  Convert label to region code
                                FIP = STR2INT( RPT_%REGNNAM )

C.................................  Check if label is a valid region code
C.................................  REGNTMP is a dummy argument at this stage
                                CALL CHECK_REGIONS( FIP, REGNTMP, IOS )

C.................................  If code is valid, store line number of
C                                   record
                                IF( IOS .EQ. 0 ) THEN
                                    N = NCNT( REG_IDX ) + 1
                                    RGNLINES( N ) = IREC
                                    NCNT( REG_IDX ) = N

C.................................  Otherwise, give warning
                                ELSE
                                    L = LEN_TRIM( RPT_%REGNNAM )
                                    WRITE( MESG,94010 ) 
     &                                'WARNING: Region label "' //
     &                                RPT_%REGNNAM( 1:L )// '" at line', 
     &                                IREC, 'does not match any groups'
     &                                // CRLF() // BLANK10 // 'or ' //
     &                                'region codes in the inventory. '
     &                                //'SELECT REGION will be ignored.'
                                    CALL M3MSG2( MESG )

                                END IF

C.............................  Label is not an integer, label is invalid
                            ELSE
                                L = LEN_TRIM( RPT_%REGNNAM )
                                WRITE( MESG,94010 ) 
     &                             'WARNING: Region label "' //
     &                             RPT_%REGNNAM( 1:L ) // '" at line', 
     &                             IREC, 'does not match any groups.'
     &                             // CRLF() // BLANK10 //
     &                             'SELECT REGION will be ignored.'
                                CALL M3MSG2( MESG )

                            END IF  ! If label is an integer or not

                        END IF      ! If region label not found in groups list
                    
C.........................  Store current region label for use in next iteration
                        PREGNNAM = RPT_%REGNNAM
                    
                    END IF         ! If region selected
                                              
C.....................  A subgrid is being selected
                    IF( LSUBGRID .AND. RPT_%SUBGNAM .NE. PSUBGNAM ) THEN

C.........................  Search for subgrid names in defined subgrids
                        J = INDEX1( RPT_%SUBGNAM, NSBGRAW, SUBGLABL )
                    
C.........................  If subgrid name is not found, then try to interpret
C                           entry as a subgrid definition.
                        IF( J .LE. 0 ) THEN

C............................. Check if subgrid is defined in-line
                            LCELSTAT = .FALSE.   ! array
                            CALL PARSE_SUBGRID( RPT_%SUBGNAM, NGRID, 
     &                                          LCELSTAT, NINCL )

C.............................  If subgrid is valid, store line number of record
                            IF( NINCL .GT. 0 ) THEN
                                N = NCNT( SBG_IDX ) + 1
                                SBGLINES( N ) = IREC
                                NCNT( SBG_IDX ) = N

C.............................  If subgrid is invalid
                            ELSE
                                L = LEN_TRIM( RPT_%SUBGNAM )
                                WRITE( MESG,94010 ) 
     &                             'WARNING: Subgrid definition "' //
     &                             RPT_%SUBGNAM( 1:L ) // '" at line', 
     &                             IREC, 'is not defined and cannot' //
     &                             // CRLF() // BLANK10 //
     &                             'be interpreted. SELECT SUBGRID ' //
     &                             'will be ignored.'
                                CALL M3MSG2( MESG )

                            END IF
                    

                        END IF        ! If subgrid name not found in list

C.........................  Store current subgrid label for use in next 
C                           iteration
                        PSUBGNAM = RPT_%SUBGNAM

                    END IF            ! If subgrid selected on current line

C.................  Store the raw information for the groups
                CASE( 'STORE' )

C.....................  If line is a group entry
                    IF( LIN_GROUP ) THEN

C.........................  Get number of packet type 
                        RCNT = PKT_COUNT( PKT_IDX )

C.........................  Get records count for current packet
                        N = GRPNRECS

C.........................  Choose group type
                        SELECT CASE( PKT_IDX )

C.........................  Store region 
                        CASE( REG_IDX )

C.............................  Store group label
                            IF( LIN_DEFGRP ) REGNLABL( RCNT )= GRP_LABEL

C.............................  Make sure region code is an integer
                            IF( CHKINT( SEGMENT( 1 ) ) ) THEN
                                REGRAW ( N,RCNT )= STR2INT( SEGMENT(1) )
                                REGSTAT( N,TCNT )= GRP_INCLSTAT

C.............................  Give an error if code not an integer
                            ELSE
                                EFLAG = .TRUE.
                                WRITE( MESG,94010 )
     &                                 'ERROR: Region code not an ' //
     &                                 'integer in group definition ' //
     &                                 'at line', IREC
                                CALL M3MESG( MESG )

                            END IF
                            
C.........................  Store subgrid label in local array
                        CASE( SBG_IDX )

C.............................  Store group label
                            IF( LIN_DEFGRP ) SUBGLABL( RCNT )= GRP_LABEL

C.............................  Make sure subgrid entries are specified properly
                            NINCL = 0
                            CALL PARSE_SUBGRID( BUFFER, NGRID, 
     &                                          LCELSTAT, NINCL )

C.............................  If the string could be parse as a subgrid, then
C                               store unparsed string
                            IF( NINCL .GT. 0 ) THEN

                                SBGRAW ( N,RCNT ) = BUFFER
                                SBGSTAT( N,RCNT ) = GRP_INCLSTAT

                            END IF

                        END SELECT

C.....................  If not in group, is this line a valid Select-specified
C                       region?
                    ELSE IF( K1 .GT. 0 ) THEN
c note: STOPPED HERE
c note: Need to define K1 and K2 using stored line numbers.

C.....................  If not in group, is this line a valid Select-specified
C                       subgrid?
                    ELSE IF( K2 .GT. 0 ) THEN

                    END IF  ! If not group entry or SELECT-specified

C.................  No default case for calling internal subprogram
                DEFAULT
                    L = LEN_TRIM( READTYPE )
                    MESG = 'INTERNAL ERROR: Can not call READ_GROUPS '//
     &                     'subprogram with type "' // READTYPE( 1:L )//
     &                     '"'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            
                END SELECT
                  
            END DO

C.............  Successful completion of subprogram
            RETURN

C.............  Problem(s) reading input file...
999         WRITE( MESG,94010 ) 'INTERNAL ERROR: Unexpected end of ' //
     &             'file at line', IREC
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
 
C......................  FORMAT  STATEMENTS   ..........................

C...........   Formatted file I/O formats............ 93xxx
93000       FORMAT( A )

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I10, :, 1X ) )

            END SUBROUTINE READ_GROUPS

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This subprogram compares a region code to the valid country,
C               state, and county codes and sets corresponding entries in a 
C               logical array aligning with the CNTYCOD array
            SUBROUTINE CHECK_REGIONS( REGN, LREGSTAT, STATUS )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: REGN                ! region code
            INTEGER, INTENT (OUT):: LREGSTAT( NCOUNTY ) ! true: region selectd
            INTEGER, INTENT (OUT):: STATUS              ! exit status

C.............  Local variables
            INTEGER  K, L, N     ! counters and indices

            INTEGER  LEVEL       ! sub-region level code
            INTEGER  RCHK        ! region code for comparison

            CHARACTER*300 MESG   ! mesg buffer
C----------------------------------------------------------------------

C.............  Initialize exit status
            STATUS = 0
c note: The country codes read in my rdstcy are not in 6-digit format, but the
c    n: state and county codes are.  This should be corrected to be consistent, 
c    n: and the places where the country codes are used should be updated
c    n: accordingly.    
C.............  Find in country list                      
            IF( MOD( REGN,100000 ) .EQ. 0 ) THEN
                K = FIND1( REGN, NCOUNTRY, CTRYCOD )
                LEVEL = 1

C.............  Find in state list
            ELSE IF( MOD( FIP,1000 ) .EQ. 0 ) THEN
                K = FIND1( REGN, NSTATE, STATCOD )
                LEVEL = 2

C.............  Find in county list                      
            ELSE
                K = FIND1( REGN, NCOUNTY, CNTYCOD )
                LEVEL = 3

            END IF
                            
C.............  If cy/st/co code matches inventory, set status of codes in
C               input array
            IF( K .GT. 0 ) THEN

                DO N = 1, NCOUNTY                    

                    SELECT CASE( LEVEL )
                    CASE( 1 )
                        RCHK = ( CNTYCOD( N ) / 100000 ) * 100000
                        
                    CASE( 2 )
                        RCHK = ( CNTYCOD( N ) / 1000 ) * 1000

                    CASE( 3 )
                        RCHK = CNTYCOD( N )

                    END SELECT

                    LREGSTAT( N ) = ( REGN .EQ. RCHK )

                END DO

C.............  If cy/st/co code does not match inventory...
            ELSE

                STATUS = 1

            END IF

            RETURN
                               
            END SUBROUTINE CHECK_REGIONS

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This subprogram compares interprets a subgrid definition
C               and sets a logical array with one record per grid cell
C               to true for the records in the subgrid.
C.............  NINCL is > 0 when subgrid definition is valid
C.............  LCELSTAT must be initialized by caller    
            SUBROUTINE PARSE_SUBGRID( SGRANGE, NGRID, LCELSTAT, NINCL )

C.............  External subroutines
            LOGICAL   CHKINT
            INTEGER   STR2INT

            EXTERNAL  CHKINT, STR2INT

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: SGRANGE           ! ASCII cell range
            INTEGER     , INTENT (IN) :: NGRID             ! no. cells
            INTEGER     , INTENT (OUT):: LCELSTAT( NGRID ) ! true: cell selected
            INTEGER     , INTENT (OUT):: NINCL             ! no. cells selected

C.............  Local variables
            INTEGER  C, C1, I, J, L, L1, L2     ! counters and indices

            INTEGER  X1, Y1        ! starting coordinate
            INTEGER  X2, Y2        ! ending corrdinate

            CHARACTER*100 BUFFER   ! input buffer
            CHARACTER*100 XBUF     ! tmp x-cell buffer
            CHARACTER*100 YBUF     ! tmp y-cell buffer
            CHARACTER*300 MESG     ! mesg buffer

C----------------------------------------------------------------------

C.............  Initialize number of valid cells
            NINCL = 0

C.............  Transfer range to upper case
            BUFFER = SGRANGE
            CALL UPCASE( BUFFER )

C.............  Find the first set of cell numbers
            L1 = INDEX( BUFFER, '(' )
            L2 = INDEX( BUFFER, ')' )
            C1 = INDEX( BUFFER, ',' )

C.............  Extract the cell positions and compare to the valid ranges
C.............  Give errors if bad entries for first corrdinate
            IF( L1 .GT. 0 .AND. 
     &          L2 .GT. 0 .AND.
     &          C1 .GT. 0       ) THEN

                XBUF = BUFFER( L1+1:C1-1 )
                YBUF = BUFFER( C1+1:L2-1 )

C.................  Ensure x-coordinate buffer is an integer
                IF( CHKINT( XBUF ) ) THEN
                    X1 = STR2INT( XBUF )

c NOTE: routine would need to give a WARNING when cell ranges are inconsistent 
c    n:   with grid definition, but only the first time!
C.....................  Ensure x-coordinate range is valid
C.....................  Check minimum value and reset of out of range
                    IF( X1 .LT. 1 ) THEN
                        L = LEN_TRIM( XBUF )
                        WRITE( MESG,94010 )
     &                         'WARNING: resetting x-coordinate "' //
     &                         XBUF( 1:L ) // '" at input line', IREC, 
     &                         'to minimum value of 1.'
                        CALL M3MESG( MESG ) 
                        X1 = 1

C.....................  Check maximum value and reset of out of range
                    ELSE IF( X1 .GT. NCOLS ) THEN
                        L = LEN_TRIM( XBUF )
                        WRITE( MESG,94010 )
     &                         'WARNING: resetting x-coordinate "' // 
     &                         XBUF( 1:L ) // '" at input line', IREC,
     &                         'to grid maximum of', NCOLS
                        CALL M3MESG( MESG ) 
                        X1 = NCOLS
                        
                    END IF

C.................  Give error if value is not an integer
                ELSE
                    EFLAG = .TRUE.
                    L = LEN_TRIM( XBUF )
                    WRITE( MESG,94010 )
     &                     'ERROR: Bad x-coordinate "' // XBUF( 1:L ) //
     &                     '" in subgrid definition at line', IREC
                    CALL M3MESG( MESG )

                ENDIF

C.................  Ensure y-coordinate buffer is an integer
                IF( CHKINT( YBUF ) ) THEN
                    Y1 = STR2INT( YBUF )

C.....................  Ensure y-coordinate range is valid
C.....................  Check minimum value and reset of out of range
                    IF( Y1 .LT. 1 ) THEN
                        L = LEN_TRIM( YBUF )
                        WRITE( MESG,94010 )
     &                         'WARNING: resetting y-coordinate "' //
     &                         YBUF( 1:L ) // '" at input line', IREC, 
     &                         'to minimum value of 1.'
                        CALL M3MESG( MESG ) 
                        Y1 = 1

C.....................  Check maximum value and reset of out of range
                    ELSE IF( Y1 .GT. NROWS ) THEN
                        L = LEN_TRIM( YBUF )
                        WRITE( MESG,94010 )
     &                         'WARNING: resetting y-coordinate "' // 
     &                         YBUF( 1:L ) // '" at input line', IREC,
     &                         'to grid maximum of', NROWS
                        CALL M3MESG( MESG ) 
                        Y1 = NROWS
                        
                    END IF

C.................  Give error if value is not an integer
                ELSE
                    EFLAG = .TRUE.
                    L = LEN_TRIM( YBUF )
                    WRITE( MESG,94010 )
     &                     'ERROR: Bad y-coordinate "' // YBUF( 1:L ) //
     &                     '" in subgrid definition at line', IREC
                    CALL M3MESG( MESG )

                END IF

C.............  Otherwise, bad entry should be ignored 
            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Could not find starting ' //
     &                 'coordinate of subgrid at line', IREC
                CALL M3MESG( MESG )

            END IF

C.............  Find the "TO" divider
            L1 = INDEX( BUFFER, 'TO' )

C.............  If "TO" is not found, bad entry should be ignored 
            IF( L1 .LE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Could not find "TO" ' //
     &                 'separator of subgrid definition at line', IREC
                CALL M3MESG( MESG )

            END IF

C.............  Remove the first cells and find the second ones
            L2 = LEN_TRIM( BUFFER )
            BUFFER = BUFFER( L1+2: L2 )

            L1 = INDEX( BUFFER, '(' )
            L2 = INDEX( BUFFER, ')' ) 
            C1 = INDEX( BUFFER, ',' )
            IF( L1 .GT. 0 .AND. 
     &          L2 .GT. 0 .AND.
     &          C1 .GT. 0       ) THEN

                XBUF = BUFFER( L1+1:C1-1 )
                YBUF = BUFFER( C1+1:L2-1 )

C.................  Ensure x-coordinate buffer is an integer
                IF( CHKINT( XBUF ) ) THEN
                    X2 = STR2INT( XBUF )

C.....................  Ensure x-coordinate range is valid
C.....................  Check minimum value and reset of out of range
                    IF( X2 .LT. 1 ) THEN
                        L = LEN_TRIM( XBUF )
                        WRITE( MESG,94010 )
     &                         'WARNING: resetting x-coordinate "' //
     &                         XBUF( 1:L ) // '" at input line', IREC, 
     &                         'to minimum value of 1.'
                        CALL M3MESG( MESG ) 
                        X2 = 1

C.....................  Check maximum value and reset of out of range
                    ELSE IF( X2 .GT. NCOLS ) THEN
                        L = LEN_TRIM( XBUF )
                        WRITE( MESG,94010 )
     &                         'WARNING: resetting x-coordinate "' // 
     &                         XBUF( 1:L ) // '" at input line', IREC,
     &                         'to grid maximum of', NCOLS
                        CALL M3MESG( MESG ) 
                        X2 = NCOLS
                        
                    END IF

C.................  Give error if value is not an integer
                ELSE
                    EFLAG = .TRUE.
                    L = LEN_TRIM( XBUF )
                    WRITE( MESG,94010 )
     &                     'ERROR: Bad x-coordinate "' // XBUF( 1:L ) //
     &                     '" in subgrid definition at line', IREC
                    CALL M3MESG( MESG )

                END IF

C.................  Ensure y-coordinate buffer is an integer
                IF( CHKINT( YBUF ) ) THEN
                    Y2 = STR2INT( YBUF )

C.....................  Ensure y-coordinate range is valid
C.....................  Check minimum value and reset of out of range
                    IF( Y2 .LT. 1 ) THEN
                        L = LEN_TRIM( YBUF )
                        WRITE( MESG,94010 )
     &                         'WARNING: resetting y-coordinate "' //
     &                         YBUF( 1:L ) // '" at input line', IREC, 
     &                         'to minimum value of 1.'
                        CALL M3MESG( MESG ) 
                        Y2 = 1

C.....................  Check maximum value and reset of out of range
                    ELSE IF( Y2 .GT. NROWS ) THEN
                        L = LEN_TRIM( YBUF )
                        WRITE( MESG,94010 )
     &                         'WARNING: resetting y-coordinate "' // 
     &                         YBUF( 1:L ) // '" at input line', IREC,
     &                         'to grid maximum of', NROWS
                        CALL M3MESG( MESG ) 
                        Y2 = NROWS
                        
                    END IF

                ELSE
                    EFLAG = .TRUE.
                    L = LEN_TRIM( YBUF )
                    WRITE( MESG,94010 )
     &                     'ERROR: Bad y-coordinate "' // YBUF( 1:L ) //
     &                     '" in subgrid definition at line', IREC
                    CALL M3MESG( MESG )

                END IF

C.............  If coordinate badly formed, bad entry should be ignored 
            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Could not find end ' //
     &                 'coordinate of subgrid at line', IREC
                CALL M3MESG( MESG )

            END IF

C.............  If error found, return with NINCL = 0
            IF( EFLAG ) THEN
                RETURN
            END IF

C.............  Now loop through the x-cell range and y-cell range and set
C               the status of the grid cells accordingly
            DO J = Y1, Y2
                DO I = X1, X2

                    C = ( J - 1 ) * NCOLS + X1
                    CELSTAT( C ) = .TRUE.
                    NINCL = NINCL + 1

                END DO
            END DO

            RETURN
                               
            END SUBROUTINE PARSE_SUBGRID

C**********************  FORMAT  STATEMENTS   **************************

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE RDGRPS

