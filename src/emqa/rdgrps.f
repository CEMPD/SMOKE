
        SUBROUTINE RDGRPS( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     The RDGRPS routine reads the groups from the REPCONFIG file, interprets
C     the entries, and stores the fully detailed group information in a table
C     for each type of group.
C      - Subgrid will store cell numbers that are included
C      - Region will store Regions that are excluded.
CC
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

C.........  MODULES for public variables
C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: NALLPCKT, NLINE_RC, REGNNAM, SUBGNAM,
     &                      NREGNGRP, NREGRAW, REG_IDX, NSUBGRID,
     &                      NSBGRAW, SBG_IDX, MXGRPREC, EXCLDRGN,
     &                      VALIDCEL, LENLAB3, PKT_IDX, LIN_DEFGRP,
     &                      GRP_LABEL, INREPORT, LIN_SUBREGN,
     &                      LIN_GROUP, GRPNRECS, GRP_INCLSTAT, RPT_IDX,
     &                      ALLRPT, RPT_, NREGREC, PKTCOUNT, LIN_SUBGRID

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, NCOLS, NROWS

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTRY, NSTATE, NCOUNTY,
     &                     CTRYCOD, STATCOD, CNTYCOD

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL        BLKORCMT
        CHARACTER(2)   CRLF

        EXTERNAL   BLKORCMT, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV       ! output file unit number

C...........   Local allocatable arrays... 

C...........   Region group input allocatable arrays
        INTEGER, ALLOCATABLE :: REGNREC ( : )    ! no. records per group
        CHARACTER(FIPLEN3), ALLOCATABLE :: REGRAW  ( :,: )  ! raw codes from input file
        INTEGER, ALLOCATABLE :: REGTYPE ( :,: )  ! 0=bad,1=country,2=st,3=county
        LOGICAL, ALLOCATABLE :: REGSTAT ( :,: )  ! raw status (true=include)

C...........   Subgrid input allocatable arrays
        INTEGER,        ALLOCATABLE :: SBGNREC ( : )   ! no. records per subgrid

        CHARACTER(100), ALLOCATABLE :: SBGRAW  ( :,: ) ! raw info from input file
        LOGICAL,        ALLOCATABLE :: SBGSTAT ( :,: ) ! raw status (true=incld)

C...........   Per grid-cell local allocatable arrays
        LOGICAL, ALLOCATABLE :: LCELSTAT( : )
        LOGICAL, ALLOCATABLE :: LCEL    ( : )

C...........   Per county allocatable arrays
        LOGICAL, ALLOCATABLE :: LRGN( : )     ! true: county included

C...........   Per input line local allocatable arrays
        INTEGER, ALLOCATABLE :: LINECODE( : ) ! 1= in-line region; 2= in-line subgrid

C...........   Other local arrays
        INTEGER   NCNT( NALLPCKT )   ! Count of groups defined by SELECT 
        
C...........   Other local variables
        INTEGER         C, I, J, K, N          ! counters and indices

        INTEGER         IC                  ! tmp partial region code
        INTEGER         IOS                 ! i/o status
        INTEGER         IREC                ! line number
        INTEGER         LEVEL               ! match level for region groups
        INTEGER      :: NINCL = 0           ! tmp number of included cells

        LOGICAL      :: EFLAG = .FALSE.     ! true: error found

        CHARACTER(200)  BUFFER   ! tmp label buffer
        CHARACTER(300)  MESG     ! tmp message buffer
        CHARACTER(FIPLEN3) CREGN ! tmp partial region code

        CHARACTER(16) :: PROGNAME = 'RDGRPS' ! program name

C***********************************************************************
C   begin body of subroutine RDGRPS

C.........  Allocate memory for reading and storing raw information from 
C           the input file.
        ALLOCATE( REGNNAM( NREGRAW ), STAT=IOS )
        CALL CHECKMEM( IOS, 'REGNNAM', PROGNAME )
        ALLOCATE( SUBGNAM( NSBGRAW ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUBGNAM', PROGNAME )
        ALLOCATE( LINECODE( NLINE_RC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LINECODE', PROGNAME )
        ALLOCATE( LCELSTAT( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LCELSTAT', PROGNAME )
        ALLOCATE( LCEL( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LCEL', PROGNAME )
        ALLOCATE( LRGN( NCOUNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LRGN', PROGNAME )
        
        REGNNAM = ' '     ! array
        SUBGNAM = ' '     ! array
        LINECODE = 0       ! array
        LCELSTAT = .FALSE. ! array
        LCEL     = .FALSE. ! array
        LRGN     = .FALSE. ! array

C.........  Read in defined group labels
        CALL READ_GROUPS( FDEV, NGRID, 'DEFINED LABELS', LCELSTAT )

C.........  Initialize number of groups of each type defined with SELECT
        NCNT = 0   ! array
        
C.........  Read in inline group labels and compare to defined groups.  If not
C           defined, try to match as country, state, or county and store
C           additional names.
        CALL READ_GROUPS( FDEV, NGRID, 'SELECT LABELS', LCELSTAT )

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
            DEALLOCATE( REGNNAM )
            ALLOCATE( REGNNAM( NREGNGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REGNNAM', PROGNAME )
            ALLOCATE( REGNREC( NREGNGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REGNREC', PROGNAME )
            ALLOCATE( REGRAW( MXGRPREC,NREGNGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REGRAW', PROGNAME )
            ALLOCATE( REGSTAT( MXGRPREC,NREGNGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REGSTAT', PROGNAME )
            ALLOCATE( REGTYPE( MXGRPREC,NREGNGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REGTYPE', PROGNAME )

            ALLOCATE( NREGREC( NREGNGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NREGREC', PROGNAME )
            ALLOCATE( EXCLDRGN( NCOUNTY,NREGNGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EXCLDRGN', PROGNAME )

            REGNNAM  = ' '     ! array
            REGNREC  = 0       ! array
            REGRAW   = ' '     ! array
            REGSTAT  = .TRUE.  ! array (default is include)
            REGTYPE  = 0       ! array (default is record invalid)
            EXCLDRGN = ' '     ! array 

        END IF

C.........  Same notes as above, but do for subgrids.
        IF( NSUBGRID .GT. 0 ) THEN
            DEALLOCATE( SUBGNAM )
            ALLOCATE( SUBGNAM( NSUBGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SUBGNAM', PROGNAME )
            ALLOCATE( SBGNREC( NSUBGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SBGNREC', PROGNAME )
            ALLOCATE( SBGRAW( MXGRPREC,NSUBGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SBGRAW', PROGNAME )
            ALLOCATE( SBGSTAT( MXGRPREC,NSUBGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SBGSTAT', PROGNAME )
            ALLOCATE( VALIDCEL( NGRID,NSUBGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VALIDCEL', PROGNAME )

            SUBGNAM  = ' '     ! array
            SBGNREC  = 0       ! array
            SBGRAW   = ' '     ! 2-d array
            SBGSTAT  = .TRUE.  ! 2-d array (initialize to "include")
            VALIDCEL = 0       ! 2-d array

        END IF

C.........  Read raw group information.  "Raw" means the entries as they
C           appear in the input file, but not converted to the data structures
C           needed for further processing.
        CALL READ_GROUPS( FDEV, NGRID, 'STORE', LCELSTAT )

C.........  Convert region group entries to data structures needed for
C           further processing
        IF( NREGNGRP .GT. 0 ) THEN

C.............  Loop through different regions
            DO N = 1, PKTCOUNT( REG_IDX )

C.................  Initialize region list indicator depending on first subgrid
C                   entry
                IF( REGSTAT( 1,N ) ) THEN   ! include
                    LRGN = .FALSE.          ! array
                ELSE                        ! exclude
                    LRGN = .TRUE.           ! array
                END IF

C.................  Loop through records in each packet or in-line region group
                DO I = 1, REGNREC( N )

C.....................  Depending on type of the record (country, state, or
C                       county), set comparison length
                    SELECT CASE( REGTYPE( I,N ) )
                    CASE( 1 )                     ! country
                        IC = FIPEXPLEN3 + 1
                    CASE( 2 )                     ! state
                        IC = STALEN3
                    CASE( 3 )                     ! county
                        IC = FIPLEN3
                    CASE DEFAULT
                        CYCLE

                    END SELECT
                    
                    CREGN = REGRAW( I,N )( 1:IC )

C..................... Loop through counties and determine which ones 
                    DO J = 1, NCOUNTY
                    
                        IF( CNTYCOD( J )( 1:IC ) == CREGN ) THEN
                            LRGN( J ) = REGSTAT( I,N )
                        END IF

                    END DO  ! End loop on counties

                END DO      ! End loop on entries in region group

C.................  Create list of counties to exclude from group
                K = 0
                DO J = 1, NCOUNTY

                    IF( .NOT. LRGN( J ) ) THEN
                        K = K + 1
                        EXCLDRGN( K,N ) = CNTYCOD( J )
                    END IF
                    NREGREC( N ) = K

                END DO      ! End loop on counties

            END DO          ! End loop on region groups

        END IF

C.........  Convert subgrid entries to data structures needed for further 
C           processing.
        IF( NSUBGRID .GT. 0 ) THEN

C.............  Loop through different subgrids
            DO N = 1, PKTCOUNT( SBG_IDX )

C.................  Initialize cell list indicator depending on first subgrid
C                   entry
                IF( SBGSTAT( 1,N ) ) THEN   ! include
                    LCEL = .FALSE. 
                ELSE                        ! exclude
                    LCEL = .TRUE. 
                END IF

C.................  Loop through records in each packet or in-line subgrid
                DO I = 1, SBGNREC( N )

C.....................  Determine which cells current entry applies
                    BUFFER = SBGRAW( I,N )
                    NINCL = 0
                    LCELSTAT = .FALSE.
                    CALL PARSE_SUBGRID( BUFFER, NGRID, LCELSTAT, NINCL )

C.....................  If current entry is an include, then include cells
                    IF( SBGSTAT( I,N ) ) THEN  ! include

                        DO C = 1, NGRID
                            IF( LCELSTAT( C ) ) LCEL( C ) = .TRUE.
                        END DO

C.....................  If current entry is an exclude, then exclude cells
                    ELSE                       ! exclude

                        DO C = 1, NGRID
                            IF( LCELSTAT( C ) ) LCEL( C ) = .FALSE.
                        END DO

                    END IF

                END DO

C.................  Based on global cell status, create list of valid cells
                J = 0
                DO C = 1, NGRID
                    IF( LCEL( C ) ) THEN
                        J = J + 1
                        VALIDCEL( J,N ) = C
                    END IF
                END DO

            END DO

        END IF 
            
C.........  Deallocate raw group information
        IF( ALLOCATED( REGRAW ) ) DEALLOCATE( REGNREC, REGRAW, 
     &                                        REGSTAT, REGTYPE )
        IF( ALLOCATED( SBGRAW ) ) DEALLOCATE( SBGRAW, SBGSTAT )
        IF( ALLOCATED( LINECODE ) ) DEALLOCATE( LINECODE )
        IF( ALLOCATED( LCELSTAT ) ) DEALLOCATE( LCELSTAT, LCEL )
        IF( ALLOCATED( LRGN ) ) DEALLOCATE( LRGN )

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
            SUBROUTINE READ_GROUPS( FDEV, NGRID, READTYPE, LGRDSTAT )

C.............  External functions
            LOGICAL      CHKINT
            CHARACTER(2) CRLF
            INTEGER      GETNLIST
            INTEGER      INDEX1
            INTEGER      STR2INT

            EXTERNAL     CHKINT, CRLF, GETNLIST, INDEX1, STR2INT

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: FDEV       ! input file unit
            INTEGER     , INTENT (IN) :: NGRID      ! no. grid cells
            CHARACTER(*), INTENT (IN) :: READTYPE   ! Reading type
            LOGICAL     , INTENT (IN) :: LGRDSTAT( NGRID ) ! true: report cell

C.............   Local parameters   
            INTEGER, PARAMETER :: MXSEG = 100

C.............  Subprogram local arrays
            CHARACTER(256) SEGMENT( MXSEG )

C.............  Local variables
            INTEGER       I, J, L, N       ! counters and indices

            INTEGER       IOS      ! i/o status
            INTEGER    :: NS = 1   ! no. segments in line
            INTEGER       RCNT     ! record count

            CHARACTER(FIPLEN3) CFIP ! tmp region code
            CHARACTER(300) BUFFER   ! tmp line buffer as uppercase
            CHARACTER(300) LINE     ! tmp line buffer
            CHARACTER(300) MESG     ! mesg buffer

            CHARACTER(LENLAB3) :: PREGNNAM   ! previous region name
            CHARACTER(LENLAB3) :: PSUBGNAM   ! previous subgrid name

C----------------------------------------------------------------------

C.............  Rewind input file
            REWIND( FDEV )

C.............  Loop though file to store local array of labeled group names
            SEGMENT  = ' '     ! array
            PREGNNAM = ' '
            PSUBGNAM = ' '
            IREC = 0
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
                    RCNT = PKTCOUNT( PKT_IDX )
                
                    SELECT CASE( PKT_IDX )

C.....................  Store region group label
                    CASE( REG_IDX )
                        IF( LIN_DEFGRP ) REGNNAM( RCNT ) = GRP_LABEL
                
C.....................  Store subgrid label in local array
                    CASE( SBG_IDX )
                        IF( LIN_DEFGRP ) SUBGNAM( RCNT ) = GRP_LABEL
                
                    END SELECT

C.................  Count SELECT statements that do not use defined 
C                   regions.  Determine those that coorespond to valid 
C                   entries and those that should be ignored.
                CASE( 'SELECT LABELS' )
            
C.....................  Get current count of current packet 
                    RCNT = PKTCOUNT( PKT_IDX )
                
                    SELECT CASE( PKT_IDX )

C.....................  Store region group label
                    CASE( REG_IDX )
                        IF( LIN_DEFGRP ) REGNNAM( RCNT ) = GRP_LABEL
                
C.....................  Store subgrid label in local array
                    CASE( SBG_IDX )
                        IF( LIN_DEFGRP ) SUBGNAM( RCNT ) = GRP_LABEL
                
                    END SELECT

C.....................  Skip if report section not started yet.
                    IF( .NOT. INREPORT ) CYCLE

C.....................  A region is being selected
                    IF( LIN_SUBREGN ) THEN

C.........................  Search for region name in defined regions
                        J = INDEX1( RPT_%REGNNAM, NREGRAW, REGNNAM )
                    
C.........................  If region name is not found, then try to convert
C                           to a region code and compare with the inventory.
                        IF( J .LE. 0 ) THEN

C.............................  Check if label is an integer
                            IF( .TRUE. .OR. CHKINT( RPT_%REGNNAM ) ) THEN

C.................................  Convert label to region code
                                CFIP = RPT_%REGNNAM

C.................................  Check if label is a valid region code
C.................................  REGNTMP is a dummy argument at this stage
                                CALL CHECK_REGIONS( CFIP, LEVEL, IOS )
                            
C.................................  If code is valid, store line number of
C                                   record
                                IF( IOS .EQ. 0 ) THEN
                                    N = NCNT( REG_IDX ) + 1
                                    NCNT( REG_IDX ) = N
                                    LINECODE( IREC ) = 1

C.................................  Otherwise, give warning because the label
C                                   does not match and it's not an inline code
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
                    IF( LIN_SUBGRID ) THEN

C.........................  Search for subgrid names in defined subgrids
                        J = 0
                        IF( NSBGRAW .GT. 0 ) 
     &                      J = INDEX1( RPT_%SUBGNAM, NSBGRAW, SUBGNAM )
                    
C.........................  If subgrid name is not found, then try to interpret
C                           entry as a subgrid definition.
                        IF( J .LE. 0 ) THEN

C............................. Check if subgrid is defined in-line
                            NINCL    = 0
                            LCELSTAT = .FALSE.   ! array 
                            CALL PARSE_SUBGRID( RPT_%SUBGNAM, NGRID, 
     &                                          LCELSTAT, NINCL )

C.............................  If subgrid is valid, increase count and
C                               flag line as a in-line subgrid
                            IF( NINCL .GT. 0 ) THEN
                                N = NCNT( SBG_IDX ) + 1
                                NCNT( SBG_IDX ) = N
                                LINECODE( IREC ) = 2

C.............................  If subgrid is invalid
                            ELSE
                                L = LEN_TRIM( RPT_%SUBGNAM )
                                WRITE( MESG,94010 ) 
     &                             'WARNING: Subgrid definition "' //
     &                             RPT_%SUBGNAM( 1:L ) // '" at line', 
     &                             IREC, 'is not defined and cannot' //
     &                             CRLF() // BLANK10 //
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

C.....................  If line has the group label
                    IF( LIN_DEFGRP ) THEN

C.........................  Get number of packet type 
                        RCNT = PKTCOUNT( PKT_IDX )

                        SELECT CASE( PKT_IDX )

C.........................  Store region group label
                        CASE( REG_IDX )
                            IF( LIN_DEFGRP ) REGNNAM( RCNT ) = GRP_LABEL
                
C.........................  Store subgrid label in local array
                        CASE( SBG_IDX )
                            IF( LIN_DEFGRP ) SUBGNAM( RCNT ) = GRP_LABEL
                
                        END SELECT

C.....................  If line is inside a group packet
                    ELSE IF( LIN_GROUP ) THEN

C.........................  Get number of packet type 
                        RCNT = PKTCOUNT( PKT_IDX )

C.........................  Get records count for current packet
                        N = GRPNRECS

C.........................  Choose group type
                        SELECT CASE( PKT_IDX )

C.........................  Store region 
                        CASE( REG_IDX )

C.............................  Make sure region code is an integer
                            IF( .TRUE. .OR. CHKINT( SEGMENT( 1 ) ) ) THEN
                                CFIP = SEGMENT(1)

                                CALL CHECK_REGIONS( CFIP, LEVEL, IOS )

                                IF( IOS .EQ. 0 ) THEN

                                    REGNREC( RCNT )  = N
                                    REGRAW ( N,RCNT )= CFIP
                                    REGSTAT( N,RCNT )= GRP_INCLSTAT
                                    REGTYPE( N,RCNT )= LEVEL

                                END IF

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

C.............................  Make sure subgrid entries are specified properly
                            NINCL = 0
                            LCELSTAT = .FALSE.   ! array 
                            CALL PARSE_SUBGRID( BUFFER, NGRID, 
     &                                          LCELSTAT, NINCL )

C.............................  If the string could be parse as a subgrid, then
C                               store unparsed string
                            IF( NINCL .GT. 0 ) THEN

                                SBGNREC( RCNT )   = GRPNRECS
                                SBGRAW ( N,RCNT ) = BUFFER
                                SBGSTAT( N,RCNT ) = GRP_INCLSTAT

                            END IF

                        END SELECT

C.....................  If not in group, is this line a valid Select-specified
C                       region?  If so, make sure not a duplicate in-line and
C                       store it.
                    ELSE IF( LINECODE( IREC ) .EQ. 1 ) THEN

C.........................  Convert in-line code to region code and rename it.
                        CFIP = RPT_%REGNNAM
                        CALL CHECK_REGIONS( CFIP, LEVEL, IOS )
                        CALL RENAME_REGION( RPT_%REGNNAM, CFIP )

C.........................  Get report number
                        N = PKTCOUNT( RPT_IDX )
                        ALLRPT( N )%REGNNAM = RPT_%REGNNAM

C.........................  See if this name is already stored and if not,
C                           store it.
                        RCNT = PKTCOUNT( REG_IDX )
                        J = INDEX1( RPT_%REGNNAM, RCNT, REGNNAM )

C.........................  Store without checking status because LINECODE = 1
C                           only if code has already been through CHECK_REGIONS
                        IF( J .LE. 0 ) THEN
                            RCNT = RCNT + 1
                            REGNNAM( RCNT ) = RPT_%REGNNAM
                            PKTCOUNT( REG_IDX ) = RCNT

                            REGNREC( RCNT )   = 1
                            REGRAW ( 1,RCNT ) = CFIP
                            REGSTAT( 1,RCNT ) = .TRUE.
                            REGTYPE( 1,RCNT ) = LEVEL
                        END IF

C.....................  If not in group, is this line a valid Select-specified
C                       subgrid?  If so, make sure not a duplicate in-line and
C                       store it.
                    ELSE IF( LINECODE( IREC ) .EQ. 2 ) THEN

C.........................  Rename subgrid name
                        CALL RENAME_SUBGRID( RPT_%SUBGNAM )

C.........................  Get report number
                        N = PKTCOUNT( RPT_IDX )
                        ALLRPT( N )%SUBGNAM = RPT_%SUBGNAM

C.........................  See if this name is already stored and if not,
C                           store it.
                        RCNT = PKTCOUNT( SBG_IDX )
                        J = INDEX1( RPT_%SUBGNAM, RCNT, SUBGNAM )

                        IF ( J .LE. 0 ) THEN
                            RCNT = RCNT + 1
                            SUBGNAM( RCNT ) = RPT_%SUBGNAM 
                            PKTCOUNT( SBG_IDX ) = RCNT

                            SBGNREC( RCNT )   = 1 
                            SBGRAW ( 1,RCNT ) = BUFFER
                            SBGSTAT( 1,RCNT ) = .TRUE.
                        END IF

                    END IF  ! If not group entry or SELECT-specified

C.................  No default case for calling internal subprogram
                CASE DEFAULT
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
            SUBROUTINE CHECK_REGIONS( REGN, LEVEL, STATUS )

C.............  Exernal subroutines
            INTEGER    FINDC
            EXTERNAL   FINDC

C.............  Subprogram arguments
            CHARACTER(FIPLEN3), INTENT (IN) :: REGN        ! region code
            INTEGER, INTENT (OUT):: LEVEL       ! sub-region level code
            INTEGER, INTENT (OUT):: STATUS      ! exit status

C.............  Local variables
            INTEGER  K, L, N     ! counters and indices

            INTEGER  FIP         ! tmp country/state/county code
            INTEGER  RCHK        ! region code for comparison

            CHARACTER(300) MESG  ! mesg buffer
C----------------------------------------------------------------------

C.............  Initialize output variables
            STATUS = 1
            LEVEL  = 0
            CALL PADZERO( REGN )

C.............  Find in country list                      
            IF( REGN( FIPEXPLEN3+2:FIPLEN3 ) == '00000' ) THEN
                K = FINDC( REGN, NCOUNTRY, CTRYCOD )
                STATUS = 0
                LEVEL = 1

C.............  Find in state list
            ELSE IF( REGN( STALEN3+1:FIPLEN3 ) == '000' ) THEN
                K = FINDC( REGN, NSTATE, STATCOD )
                STATUS = 0
                LEVEL = 2

C.............  Find in county list                      
            ELSE
                K = FINDC( REGN, NCOUNTY, CNTYCOD )
                STATUS = 0
                LEVEL = 3
                            
            END IF

            RETURN
                               
C......................  FORMAT  STATEMENTS   ..........................

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I10, :, 1X ) )

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
            LOGICAL     , INTENT (OUT):: LCELSTAT( NGRID ) ! true: cell selected
            INTEGER     , INTENT (OUT):: NINCL             ! no. cells selected

C.............  Local variables
            INTEGER  C, C1, I, J, L, L1, L2     ! counters and indices

            INTEGER  X1, Y1        ! starting coordinate
            INTEGER  X2, Y2        ! ending corrdinate

            CHARACTER(100) BUFFER  ! input buffer
            CHARACTER(100) XBUF    ! tmp x-cell buffer
            CHARACTER(100) YBUF    ! tmp y-cell buffer
            CHARACTER(300) MESG    ! mesg buffer

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

c note:  This could be updated to only give the warning the first time
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

                END IF

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

                    C = ( J - 1 ) * NCOLS + I
                    LCELSTAT( C ) = .TRUE.
                    NINCL = NINCL + 1

                END DO
            END DO

            RETURN
                               
C......................  FORMAT  STATEMENTS   ..........................

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I10, :, 1X ) )

            END SUBROUTINE PARSE_SUBGRID

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This subprogram renames internal region names
            SUBROUTINE RENAME_REGION( NAM, CFIP )

            CHARACTER(*), INTENT (IN OUT) :: NAM
            CHARACTER(FIPLEN3), INTENT (IN    ) :: CFIP

C----------------------------------------------------------------------

            WRITE( NAM, '(A,I12.12)' ) 'In-line region ', CFIP

            END SUBROUTINE RENAME_REGION

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This subprogram renames internal subgrid names
            SUBROUTINE RENAME_SUBGRID( NAM )

            CHARACTER(*), INTENT (IN OUT) :: NAM
            CHARACTER(LENLAB3) TMPNAM

            INTEGER L

C----------------------------------------------------------------------
            L = LEN_TRIM( NAM )
            TMPNAM = 'In-line region ' // NAM( 1:L )
            NAM = TMPNAM

            END SUBROUTINE RENAME_SUBGRID

        END SUBROUTINE RDGRPS

