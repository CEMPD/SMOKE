

        SUBROUTINE PRCLINRC( IREC, NSEGS, LINE, SEGMENT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C    The PRCLINRC processes the REPCONFIG instructions by setting variables
C    that can be used by the calling routine(s) for futher processing.
C    This routine is used for several cases.  It is be able to interact
C    with various other routines at various stages of processing. It resolves 
C    any implies report settings that affect input file requirements (e.g., 
C    BY CELL implies GRIDDING). Settings are passed to other routines with 
C    modules.
C
C  PRECONDITIONS REQUIRED:
C    REPCONFIG file is opened
C    LINE read, separated into SEGMENTS, and left-justified
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
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2   CRLF
        INTEGER       INDEX1
        INTEGER       STR2INT

        EXTERNAL   CRLF, INDEX1, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: IREC        ! line counter
        INTEGER     , INTENT (IN) :: NSEGS       ! no. line segments
        CHARACTER(*), INTENT (IN) :: LINE        ! full input record less commnt
        CHARACTER(*), INTENT (IN) :: SEGMENT( * )! parsed input record

C...........   Other local variables
        INTEGER          I, J, L, L2            ! counters and indices
        INTEGER          IOS                 ! i/o status

        INTEGER, SAVE :: IREC_MIN = 99999999 ! minimum record number
        INTEGER, SAVE :: PREC     = -9       ! previous call record number

        LOGICAL, SAVE :: FIRSTIME  = .TRUE.   ! true: first time routine called
        LOGICAL, SAVE :: FIRSTLOOP = .TRUE.   ! true: file already read once
        LOGICAL, SAVE :: LCATSET   = .FALSE.  ! true: source category set

        CHARACTER*1      O3SYN             !  Y or N for ozone season
        CHARACTER*6      FILNUM            !  tmp file number string
        CHARACTER*300    MESG              !  message buffer

        CHARACTER*16 :: PROGNAME = 'PRCLINRC' ! program name

C***********************************************************************
C   begin body of subroutine PRCLINRC

C.........  Initialize for start of file
        IF( FIRSTIME .OR. IREC .EQ. IREC_MIN ) THEN

            GRPNRECS  = 0 
            INPACKET  = .FALSE.
            INGROUP   = .FALSE.
            INREPORT  = .FALSE.
            INSPCIFY  = .FALSE.

            PKTCOUNT  = 0         ! array
            PKTSTATUS = .FALSE.   ! array

            RPT_%O3SEASON = .FALSE.  ! default to not use ozone season data
            RPT_%OUTTIME  = 230000   ! default to output at 2300 hours
            RPT_%DELIM    = ';'      ! default to semi-colon

            FIRSTIME = .FALSE.

        END IF

C.........  Section for resetting for start of file (only access for read loops
C           2 and higher)
        IF( IREC .LT. PREC ) THEN
            FIRSTLOOP = .FALSE.

C.............  As a double check, make sure that the source category has been
C               set
            IF( .NOT. LCATSET ) THEN
                RC_ERROR = .TRUE.
                MESG = 'ERROR: Source category not set in reports ' //
     &                 'configuration file. Use SMK_SOURCE instruction.'
                CALL M3MSG2( MESG )
            END IF

        END IF     

C.........  Section for first loop
        IF( FIRSTLOOP ) THEN

C.............  Check for source category
C.............  Warning if it is reset to something else
            IF( LCATSET .AND. SEGMENT( 1 ) .EQ. 'SMK_SOURCE' ) THEN

                IF( SEGMENT( 2 ) .NE. CATEGORY ) THEN
                    WRITE( MESG,94010 ) 'Source category reset ' //
     &                'attempted at line', IREC, ', but was ignored.'
                    CALL M3MSG2( MESG )

                END IF

                GO TO 999    ! to end of routine

C.............  Store source category setting
            ELSE IF( SEGMENT( 1 ) .EQ. 'SMK_SOURCE' ) THEN

                LCATSET = .TRUE.

                IF( SEGMENT( 2 ) .EQ. 'A' ) THEN        ! area
		    MINC     = 3
                    CRL      = 'A'
                    CATEGORY = 'AREA'
                    CATDESC  = 'Area'

                ELSE IF( SEGMENT( 2 ) .EQ. 'M' ) THEN   ! mobile
		    MINC     = 2
                    CRL      = 'M'
                    CATEGORY = 'MOBILE'
                    CATDESC  = 'Mobile'

                ELSE IF( SEGMENT( 2 ) .EQ. 'P' ) THEN   ! point
		    MINC     = 2
                    CRL      = 'P'
                    CATEGORY = 'POINT'
                    CATDESC  = 'Point'

                ELSE                                    ! error
                    RC_ERROR = .TRUE.
                    LCATSET  = .FALSE.
                    L = LEN_TRIM( SEGMENT( 2 ) )
                    WRITE( MESG,94010 ) 'ERROR: Source category ' //
     &                'setting "' // SEGMENT( 2 )( 1:L ) // 
     &                '" at line', IREC, 'is invalid. Must be ' //
     &                '"A", "M", or "P".'
                    CALL M3MSG2( MESG )

                END IF

                GO TO 999    ! to end of routine

            END IF     ! End setting of source category

        END IF         ! End section for first loop through input file

C.........  Initializations needed for every line
        LIN_DEFGRP  = .FALSE.
        LIN_GROUP   = .FALSE.
        LIN_SUBDATA = .FALSE.
        LIN_SUBGRID = .FALSE.
        LIN_SUBREGN = .FALSE.
        LIN_TITLE   = .FALSE.
        LIN_UNIT    = .FALSE.
        LIN_SPCIFY  = .FALSE.

C.........  Set length of full line
        L2 = LEN_TRIM( LINE )

C.........  If a packet is not currently being read, determine what packet
C           starts on the current line
C.........  Will skip line if no packet has been found and the line does not
C           start with "/"
        IF( .NOT. INPACKET ) THEN

            J = INDEX( LINE( 2:L2 ), '/' )

C.............  Error if packet started but does not end with a slash
            IF( LINE( 1:1 ) .EQ. '/' .AND.
     &                J     .LE. 0         ) THEN

                RC_ERROR = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Packet at line', IREC,
     &                 'started but not finished.'
                CALL M3MSG2( MESG )

C.............  Packet found. Interpret it...
            ELSE IF( J .GT. 0 ) THEN
                J = J + 1

C.................  Set packet name and extract group label
                PCKTNAM = LINE( 1:J )

C.................  Find packet name in list of valid names
                PKT_IDX = INDEX1( PCKTNAM, NALLPCKT, ALLPCKTS )

C.................  Warning if packet is not recognized
                IF( PKT_IDX .LE. 0 ) THEN
                    INPACKET = .TRUE.
                    WRITE( MESG,94010 ) 'WARNING: Packet "' //
     &                     PCKTNAM( 1:J ) // '" in report configuration'
     &                     // CRLF() // BLANK10 //' file at line', IREC,
     &                     'is not recognized. Packet will be ignored.'
                    CALL M3MESG( MESG )

C.................  Packet recognized
                ELSE

C.....................  Count the number of groups and the number of reports
                    PKTCOUNT ( PKT_IDX ) = PKTCOUNT( PKT_IDX ) + 1
                    PKTSTATUS( PKT_IDX ) = .TRUE.

                    INPACKET = .TRUE.

C.....................  Set first line number of packet
                    PKTSTART = IREC

C.....................  Packet-specific processing
                    SELECT CASE( PKT_IDX )

C.....................  A reporting time packet
                    CASE( TIM_IDX )

                        IF( J+1 .GE. L2 ) THEN
C.........................  Ensure reporting time is set
                            CALL NO_SETTING_FOUND( IREC, PKT_IDX )
                        ELSE                            
                            RPT_%OUTTIME = STR2INT( LINE( J+1:L2 ) )
                        END IF

                        PKTEND   = IREC
                        INPACKET = .FALSE.             ! end implied

C.....................  A ozone season data usage packet
                    CASE( O3S_IDX )

C.........................  Ensure ozone season setting made
                        IF( J+1 .GE. L2 ) THEN
                            CALL NO_SETTING_FOUND( IREC, PKT_IDX )

                        ELSE                            
                            O3SYN = ADJUSTL( LINE( J+1:L2 ) )
                            CALL UPCASE( O3SYN )
                            IF( O3SYN .EQ. 'Y' ) THEN
                                RPT_%O3SEASON = .TRUE.

                            ELSE IF( O3SYN .EQ. 'N' ) THEN
                                RPT_%O3SEASON = .FALSE.

                            ELSE
                                L = LEN_TRIM( O3SYN )
                                WRITE( MESG,94010 ) 
     &                            'WARNING: Unrecognized /O3SEASON/ ' //
     &                            'settting "' // O3SYN( 1:L ) // 
     &                            '" at line', IREC, '. Setting to ' //
     &                            'default of FALSE.'  
                                CALL M3MSG2( MESG )
                                RPT_%O3SEASON = .FALSE.

                            END IF

                        END IF

                        PKTEND   = IREC
                        INPACKET = .FALSE.             ! end implied

C.....................  A new file packet
                    CASE( FIL_IDX )

C.........................  Set default logical/physical file name
                        WRITE( FILNUM,'(I6)' ) PKTCOUNT ( FIL_IDX )
                        FILNUM = ADJUSTL( FILNUM )
                        FIL_ONAME = 'REPORT' // FILNUM

C.........................  Extract file name (logical or physical, if any)
                        IF( J+1 .LE. L2 ) THEN
                             FIL_ONAME = ADJUSTL( LINE( J+1:L2 ) )
                        END IF

                        PKTEND   = IREC
                        INPACKET = .FALSE.             ! end implied

C.....................  A new delimiter packet
                    CASE( DEL_IDX )

C.........................  Ensure ozone season setting made
                        IF( J+1 .GE. L2 ) THEN
                            CALL NO_SETTING_FOUND( IREC, PKT_IDX )
                        ELSE                            
                            RPT_%DELIM = ADJUSTL( LINE( J+1:L2 ) )
                        END IF

                        PKTEND   = IREC
                        INPACKET = .FALSE.             ! end implied

C.....................  A region group or subgrid definition packet
C.....................  Reset group specific to defaults
                    CASE( REG_IDX, SBG_IDX )
                        INGROUP      = .TRUE.
                        GRP_INCLSTAT = .TRUE.
                        LIN_DEFGRP   = .TRUE.                        

                        IF( J+1 .GE. L2 ) THEN
                            CALL NO_SETTING_FOUND( IREC, PKT_IDX )
                        ELSE
                            GRP_LABEL = ADJUSTL( LINE( J+1:L2 ) )
                        END IF

C.....................  Specification of elevated groups, PinG sources, or
C                       elevated sources
                    CASE( ELG_IDX, PNG_IDX, ELV_IDX )
                        INSPCIFY    = .TRUE.
                        LIN_SPCIFY  = .TRUE.
                        SPCF_NOR    = 0
                        SPCF_NAND   = 0

C.....................  A report packet
                    CASE( RPT_IDX )

C.........................  If have not found a NEWFILE packet, then
C                           set the first output file name
                        IF( PKTCOUNT( FIL_IDX ) .EQ. 0 ) THEN

                            PKTCOUNT( FIL_IDX ) = 1
                            WRITE( FILNUM,'(I6)' ) PKTCOUNT ( FIL_IDX )
                            FILNUM = ADJUSTL( FILNUM )
                            FIL_ONAME = 'REPORT' // FILNUM

                        END IF

C.........................  Reset report settings to defaults
                        INREPORT      = .TRUE.
                        RPT_%BYCELL   = .FALSE.
                        RPT_%BYCNRY   = .FALSE.
                        RPT_%BYDATE   = .FALSE.
                        RPT_%BYDIU    = .FALSE.
                        RPT_%BYSTAT   = .FALSE.
                        RPT_%BYCNTY   = .FALSE.
                        RPT_%BYELEV   = .FALSE.
                        RPT_%BYHOUR   = .FALSE.
                        RPT_%BYLAYER  = .FALSE.
                        RPT_%BYMON    = .FALSE.
                        RPT_%BYSCC    = .FALSE.
                        RPT_%BYSPC    = .FALSE.
                        RPT_%BYSRC    = .FALSE.
                        RPT_%BYSRG    = .FALSE.
                        RPT_%BYCONAM  = .FALSE.
                        RPT_%BYSTNAM  = .FALSE.
                        RPT_%BYCYNAM  = .FALSE.
                        RPT_%BYRCL    = .FALSE.
                        RPT_%BYWEK    = .FALSE.
                        RPT_%LAYFRAC  = .FALSE.
                        RPT_%NORMCELL = .FALSE.
                        RPT_%NORMPOP  = .FALSE.
                        RPT_%SCCNAM   = .FALSE.
                        RPT_%SRCNAM   = .FALSE.
                        RPT_%STKPARM  = .FALSE.
                        RPT_%USEGMAT  = .FALSE.
                        RPT_%USESSMAT = .FALSE.
                        RPT_%USESLMAT = .FALSE.
                        RPT_%USEHOUR  = .FALSE.
                        LREGION       = .FALSE.
                        LSUBGRID      = .FALSE.

                        RPT_%ELEVSTAT = 0
                        RPT_%NUMDATA  = -9   ! zero is legitimate
                        RPT_%NUMTITLE = 0
                        RPT_%RENDLIN  = 0    ! init for consistency
                        RPT_%RSTARTLIN= 0    ! init for consistency
                        RPT_%SCCRES   = 10
                        RPT_%SRGRES   = 0

                        RPT_%DATAFMT  = 'E8.3'
                        RPT_%OFILENAM = ' '    ! init for consistency
                        RPT_%REGNNAM  = ' '
                        RPT_%SUBGNAM  = ' '
                        RPT_%SPCPOL   = ' '
                        TITLE         = ' '

                    END SELECT

                END IF
 
            END IF

            GO TO 999    ! to end of routine

C........  End of packet...
       ELSE IF( SEGMENT( 1 ) .EQ. '/END/' ) THEN

C.............  Resolve status of settings that cannot be determined
C               definitively until the end of the packet...

C.............  If packet was a report packet and no data have
C               been selected, then set input file global setting
            IF( INREPORT .AND. RPT_%NUMDATA .LE. 0 ) THEN
                DATAMISS = .TRUE.
            END IF

C.............  If using elevated settings and no layer fractions, make sure
C               PELV file will be read in 
            IF( .NOT. VFLAG  .AND. .NOT. RPT_%LAYFRAC .AND.
     &        ( RPT_%ELEVSTAT .GT. 0 .OR. RPT_%BYELEV )      ) THEN
                VFLAG = .TRUE.
            END IF

C.............  In in a group and no records read, give warning
            IF( INGROUP .AND. GRPNRECS .EQ. 0 .AND. FIRSTLOOP ) THEN
                L = LEN_TRIM( GRP_LABEL )
                MESG = 'WARNING: No valid records found for group "'//
     &                 GRP_LABEL( 1:L ) // '"'
                CALL M3MSG2( MESG )
            END IF

C.............  Reset all packet-specific settings
            GRPNRECS   = 0
            PKTEND     = IREC
            INPACKET   = .FALSE.
            INGROUP    = .FALSE.
            INREPORT   = .FALSE.
            INSPCIFY   = .FALSE.

            GRP_LABEL  = ' '

            GO TO 999          ! to end of routine

C........  Missing end of packet and already in a packet
       ELSE IF( LINE( 1:1 ) .EQ. '/' ) THEN
            RC_ERROR = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Missing end of packet ' //
     &             'at line', IREC
            CALL M3MSG2( MESG )

            GO TO 999          ! to end of routine

        END IF

C.........  General group processing
        IF( INGROUP ) THEN

C.............  Check the status of INCLUDE versus EXCLUDE
            IF( SEGMENT( 1 ) .EQ. 'INCLUDE' ) THEN
                IF( SEGMENT( 2 ) .NE. ' ' .AND. FIRSTLOOP ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Other data ' //
     &                'appended to INCLUDE instruction at line', IREC
                    CALL M3MSG2( MESG )
                END IF
                GRP_INCLSTAT = .TRUE.                

            ELSE IF( SEGMENT( 1 ) .EQ. 'EXCLUDE' ) THEN
                IF( SEGMENT( 2 ) .NE. ' ' .AND. FIRSTLOOP ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Other data ' //
     &                'appended to EXCLUDE instruction at line', IREC
                    CALL M3MSG2( MESG )
                END IF
                GRP_INCLSTAT = .FALSE.

C.............  If neither, count entries in the current group definition
            ELSE
                LIN_GROUP = .TRUE.
                GRPNRECS = GRPNRECS + 1

            END IF
            
        END IF      ! End group packet

C.............  Specification processing
        IF( INSPCIFY ) THEN

C.............  Determine how many conditions appear for this line
            J = NSEGS
            I = 0
            DO
                I = I + 1
                J = J - 4              ! three fields plus AND
                IF ( J .LE. 0 ) EXIT
            END DO
            SPCF_NAND = I

C.............  Increment the OR field count (one per row)
            SPCF_NOR  = SPCF_NOR + 1

        END IF

C.............  Check current line for report settings
C.............  Set global and line-by-line logical variables, shared through
C               MODREPRT
        IF( INREPORT ) THEN

            SELECT CASE( SEGMENT( 1 ) )

C.............  Gridding used for report
            CASE( 'GRIDDING' )
                GFLAG      = .TRUE.
                RPT_%USEGMAT  = .TRUE.

C.............  Speciation used for report
            CASE( 'SPECIATION' )
                IF( SEGMENT( 2 ) .EQ. 'MASS' ) THEN
                    SSFLAG      = .TRUE.
                    RPT_%USESSMAT  = .TRUE.

                ELSE IF( SEGMENT( 2 ) .EQ. 'MOLE' ) THEN
                    SLFLAG      = .TRUE.
                    RPT_%USESLMAT  = .TRUE.

                ELSE
                    IF( FIRSTLOOP ) THEN
                        L = LEN_TRIM( SEGMENT( 2 ) )
                	WRITE( MESG,94010 ) 
     &                    'WARNING: Speciation type "'// 
     &                    SEGMENT( 2 )( 1:L )// '" at line', IREC, 
     &                    'is not known.' // CRLF()// BLANK10 // 
     &                    'Will assume mass-based speciation.'
                        CALL M3MSG2( MESG )
                    END IF

                    SSFLAG = .TRUE.

                END IF

C.............  Temporal allocated emission used for report
            CASE( 'TEMPORAL' )
                TFLAG     = .TRUE.
                RPT_%USEHOUR = .TRUE.
                RPT_%BYDATE = .TRUE.

C.............  BY options affecting inputs needed
            CASE( 'BY' )

                SELECT CASE( SEGMENT( 2 ) )

                CASE( 'CELL' )
                    GFLAG      = .TRUE.
                    RPT_%USEGMAT  = .TRUE.
                    RPT_%BYCELL = .TRUE.

                CASE( 'COUNTRY' )
                    YFLAG      = .TRUE.
                    RPT_%BYCNRY = .TRUE.
                    IF( SEGMENT( 3 ) .EQ. 'NAME' ) RPT_%BYCONAM = .TRUE.

                CASE( 'STATE' )
                    YFLAG = .TRUE.
                    RPT_%BYSTAT = .TRUE.
                    IF( SEGMENT( 3 ) .EQ. 'NAME' ) RPT_%BYSTNAM = .TRUE.

                CASE( 'COUNTY' )
                    YFLAG = .TRUE.
                    RPT_%BYCNTY = .TRUE.
                    IF( SEGMENT( 3 ) .EQ. 'NAME' ) RPT_%BYCYNAM = .TRUE.

                CASE( 'ELEVSTAT' )
                    IF( CATEGORY .EQ. 'POINT' ) THEN
                	RPT_%BYELEV = .TRUE.

                    ELSE IF( FIRSTLOOP ) THEN

                        CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )

                    END IF

                CASE( 'HOUR' )
                    TFLAG      = .TRUE.           ! Implies temporal allocation
                    RPT_%USEHOUR  = .TRUE.
                    RPT_%BYHOUR = .TRUE.

                CASE( 'LAYER' )
                    IF( CATEGORY .EQ. 'POINT' ) THEN
                        LFLAG        = .TRUE.     ! Implies layer fractions file
                        TFLAG        = .TRUE.     ! Implies temporal allocation
                        RPT_%BYLAYER = .TRUE.
                        RPT_%LAYFRAC = .TRUE.
                        RPT_%USEHOUR = .TRUE.     ! Implies temporal allocation
                        RPT_%BYHOUR  = .TRUE.

                    ELSE IF( FIRSTLOOP ) THEN

                        CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )

                    END IF

                CASE( 'ROADCLASS' )
                    IF( CATEGORY .EQ. 'MOBILE' ) THEN
                	RPT_%BYRCL = .TRUE.

                    ELSE IF( FIRSTLOOP ) THEN

                        CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )

                    END IF

                CASE( 'SCC10' )
                    RPT_%BYSCC  = .TRUE.
                    RPT_%SCCRES = 10
                    IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                        NFLAG = .TRUE.
                        RPT_%SCCNAM = .TRUE.
                    END IF

                CASE( 'SOURCE' )
                    RPT_%BYSRC = .TRUE.
                    RPT_%BYCNTY = .TRUE.
                    RPT_%BYSCC  = .TRUE.
                    RPT_%SCCRES = 10
                    IF( SEGMENT( 3 ) .EQ. 'NAME' .OR.
     &                  SEGMENT( 4 ) .EQ. 'NAME' ) RPT_%SRCNAM = .TRUE.

                    IF( CATEGORY     .EQ. 'POINT'    .AND.
     &                ( SEGMENT( 3 ) .EQ. 'STACKPARM' .OR.
     &                  SEGMENT( 4 ) .EQ. 'STACKPARM'     ) )
     &                  RPT_%STKPARM = .TRUE.

                CASE( 'SPCCODE' )
                    PSFLAG = .TRUE.
                    RPT_%BYSPC = .TRUE.
                    RPT_%SPCPOL = SEGMENT( 3 )

                CASE( 'SRGCODE' )
                    IF( CATEGORY .NE. 'POINT' ) THEN
                        GSFLAG = .TRUE.
                        RPT_%BYSRG = .TRUE.
                        RPT_%SRGRES = 1

                        IF( SEGMENT(3) .EQ. 'FALLBACK' ) RPT_%SRGRES = 2

                    ELSE
                        CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )

                    END IF

                CASE( 'MONCODE' )
                    TSFLAG = .TRUE.
                    RPT_%BYMON = .TRUE.

                CASE( 'WEKCODE' )
                    TSFLAG = .TRUE.
                    RPT_%BYWEK = .TRUE.

                CASE( 'DIUCODE' )
                    TSFLAG = .TRUE.
                    RPT_%BYDIU = .TRUE.

                CASE DEFAULT
                    IF( FIRSTLOOP ) CALL WRITE_IGNORE_MESSAGE

                END SELECT

C.............  Setting for the use of layer fractions
            CASE( 'LAYFRAC' )

                IF( CATEGORY .EQ. 'POINT' ) THEN
                    LFLAG = .TRUE.
                    RPT_%LAYFRAC = .TRUE.

                ELSE IF( FIRSTLOOP ) THEN

                    CALL WRONG_SOURCE_CATEGORY( SEGMENT( 1 ) )

                END IF

C.............  Setting for the normalize instruction
            CASE( 'NORMALIZE' )

                SELECT CASE( SEGMENT( 2 ) )
                CASE( 'CELLAREA' ) 
                    GFLAG = .TRUE.                     ! Implies gridding
                    RPT_%NORMCELL = .TRUE.
                    RPT_%USEGMAT  = .TRUE.

                CASE( 'POPULATION' )
                    YFLAG = .TRUE.                     ! read cy/st/cy 
                    RPT_%NORMPOP = .TRUE.

                CASE DEFAULT
                    IF( FIRSTLOOP ) CALL WRITE_IGNORE_MESSAGE

                END SELECT

C.............  Settings the output data format
            CASE( 'NUMBER' )
                RPT_%DATAFMT = SEGMENT( 2 )

C.............  SELECT options affecting inputs needed
            CASE( 'SELECT' )

                SELECT CASE( SEGMENT( 2 ) )

                CASE( 'DATA' )
                    LIN_SUBDATA = .TRUE.
                    RPT_%NUMDATA = NSEGS - 2

                CASE( 'ELEVATED' )

                    IF( CATEGORY .EQ. 'POINT' ) THEN
                	RPT_%ELEVSTAT = ELEVOUT3
                	RPT_%BYELEV = .TRUE.

                	IF( SEGMENT( 3 ) .EQ. 'PING' ) THEN
                            VFLAG = .TRUE.
                            RPT_%ELEVSTAT = PINGOUT3
                        END IF

                    ELSE IF( FIRSTLOOP ) THEN

                        CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )

                    END IF

                CASE( 'NOELEVATED' )
                    IF( CATEGORY .EQ. 'POINT' ) THEN
                	RPT_%BYELEV = .TRUE.
                	RPT_%ELEVSTAT = NOELOUT3

                    ELSE IF( FIRSTLOOP ) THEN

                        CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )

                    END IF

                CASE( 'REGION' )

C.....................  Warning if used more than once
                    IF( FIRSTLOOP .AND. LREGION ) THEN

                        WRITE( MESG,94010 ) 'WARNING: Multiple ' //
     &                         'REGION selections made in same ' //
     &                         'report at line', IREC, CRLF() // 
     &                         BLANK10 // 'Final label will be applied.'
                        CALL M3MSG2( MESG )

                    END IF

                    LREGION = .TRUE.
                    LIN_SUBREGN = .TRUE.
                    CALL EXTRACT_LABEL( IREC, 'REGION', LINE, 
     &                                  RPT_%REGNNAM )

                CASE( 'SUBGRID' )

C.....................  Warning if used more than once
                    IF( FIRSTLOOP .AND. LSUBGRID ) THEN

                        WRITE( MESG,94010 ) 'WARNING: Multiple ' //
     &                         'SUBGRID selections made in same ' //
     &                         'report at line', IREC, CRLF() // 
     &                         BLANK10 // 'Final label will be applied.'
                        CALL M3MSG2( MESG )

                    END IF

                    GFLAG        = .TRUE.
                    RPT_%USEGMAT = .TRUE.
                    LSUBGRID = .TRUE.
                    LIN_SUBGRID = .TRUE.
                    CALL EXTRACT_LABEL( IREC, 'SUBGRID', LINE, 
     &                                  RPT_%SUBGNAM )

                CASE DEFAULT
                    IF( FIRSTLOOP ) CALL WRITE_IGNORE_MESSAGE

                END SELECT   ! on SELECT options                

            CASE( 'TITLE' )
                LIN_TITLE = .TRUE.
                RPT_%NUMTITLE = RPT_%NUMTITLE + 1
           
C.................  Keep any leading spaces in title with WRITE
                TITLE = ADJUSTL( LINE( 7:L2 ) )

            CASE( 'UNITS' )

C.................  Get units to use for all output data
                IF( SEGMENT( 2 ) .EQ. 'ALL' ) THEN

                    L = INDEX( LINE, 'ALL' )
                    UNITSET = ADJUSTL( LINE( L+3:L2 ) )
                    LIN_UNIT = .TRUE.

C.................  In future, support data-specific unit setting
                ELSE

                    IF( FIRSTLOOP ) CALL WRITE_IGNORE_MESSAGE

                END IF

            CASE DEFAULT
                IF( FIRSTLOOP ) CALL WRITE_IGNORE_MESSAGE

            END SELECT      ! on report instructions

        END IF              ! Report packet

C.........  Save variables from current call and exit from subroutine
999     CONTINUE

        PREC = IREC
        IREC_MIN  = MIN( IREC, IREC_MIN )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram prints a warning that lines of the
C               report configuration file are being ignored
            SUBROUTINE WRITE_IGNORE_MESSAGE

C.............  Local variables
            CHARACTER*300  MESG

C----------------------------------------------------------------------

            WRITE( MESG,94010 ) 'WARNING: Instructions not recognized '
     &             // 'at line', IREC, '. Line is ignored.'
            CALL M3MSG2( MESG )

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE WRITE_IGNORE_MESSAGE

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram prints a warning that instruction
C               is not valid for the current source category
            SUBROUTINE WRONG_SOURCE_CATEGORY( COMMAND )

C.............  Subroutine arguments
            CHARACTER(*), INTENT (IN) :: COMMAND

C.............  Local variables
            INTEGER        L

            CHARACTER*300  MESG

C----------------------------------------------------------------------

            L = LEN_TRIM( COMMAND )
            WRITE( MESG,94010 ) 'WARNING: Instruction "' // 
     &             COMMAND( 1:L ) // '" at line', IREC, 
     &             'not valid for ' // CATDESC // ' sources.' //
     &             CRLF() // BLANK10 // 'Line is ignored.'
            CALL M3MSG2( MESG )

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE WRONG_SOURCE_CATEGORY

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram prints a warning that no setting
C               was made by a packet that could have, and the previous
C               setting will be used.
            SUBROUTINE NO_SETTING_FOUND( IREC, IDX )

C.............  Subroutine arguments
            INTEGER     , INTENT (IN) :: IREC
            INTEGER     , INTENT (IN) :: IDX

C.............  Local variables
            INTEGER        L
            CHARACTER*300  MESG

C----------------------------------------------------------------------

            L = LEN_TRIM( ALLPCKTS( IDX ) )
            WRITE( MESG,94010 ) 'WARNING: No value listed with "' //
     &             ALLPCKTS( IDX )( 1:L ) // '" packet at' // CRLF() //
     &             BLANK10 // 'line', IREC, '. Previous value ' //
     &             '(or default) will continue to apply.'
            CALL M3MSG2( MESG )

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE NO_SETTING_FOUND

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram extracts the label of the group or
C               subgrid from the SELECT instruction
            SUBROUTINE EXTRACT_LABEL( IREC, KEYWORD, LINE, LABEL )

C.............  Subroutine arguments
            INTEGER     , INTENT (IN) :: IREC
            CHARACTER(*), INTENT (IN) :: KEYWORD
            CHARACTER(*), INTENT (IN) :: LINE
            CHARACTER(*), INTENT(OUT) :: LABEL

C.............  Local variables
            INTEGER        K        ! counters and indices
            INTEGER        LK, LL   ! string lengths

            CHARACTER*300  MESG

C----------------------------------------------------------------------

C.............  Determine lengths of input strings
            LK = LEN_TRIM( KEYWORD )
            LL = LEN_TRIM( LINE )

C.............  Determine first position after keyword
            K = INDEX( LINE, KEYWORD ) + LK + 1

C.............  Make sure a label is there. If not, error.
            IF( LL .LE. K ) THEN

                RC_ERROR = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: No label found for "' //
     &                 'SELECT ' // KEYWORD // '" instruction at ' //
     &                 'line', IREC, CRLF() // BLANK10 // 
     &                 'in report configuration file.'
                CALL M3MSG2( MESG )

C.............  If so, store it in output argument
            ELSE

                LABEL = LINE( K:LL )
                CALL UPCASE( LABEL )

            END IF

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE EXTRACT_LABEL

        END SUBROUTINE PRCLINRC

