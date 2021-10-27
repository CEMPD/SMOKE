
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
C     Created 7/2000 by M. Houyoux
C     Revised 7/2003 by A. Holland
C
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C  
C COPYRIGHT (C) 2005, Environmental Modeling for Policy Development
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
        USE MODREPRT, ONLY: GRPNRECS, INPACKET, INGROUP, INREPORT,
     &                      INSPCIFY, PKTCOUNT, PKTSTATUS, RPT_,
     &                      RC_ERROR, MINC, LIN_DEFGRP, LIN_GROUP,
     &                      LIN_SUBDATA, LIN_SUBGRID, LIN_SUBREGN,
     &                      LIN_TITLE, LIN_UNIT, LIN_SPCIFY, PCKTNAM,
     &                      PKT_IDX, NALLPCKT, ALLPCKTS, PKTSTART, 
     &                      TIM_IDX, PKTEND, ADY_IDX, FIL_IDX,
     &                      DEL_IDX, REG_IDX, SBG_IDX, GRP_INCLSTAT,
     &                      GRP_LABEL, YFLAG, ELG_IDX, PNG_IDX, ELV_IDX,
     &                      SPCF_NOR, SPCF_NAND, RPT_IDX, LREGION,
     &                      LSUBGRID, TITLE, DATAMISS, VFLAG, AFLAG,
     &                      PRFLAG, PRRPTFLG, CUFLAG, GFLAG, CRFLAG,
     &                      SSFLAG, SLFLAG, TFLAG, LFLAG, NFLAG, PSFLAG,
     &                      GSFLAG, TSFLAG, UNITSET, MXRPTNVAR,
     &                      ELEVOUT3, PINGOUT3, NOELOUT3, FIL_ONAME,
     &                      NIFLAG, NMFLAG, NNFLAG, NOFLAG, SDFLAG,
     &                      LAB_IDX, LENLAB3,
     &                      DLFLAG, MATFLAG, NFDFLAG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CRL, CATDESC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)  CRLF
        INTEGER       INDEX1
        INTEGER       STR2INT

        EXTERNAL   CRLF, INDEX1, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: IREC        ! line counter
        INTEGER     , INTENT (IN) :: NSEGS       ! no. line segments
        CHARACTER(*), INTENT (IN) :: LINE        ! full input record less commnt
        CHARACTER(*), INTENT (IN) :: SEGMENT( * )! parsed input record

C...........   Local parameters
        CHARACTER(1), PARAMETER :: SCCLEV( NSCCLV3 ) =
     &                             (/ '1', '2', '3', '4', '5' /)

C...........   Other local variables
        INTEGER          I, J, K, L, L2      ! counters and indices
        INTEGER          IOS                 ! i/o status

        INTEGER, SAVE :: IREC_MIN = 99999999 ! minimum record number
        INTEGER, SAVE :: PREC     = -9       ! previous call record number

        LOGICAL, SAVE :: ASCFLAG   = .FALSE.  ! true: using ascii file
        LOGICAL, SAVE :: FIRSTIME  = .TRUE.   ! true: first time routine called
        LOGICAL, SAVE :: FIRSTLOOP = .TRUE.   ! true: file already read once
        LOGICAL, SAVE :: LCATSET   = .FALSE.  ! true: source category set
        LOGICAL, SAVE :: LDELIM    = .FALSE.  ! true: manual delimeter set
        LOGICAL, SAVE :: LLABEL    = .FALSE.  ! true: user-defined label set
        LOGICAL, SAVE :: HHFLAG    = .FALSE.  ! true: indicator of BY HOUR set for BY LAYER

        CHARACTER        DAYYN             !  Y or N for average day
        CHARACTER(6)     FILNUM            !  tmp file number string
        CHARACTER(300)   MESG              !  message buffer

        CHARACTER(LENLAB3) ENVBUF              ! environment variable buffer
        CHARACTER(LENLAB3), SAVE :: TMPLABEL   ! tmp label name
        CHARACTER(LENLAB3) TMPLABEL2           ! another tmp label name

        CHARACTER(16) :: PROGNAME = 'PRCLINRC' ! program name

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

            RPT_%AVEDAY   = .FALSE.  ! default to not use average day data
            RPT_%OUTTIME  = 230000   ! default to output at 2300 hours
            RPT_%DELIM    = ';'      ! default to semi-colon
            RPT_%USELABEL = .FALSE.  ! default to not use a label in report

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

C.....................  An average day data usage packet
                    CASE( ADY_IDX )

C.........................  Ensure average day setting made
                        IF( J+1 .GE. L2 ) THEN
                            CALL NO_SETTING_FOUND( IREC, PKT_IDX )

                        ELSE                            
                            DAYYN = ADJUSTL( LINE( J+1:L2 ) )
                            CALL UPCASE( DAYYN )
                            IF( DAYYN .EQ. 'Y' ) THEN
                                RPT_%AVEDAY = .TRUE.

                            ELSE IF( DAYYN .EQ. 'N' ) THEN
                                RPT_%AVEDAY = .FALSE.

                            ELSE
                                L = LEN_TRIM( DAYYN )
                                WRITE( MESG,94010 ) 
     &                            'WARNING: Unrecognized /AVEDAY/ ' //
     &                            'settting "' // DAYYN( 1:L ) // 
     &                            '" at line', IREC, '. Setting to ' //
     &                            'default of FALSE.'  
                                CALL M3MSG2( MESG )
                                RPT_%AVEDAY = .FALSE.

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

C.........................  Ensure delimeter setting made
                        IF( J+1 .GE. L2 ) THEN
                            CALL NO_SETTING_FOUND( IREC, PKT_IDX )
                        ELSE                            
                            RPT_%DELIM = ADJUSTL( LINE( J+1:L2 ) )
                            LDELIM = .TRUE.
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

C.......................  If a region packet, read state/county file
                        IF( PKT_IDX .EQ. REG_IDX ) YFLAG = .TRUE.

C.....................  Specification of elevated groups, PinG sources, or
C                       elevated sources
                    CASE( ELG_IDX, PNG_IDX, ELV_IDX )
                        INSPCIFY    = .TRUE.
                        LIN_SPCIFY  = .TRUE.
                        SPCF_NOR    = 0
                        SPCF_NAND   = 0

C.....................  Set a label - figure out what the label is...
                    CASE( LAB_IDX )

C.........................  Ensure a label is provided
                        IF( J+1 >= L2 ) THEN
                            CALL NO_SETTING_FOUND( IREC, PKT_IDX )
                        ELSE
                            TMPLABEL2 = ADJUSTL( LINE( J+1:L2 ) )
                        END IF

C.........................  If label is too long, it will be truncated
                        IF( L2 - J > LENLAB3 ) THEN
                            WRITE( MESG,94010 ) 
     &                            'WARNING: Label longer than',
     &                            LENLAB3, '-character maximum at line',
     &                            IREC, 'and will be truncated.'

                            CALL M3MSG2( MESG )
                        END IF

C.........................  Check if label is an environment variable
                        IF( TMPLABEL2( 1:1 ) == '$' ) THEN

C.............................  Try to evaluate environment variable
                            ENVBUF = TMPLABEL2( 2:LENLAB3 )
                            CALL ENVSTR( ENVBUF,' ',' ',TMPLABEL2,IOS )

C.............................  If E.V. can't be evaluated, write message
                            IF( IOS /= 0 ) THEN
                                WRITE( MESG,94010 )
     &                            'WARNING: Environment variable "' //  
     &                            TRIM( ENVBUF ) // '" used in label '//
     &                            'at line', IREC, 'is misformatted, '//
     &                            'not defined, or defined as blank.'
                                CALL M3MSG2( MESG )
                                LLABEL = .FALSE.

C.............................  Otherwise, store label
                            ELSE
                                LLABEL = .TRUE.
                                TMPLABEL = TMPLABEL2

                            END IF

C.........................  If label is set to OFF, then turn it off
                        ELSE IF( TMPLABEL2( 1:3 ) == 'OFF' ) THEN
                            LLABEL = .FALSE.

C.........................  If not environment variable, then set label for report
                        ELSE
                            LLABEL = .TRUE.
                            TMPLABEL = TMPLABEL2
                        END IF

                        PKTEND   = IREC
                        INPACKET = .FALSE.            ! end implied 

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
                        INREPORT       = .TRUE.
                        RPT_%BYCELL    = .FALSE.
                        RPT_%BYCNRY    = .FALSE.
                        RPT_%BYCNTY    = .FALSE.
                        RPT_%BYCONAM   = .FALSE.
                        RPT_%BYCYNAM   = .FALSE.
                        RPT_%BYDATE    = .FALSE.
                        RPT_%BYELEV    = .FALSE.
                        RPT_%BYERPTYP  = .FALSE.
                        RPT_%BYGEO1    = .FALSE.
                        RPT_%BYGEO1NAM = .FALSE.
                        RPT_%BYHOUR    = .FALSE.
                        RPT_%BYLAYER   = .FALSE.
                        RPT_%BYLATLON  = .FALSE.
                        RPT_%BYMON     = .FALSE.
                        RPT_%BYDOM     = .FALSE.
                        RPT_%BYWEK     = .FALSE.
                        RPT_%BYMND     = .FALSE.
                        RPT_%BYTUE     = .FALSE.
                        RPT_%BYWED     = .FALSE.
                        RPT_%BYTHU     = .FALSE.
                        RPT_%BYFRI     = .FALSE.
                        RPT_%BYSAT     = .FALSE.
                        RPT_%BYSUN     = .FALSE.
                        RPT_%BYMET     = .FALSE.
                        RPT_%BYPLANT   = .FALSE.
                        RPT_%BYFACILITY = .FALSE.
                        RPT_%BYRCL     = .FALSE.
                        RPT_%BYSIC     = .FALSE.
                        RPT_%BYSCC     = .FALSE.
                        RPT_%BYINTGR   = .FALSE.
                        RPT_%BYMACT    = .FALSE.
                        RPT_%BYNAICS   = .FALSE.
                        RPT_%BYORIS    = .FALSE.
                        RPT_%BYSRCTYP  = .FALSE.
                        RPT_%MACTNAM   = .FALSE.
                        RPT_%NAICSNAM  = .FALSE.
                        RPT_%ORISNAM   = .FALSE.
                        RPT_%BYBOILER  = .FALSE.
                        RPT_%BYSPC     = .FALSE.
                        RPT_%BYSRC     = .FALSE.
                        RPT_%BYSRG     = .FALSE.
                        RPT_%BYSTACK   = .FALSE.
                        RPT_%BYSTKPARM = .FALSE.
                        RPT_%BYSTAT    = .FALSE.
                        RPT_%BYSTNAM   = .FALSE.
                        RPT_%CARB      = .FALSE.
                        RPT_%CHKPROJ   = .FALSE.
                        RPT_%CHKCNTL   = .FALSE.
                        RPT_%ELVSTKGRP = .FALSE.
                        RPT_%LATLON    = .FALSE.
                        RPT_%GRDCOR    = .FALSE.
                        RPT_%GRDPNT    = .FALSE.
                        RPT_%LAYFRAC   = .FALSE.
                        RPT_%NORMCELL  = .FALSE.
                        RPT_%NORMPOP   = .FALSE.
                        RPT_%GSPRONAM  = .FALSE.
                        RPT_%SICNAM    = .FALSE.
                        RPT_%SCCNAM    = .FALSE.
                        RPT_%SRCNAM    = .FALSE.
                        RPT_%STKPARM   = .FALSE.
                        RPT_%FUGPARM   = .FALSE.
                        RPT_%USEASCELEV= .FALSE.
                        RPT_%USECRMAT  = .FALSE.
                        RPT_%USECUMAT  = .FALSE.
                        RPT_%USEGMAT   = .FALSE.
                        RPT_%USEHOUR   = .FALSE.
                        RPT_%USELABEL  = LLABEL 
                        RPT_%USEPRMAT  = .FALSE.
                        RPT_%USESLMAT  = .FALSE.
                        RPT_%USESSMAT  = .FALSE.
                        RPT_%SRCMAP    = .FALSE.
                        LREGION        = .FALSE.
                        LSUBGRID       = .FALSE.

                        RPT_%ELEVSTAT  = 0
                        RPT_%NUMDATA   = -9   ! zero is legitimate
                        RPT_%NUMFILES  = 0
                        RPT_%NUMSECT   = 0
                        RPT_%NUMTITLE  = 0
                        RPT_%RENDLIN   = 0    ! init for consistency
                        RPT_%RPTMODE   = 0
                        RPT_%RPTNVAR   = 0
                        RPT_%RSTARTLIN = 0    ! init for consistency
                        RPT_%SCCRES    = 4
                        RPT_%SRGRES    = 0

                        RPT_%DATAFMT   = 'E8.3'
                        RPT_%OFILENAM  = ' '    ! init for consistency
                        RPT_%REGNNAM   = ' '
                        RPT_%SUBGNAM   = ' '
                        RPT_%SPCPOL    = ' '
                        TITLE          = ' '
                        IF( .NOT. LDELIM )
     &                     RPT_%DELIM    = ';'      ! default to semi-colon
                        IF( LLABEL )
     &                     RPT_%LABEL = TMPLABEL

                    END SELECT

                END IF
 
            END IF

            GO TO 999    ! to end of routine

C........  End of packet...
       ELSE IF( SEGMENT( 1 ) .EQ. '/END/' ) THEN

C.............  Resolve status of settings that cannot be determined
C               definitively until the end of the packet...

C.............  If one report contains the ASCIIELEV instruction,
C               all reports must contain it.
               IF( PKTCOUNT( RPT_IDX ) .GT. 1 .AND.
     &                                 ASCFLAG ) THEN
                   IF( .NOT. RPT_%USEASCELEV ) THEN
                       MESG = 'ERROR: If one report contains ' //
     &                        'the ASCIIELEV instruction, '//
     &                        'all reports must contain it.'
                       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                   END IF
               END IF

C.............  If packet was a report packet and no data have
C               been selected, then set input file global setting
            IF( INREPORT .AND. RPT_%NUMDATA .LE. 0 ) THEN
                DATAMISS = .TRUE.
            END IF

C.............  If using elevated settings and no layer fractions, make sure
C               PELV file will be read in 
            IF( .NOT. VFLAG  .AND. .NOT. RPT_%LAYFRAC .AND.
     &        ( RPT_%ELEVSTAT .GT. 0 .OR. RPT_%BYELEV )
     &          .AND. .NOT. AFLAG      ) THEN
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

C.............  Report arrangement
            CASE( 'ARRANGE' )

                IF( SEGMENT( 2 ) .EQ. 'MULTIFILE' ) THEN
                    RPT_%RPTMODE = 1
                    READ( SEGMENT( 3 ), * ) RPT_%RPTNVAR

C...............  If number of variables per report is greater than the maximum,
C                 then reset to maximum value
                    IF( RPT_%RPTNVAR .GT. MXRPTNVAR ) THEN
                        RPT_%RPTNVAR = MXRPTNVAR

                        IF( FIRSTLOOP ) THEN
                          WRITE( MESG, 94010 )
     &                     'WARNING: Number of variables per report ' //
     &                     'at line', IREC, ' is greater than the ' //
     &                     'maximum allowed. Resetting to maximum ' //
     &                     'value.'
                          CALL M3MSG2( MESG )
                        END IF
                    END IF
                
                ELSE IF( SEGMENT( 2 ) .EQ. 'ONEFILE' ) THEN
                    RPT_%RPTMODE = 2
                    READ( SEGMENT( 3 ), * ) RPT_%RPTNVAR

C...............  If number of variables per report is greater than the maximum,
C                 then reset to maximum value
                    IF( RPT_%RPTNVAR .GT. MXRPTNVAR ) THEN
                        RPT_%RPTNVAR = MXRPTNVAR

                        IF( FIRSTLOOP ) THEN
                          WRITE( MESG, 94010 )
     &                     'WARNING: Number of variables per report ' //
     &                     'at line', IREC, ' is greater than the ' //
     &                     'maximum allowed. Resetting to maximum ' //
     &                     'value.'
                          CALL M3MSG2( MESG )
                        END IF
                    END IF

                ELSE IF( SEGMENT( 2 ) .EQ. 'DATABASE' ) THEN
                    RPT_%RPTMODE = 3

                ELSE
                    IF( FIRSTLOOP ) THEN
                        L = LEN_TRIM( SEGMENT( 2 ) )
                        WRITE( MESG, 94010 )
     &                      'WARNING: Arrangement type "'//
     &                      SEGMENT( 2 )( 1:L )// '" at line', IREC,
     &                      ' is not known.' // CRLF()// BLANK10 //
     &                      'Will use no arrangement.'
                        CALL M3MSG2( MESG )
                    END IF

                END IF

C.............  ASCII elevated sources file used for report
            CASE( 'ASCIIELEV' )

                IF( PKTCOUNT( RPT_IDX ) .EQ. 1 ) THEN
                    ASCFLAG = .TRUE.
                ELSE IF( PKTCOUNT( RPT_IDX ) .GT. 1 .AND.
     &                   .NOT. ASCFLAG ) THEN
                    MESG = 'ERROR: If one report contains the ' //
     &                     'ASCIIELEV instruction, all reports ' //
     &                     'must contain it.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                RPT_%USEASCELEV = .TRUE.
                AFLAG = .TRUE.

C.............  Check control or projection matrix versus report
            CASE( 'CHECK' )

                SELECT CASE( SEGMENT( 2 ) )

                CASE( 'PROJECTION' )
                    PRFLAG        = .TRUE.
                    PRRPTFLG      = .TRUE.
                    RPT_%CHKPROJ  = .TRUE.
                    RPT_%USEPRMAT = .TRUE.
                    
C             CASE( 'CONTROL' )
C             CASE( 'REACTIVITY' )

                END SELECT

C.............  Multiplicative controls used for report
            CASE( 'CONTROL' )
                CUFLAG        = .TRUE.
                RPT_%USECUMAT = .TRUE.

C.............  Gridding used for report
            CASE( 'GRIDDING' )
                IF( .NOT. AFLAG ) THEN
                    GFLAG      = .TRUE.
                    RPT_%USEGMAT  = .TRUE.
                END IF

C.............  Source mapping for report
            CASE( 'CROSSWALK' )
                IF( CATEGORY .EQ. 'POINT' ) RPT_%SRCMAP = .TRUE.

C.............  Projection used for report
            CASE( 'PROJECTION' )
                PRFLAG        = .TRUE.
                RPT_%USEPRMAT = .TRUE.

C.............  Reactivity controls used for report
            CASE( 'REACTIVITY' )
                CRFLAG        = .TRUE.
                RPT_%USECRMAT = .TRUE.

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
                    RPT_%USESSMAT  = .TRUE.

                END IF

C.............  Temporal allocated emission used for report
            CASE( 'TEMPORAL' )
                TFLAG     = .TRUE.
                RPT_%USEHOUR = .TRUE.
                RPT_%BYDATE = .TRUE.

C.............  AERMOD support report
            CASE( 'AERMOD' )

                GFLAG           = .TRUE.
                RPT_%USEGMAT    = .TRUE.
                RPT_%BYCELL     = .TRUE.    ! reporty by cell
                YFLAG           = .TRUE.    ! read costcy input file
                RPT_%USEASCELEV = .FALSE.

                IF( SEGMENT( 2 ) .EQ. 'POINT' ) THEN

                    RPT_%BYSRCTYP  = .TRUE.      ! By facility source type
                    RPT_%BYSTAT    = .TRUE.      ! report by state
                    RPT_%GRDCOR    = .TRUE.      ! calculate grid lambert_x-y and utm x_y with zone
                    RPT_%SRCMAP    = .TRUE.      ! output source mapping output file

                    RPT_%BYPLANT   = .TRUE.      ! By Plant ID
                    RPT_%SRCNAM    = .TRUE.      ! By Plant Name

                    TSFLAG         = .TRUE.      ! By TSUP file
                    RPT_%BYMON     = .TRUE.      ! By monthly profile ID
                    RPT_%BYDOM     = .TRUE.      ! By day-month profile ID
                    RPT_%BYWEK     = .TRUE.      ! By weekly profile ID
                    RPT_%BYMND     = .TRUE.
                    RPT_%BYTUE     = .TRUE.
                    RPT_%BYWED     = .TRUE.
                    RPT_%BYTHU     = .TRUE.
                    RPT_%BYFRI     = .TRUE.
                    RPT_%BYSAT     = .TRUE.
                    RPT_%BYSUN     = .TRUE.

                    RPT_%LATLON    = .TRUE.      ! include lat-lon coords in report

                    IF( SEGMENT( 3 ) .EQ. 'PTNONIPM' .OR.
     &                  SEGMENT( 3 ) .EQ. 'PTEGU'         ) THEN

                        RPT_%BYERPTYP  = .TRUE.   ! By release point type
                        RPT_%BYSTKPARM = .TRUE.   ! By stack parameters
                        RPT_%STKPARM   = .TRUE.   ! stack parameters
                        RPT_%FUGPARM   = .TRUE.   ! By fugitive parameters
                        RPT_%BYLATLON  = .TRUE.   ! By lat and long coordinates

                        IF( SEGMENT( 4 ) .EQ. 'EMIS' ) THEN
                            TFLAG          = .TRUE.    ! read PTMP file
                            RPT_%USEHOUR   = .TRUE.    ! By hourly emissions
                        END IF

                    END IF

                ELSE IF( SEGMENT( 2 ) .EQ. 'NONPOINT' ) THEN

                    RPT_%BYCNTY = .TRUE.   ! By county
                    RPT_%BYSCC  = .TRUE.   ! By SCC
                    RPT_%GRDPNT = .TRUE.   ! report grid corner coordinates

                    TSFLAG      = .TRUE.   ! By TSUP file
                    RPT_%BYMON  = .TRUE.   ! By monthly profile ID
                    RPT_%BYDOM  = .TRUE.   ! By day-month profile ID
                    RPT_%BYWEK  = .TRUE.   ! By weekly profile ID
                    RPT_%BYMND  = .TRUE.
                    RPT_%BYTUE  = .TRUE.
                    RPT_%BYWED  = .TRUE.
                    RPT_%BYTHU  = .TRUE.
                    RPT_%BYFRI  = .TRUE.
                    RPT_%BYSAT  = .TRUE.
                    RPT_%BYSUN  = .TRUE.

                ELSE IF( SEGMENT( 2 ) .EQ. 'ANNUAL' ) THEN

                    RPT_%BYCNTY = .TRUE.   ! by county level report
                    RPT_%BYSCC  = .TRUE.   ! by SCC
                    RPT_%BYSRC  = .TRUE.   ! by source level report
                    RPT_%GRDPNT = .TRUE.   ! report grid corner coordinates


                ELSE

                    IF( FIRSTLOOP ) THEN
                        L = LEN_TRIM( SEGMENT( 2 ) )
                        WRITE( MESG,94010 ) 
     &                    'WARNING: AERMOD sector type "'// 
     &                    SEGMENT( 2 )( 1:L )// '" at line', IREC, 
     &                    'is not currently supported.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                END IF

C.............  BY options affecting inputs needed
            CASE( 'BY' )

                SELECT CASE( SEGMENT( 2 ) )

                CASE( 'CELL' )
                    IF( .NOT. AFLAG ) THEN
                        GFLAG      = .TRUE.
                        RPT_%USEGMAT  = .TRUE.
                    END IF

                    IF( RPT_%CARB ) THEN
                        RPT_%BYCNTY = .FALSE. ! disable county-level report
                    END IF

                    RPT_%BYCELL = .TRUE.

                CASE( 'GEOCODE1' )
                    YFLAG = .TRUE.
                    RPT_%BYGEO1 = .TRUE.
                    IF( SEGMENT( 3 ) .EQ. 'NAME' ) RPT_%BYGEO1NAM = .TRUE.

                CASE( 'COUNTRY', 'GEOCODE2' )
                    YFLAG      = .TRUE.
                    RPT_%BYCNRY = .TRUE.
                    IF( SEGMENT( 3 ) .EQ. 'NAME' ) RPT_%BYCONAM = .TRUE.

                CASE( 'STATE', 'GEOCODE3' )
                    YFLAG = .TRUE.
                    RPT_%BYSTAT = .TRUE.
                    IF( SEGMENT( 3 ) .EQ. 'NAME' ) RPT_%BYSTNAM = .TRUE.

                CASE( 'COUNTY', 'GEOCODE4' )
                    YFLAG = .TRUE.
                    RPT_%BYCNTY = .TRUE.
                    IF( SEGMENT( 3 ) .EQ. 'NAME' ) RPT_%BYCYNAM = .TRUE.

                CASE( 'ELEVSTAT' )
                    IF( CATEGORY .EQ. 'POINT' ) THEN
                        RPT_%BYELEV = .TRUE.
                        IF( SEGMENT( 3 ) .EQ. 'STACKGROUP' ) THEN
                            RPT_%ELVSTKGRP = .TRUE.
                        END IF

                    ELSE IF( FIRSTLOOP ) THEN

                        CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )

                    END IF

                CASE( 'ERPTYPE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        IF( CATEGORY .EQ. 'POINT' ) THEN
                            RPT_%BYERPTYP = .TRUE.
                        ELSE IF( FIRSTLOOP ) THEN
                            CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )
                        END IF
                    END IF

                CASE( 'HOUR' )
                    IF( .NOT. AFLAG ) THEN
                        TFLAG      = .TRUE.           ! Implies temporal allocation
                        RPT_%USEHOUR  = .TRUE.
                    END IF

                    RPT_%BYHOUR = .TRUE.
                    HHFLAG      = .TRUE.          ! indicator flag for BY LAYER instruction only

                CASE( 'LAYER' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                      IF( CATEGORY .EQ. 'POINT' ) THEN
                        LFLAG        = .TRUE.     ! Implies layer fractions file
                        TFLAG        = .TRUE.     ! Implies temporal allocation
                        RPT_%BYLAYER = .TRUE.
                        RPT_%LAYFRAC = .TRUE.
                        RPT_%USEHOUR = .TRUE.     ! Implies temporal allocation
                        DLFLAG       = .FALSE.    ! initializing every report

C.........................  Daily layered emission is set to Y if BYHOUR is not set at this point
                        IF( .NOT. HHFLAG ) DLFLAG = .TRUE.

                        RPT_%BYHOUR  = .TRUE.

                      ELSE IF( FIRSTLOOP ) THEN

                        CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )

                      END IF

                    END IF

                CASE( 'LATLON' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        IF( CATEGORY .EQ. 'POINT' ) THEN
                            RPT_%BYLATLON = .TRUE.
                            RPT_%LATLON = .TRUE.
                        ELSE IF( FIRSTLOOP ) THEN
                            CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )
                        END IF
                    END IF

                CASE( 'ORIS' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        RPT_%BYORIS = .TRUE.
                        IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                            NOFLAG = .TRUE.
                            RPT_%ORISNAM = .TRUE.
                        END IF
                    END IF

                 CASE( 'BOILER' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        RPT_%BYBOILER = .TRUE.
                        IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                            NOFLAG = .TRUE.
                            RPT_%ORISNAM = .TRUE.
                        END IF
                    END IF

                CASE( 'PLANT' )
                    RPT_%BYPLANT = .TRUE.
                    IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                        RPT_%SRCNAM = .TRUE.
                    END IF

                CASE( 'FACILITY' )
                    RPT_%BYFACILITY = .TRUE.
                    IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                        RPT_%SRCNAM = .TRUE.
                    END IF

                CASE( 'ROADCLASS' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        IF( CATEGORY .EQ. 'MOBILE' ) THEN
                            RPT_%BYRCL = .TRUE.

                        ELSE IF( FIRSTLOOP ) THEN

                            CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )

                        END IF
                    END IF

                CASE( 'SIC' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        RPT_%BYSIC  = .TRUE.
                        IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                            NIFLAG = .TRUE.
                            RPT_%SICNAM = .TRUE.
                            IF( .NOT. LDELIM ) RPT_%DELIM  = '|'
                        END IF
                    END IF

                CASE( 'MATBURNED' )    ! using SIC as an alias for MATBURNED in wildfire
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        RPT_%BYSIC  = .TRUE.
                        MATFLAG     = .TRUE.
                        IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                            NIFLAG = .TRUE.
                            RPT_%SICNAM = .TRUE.
                            IF( .NOT. LDELIM ) RPT_%DELIM  = '|'
                        END IF
                    END IF

                CASE( 'INTEGRATE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        RPT_%BYINTGR  = .TRUE.
                    END IF

                CASE( 'MACT' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        RPT_%BYMACT  = .TRUE.
                        IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                            NMFLAG = .TRUE.
                            RPT_%MACTNAM = .TRUE.
                            IF( .NOT. LDELIM ) RPT_%DELIM  = '|'
                        END IF
                    END IF

                CASE( 'NFDRSCODE' )    ! using MACT as an alias for NFDRSCODE in wildfire
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        RPT_%BYMACT  = .TRUE.
                        NFDFLAG      = .TRUE.
                        IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                            NMFLAG = .TRUE.
                            RPT_%MACTNAM = .TRUE.
                            IF( .NOT. LDELIM ) RPT_%DELIM  = '|'
                        END IF
                    END IF

                CASE( 'NAICS' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        RPT_%BYNAICS  = .TRUE.
                        IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                            NNFLAG = .TRUE.
                            RPT_%NAICSNAM = .TRUE.
                            IF( .NOT. LDELIM ) RPT_%DELIM  = '|'
                        END IF
                    END IF

                CASE( 'SCC10' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        RPT_%BYSCC  = .TRUE.
                        RPT_%SCCRES = 4
                        IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                            NFLAG = .TRUE.
                            RPT_%SCCNAM = .TRUE.
                            IF( .NOT. LDELIM ) RPT_%DELIM  = '|'
                        END IF
                    END IF

                CASE( 'SRCTYPE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        RPT_%BYSRCTYP  = .TRUE.
                    END IF

                CASE( 'SCC' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        RPT_%BYSCC = .TRUE.
                        K = INDEX1( SEGMENT(3)(1:1), NSCCLV3, SCCLEV )

                        IF( K > 0 ) THEN
                            RPT_%SCCRES = STR2INT( SEGMENT( 3 ) )

                        ELSE
                            WRITE( MESG,94010 )
     &                        'WARNING: BY SCC instruction at ' //
     &                        'line', IREC, 'does not include proper '//
     &                        'SCC aggregation level 1-4. Assuming ' //
     &                        'full SCC.'
                            CALL M3MSG2( MESG )
                            RPT_%SCCRES = 4
                        END IF

                        IF( SEGMENT( 3 ) == 'NAME' .OR.
     &                      SEGMENT( 4 ) == 'NAME'      ) THEN
                            NFLAG = .TRUE.
                            RPT_%SCCNAM = .TRUE.
                            IF( .NOT. LDELIM ) RPT_%DELIM = '|'
                        END IF
                    END IF
                        
                CASE( 'SOURCE' )
                    RPT_%BYSRC   = .TRUE.
                    RPT_%BYFACILITY = .FALSE. ! would be a duplicate
                    RPT_%BYPLANT = .FALSE.  ! would be a duplicate
                    RPT_%BYCNTY  = .TRUE.
                    IF( .NOT. AFLAG ) THEN
                        RPT_%BYSCC   = .TRUE.
                        RPT_%SCCRES  = 4
                        IF ( CATEGORY .EQ. 'POINT' ) RPT_%BYSIC = .TRUE.
                    END IF
                    IF( SEGMENT( 3 ) .EQ. 'NAME' .OR.
     &                  SEGMENT( 4 ) .EQ. 'NAME' .OR.
     &                  SEGMENT( 5 ) .EQ. 'NAME' .OR.
     &                  SEGMENT( 6 ) .EQ. 'NAME'      ) THEN
                        IF( CATEGORY .EQ. 'POINT' ) THEN
                            RPT_%SRCNAM = .TRUE.
                        ELSE
                            NFLAG = .TRUE.
                            RPT_%SCCNAM = .TRUE.
                            IF( .NOT. LDELIM ) RPT_%DELIM  = '|'
                        END IF
                    END IF

                    IF( CATEGORY     .EQ. 'POINT'    .AND.
     &                ( SEGMENT( 3 ) .EQ. 'STACKPARM' .OR.
     &                  SEGMENT( 4 ) .EQ. 'STACKPARM' .OR.
     &                  SEGMENT( 5 ) .EQ. 'STACKPARM' .OR.
     &                  SEGMENT( 6 ) .EQ. 'STACKPARM'      ) )
     &                  RPT_%STKPARM = .TRUE.

                    IF( CATEGORY     .EQ. 'POINT'    .AND.
     &                ( SEGMENT( 3 ) .EQ. 'FUGPARM' .OR.
     &                  SEGMENT( 4 ) .EQ. 'FUGPARM' .OR.
     &                  SEGMENT( 5 ) .EQ. 'FUGPARM' .OR.
     &                  SEGMENT( 6 ) .EQ. 'FUGPARM'      ) )
     &                  RPT_%FUGPARM = .TRUE.

                    IF( SEGMENT( 3 ) .EQ. 'LATLON' .OR.
     &                  SEGMENT( 4 ) .EQ. 'LATLON' .OR.
     &                  SEGMENT( 5 ) .EQ. 'LATLON' .OR.
     &                  SEGMENT( 6 ) .EQ. 'LATLON'      )
     &                  RPT_%LATLON = .TRUE.

                CASE( 'SPCCODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        PSFLAG = .TRUE.
                        SSFLAG = .TRUE.    ! open SMAT intermed files
                        RPT_%USESSMAT  = .TRUE.  ! use SMAT intermed files
                        RPT_%BYSPC = .TRUE.
                        RPT_%SPCPOL = SEGMENT( 3 )
                        IF( GFLAG .OR. RPT_%USEGMAT ) THEN
                            MESG = 'CRITICAL: GRIDDING command is not '
     &                          // 'allowed in BY SPCCODE report'
                            CALL M3MSG2( MESG)
                            GFLAG      = .FALSE.
                            RPT_%USEGMAT  = .FALSE.
                        END IF
                        IF( .NOT. LDELIM ) RPT_%DELIM  = '|'
                        IF( SEGMENT( 3 ) .EQ. 'NAME' ) THEN
                            RPT_%SPCPOL = SEGMENT( 4 )
                            SDFLAG = .TRUE.                    ! read GSPRODESC input file
                            RPT_%GSPRONAM = .TRUE.
                        END IF

                    END IF

                CASE( 'SRGCODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        IF( CATEGORY .NE. 'POINT' ) THEN
                            GSFLAG = .TRUE.
                            RPT_%BYSRG = .TRUE.
                            RPT_%SRGRES = 1

                            IF( SEGMENT(3) .EQ. 'FALLBACK' )
     &                                      RPT_%SRGRES = 2

                        ELSE
                            CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )

                        END IF
                    END IF

                CASE( 'STACK' )
                    RPT_%BYSTACK = .TRUE.
                    RPT_%BYPLANT = .TRUE.
                    RPT_%BYFACILITY = .TRUE.
                    IF( SEGMENT( 3 ) .EQ. 'STACKPARM' .OR.
     &                  SEGMENT( 4 ) .EQ. 'STACKPARM'      )
     &                  RPT_%STKPARM = .TRUE.

                    IF( SEGMENT( 3 ) .EQ. 'LATLON' .OR.
     &                  SEGMENT( 4 ) .EQ. 'LATLON'      )
     &                  RPT_%LATLON = .TRUE.

                CASE( 'STACKPARM' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                        IF( CATEGORY .EQ. 'POINT' ) THEN
                            RPT_%BYSTKPARM = .TRUE.
                            RPT_%STKPARM = .TRUE.
                            RPT_%FUGPARM = .TRUE.
                        ELSE IF( FIRSTLOOP ) THEN
                            CALL WRONG_SOURCE_CATEGORY( SEGMENT( 2 ) )
                        END IF
                    END IF

                CASE( 'MONCODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                       TSFLAG = .TRUE.
                       RPT_%BYMON = .TRUE.   ! monthly profile
                   END IF

                CASE( 'WEKCODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                       TSFLAG = .TRUE.
                       RPT_%BYWEK = .TRUE.   ! weekly profile
                   END IF

                CASE( 'DOMCODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                       TSFLAG = .TRUE.
                       RPT_%BYDOM = .TRUE.   ! day of month profile
                   END IF

                CASE( 'MNDCODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                       TSFLAG = .TRUE.
                       RPT_%BYMND = .TRUE.   ! monday hourly profile
                   END IF

                CASE( 'TUECODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                       TSFLAG = .TRUE.
                       RPT_%BYTUE = .TRUE.   ! tuesday hourly profile
                   END IF

                CASE( 'WEDCODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                       TSFLAG = .TRUE.
                       RPT_%BYWED = .TRUE.   ! Wednesday hourly profile
                   END IF

                CASE( 'THUCODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                       TSFLAG = .TRUE.
                       RPT_%BYTHU = .TRUE.   ! thursday hourly profile
                   END IF

                CASE( 'FRICODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                       TSFLAG = .TRUE.
                       RPT_%BYFRI = .TRUE.   ! Friday hourly profile
                   END IF

                CASE( 'SATCODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                       TSFLAG = .TRUE.
                       RPT_%BYSAT = .TRUE.   ! Saturday hourly profile
                   END IF

                CASE( 'SUNCODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                       TSFLAG = .TRUE.
                       RPT_%BYSUN = .TRUE.   ! Sunday hourly profile
                   END IF

                CASE( 'METCODE' )
                    IF( NOT_ASCIIELEV( 'BY ' // SEGMENT( 2 ) ) ) THEN
                       TSFLAG = .TRUE.
                       RPT_%BYMET = .TRUE.   ! Met-based hourly profile
                   END IF

                CASE DEFAULT
                    IF( FIRSTLOOP ) CALL WRITE_IGNORE_MESSAGE

                END SELECT

C.............  Setting for the QA extract data format for CARB
C               Fixed setting coditions
C               BY EIC and BY COABDIST (if not gridding)
C               BY DAY is default and BY HOUR (optional)
C
            CASE( 'CARB_QADEF' )

                RPT_%CARB       = .TRUE.  ! true: generating CARB QA extract table
                RPT_%RPTMODE    = 3       ! ARRANGE DATABASE output option. ORL-CARB format
                RPT_%USEASCELEV = .FALSE.

                RPT_%BYSCC  = .TRUE.      ! 20-digit EIC code report
                RPT_%SCCRES = 5           ! EIC3 support

                IF( RPT_%USEGMAT .OR. RPT_%BYCELL ) THEN

                    RPT_%USEGMAT = .TRUE.
                    RPT_%BYCELL  = .TRUE.      ! by cell is default when gridding is applied
                    RPT_%BYCNTY  = .FALSE.     ! default is full 12-digit coabdist code (county-level)
                    RPT_%USEHOUR = .FALSE.     ! hourly report is disabled when gridding is applied

                ELSE

                    RPT_%USEGMAT = .FALSE.
                    RPT_%BYCELL  = .FALSE.     ! by cell is default when gridding is applied
                    YFLAG        = .TRUE.      ! geocode level summary report
                    RPT_%BYCNTY  = .TRUE.      ! default is full 12-digit coabdist code (county-level)

                END IF

C.............  Setting for the use of layer fractions
            CASE( 'LAYFRAC' )

                IF( NOT_ASCIIELEV( SEGMENT( 1 ) ) ) THEN

                    IF( CATEGORY .EQ. 'POINT' ) THEN
                        LFLAG = .TRUE.
                        RPT_%LAYFRAC = .TRUE.

                    ELSE IF( FIRSTLOOP ) THEN

                        CALL WRONG_SOURCE_CATEGORY( SEGMENT( 1 ) )

                    END IF

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
            CHARACTER(300)  MESG

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

            CHARACTER(300)  MESG

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
            CHARACTER(300) MESG

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

            CHARACTER(300) MESG

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

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal function returns true if the ASCIIELEV instruction
C               is not in use; otherwise it prints a warning and returns false
            LOGICAL FUNCTION NOT_ASCIIELEV( COMMAND )

C.............  Subroutine arguments
            CHARACTER(*), INTENT (IN) :: COMMAND

C.............  Local variables
            INTEGER        L

            CHARACTER(300) MESG

C----------------------------------------------------------------------

            NOT_ASCIIELEV = .TRUE.

            IF( RPT_%USEASCELEV ) THEN
                NOT_ASCIIELEV = .FALSE.

                L = LEN_TRIM( COMMAND )
                WRITE( MESG, 94010 )
     &             'WARNING: ' // COMMAND( 1:L ) // ' instruction at ' //
     &             'line', IREC, 'is not allowed with ' //
     &             'the ASCIIELEV instruction.'
                CALL M3MSG2( MESG )
            END IF

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END FUNCTION NOT_ASCIIELEV

        END SUBROUTINE PRCLINRC

