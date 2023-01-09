
        SUBROUTINE RDRPRTS( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      The RDRPRTS routine reads the report information from the REPCONFIG 
C      and sets report arrays from the MODRPRT module
C
C  PRECONDITIONS REQUIRED:
C    REPCONFIG file is opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 7/2000 by M Houyoux
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
        USE MODREPRT, ONLY: ALLRPT, ALLOUTHR, ALLUSET, INDNAM,
     &                      TITLES, NLINE_RC, INREPORT, RPT_IDX,
     &                      LIN_TITLE, LIN_SUBDATA, LIN_UNIT, PSFLAG,
     &                      SPCPOL, LSPCPOL, NREPORT, NSPCPOL, RC_ERROR,
     &                      MXINDAT, MXTITLE, RPT_, PKTEND, PKTSTART,
     &                      FIL_ONAME, UNITSET, TITLE, PKTCOUNT

C...........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPOL, EINAM

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)  CRLF
        LOGICAL       BLKORCMT
        INTEGER       GETNLIST
        INTEGER       INDEX1

        EXTERNAL   CRLF, BLKORCMT, GETNLIST, INDEX1

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV    ! File unit number

C...........   Line parsing array
        CHARACTER(300) SEGMENT( 100 )
 
C...........   Other local variables
        INTEGER          H, I, J, L, N    ! counters and indices

        INTEGER          IOS     !  i/o status
        INTEGER          IREC    !  record counter
        INTEGER       :: NS = 1  !  no. segments in line
        INTEGER          NUNIT   !  tmp number of units

        LOGICAL       :: EFLAG   = .FALSE. !  true: error found

        CHARACTER(1024)  BUFFER            !  line work buffer
        CHARACTER(1024)  LINE              !  line input buffer
        CHARACTER(300)   MESG              !  message buffer

        CHARACTER(IOVLEN3) :: CBUF     !  tmp pollutant buffer

        CHARACTER(16) :: PROGNAME = 'RDRPRTS' ! program name

C***********************************************************************
C   begin body of subroutine RDRPRTS

C.........  Allocate memory for report arrays based on previous read of file 
C           and previously determined settings...

C.........  Allocate and initialize report arrays
        ALLOCATE( ALLRPT( NREPORT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLRPT', PROGNAME )
        ALLOCATE( ALLOUTHR( 24,NREPORT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLOUTHR', PROGNAME )
        ALLOCATE( ALLUSET( MXINDAT, NREPORT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLUSET', PROGNAME )
        ALLOCATE( INDNAM( MXINDAT, NREPORT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDNAM', PROGNAME )
        ALLOCATE( TITLES( MXTITLE, NREPORT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TITLES', PROGNAME )

        ALLRPT%BYCELL     = .FALSE.
        ALLRPT%BYCNRY     = .FALSE.
        ALLRPT%BYCNTY     = .FALSE.
        ALLRPT%BYCONAM    = .FALSE.
        ALLRPT%BYCYNAM    = .FALSE.
        ALLRPT%BYDATE     = .FALSE.
        ALLRPT%BYELEV     = .FALSE.
        ALLRPT%BYERPTYP   = .FALSE.
        ALLRPT%ELVSTKGRP  = .FALSE.
        ALLRPT%BYGEO1     = .FALSE.
        ALLRPT%BYGEO1NAM  = .FALSE.
        ALLRPT%BYHOUR     = .FALSE.
        ALLRPT%BYLAYER    = .FALSE.
        ALLRPT%BYPLANT    = .FALSE.
        ALLRPT%BYRCL      = .FALSE.
        ALLRPT%BYSCC      = .FALSE.
        ALLRPT%BYSIC      = .FALSE.
        ALLRPT%BYINTGR    = .FALSE.
        ALLRPT%BYMACT     = .FALSE.
        ALLRPT%BYNAICS    = .FALSE.
        ALLRPT%BYORIS     = .FALSE.
        ALLRPT%BYSRCTYP   = .FALSE.
        ALLRPT%MACTNAM    = .FALSE.
        ALLRPT%NAICSNAM   = .FALSE.
        ALLRPT%ORISNAM    = .FALSE.
        ALLRPT%BYSPC      = .FALSE.
        ALLRPT%BYSRC      = .FALSE.
        ALLRPT%BYSRG      = .FALSE.
        ALLRPT%BYSTACK    = .FALSE.
        ALLRPT%BYSTKPARM  = .FALSE.
        ALLRPT%BYSTAT     = .FALSE.
        ALLRPT%BYSTNAM    = .FALSE.
        ALLRPT%BYMON      = .FALSE.
        ALLRPT%BYWEK      = .FALSE.
        ALLRPT%BYDOM      = .FALSE.
        ALLRPT%BYMND      = .FALSE.
        ALLRPT%BYTUE      = .FALSE.
        ALLRPT%BYWED      = .FALSE.
        ALLRPT%BYTHU      = .FALSE.
        ALLRPT%BYFRI      = .FALSE.
        ALLRPT%BYSAT      = .FALSE.
        ALLRPT%BYSUN      = .FALSE.
        ALLRPT%BYMET      = .FALSE.
        ALLRPT%CHKPROJ    = .FALSE.
        ALLRPT%CHKCNTL    = .FALSE.
        ALLRPT%LAYFRAC    = .FALSE.
        ALLRPT%NORMCELL   = .FALSE.
        ALLRPT%NORMPOP    = .FALSE.
        ALLRPT%AVEDAY     = .FALSE.
        ALLRPT%SCCNAM     = .FALSE.
        ALLRPT%SRCNAM     = .FALSE.
        ALLRPT%STKPARM    = .FALSE.
        ALLRPT%FUGPARM    = .FALSE.
        ALLRPT%USEASCELEV = .FALSE.
        ALLRPT%USECRMAT   = .FALSE.
        ALLRPT%USECUMAT   = .FALSE.
        ALLRPT%USEGMAT    = .FALSE.
        ALLRPT%USEHOUR    = .FALSE.
        ALLRPT%USELABEL   = .FALSE.
        ALLRPT%USEPRMAT   = .FALSE.
        ALLRPT%USESLMAT   = .FALSE.
        ALLRPT%USESSMAT   = .FALSE.
        ALLRPT%DELIM      = ' '
        ALLRPT%DATAFMT    = ' '
        ALLRPT%LABEL      = ' '
        ALLRPT%OFILENAM   = ' '
        ALLRPT%REGNNAM    = ' '
        ALLRPT%SUBGNAM    = ' '

        ALLRPT%BEGSUMHR   = 0
        ALLRPT%ELEVSTAT   = 0
        ALLRPT%NUMDATA    = -9      ! zero is legitimate
        ALLRPT%NUMTITLE   = 0
        ALLRPT%NUMFILES   = 0
        ALLRPT%NUMSECT    = 0
        ALLRPT%OUTTIME    = 0
        ALLRPT%RENDLIN    = 0
        ALLRPT%RPTMODE    = 0
        ALLRPT%RPTNVAR    = 0
        ALLRPT%RSTARTLIN  = 0
        ALLRPT%SCCRES     = 4
        ALLRPT%SRGRES     = 0

        ALLOUTHR = .FALSE.
        ALLUSET  = ' '
        INDNAM  = ' '
        TITLES   = ' '

C.........  Read lines of file and store report characteristics
        IREC = 0
        DO I = 1, NLINE_RC

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading report configuration file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank lines and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Screen for appended comments and remove them
            CALL RMCOMMNT( '##', LINE )

C.............  Left-justify and convert line to upper case
            LINE = ADJUSTL( LINE )
            BUFFER = LINE
            CALL UPCASE( BUFFER )

C.............  Initialize segment from previous iteration
            SEGMENT( 1:NS ) = ' '

C.............  Parse line into segments
            L = LEN_TRIM( BUFFER )
            NS = GETNLIST( L, BUFFER )
            CALL PARSLINE( BUFFER, NS, SEGMENT )

C.............  Interpret line of code.  Set global variables in MODREPRT.
            CALL PRCLINRC( IREC, NS, LINE, SEGMENT )

C.............  Skip if report section not started yet.
            IF( .NOT. INREPORT ) CYCLE

C.............  Get count of report packets
            N = PKTCOUNT( RPT_IDX )

C.............  Store settings for current report
            ALLRPT( N ) = RPT_

            ALLRPT( N )%RENDLIN  = PKTEND
            ALLRPT( N )%RSTARTLIN= PKTSTART
            ALLRPT( N )%OFILENAM = FIL_ONAME

C.............  Conditional settings - only set if current line is type...
C.............  Title line
            IF( LIN_TITLE   ) 
     &          TITLES( RPT_%NUMTITLE , N ) = TITLE

C.............  Data subselection
            IF( LIN_SUBDATA ) 
     &          INDNAM( 1:RPT_%NUMDATA, N ) = SEGMENT( 3:NS )

C.............  Units - for now, one unit applies to all
C               note: Must edit here to permit units to be different by variable
            IF( LIN_UNIT ) THEN
                ALLUSET( :, N ) = UNITSET
            END IF

        END DO    ! End read loop of report configuration file

C.........  Rewind file
        REWIND( FDEV )

C.........  If speciation codes selected, allocate memory for pollutants list
C           and flag to indicate which of the selected pollutants were in the
C           supplemental speciation file.
        IF ( PSFLAG ) THEN
            ALLOCATE( SPCPOL( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPCPOL', PROGNAME )
            ALLOCATE( LSPCPOL( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LSPCPOL', PROGNAME )
            SPCPOL  = ' '       ! array
            LSPCPOL = .FALSE.   ! array
        END IF

C.........  Post-process reports...
        DO N = 1, NREPORT

C.............  Post-process output times and update logical array
C.............  If reporting "BY HOUR" or if temporal allocation not used for 
C               the report, set output for all "hours" to true. 
            IF( ALLRPT( N )%BYHOUR .OR. .NOT. ALLRPT( N )%USEHOUR ) THEN

                ALLOUTHR( 1:24, N ) = .TRUE.

C.............  Output for a specific hour (set as HHMMSS)
            ELSE
                H = MIN( ( ALLRPT( N )%OUTTIME/10000 ) + 1, 24 )
                ALLOUTHR( H, N ) = .TRUE.

            END IF

C.............  Check pollutants selected with SPCCODE and build a list of
C               valid pollutants
            IF ( ALLRPT( N )%BYSPC ) THEN

                CBUF = ALLRPT( N )%SPCPOL

C.................  Search for pollutant from REPCONFIG in inventory list
                J = INDEX1( CBUF, NIPOL, EINAM )

C.................  Error if not found
                IF ( J .LE. 0 ) THEN

                    L = LEN_TRIM( CBUF )
                    WRITE( MESG,94010 ) 'WARNING: Pollutant "'//
     &                     CBUF( 1:L ) // '" for report', N,
     &                     'not found in inventory ' // CRLF()//
     &                     BLANK10// 'pollutant list.  Ignoring '//
     &                     'instruction for speciation codes.'
                    CALL M3MSG2( MESG )
                    ALLRPT( N )%BYSPC = .FALSE.

C.................  Store first valid pollutant
                ELSE IF ( NSPCPOL .EQ. 0 ) THEN
                    NSPCPOL = 1
                    SPCPOL( 1 ) = CBUF
                    INDNAM( 1, N ) = CBUF      ! select pollutant
                    ALLRPT( N )%NUMDATA = 1

C.................  Only allow one pollutant in REPCONFIG file 
                ELSE IF ( CBUF .EQ. SPCPOL( 1 ) ) THEN
                    INDNAM( 1, N ) = CBUF      ! select pollutant
                    ALLRPT( N )%NUMDATA = 1

C.................  Store other valid pollutants, but ensure not already in list
                ELSE 
                    L = LEN_TRIM( CBUF )
                    WRITE( MESG,94010 ) 'WARNING: Skipping pollutant '//
     &                      '"'// CBUF( 1:L ) // '" at line', IREC, 
     &                      CRLF()// BLANK10 // 'Only allowed one '
     &                      // 'pollutant per REPCONFIG file '//
     &                      'with SPCCODE instruction.'
                    CALL M3MSG2( MESG )

                    INDNAM( 1, N ) = CBUF      ! select pollutant
                    ALLRPT( N )%NUMDATA = 1

c                    J = INDEX1( CBUF, NSPCPOL, SPCPOL )
c                    IF ( J .LE. 0 ) THEN
c                        NSPCPOL = NSPCPOL + 1
c                        SPCPOL( NSPCPOL ) = CBUF
c                    END IF

                END IF

            END IF    
                
        END DO

C.........  If there was an error reading the file
        IF( RC_ERROR ) EFLAG = .TRUE. 

C.........  If there was any error, exit 
        IF( EFLAG ) THEN
             MESG = 'Problem reading reports configuration file.'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C.........  Error message for reaching the end of file too soon
999     WRITE( MESG,94010 )
     &         'End of file reached unexpectedly at line', IREC, CRLF()
     &         //BLANK10 //'Check format of reports configuration file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE RDRPRTS
