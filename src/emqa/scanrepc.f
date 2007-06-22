
        SUBROUTINE SCANREPC( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C    The SCANREPC routine will be used to scan the REPCONFIG file and make
C    global program settings including what input files are needed and the 
C    source-category with related settings.  It determines the maximum number
C    of entries per group definition and the maximum number of titles per
C    report. It also determines the number of reports and the number of 
C    group packets of each type.
C
C  PRECONDITIONS REQUIRED:
C    REPCONFIG file is opened
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
        USE MODREPRT, ONLY: NLINE_RC, MXGRPREC, GRPNRECS, MXTITLE,
     &                      MXINDAT, NFILE, FIL_IDX, NREPORT, RPT_IDX,
     &                      NREGRAW, REG_IDX, NSBGRAW, SBG_IDX,
     &                      RC_ERROR, RPT_, PKTCOUNT

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)  CRLF
        LOGICAL       BLKORCMT
        INTEGER       GETFLINE
        INTEGER       GETNLIST

        EXTERNAL   CRLF, BLKORCMT, GETFLINE, GETNLIST

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! report configuration file unit no.

C...........   Line parsing array
        CHARACTER(300) SEGMENT( 100 )
 
C...........   Other local variables
        INTEGER         I, L    ! counters and indices

        INTEGER      :: ICNT = 0    !  no. non-blank/non-comment lines
        INTEGER          IOS     !  i/o status
        INTEGER          IREC    !  record counter
        INTEGER       :: N = 1   !  no. segments in line

        LOGICAL       :: EFLAG   = .FALSE. !  true: error found
        LOGICAL       :: LCATSET = .FALSE. !  true: source category set

        CHARACTER(1024)  BUFFER            !  line work buffer
        CHARACTER(1024)  LINE              !  line input buffer
        CHARACTER(300)   MESG              !  message buffer

        CHARACTER(16) :: PROGNAME = 'SCANREPC' ! program name

C***********************************************************************
C   begin body of subroutine SCANREPC
        
C.........  Get the number of lines in the file
        NLINE_RC = GETFLINE( FDEV, 'Report configuration file' )

C.........  Read lines of file to extract settings
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
            SEGMENT( 1:N ) = ' '

C.............  Parse line into segments
            L = LEN_TRIM( BUFFER )
            N = GETNLIST( L, BUFFER )
            CALL PARSLINE( BUFFER, N, SEGMENT )

C.............  Interpret line of code.  Set global variables in MODREPRT.
            CALL PRCLINRC( IREC, N, LINE, SEGMENT )

C.............  Maximum number of records per group
            MXGRPREC = MAX( MXGRPREC, GRPNRECS, 1 )

C.............  Maximum number of titles per report
            MXTITLE  = MAX( MXTITLE, RPT_%NUMTITLE, 1 )

C.............  Maximum number of selected input data variables per report
            MXINDAT = MAX( MXINDAT, RPT_%NUMDATA )

C.............  Keep count of non-blank/non-comment lines
            ICNT = ICNT + 1

        END DO    ! End read loop

C.........  Rewind file
        REWIND( FDEV )

C.........  Store the number of packets of each type
        NFILE   = PKTCOUNT( FIL_IDX )   ! Number of output files
        NREPORT = PKTCOUNT( RPT_IDX )   ! Number of reports
        NREGRAW = PKTCOUNT( REG_IDX )   ! Number of region groups
        NSBGRAW = PKTCOUNT( SBG_IDX )   ! Number of subgrids

C.........  If there were no non-blank and non-comment lines, error
        IF( ICNT .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: No non-blank or non-comments lines found ' //
     &             'in reports configuration file.'
            CALL M3MSG2( MESG )
        END IF

C.........  If there was an error reading the file
        IF( RC_ERROR ) EFLAG = .TRUE. 

C.........  If the source category has not been set, then error
        IF( CATEGORY .EQ. ' ' ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Source category not set in reports ' //
     &             'configuration file. Use SMK_SOURCE instruction.'
            CALL M3MSG2( MESG )
        END IF

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

        END SUBROUTINE SCANREPC
