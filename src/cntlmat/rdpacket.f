
        SUBROUTINE RDPACKET( FDEV, PKTTYP, USEPOL, IREC, 
     &                       PKTINFO, CFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads one line of a control packet, with different
C      read formats for different packet types
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     
C
C  REVISION  HISTORY:
C      Started 3/99 by M. Houyoux
C      09/2025 by HT UNC-IE:  Use M3UTILIO
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
C***************************************************************************
        USE M3UTILIO

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, EANAM

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
         INCLUDE 'CPKTDAT.EXT'   !  control packet contents

C...........   EXTERNAL FUNCTIONS:
c       LOGICAL       CHKREAL, BLKORCMT
c       CHARACTER(2)  CRLF
c       INTEGER       INDEX1
c       INTEGER       STR2INT
c       REAL          STR2REAL

c       EXTERNAL      BLKORCMT, CHKREAL, CRLF, INDEX1, STR2INT, STR2REAL
        LOGICAL, EXTERNAL :: BLKORCMT

C...........   SUBROUTINE ARGUMENTS:
        INTEGER        , INTENT (IN) :: FDEV      ! in file unit number
        CHARACTER(*)   , INTENT (IN) :: PKTTYP    ! packet type
        LOGICAL     , INTENT(IN OUT) :: USEPOL( NIPPA ) ! true: use pollutant
        INTEGER     , INTENT(IN OUT) :: IREC      ! file line number
        TYPE( CPACKET ), INTENT(OUT) :: PKTINFO   ! packet information
        LOGICAL        , INTENT(OUT) :: CFLAG     ! true: line is a comment
        LOGICAL        , INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters
        INTEGER, PARAMETER :: MXSEG = 17   ! number of potential line segments

C...........   Other arrays
        CHARACTER(32) SEGMENT( MXSEG )     ! Segments of parsed packet lines

C...........   Other local variables
        INTEGER         K          ! index
        INTEGER         LC, LN     ! line positions for using comments
        INTEGER         IOS        ! i/o error status
        INTEGER         CYID       ! tmp county ID
        INTEGER         FIP        ! tmp state/county FIPS code
        INTEGER         STID       ! tmp state ID

        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called

        CHARACTER(300)        LINE       ! line reading buffer
        CHARACTER(300)        MESG       ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDPACKET' ! program name

C***********************************************************************
C   Begin body of subroutine RDPACKET

        READ( FDEV, 93000, END = 999, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF( IOS .GT. 0 ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'ERROR: I/O error', IOS,
     &             'reading control packets file at line', IREC
            CALL M3MESG( MESG )

        END IF

C.........  Check for comment and blank lines
        CFLAG = .FALSE.
        IF( BLKORCMT( LINE ) ) THEN
            CFLAG = .TRUE.
            RETURN
        END IF

C.........  Parse the line of data into segments based on the rules
C           for "list-formatted" in fortran, but not requiring 
C           quotes around the text strings
        CALL PARSLINE( LINE, MXSEG, SEGMENT )

        LC = INDEX( LINE, '!' )
        IF ( LC .LE. 0 ) LC = 0
        LN = LEN( LINE )

        PKTINFO%COMMENT = " "
        IF ( LC .GT. 0 ) PKTINFO%COMMENT = LINE( LC+1:LN )

C.........  Process the line of data, depending on packet type
        SELECT CASE( PKTTYP )

        CASE( 'CTG' )
            PKTINFO%CSIC =           ' '
            PKTINFO%CFIP =           SEGMENT( 1 )
            PKTINFO%TSCC =           SEGMENT( 2 )
            PKTINFO%CPOL =           SEGMENT( 3 )
            PKTINFO%FAC1 = STR2REAL( SEGMENT( 4 ) )
            PKTINFO%FAC2 = STR2REAL( SEGMENT( 5 ) )
            PKTINFO%FAC3 = STR2REAL( SEGMENT( 6 ) )
            PKTINFO%CMCT =           ' '
            PKTINFO%PLT  =           ' '
            PKTINFO%CHAR1=           ' '
            PKTINFO%CHAR2=           ' '
            PKTINFO%CHAR3=           ' '
            PKTINFO%CHAR4=           ' '
            PKTINFO%CHAR5=           ' '

C.........  Check to see if cutoff value is missing, and give error
C           if it is.
            IF ( PKTINFO%FAC1 .LT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94020 ) 'ERROR: CUTOFF value missing '//
     &                'from CTG packet record at line', IREC
                CALL M3MSG2( MESG )
            END IF

C.........  Check to see if last column is blank. If blank, then
C           set PKTINFO%FAC4 = -9 and issue warning.
            IF ( SEGMENT( 7 ) .EQ. ' ' ) THEN
               PKTINFO%FAC4 = -9
               WRITE( MESG, 94020 ) 'WARNING: RACT value missing from'//
     &                ' CTG packet record at line', IREC
               CALL M3MESG( MESG )
            ELSE
               PKTINFO%FAC4 = STR2REAL( SEGMENT( 7 ) )
            END IF

        CASE( 'CONTROL' )
            PKTINFO%CFIP    =           SEGMENT( 1 )
            PKTINFO%TSCC    =           SEGMENT( 2 )
            PKTINFO%CPOL    =           SEGMENT( 3 )
            PKTINFO%FAC1    = STR2INT ( SEGMENT( 4 ) )
            PKTINFO%FAC2    = STR2REAL( SEGMENT( 5 ) )
            PKTINFO%FAC3    = STR2REAL( SEGMENT( 6 ) )
            PKTINFO%FAC4    = STR2REAL( SEGMENT( 7 ) )
            PKTINFO%CSIC    =           SEGMENT( 8 )
            PKTINFO%CMCT    =           SEGMENT( 9 )
            PKTINFO%APPFLAG = ADJUSTL ( SEGMENT( 10 ) )
            PKTINFO%REPFLAG = ADJUSTL ( SEGMENT( 11 ) )
            PKTINFO%PLT     =           SEGMENT( 12 )
            PKTINFO%CHAR1   =           SEGMENT( 13 )
            PKTINFO%CHAR2   =           SEGMENT( 14 )
            PKTINFO%CHAR3   =           SEGMENT( 15 )
            PKTINFO%CHAR4   =           SEGMENT( 16 )
            PKTINFO%CHAR5   =           SEGMENT( 17 )

        CASE( 'ALLOWABLE' )
            PKTINFO%CFIP =           SEGMENT( 1 )
            PKTINFO%TSCC =           SEGMENT( 2 )
            PKTINFO%CPOL =           SEGMENT( 3 )
            PKTINFO%FAC1 = STR2REAL( SEGMENT( 4 ) )
            PKTINFO%FAC2 = STR2REAL( SEGMENT( 5 ) )
            PKTINFO%FAC3 = STR2REAL( SEGMENT( 6 ) )
            PKTINFO%CSIC =           SEGMENT( 7 )
            PKTINFO%PLT  =           SEGMENT( 8 )
            PKTINFO%CHAR1=           SEGMENT( 9 )
            PKTINFO%CHAR2=           SEGMENT( 10 )
            PKTINFO%CHAR3=           SEGMENT( 11 )
            PKTINFO%CHAR4=           SEGMENT( 12 )
            PKTINFO%CHAR5=           SEGMENT( 13 )
            PKTINFO%CMCT =           ' '

C.........  Check to see if both CAP and REPLACE are missing. If so, issue
C           a warning.
            IF ( PKTINFO%FAC2 .LT. 0 .AND. PKTINFO%FAC3 .LT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94020 ) 'ERROR: Neither CAP or REPLACE '//
     &                'defined in allowable packet record at line', IREC
                CALL M3MSG2( MESG )
            END IF

        CASE( 'REACTIVITY' )
            PKTINFO%CFIP   =           SEGMENT( 1 )
            PKTINFO%TSCC   =           SEGMENT( 2 )
            PKTINFO%CPOL   =           SEGMENT( 3 )
            PKTINFO%FAC1   = STR2REAL( SEGMENT( 4 ) )
            PKTINFO%FAC2   = STR2REAL( SEGMENT( 5 ) )
            PKTINFO%NSCC   =           SEGMENT( 6 )
            PKTINFO%TMPPRF =           SEGMENT( 7 )
            PKTINFO%FAC3   = STR2REAL( SEGMENT( 8 ) )
            PKTINFO%CSIC   =           SEGMENT( 9 )
            PKTINFO%PLT    =           SEGMENT( 10 )
            PKTINFO%CHAR1  =           SEGMENT( 11 )
            PKTINFO%CHAR2  =           SEGMENT( 12 )
            PKTINFO%CHAR3  =           SEGMENT( 13 )
            PKTINFO%CHAR4  =           SEGMENT( 14 )
            PKTINFO%CHAR5  =           SEGMENT( 15 )
            PKTINFO%CMCT   =           ' '

            CALL PADZERO( PKTINFO%NSCC )

        CASE( 'PROJECTION' )
            PKTINFO%CFIP  =           SEGMENT( 1 )
            PKTINFO%TSCC  =           SEGMENT( 2 )
            PKTINFO%FAC1  = STR2REAL( SEGMENT( 3 ) )
            PKTINFO%CPOL  =           SEGMENT( 4 )
            PKTINFO%CSIC  =           SEGMENT( 5 )
            PKTINFO%CMCT  =           SEGMENT( 6 )
            PKTINFO%PLT   =           SEGMENT( 7 )
            PKTINFO%CHAR1 =           SEGMENT( 8 )
            PKTINFO%CHAR2 =           SEGMENT( 9 )
            PKTINFO%CHAR3 =           SEGMENT( 10)
            PKTINFO%CHAR4 =           SEGMENT( 11)
            PKTINFO%CHAR5 =           SEGMENT( 12)
            
        CASE( 'MACT' )
            PKTINFO%CFIP    = ' '
            PKTINFO%CSIC    = ' '
            PKTINFO%CMCT    =           SEGMENT( 1 )
            PKTINFO%TSCC    =           SEGMENT( 2 )
            PKTINFO%CSTYP   =           SEGMENT( 3 )
            PKTINFO%APPFLAG = ADJUSTL ( SEGMENT( 4 ) )
            PKTINFO%CPOL    =           SEGMENT( 5 )
            PKTINFO%FAC1    = STR2REAL( SEGMENT( 6 ) )
            PKTINFO%FAC2    = STR2REAL( SEGMENT( 7 ) )
            PKTINFO%FAC3    = STR2REAL( SEGMENT( 8 ) )

        CASE DEFAULT
            MESG = 'INTERNAL ERROR: Packet type ' // PKTTYP // 
     &             ' not recognized in program ' // PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

C.........  Set status of pollutants for current packet
        K = INDEX1( PKTINFO%CPOL, NIPPA, EANAM )
        IF( K .GT. 0 ) THEN
            USEPOL( K ) = .TRUE.

        ELSE IF( PKTINFO%CPOL .EQ. '-9' .OR.
     &           PKTINFO%CPOL .EQ. ' '       ) THEN
            USEPOL = .TRUE.   ! all pollutants

        END IF

        RETURN

C.........  End of file reached unexpectedly
999     MESG = 'End of control packets file reached unexpectedly!'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( A, I8 )
94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE RDPACKET

