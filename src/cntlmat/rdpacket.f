
        SUBROUTINE RDPACKET( FDEV, PKTTYP, FIXEDFMT, USEPOL, IREC, 
     &                       PKTINFO, EFLAG )

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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
         INCLUDE 'CPKTDAT.EXT'   !  control packet contents

C...........   EXTERNAL FUNCTIONS:
        LOGICAL       CHKREAL
        CHARACTER*2   CRLF
        INTEGER       INDEX1
	INTEGER       STR2INT
        REAL          STR2REAL

        EXTERNAL      CHKREAL, CRLF, INDEX1, STR2INT, STR2REAL

C...........   SUBROUTINE ARGUMENTS:
        INTEGER        , INTENT (IN) :: FDEV      ! in file unit number
        CHARACTER(*)   , INTENT (IN) :: PKTTYP    ! packet type
        LOGICAL        , INTENT (IN) :: FIXEDFMT  ! true: packet has fixed fmt
        LOGICAL     , INTENT(IN OUT) :: USEPOL( NIPPA ) ! true: use pollutant
        INTEGER     , INTENT(IN OUT) :: IREC      ! file line number
        TYPE( CPACKET ), INTENT(OUT) :: PKTINFO   ! packet information
        LOGICAL        , INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters
        INTEGER, PARAMETER :: MXSEG = 16   ! number of potential line segments

C...........   Other arrays
        CHARACTER*20 SEGMENT( MXSEG )      ! Segments of parsed packet lines

C...........   Other local variables
        INTEGER         K          ! index
        INTEGER         IOS        ! i/o error status
        INTEGER         CYID       ! tmp county ID
        INTEGER         FIP        ! tmp state/county FIPS code
        INTEGER         STID       ! tmp state ID

        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called

        CHARACTER*8, SAVE :: FMTFIP     ! format for writing FIPS codeSS
        CHARACTER*300        LINE       ! line reading buffer
        CHARACTER*300        MESG       ! message buffer

        CHARACTER*16 :: PROGNAME = 'RDPACKET' ! program name

C***********************************************************************
C   Begin body of subroutine RDPACKET

        IF( FIRSTIME ) THEN
            WRITE( FMTFIP, 94300 ) '(I', FIPLEN3, '.', FIPLEN3, ')'
            FIRSTIME = .FALSE.
        END IF

        READ( FDEV, 93000, END = 999, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF( IOS .GT. 0 ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'ERROR: I/O error', IOS,
     &             'reading control packets file at line', IREC
            CALL M3MESG( MESG )

        END IF

C.........  When packet has a free format...
        IF( .NOT. FIXEDFMT ) THEN

C.............  Parse the line of data into segments based on the rules
C               for "list-formatted" in fortran, but not requiring 
C               quotes around the text strings
            CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.............  Store country/state/county code for any packet format
            PKTINFO%CFIP = SEGMENT( 1 )

        END IF

C.........  Process the line of data, depending on packet type
        SELECT CASE( PKTTYP )

        CASE( 'CTG' )
            PKTINFO%CSIC =           ' '
            PKTINFO%TSCC =           SEGMENT( 2 )
            PKTINFO%CPOL =           SEGMENT( 3 )
            PKTINFO%FAC1 = STR2REAL( SEGMENT( 4 ) )
            PKTINFO%FAC2 = STR2REAL( SEGMENT( 5 ) )
            PKTINFO%FAC3 = STR2REAL( SEGMENT( 6 ) )

C.........  Check to see if last column is blank. If blank, then
C           set PKTINFO%FAC4 = -9 and issue warning.
            IF ( SEGMENT( 7 ) .EQ. ' ' ) THEN
               PKTINFO%FAC4 = -9
               WRITE( MESG, 94020 ) 'RACT value missing from CTG' 
     &         // ' packet record at line', IREC
               CALL M3WARN( PROGNAME, 0, 0, MESG )
            ELSE
               PKTINFO%FAC4 = STR2REAL( SEGMENT( 7 ) )
            END IF

        CASE( 'CONTROL' )
            PKTINFO%TSCC =           SEGMENT( 2 )
            PKTINFO%CPOL =           SEGMENT( 3 )
            PKTINFO%FAC1 = STR2INT ( SEGMENT( 4 ) )
            PKTINFO%FAC2 = STR2REAL( SEGMENT( 5 ) )
            PKTINFO%FAC3 = STR2REAL( SEGMENT( 6 ) )
            PKTINFO%FAC4 = STR2REAL( SEGMENT( 7 ) )
            PKTINFO%CSIC =           SEGMENT( 8 )
            PKTINFO%PLT  =           SEGMENT( 9 )
            PKTINFO%CHAR1=           SEGMENT( 10 )
            PKTINFO%CHAR2=           SEGMENT( 11 )
            PKTINFO%CHAR3=           SEGMENT( 12 )
            PKTINFO%CHAR4=           SEGMENT( 13 )
            PKTINFO%CHAR5=           SEGMENT( 14 )

        CASE( 'ALLOWABLE' )
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

C.........  Check to see if both CAP and REPLACE are missing. If so, issue
C           a warning.
            IF ( PKTINFO%FAC2 .LT. 0 .AND. PKTINFO%FAC3 .LT. 0) THEN
               WRITE( MESG, 94020 ) 'Neither CAP or REPLACE defined'
     &         // ' in allowable packet record at line', IREC
               CALL M3WARN( PROGNAME, 0, 0, MESG )
            END IF

        CASE( 'REACTIVITY' )
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

            CALL PADZERO( PKTINFO%NSCC )

        CASE( 'PROJECTION' )
            PKTINFO%TSCC  =           SEGMENT( 2 )
            PKTINFO%FAC1  = STR2REAL( SEGMENT( 3 ) ) 
            PKTINFO%CPOL  =            ' '          !     SEGMENT( 4 )
            PKTINFO%CSIC  =           SEGMENT( 4 )  !     SEGMENT( 5 )
            PKTINFO%PLT   = ' '
            PKTINFO%CHAR1 = ' '
            PKTINFO%CHAR2 = ' '
            PKTINFO%CHAR3 = ' '
            PKTINFO%CHAR4 = ' '
            PKTINFO%CHAR5 = ' '

        CASE( 'EMS_CONTROL' )
            STID = STR2INT( LINE( 1:2 ) )
            CYID = STR2INT( LINE( 3:5 ) )
            FIP  = STID * 1000 + CYID
            WRITE( PKTINFO%CFIP, FMTFIP ) FIP
            PKTINFO%TSCC = ADJUSTL ( LINE(  10:17  ) )
            PKTINFO%CPOL = ADJUSTL ( LINE(  18:22  ) )
            PKTINFO%FAC1 = STR2REAL( LINE(  82:85  ) )
            PKTINFO%FAC2 = STR2REAL( LINE(  86:89  ) )
            PKTINFO%FAC3 = STR2REAL( LINE(  90:93  ) )
            PKTINFO%FAC4 = STR2REAL( LINE(  94:97  ) ) * 100
            PKTINFO%FAC5 = STR2REAL( LINE(  98:101 ) ) * 100
            PKTINFO%FAC6 = STR2REAL( LINE( 102:105 ) ) * 100
            PKTINFO%FAC7 = STR2REAL( LINE( 106:110 ) )
            PKTINFO%FAC8 = STR2REAL( LINE( 151:154 ) )
            PKTINFO%CSIC = ADJUSTL ( LINE(   6:9   ) )
            PKTINFO%PLT  = ADJUSTL ( LINE(  23:37  ) )
            PKTINFO%CHAR1= ADJUSTL ( LINE(  50:61  ) )
            PKTINFO%CHAR2= ADJUSTL ( LINE(  38:49  ) )
            PKTINFO%CHAR3= ADJUSTL ( LINE(  62:73  ) )
            PKTINFO%CHAR4= ' '
            PKTINFO%CHAR5= ' '
            PKTINFO%NSCC = ADJUSTL ( LINE(  74:81  ) ) 

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

