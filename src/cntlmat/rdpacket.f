
        SUBROUTINE RDPACKET( FDEV, PKTTYP, IREC, PKTINFO, EFLAG )

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
C****************************************************************************/
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
         INCLUDE 'CPKTDAT.EXT'   !  control packet contents

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2   CRLF
        REAL          STR2REAL

        EXTERNAL      CRLF, STR2REAL

C...........   SUBROUTINE ARGUMENTS:
        INTEGER        , INTENT (IN) :: FDEV      ! in file unit number
        CHARACTER(*)   , INTENT (IN) :: PKTTYP    ! packet type
        INTEGER     , INTENT(IN OUT) :: IREC      ! file line number
        TYPE( CPACKET ), INTENT(OUT) :: PKTINFO   ! packet information
        LOGICAL        , INTENT(OUT) :: EFLAG     ! error flag

C...........   Local parameters
        INTEGER, PARAMETER :: MXSEG = 16   ! number of potential line segments

C...........   Other arrays
        CHARACTER*20 SEGMENT( MXSEG )      ! Segments of parsed packet lines

C...........   Other local variables
        INTEGER         IOS        ! i/o error status

        CHARACTER*300   LINE       ! line reading buffer
        CHARACTER*300   MESG       ! message buffer

        CHARACTER*16 :: PROGNAME = 'RDPACKET' ! program name

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

C.........  Parse the line of data into segments based on the rules
C           for "list-formatted" in fortran, but not requiring 
C           quotes around the text strings
        CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.........  Store country/state/county code for any packet format
        PKTINFO%CFIP = SEGMENT( 1 )

C.........  Process the line of data, depending on packet type
        SELECT CASE( PKTTYP )

        CASE( 'CTG' )
            PKTINFO%CSIC =           ' '
            PKTINFO%TSCC =           SEGMENT( 2 )
            PKTINFO%CPOL =           SEGMENT( 3 )
            PKTINFO%FAC1 = STR2REAL( SEGMENT( 4 ) )
            PKTINFO%FAC2 = STR2REAL( SEGMENT( 5 ) )
            PKTINFO%FAC3 = STR2REAL( SEGMENT( 6 ) )
            PKTINFO%FAC4 = STR2REAL( SEGMENT( 7 ) )

        CASE( 'CONTROL' )
            PKTINFO%CSIC =           SEGMENT( 2 )
            PKTINFO%TSCC =           SEGMENT( 3 )
            PKTINFO%CPOL =           SEGMENT( 4 )
            PKTINFO%FAC1 = STR2REAL( SEGMENT( 5 ) )
            PKTINFO%FAC2 = STR2REAL( SEGMENT( 6 ) )
            PKTINFO%FAC3 = STR2REAL( SEGMENT( 7 ) )
            PKTINFO%FAC4 = STR2REAL( SEGMENT( 8 ) )

        CASE( 'ALLOWABLE' )
            PKTINFO%TSCC =           SEGMENT( 2 )
            PKTINFO%CPOL =           SEGMENT( 3 )
            PKTINFO%FAC1 = STR2REAL( SEGMENT( 4 ) )
            PKTINFO%FAC2 = STR2REAL( SEGMENT( 5 ) )
            PKTINFO%FAC3 = STR2REAL( SEGMENT( 6 ) )
            PKTINFO%CSIC =           SEGMENT( 7 )
            PKTINFO%PLT  =           SEGMENT( 8 )
            PKTINFO%CHAR1=           SEGMENT( 8 )
            PKTINFO%CHAR2=           SEGMENT( 9 )
            PKTINFO%CHAR3=           SEGMENT( 10 )
            PKTINFO%CHAR4=           SEGMENT( 11 )
            PKTINFO%CHAR5=           SEGMENT( 12 )

        CASE( 'ADD' )
            PKTINFO%TSCC =           SEGMENT( 2 )
            PKTINFO%CPOL =           SEGMENT( 3 )
            PKTINFO%FAC1 = STR2REAL( SEGMENT( 4 ) )
            PKTINFO%CSIC =           SEGMENT( 5 )
            PKTINFO%PLT  =           SEGMENT( 6 )
            PKTINFO%CHAR1=           SEGMENT( 7 )
            PKTINFO%CHAR2=           SEGMENT( 8 )
            PKTINFO%CHAR3=           SEGMENT( 9 )
            PKTINFO%CHAR4=           SEGMENT( 10 )
            PKTINFO%CHAR5=           SEGMENT( 11 )

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
            PKTINFO%CSIC  =           SEGMENT( 2 ) 
            PKTINFO%CPOL  =           ' '
            PKTINFO%TSCC  =           SEGMENT( 3 )
            PKTINFO%FAC1  = STR2REAL( SEGMENT( 4 ) )
            PKTINFO%PLT   = ' '
            PKTINFO%CHAR1 = ' '
            PKTINFO%CHAR2 = ' '
            PKTINFO%CHAR3 = ' '
            PKTINFO%CHAR4 = ' '
            PKTINFO%CHAR5 = ' '

        CASE DEFAULT
            MESG = 'INTERNAL ERROR: Packet type ' // PKTTYP // 
     &             ' not recognized in program ' // PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

        RETURN

C.........  End of file reached unexpectedly
999     MESG = 'End of control packets file reached unexpectedly!'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDPACKET

