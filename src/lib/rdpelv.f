
        SUBROUTINE RDPELV( FDEV, NPSRC, MXPELV, CSOURC, NPELV, INDXE )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C	Reads 
C
C  PRECONDITIONS REQUIRED:
C
C  REVISION  HISTORY:
C	Written  1/99 by M. Houyoux
C
C***********************************************************************
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
C****************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'     ! emissions constant parameters
        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2            CRLF
        INTEGER                FINDC

        EXTERNAL    CRLF, FINDC

C...........   ARGUMENTS and their descriptions: actually-occurring ASC table

        INTEGER     FDEV                !  unit number for elev srcs file (in)
        INTEGER     NPSRC               !  no. of point sources (in)
        INTEGER     MXPELV              !  max valid entires in FDEV (in)
        INTEGER     CSOURC( NPSRC )     !  source characteristics (int)
        INTEGER     NPELV               !  actual usable records in FDEV (out)
        INTEGER     INDXE ( MXPELV )    !  elevated sources index (out)

C...........   Variables for source definition input and manipulation
        INTEGER         NCHARS          ! number of source characteristics

        LOGICAL      :: TRUEARR( 7 ) =  ! array of trues
     &                             ( / .TRUE., .TRUE., .TRUE., .TRUE., 
     &                                 .TRUE., .TRUE., .TRUE.        / )

        CHARACTER*300   CHARS( 7 )      ! specific source characteristics

C...........   OTHER LOCAL VARIABLES and their descriptions:

        INTEGER         I, J, K, L, L2   !  counters and indices
        INTEGER         IOS              !  I/O Status
        INTEGER         IREC             !  input line counter

        LOGICAL      :: EFLAG = .FALSE.  !  error flag

        CHARACTER*300   BUFFER           !  buffer for formatted source chars
        CHARACTER*300   MESG             !  message buffer
        CHARACTER(LEN=SRCLEN3) CSRC      !  source characteristics buffer
        CHARACTER(LEN=ALLLEN3) CSRCALL   !  source characteristics buffer // pol

        CHARACTER*16 :: PROGNAME = 'RDPELV' ! program name

C***********************************************************************
C   begin body of subroutine RDPELV

C.........  Read in lines knowing that they are formatted in the CSOURC
C           spacing
        IREC = 0
        J    = 0
        DO           !  head of the FDEV-read loop

            READ( FDEV, 93000, END=23, IOSTAT=IOS ) CSRC
            IREC = IREC + 1

            IF( IOS .GT. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'ERROR: I/O error', IOS,
     &              'reading elevated sources file at line', IREC
                CALL M3MESG( MESG )

            END IF

            CHARS = ' '  ! Array

C.............  Separare CSRC array into components
            CALL PARSCSRC( CSRC, TRUEARR, CHARS, NCHARS )

            CALL BLDCSRC( CHARS(1), CHARS(2), CHARS(3), CHARS(4),
     &                    CHARS(5), CHARS(6), CHARS(7), ' ', CSRCALL )

            L = LEN_TRIM( CSRCALL )
            I = FINDC( CSRCALL( 1:L ), NPSRC, CSOURC )

            IF ( I .LE. 0 ) THEN

                CALL FMTCSRC ( CSRCALL, 7, BUFFER, L2 )

                MESG = 'Input elevated source not found in inventory:'//
     &                  CRLF() // BLANK10 // BUFFER( 1:L2 )

                CALL M3MESG( MESG )

C.............  Ensure there is no overflow of INDXE
            ELSEIF( J .LE. MXPELV ) THEN
                J = J + 1
                INDXE( J ) = I

            END IF 

        ENDDO       !  to head of elevated sources loop

23      CONTINUE    !  end of read loop

        NPELV = J

C.........  Report warning for skipped sources from elevated sources file
        IF( NPELV .LT. MXPELV ) THEN

            J = MXPELV - NPELV
            WRITE( MESG,94010 ) 'Program is ignoring', J, 'records ' //
     &                          'in the input elevated sources file.'
            CALL M3WARN( PROGNAME, 0, 0, MESG )

        ENDIF

C.........  Check for overflow...
        IF( J .GT. MXPELV ) THEN 

            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Memory allocation insufficient ' //
     &             'for elevated point sources file'
            CALL M3MSG2( MESG )

        ENDIF

        IF( EFLAG ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93020   FORMAT( I7, I3 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )


        END SUBROUTINE RDPELV

