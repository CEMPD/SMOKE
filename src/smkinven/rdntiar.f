
        SUBROUTINE RDNTIAR( FDEV )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine reads a fake National Toxics inventory format
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      created by M. Houyoux (03/2002)
C
C**************************************************************************
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

C...........   MODULES for public variables
C...........   This module is the point source inventory arrays
c        USE MODSOURC

C.........  This module contains the lists of unique inventory information
c        USE MODLISTS

C.........  This module contains the arrays for state and county summaries
c        USE MODSTCY

C.........  This module contains the information about the source category
c        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
         INCLUDE 'CONST3.EXT'    !  physical and mathematical constants
         INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
      	LOGICAL                CHKINT
        LOGICAL                CHKREAL
        CHARACTER*2            CRLF
        INTEGER                ENVINT
        LOGICAL                ENVYN
        INTEGER                GETFLINE
        INTEGER                GETNLIST
        INTEGER                INDEX1
        INTEGER                STR2INT
        REAL                   STR2REAL
        REAL                   YR2DAY 

        EXTERNAL    CRLF, ENVINT, ENVYN, GETFLINE, GETNLIST, INDEX1, 
     &              STR2INT, STR2REAL, YR2DAY

C...........   SUBROUTINE ARGUMENTS
C...........   NOTE that NDROP and EDROP are not used at present
        INTEGER     , INTENT (IN) :: FDEV   ! unit number of input file

C...........   Local parameters, indpendent
c        INTEGER, PARAMETER :: MXPOLFIL = 63  ! maximum pollutants in file
c        INTEGER, PARAMETER :: AROTWIDE = 47  ! total width of all pol fields
c        INTEGER, PARAMETER :: ARNONPWD = 15  ! width of non-pol fields

C...........   Local parameters, dependent
c        INTEGER, PARAMETER :: LINSIZ  = ARNONPWD + MXPOLFIL * AROTWIDE

C...........   Local parameter arrays...
C...........   Start and end positions in the file format of the first set
C              of pollutant fields.
c        INTEGER, PARAMETER :: ISINIT( NARPPOL3 ) = 
c     &                              ( / 16,26,36,47,54,57 / )

c        INTEGER, PARAMETER :: IEINIT( NARPPOL3 ) = 
c     &                              ( / 25,35,46,53,56,62 / )

C...........   Local arrays
c        INTEGER          IS( NARPPOL3 )  ! start position for each pol char
c        INTEGER          IE( NARPPOL3 )  ! end position for each pol char

        CHARACTER*15 :: SEGMENT( 3 )

C...........   Counters of total number of input records
        INTEGER, SAVE :: NSRCSAV = 0 ! cumulative source count
        INTEGER, SAVE :: NSRCPOL = 0 ! cumulative source x pollutants count

C...........   Other local variables
        INTEGER         I, J, K, L, N, V  ! counters and indices

        INTEGER         ES      !  counter for source x pollutants
        INTEGER         FIP     !  tmp FIPS code
        INTEGER         PFIP     !  tmp FIPS code
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  line counter
        INTEGER         NPOL    !  number of pollutants per source
        INTEGER         NSRC    !  number of sources
        INTEGER         MNPOL   !  min pollutant per source
        INTEGER         MXPOL   ! max pollutant per source
        INTEGER         POLID   ! tmp pollutant ID

        LOGICAL, SAVE:: EFLAG   = .FALSE.
        LOGICAL, SAVE:: FFLAG    = .FALSE. ! true: fill in 0. annual with seasonal
        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called

        CHARACTER*300   MESG    !  message buffer

        CHARACTER(LEN=POLLEN3) CCOD  ! character pollutant index to INVDNAM
        CHARACTER(LEN=FIPLEN3) CFIP  ! character FIP code
        CHARACTER(LEN=IOVLEN3) CBUF  ! tmp pollutant code
        CHARACTER(LEN=256)  LINE  ! input line from inventory file
        CHARACTER(LEN=SCCLEN3) TSCC  ! tmp scc
        CHARACTER(LEN=SCCLEN3) PSCC  ! previous scc

        CHARACTER*16 :: PROGNAME = 'RDNTIAR' ! Program name

C***********************************************************************
C   begin body of subroutine RDNTIAR

C.........  Reinitialize for multiple subroutine calls
        EFLAG = .FALSE.

C........................................................................
C.............  Head of the main read loop  .............................
C........................................................................

        IREC = 0
        NSRC = 0
        MXPOL = 0
        MNPOL = 200
        NPOL = 0
        PFIP = 0
        PSCC = '-9'
        DO

C.............  Read a line of IDA file as a character string
            READ( FDEV, 93000, END=199, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading inventory file at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            END IF

            L = LEN_TRIM( LINE )  ! store width of line and check

C.............  Skip blank lines and first line
            IF( IREC .EQ. 1 .OR. L .EQ. 0 ) CYCLE

C.............  Scan for header lines and check to ensure all are set 
C               properly
c            CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .TRUE., 
c     &                   LINE, ICC, INY, NPOL, IOS )

C.............  Interpret error status
c            IF( IOS .EQ. 4 ) THEN
c                WRITE( MESG,94010 ) 
c     &                 'Maximum allowed data variables ' //
c     &                 '(MXPOLFIL=', MXPOLFIL, CRLF() // BLANK10 //
c     &                 ') exceeded in input file'
c                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

c            ELSE IF( IOS .GT. 0 ) THEN
c                EFLAG = .TRUE.

c            END IF

C.............  If a header line was encountered, go to next line
c            IF( IOS .GE. 0 ) CYCLE

C.............  Parse line
C.............  NOTE: This is assuming the version that I created with
C               n: just FIPS, SCC, and POLID for computing stats, which
C               n: has already been sorted.
            CALL PARSLINE( LINE, 3, SEGMENT )

C.............  Convert information as needed
            FIP   = STR2INT( SEGMENT( 1 ) )
            TSCC  = SEGMENT( 2 )
            POLID = STR2INT( SEGMENT( 3 ) )

C.............  If source is new
            IF( FIP .NE. PFIP  .OR. TSCC .NE. PSCC ) THEN

C.................  Write out previous source number of pollutants
C.................  and compute stats
                IF( NSRC .GT. 0 ) THEN
                    WRITE( 70, ('(I8)') ) NPOL
                    MXPOL = MAX( MXPOL, NPOL )
                    MNPOL = MIN( MNPOL, NPOL )
                END IF

                NSRC = NSRC + 1
                NPOL = 1

C.............  If same source, count pollutants
            ELSE
                NPOL = NPOL + 1
            END IF

            PFIP = FIP
            PSCC = TSCC

        END DO          !  to head of FDEV-read loop

199     CONTINUE        !  exit from the FDEV-read loop

        WRITE( 70,('(I8)') ) NPOL

        print *, 'min pol per source=', MNPOL
        print *, 'max pol per source=', MXPOL

        CALL M3EXIT( PROGNAME, 0, 0, 'NTI read stop', 0)

        CLOSE( FDEV )

C.........  Make sure routine knows it's been called already
        FIRSTIME = .FALSE.

C.........  Return from subroutine 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

        END SUBROUTINE RDNTIAR
