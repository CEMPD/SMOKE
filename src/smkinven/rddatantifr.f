
        SUBROUTINE RDDATAORLFR( LINE, READDATA, READPOL, IYEAR, DESC,
     &                          ERPTYP, SRCTYP, SIC, MACT, NAICS, CTYPE,
     &                          LAT, LON, UTMZ, CORS, BLID, HDRFLAG,
     &                          EFLAG, IREC, NF, NLINEFR, NFRPOL,
     &                          FIREPOL )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an ORL FIRE format wildfire
C      point-source inventory file and returns the unique source characteristics.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by B.H. Baek (02/06) based on rddatantipt.f
C
C**************************************************************************
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

C...........   MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NEM, NDY, NEF, NCE, NRE, NC1, NC2

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)           CRLF
        INTEGER                FINDC
        
        EXTERNAL   CRLF, FINDC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT  (IN) :: LINE                  ! input line
        CHARACTER(*),       INTENT (OUT) :: READDATA( 1,NPTPPOL3 )! array of data values
        CHARACTER(IOVLEN3), INTENT (OUT) :: READPOL( 1 )          ! array of pollutant names
        INTEGER,            INTENT (OUT) :: IYEAR                 ! inventory year
        CHARACTER(40),      INTENT (OUT) :: DESC                  ! plant description
        CHARACTER(ERPLEN3), INTENT (OUT) :: ERPTYP                ! emissions release point type
        CHARACTER(STPLEN3), INTENT (OUT) :: SRCTYP                ! source type code
        CHARACTER(SICLEN3), INTENT (OUT) :: SIC                   ! SIC
        CHARACTER(MACLEN3), INTENT (OUT) :: MACT                  ! MACT code
        CHARACTER(NAILEN3), INTENT (OUT) :: NAICS                 ! NAICS code
        CHARACTER,          INTENT (OUT) :: CTYPE                 ! coordinate type
        CHARACTER(9),       INTENT (OUT) :: LAT                   ! stack latitude
        CHARACTER(9),       INTENT (OUT) :: LON                   ! stack longitude
        CHARACTER(2),       INTENT (OUT) :: UTMZ                  ! UTM zone
        CHARACTER(ORSLEN3), INTENT (OUT) :: CORS                  ! DOE plant ID
        CHARACTER(BLRLEN3), INTENT (OUT) :: BLID                  ! boiler ID
        LOGICAL,            INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,            INTENT (OUT) :: EFLAG                 ! error flag
        INTEGER,            INTENT  (IN) :: IREC                  ! no of current record line buffer
        INTEGER,            INTENT(INOUT):: NF                    ! no of current FIREPOL
        INTEGER,            INTENT  (IN) :: NLINEFR               ! no of lines of current opened file
        INTEGER,            INTENT  (IN) :: NFRPOL                ! no of FIREPOL list
        CHARACTER(CHRLEN3), INTENT  (IN) :: FIREPOL( NFRPOL )     ! extra pol names in wildfire

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 60  ! arbitrary maximum pollutants in file
        INTEGER, PARAMETER :: NSEG = 63      ! number of segments in line

C...........   Other local variables
        INTEGER         I       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPOL    !  number of pollutants in file

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called
 
        CHARACTER(40)      SEGMENT( NSEG ) ! segments of line
        CHARACTER(300)     MESG            ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDDATAORLFR' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAORLFR

C.........  Scan for header lines and check to ensure all are set 
C           properly
        CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .FALSE., 
     &               LINE, ICC, INY, NPOL, IOS )

C.........  Interpret error status
        IF( IOS == 4 ) THEN
            WRITE( MESG,94010 ) 
     &             'Maximum allowed data variables ' //
     &             '(MXPOLFIL=', MXPOLFIL, CRLF() // BLANK10 //
     &             ') exceeded in input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE IF( IOS > 0 ) THEN
            EFLAG = .TRUE.
        END IF

C.........  If a header line was encountered, set flag and return
        IF( IOS >= 0 ) THEN
            HDRFLAG = .TRUE.
            IYEAR = INY
            RETURN
        ELSE
            HDRFLAG = .FALSE.
        END IF

C.........  Separate line into segments
        CALL PARSLINE( LINE, NSEG, SEGMENT )

C.........  Re-define line buffer with other pollutants (CO,NOX,SO2, and others)
        IF( IREC > NLINEFR .AND. IREC <= NLINEFR + NFRPOL ) THEN
            NF = NF + 1
            SEGMENT( 10 ) = ADJUSTL( FIREPOL( NF ) )  ! replace with other #DATA pol name
            SEGMENT( 11 ) = '1.000E-36'               ! replace with a blank for those pol
        END IF

C.........  Use the file format definition to parse the line into
C           the various data fields
        DESC   = ADJUSTL( SEGMENT( 5 ) )   ! plant description
        LAT    = SEGMENT( 6 )              ! stack latitude
        LON    = SEGMENT( 7 )              ! stack longitude
        MACT   = ADJUSTL( SEGMENT( 8 ) )   ! MACT code (NFDRSCODE)
        SIC    = SEGMENT( 9 )              ! SIC (MATBURNED)
        CTYPE  = ADJUSTL( 'L' )            ! fixed coordinate type for wildfire

        ERPTYP = ADJUSTL( ' ' )            ! dummy emissions release point type 
        SRCTYP = ADJUSTL( ' ' )            ! dummy source type code    
        NAICS  = ADJUSTL( ' ' )            ! dummy NAICS code
        UTMZ   = ADJUSTL( ' ' )            ! dummy UTM zone
        CORS   = ADJUSTL( ' ' )            ! dummy DOE plant ID
        BLID   = ADJUSTL( ' ' )            ! dummy boiler ID

        READPOL ( 1     ) = SEGMENT( 10 )  ! name of variable
        READDATA( 1,NEM ) = SEGMENT( 11 )  ! HEATCONTENT Data
        READDATA( 1,NDY ) = '-9'    ! dummy average-day emissions
        READDATA( 1,NEF ) = '-9'    ! dummy emission factor
        READDATA( 1,NCE ) = '-9'    ! dummy control efficiency
        READDATA( 1,NRE ) = '-9'    ! dummy rule effectiveness
        READDATA( 1,NC1 ) = '-9'    ! dummy primary control equipment code
        READDATA( 1,NC2 ) = '-9'    ! dummy secondary control equipment code

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

        END SUBROUTINE RDDATAORLFR
