
        SUBROUTINE RDDATAORLFR( LINE, READDATA, READPOL, NDATPERLN, 
     &                          IYEAR, DESC, SIC, MACT, CTYPE, 
     &                          LAT, LON, HDRFLAG, EFLAG )

C***********************************************************************
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
        USE MODINFO, ONLY: NEM, NDY, NEF, NCE, NRE, NC1, NC2, TMPNAM

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)           CRLF
        INTEGER                FINDC
        
        EXTERNAL   CRLF, FINDC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*),       INTENT  (IN) :: LINE                  ! input line
        CHARACTER(*),       INTENT (OUT) :: READDATA( NDATPERLN,NPTPPOL3 )! array of data values
        CHARACTER(IOVLEN3), INTENT (OUT) :: READPOL( NDATPERLN )          ! array of pollutant names
        INTEGER,            INTENT(INOUT):: NDATPERLN             ! number of data values per line
        INTEGER,            INTENT (OUT) :: IYEAR                 ! inventory year
        CHARACTER(40),      INTENT (OUT) :: DESC                  ! plant description
        CHARACTER(SICLEN3), INTENT (OUT) :: SIC                   ! Material burned code (stored in SIC)
        CHARACTER(MACLEN3), INTENT (OUT) :: MACT                  ! NFDRS code (stored in MACT)
        CHARACTER,          INTENT (OUT) :: CTYPE                 ! coordinate type
        CHARACTER(9),       INTENT (OUT) :: LAT                   ! stack latitude
        CHARACTER(9),       INTENT (OUT) :: LON                   ! stack longitude
        LOGICAL,            INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,            INTENT (OUT) :: EFLAG                 ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 1000  ! arbitrary maximum pollutants in file
        INTEGER, PARAMETER :: NSEG = 10      ! number of segments in line for format
        INTEGER, PARAMETER :: NEXTRA  = 4    ! number of extra non-data fields that need
                                             ! to be added as "pollutants)
        CHARACTER(IOVLEN3), PARAMETER :: FIREVNAM( NEXTRA ) = ! fire variable names
     &                      ( / 'HEATCONTENT     ',
     &                          'HFLUX           ',
     &                          'ENDHOUR         ',
     &                          'BEGHOUR         '  / )

C...........   Other local variables
        INTEGER         I, L       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NDAT = -1 !  number of data values as set by header

        LOGICAL, SAVE:: FIRSTDATA = .TRUE. ! true: first time data row in encountered
 
        CHARACTER(40)      SEGMENT( NSEG ) ! segments of line
        CHARACTER(300)     MESG            ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDDATAORLFR' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAORLFR

C.........  Scan for header lines and check to ensure all are set 
C           properly
        CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .TRUE., 
     &               LINE, ICC, INY, NDAT, IOS )

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
            IF( NDAT > 0 ) NDATPERLN = NDAT + NEXTRA  
            RETURN
        ELSE
            HDRFLAG = .FALSE.
        END IF

C.........  Give error if #DATA line has not defined the pollutants that
C           are contained in the day-specific data file. NOTE: This code
C           will not be reached until after the last header line)
        IF ( NDAT < 0 ) THEN
            WRITE( MESG, 94010 ) 'First data line reached in ORL '//
     &             'FIRE file with required #DATA header found.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, the READPOL array has been redfined and now can
C           be populated with the TMPNAM array set by the GETHDR routine
        ELSE IF ( FIRSTDATA ) THEN

            FIRSTDATA = .FALSE.
            READPOL( 1:NEXTRA ) = FIREVNAM( 1:NEXTRA ) 
            READPOL( NEXTRA+1:NDATPERLN ) = TMPNAM( 1:NDAT )  ! array

        END IF

C.........  Separate line into segments
        CALL PARSLINE( LINE, NSEG, SEGMENT )

C.........  Use the file format definition to parse the line into
C           the various data fields
        DESC   = SEGMENT( 5 )              ! fire description
        LAT    = SEGMENT( 6 )              ! fire latitude
        LON    = SEGMENT( 7 )              ! fire longitude
        MACT   = SEGMENT( 8 )              ! MACT code (NFDRSCODE field)
        SIC    = SEGMENT( 9 )              ! SIC (MATBURNED field)
        CTYPE  = 'L'                       ! lat-lon coordinate type part of format

C.........  Populate all of the data fields with dummy values
        DO I = 1, NDATPERLN
            READDATA( I,NEM ) = '0.0'   ! dummy annual emissions
            READDATA( I,NDY ) = '-9'    ! dummy average-day emissions
            READDATA( I,NEF ) = '-9'    ! dummy emission factor
            READDATA( I,NCE ) = '-9'    ! dummy control efficiency
            READDATA( I,NRE ) = '-9'    ! dummy rule effectiveness
            READDATA( I,NC1 ) = '-9'    ! dummy primary control equipment code
            READDATA( I,NC2 ) = '-9'    ! dummy secondary control equipment code
        END DO

C.........  Populate heat content value
        L = LEN( READDATA( 1, NEM ) )
        READDATA( 1, NEM ) = SEGMENT( 10 )( 1:L )
      
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
