
        SUBROUTINE RDDATAORLPT( LINE, READDATA, READPOL, IYEAR, DESC,
     &                          ERPTYP, SRCTYP, HT, DM, TK, FL, VL, SIC, 
     &                          MACT, NAICS, CTYPE, LAT, LON, UTMZ, 
     &                          NEID, CORS, BLID, FUGHT, FUGAR,
     &                          EXTORL, HDRFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an ORL format point-source inventory
C      file and returns the inventory data values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created by C. Seppanen (01/03) based on rddataidapt.f
C
C       Version June 2016 by Carlie Coats:  add fugitive-emissions properties
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
        CHARACTER(16),      INTENT (OUT) :: HT                    ! stack height
        CHARACTER(16),      INTENT (OUT) :: DM                    ! stack diameter
	CHARACTER(16),      INTENT (OUT) :: TK                    ! exit temperature
        CHARACTER(16),      INTENT (OUT) :: FL                    ! flow rate
        CHARACTER(16),      INTENT (OUT) :: VL                    ! exit velocity
        CHARACTER(SICLEN3), INTENT (OUT) :: SIC                   ! SIC
        CHARACTER(MACLEN3), INTENT (OUT) :: MACT                  ! MACT code
        CHARACTER(NAILEN3), INTENT (OUT) :: NAICS                 ! NAICS code
        CHARACTER,          INTENT (OUT) :: CTYPE                 ! coordinate type
        CHARACTER(16),      INTENT (OUT) :: LAT                   ! stack latitude
        CHARACTER(16),      INTENT (OUT) :: LON                   ! stack longitude
        CHARACTER(2),       INTENT (OUT) :: UTMZ                  ! UTM zone
        CHARACTER(NEILEN3), INTENT (OUT) :: NEID                  ! NEI unique ID
        CHARACTER(ORSLEN3), INTENT (OUT) :: CORS                  ! DOE plant ID
        CHARACTER(BLRLEN3), INTENT (OUT) :: BLID                  ! boiler ID
        CHARACTER(16),      INTENT (OUT) :: FUGHT                 ! RELEASE_HEIGHT_FUGITIVE
        CHARACTER(16),      INTENT (OUT) :: FUGAR                 ! HORIZONTAL_AREA_FUGITIVE
        CHARACTER(EXTLEN3), INTENT (OUT) :: EXTORL                ! additional ext vars
        LOGICAL,            INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,            INTENT (OUT) :: EFLAG                 ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 60  ! arbitrary maximum pollutants in file
        INTEGER, PARAMETER :: NSEG = 75      ! number of segments in line

C...........   Other local variables
        INTEGER         I, L, L1, LL       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPOL    !  number of pollutants in file

        LOGICAL, SAVE:: FIRSTIME = .TRUE.  ! true: first time routine is called
        LOGICAL      :: BLKFLAG  = .TRUE.  ! true when it is blank
 
        CHARACTER(40)      TMPSEG          ! tmp segments of line
        CHARACTER(40)      SEGMENT( NSEG ) ! segments of line
        CHARACTER(CASLEN3) TCAS            ! tmp cas number
        CHARACTER(300)     MESG            ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDDATAORLPT' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAORLPT

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

C.........  Use the file format definition to parse the line into
C           the various data fields
        DESC   = ADJUSTL( SEGMENT( 6 ) )   ! plant description
        ERPTYP = ADJUSTL( SEGMENT( 8 ) )   ! emissions release point type 
        SRCTYP = ADJUSTL( SEGMENT( 9 ) )   ! source type code    
        HT     = SEGMENT( 10 )             ! stack height
        DM     = SEGMENT( 11 )             ! stack diameter
        TK     = SEGMENT( 12 )             ! exit temperature
        FL     = SEGMENT( 13 )             ! flow rate
        VL     = SEGMENT( 14 )             ! exit velocity
        SIC    = SEGMENT( 15 )             ! SIC
        MACT   = ADJUSTL( SEGMENT( 16 ) )  ! MACT code
        NAICS  = ADJUSTL( SEGMENT( 17 ) )  ! NAICS code
        CTYPE  = ADJUSTL( SEGMENT( 18 ) )  ! coordinate type
        LON    = SEGMENT( 19 )             ! stack longitude
        LAT    = SEGMENT( 20 )             ! stack latitude
        UTMZ   = ADJUSTL( SEGMENT( 21 ) )  ! UTM zone

        READPOL ( 1     ) = SEGMENT( 22 )
        READDATA( 1,NEM ) = SEGMENT( 23 ) ! annual emissions
        READDATA( 1,NDY ) = SEGMENT( 24 ) ! average-day emissions
        READDATA( 1,NEF ) = ' '           ! emission factor
        READDATA( 1,NCE ) = SEGMENT( 25 ) ! control efficiency
        READDATA( 1,NRE ) = SEGMENT( 26 ) ! rule effectiveness
        READDATA( 1,NC1 ) = SEGMENT( 27 ) ! primary control equipment code
        READDATA( 1,NC2 ) = SEGMENT( 28 ) ! secondary control equipment code

        NEID   = ADJUSTL( SEGMENT( 29 ) )  ! NEI Unique ID
        IF( NEID == ' ' ) NEID = '-9'
        CORS   = ADJUSTL( SEGMENT( 30 ) )  ! DOE plant ID
        BLID   = ADJUSTL( SEGMENT( 31 ) )  ! boiler ID
        FUGAR  = SEGMENT( 38 )
        FUGHT  = SEGMENT( 39 )

C.........  Read extended orl variables and store it as string
        EXTORL = ' '
        DO I = 32, 70
            IF( SEGMENT( I ) == ' ' ) THEN
                TMPSEG = ','
            ELSE
                TMPSEG = ',' // TRIM( SEGMENT( I ) )
                BLKFLAG = .FALSE.
            ENDIF
            
            EXTORL = TRIM( EXTORL ) // TRIM( TMPSEG )
        END DO

        IF( BLKFLAG ) EXTORL = ' '

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

        END SUBROUTINE RDDATAORLPT
