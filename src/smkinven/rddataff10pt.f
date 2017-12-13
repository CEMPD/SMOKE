
        SUBROUTINE RDDATAFF10PT( LINE, READDATA, READPOL, IYEAR, DESC,
     &                          ERPTYP, SRCTYP, HT, DM, TK, FL, VL, SIC,
     &                          MACT, NAICS, CTYPE, LAT, LON, UTMZ, 
     &                          NEID, CORS, BLID,  
     &                          FUGHT, FUGWD, FUGLN, FUGAN,
     &                          EXTORL, HDRFLAG,
     &                          AVEFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 156
C
C  DESCRIPTION:
C      This subroutine processes a line from an FF10 format point-source inventory
C      file and returns the inventory data values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created by Dongmei Yang (Oct, 2011) based on rddatantipt.f
C       Version June 2016 by Carlie Coats:  add fugitive-emissions p
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
        USE MODINFO, ONLY: NEM, NDY, NEF, NCE, NRE, NC1, NC2, INV_MON

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: FF10INVFLAG

        IMPLICIT NONE

C...........   INCLUDES
         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        INTEGER         FINDC, STR2INT
        REAL            YR2DAY, STR2REAL
        LOGICAL         CHKINT
        
        
        EXTERNAL   CRLF, FINDC, STR2INT, STR2REAL, YR2DAY, CHKINT

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
        CHARACTER(16),      INTENT (OUT) :: FUGHT                 ! FUG_HEIGHT
        CHARACTER(16),      INTENT (OUT) :: FUGWD                 ! FUG_WIDTHT
        CHARACTER(16),      INTENT (OUT) :: FUGLN                 ! FUG_LENGTH
        CHARACTER(16),      INTENT (OUT) :: FUGAN                 ! FUG_ANGLE
        CHARACTER(EXTLEN3), INTENT (OUT) :: EXTORL                ! additional ext vars
        LOGICAL,            INTENT (OUT) :: HDRFLAG               ! true: line is a header line
        LOGICAL,            INTENT (OUT) :: AVEFLAG               ! true: Aveday inv is processed
        LOGICAL,            INTENT (OUT) :: EFLAG                 ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXPOLFIL = 60  ! arbitrary maximum pollutants in file
        INTEGER, PARAMETER :: NSEG = 80      ! number of segments in line

C...........   Other local variables
        INTEGER         I, L, L1, LL       ! counters and indices

        INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
        INTEGER, SAVE:: INY     !  inventory year
        INTEGER         IOS     !  i/o status
        INTEGER, SAVE:: NPOL    !  number of pollutants in file
        INTEGER      :: LYEAR   !  Leap year (366 days per year)
        INTEGER      :: MDAYS   !  days of modeling inventory month

        REAL         :: AVEINV  !  annual total estimate from monthly total VMT

        LOGICAL, SAVE:: FIRSTIME = .TRUE.  ! true: first time routine is called
        LOGICAL      :: BLKFLAG  = .TRUE.  ! true when it is blank
 
        CHARACTER(40)      TMPSEG          ! tmp segments of line
        CHARACTER(40)      SEGMENT( NSEG ) ! segments of line
        CHARACTER(CASLEN3) TCAS            ! tmp cas number
        CHARACTER(300)     MESG            ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDDATAFF10PT' ! Program name

C***********************************************************************
C   begin body of subroutine RDDATAFF10PT

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

C......... Return if the first line is a header line
        IF( .NOT. CHKINT( SEGMENT( 2 ) ) ) THEN
            HDRFLAG = .TRUE.
            RETURN
        END IF 

C.........  Use the file format definition to parse the line into
C           the various data fields
        DESC   = ADJUSTL( SEGMENT( 16 ) )  ! plant description
        ERPTYP = ADJUSTL( SEGMENT( 17 ) )  ! emissions release point type 
        SRCTYP = ADJUSTL( SEGMENT( 31 ) )  ! source type code    
        HT     = SEGMENT( 18 )             ! stack height
        DM     = SEGMENT( 19 )             ! stack diameter
        TK     = SEGMENT( 20 )             ! exit temperature
        FL     = SEGMENT( 21 )             ! flow rate
        VL     = SEGMENT( 22 )             ! exit velocity
        SIC    = ADJUSTL( SEGMENT( 30 ) )  ! SIC (optional in regulatory code field)
        MACT   = ""                        ! MACT (retired in FF10) 
        NAICS  = ADJUSTL( SEGMENT( 23 ) )  ! NAICS code
        CTYPE  = "L"                       ! coordinate type (default:lat/lon in FF10)
        LON    = SEGMENT( 24 )             ! stack longitude
        LAT    = SEGMENT( 25 )             ! stack latitude
        UTMZ   = ""                        ! UTM zone (n/a:Lat/Lon only in FF10)
 
        FUGHT = SEGMENT( 47 )
        FUGWD = SEGMENT( 48 )
        FUGLN = SEGMENT( 49 )
        FUGAN = SEGMENT( 50 )

        READPOL ( 1     ) = SEGMENT( 13 )
        READDATA( 1,NEM ) = SEGMENT( 14 )  ! annual emissions
        READDATA( 1,NDY ) = ''             ! average-day emissions
        READDATA( 1,NEF ) = '-9'           ! emission factor
        READDATA( 1,NCE ) = SEGMENT( 15)   ! control efficiency
        READDATA( 1,NRE ) = '100'          ! rule effectiveness
        READDATA( 1,NC1 ) = '-9'           ! primary control equipment code
        READDATA( 1,NC2 ) = '-9'           ! secondary control equipment code

        NEID   = '-9'                      ! NEI Unique ID
        CORS   = ADJUSTL( SEGMENT( 42 ) )  ! DOE plant ID
        BLID   = ADJUSTL( SEGMENT( 43 ) )  ! boiler ID

        EXTORL = ' '   ! extended orl (N/A)

C.........  Compute annual total based on monthly total
        AVEINV = 0.0
        AVEFLAG = .FALSE.
        IF( INV_MON > 0 ) THEN

            DO I = 1, 12
                IF( LEN_TRIM( SEGMENT( 52+I ) ) < 1 ) THEN
                    SEGMENT( 52+I ) = '0.0'
                END IF
               AVEINV = AVEINV + STR2REAL( SEGMENT( 52+I ) )
            END DO

            IF( AVEINV > 0.0 ) AVEFLAG = .TRUE.

            IF( AVEFLAG ) THEN

                READDATA( 1,NEM ) = '0.0'
                READDATA( 1,NDY ) = SEGMENT( 52 + INV_MON )

                MDAYS = MON_DAYS( INV_MON )    ! day of months

                LYEAR = INT( 1 / YR2DAY ( INY ) )                ! convert year to days
                IF( LYEAR > 365 .AND. INV_MON == 2 ) MDAYS = 29  ! leap year (feb = 29days)

                AVEINV = STR2REAL( READDATA(1,NDY) ) / MDAYS
                WRITE( READDATA( 1,NDY ), '( E15.10 )' ) AVEINV

                IF( AVEINV < 0.0 ) THEN
                    MESG = 'ERROR: Can not process negative value'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END IF

        END IF

C.........  Reset annual total inventory to zero for daily/hourly FF10 processing
        IF( FF10INVFLAG ) THEN
            READPOL ( 1     ) = SEGMENT( 9 )
            READDATA( 1,NEM ) = ''
            READDATA( 1,NDY ) = ''
        END IF

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

        END SUBROUTINE RDDATAFF10PT
