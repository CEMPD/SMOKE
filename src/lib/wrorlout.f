
        SUBROUTINE WRORLOUT( RDEV, DATNAM, NREC, NPVAR, SRCID, SRCDAT, 
     &                       STATUS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Write ORL output file
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux, U.S. EPA 6/8/2005
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: CIFIP, INVYR, XLOCA, YLOCA,
     &                      CORIS, STKHT, STKDM, STKTK, STKVE,
     &                      CSCC, CSOURC, CPDESC, CLINK, CBLRID,
     &                      CERPTYP, CMACT, CNAICS, CSRCTYP, CNEIUID,
     &                      CEXTORL, CISIC

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTRY, CTRYCOD, CTRYNAM

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NINVTBL, ITNAMA, ITCASA, ITFACA

C...........   This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NIPPA, NIPOL, NIACT,
     &                     NPPOL, NPACT, EANAM, EAUNIT, MXCHRS,
     &                     SC_BEGP, SC_ENDP, ACTVTY

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical constants
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS:
C       CHARACTER*2     CRLF
C       LOGICAL         ENVYN
C       INTEGER         FINDC
C       INTEGER         INDEX1
C       INTEGER         STR2INT

C        EXTERNAL        CRLF, ENVYN, FINDC, INDEX1, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER      , INTENT (IN) :: RDEV           ! emissions unit no.
        CHARACTER(*) , INTENT (IN) :: DATNAM         ! name of output pollutant
        INTEGER      , INTENT (IN) :: NREC           ! number of sources
        INTEGER      , INTENT (IN) :: NPVAR          ! no. pollutant vars per src
        INTEGER      , INTENT (IN) :: SRCID ( NREC ) ! Source IDs
        REAL         , INTENT (IN) :: SRCDAT( NREC, NPVAR ) ! emissions + attrb
        INTEGER      , INTENT(OUT) :: STATUS         ! exit status

C...........   ORL output variables (names same as ORL format description)

        INTEGER         FIP
        INTEGER         CPRI     ! primary control device code
        INTEGER         CSEC     ! secondary control device code
        INTEGER         UTMZ     ! tmp default UTM zone

        REAL            STKHGT, STKDIAM, STKTEMP, STKFLOW, STKVEL
        REAL            XLOC, YLOC, ANN_EMIS, AVD_EMIS, CEFF, REFF, RPEN

        CHARACTER*1        CTYPE          !  tmp coordinate system type
        CHARACTER(CASLEN3) CAS
        CHARACTER(MACLEN3) MACT
        CHARACTER(NAILEN3) NAICS
        CHARACTER(CHRLEN3) POINTID
        CHARACTER(PLTLEN3) PLANTID  
        CHARACTER(DSCLEN3) PLNTDESC
        CHARACTER(SCCLEN3) SCC
        CHARACTER(SICLEN3) SIC
        CHARACTER(CHRLEN3) SEGMENT
        CHARACTER(STPLEN3) SRCTYPE
        CHARACTER(CHRLEN3) STACKID
        CHARACTER(ERPLEN3) ERPTYPE
        CHARACTER(ORSLEN3) CORS   ! temporary DOE plant ID
        CHARACTER(BLRLEN3) CBLR   ! temporary boiler name
        CHARACTER(NEILEN3) CNEI   ! NEI unique ID
        CHARACTER(EXTLEN3) :: CEXT = ''   ! Extended ORL vars

C...........   Other local variables

        INTEGER         C, I, J, L, L1, L2, K, N, S  ! counters and indices

        INTEGER         COID     ! tmp country code
        INTEGER         FDEV     ! tmp unit no.
        INTEGER         IOS      ! i/o status
        INTEGER,SAVE :: LCOID = -9    ! previous country ID in loop
        INTEGER,SAVE :: LYEAR = -9    ! previous year in loop
        INTEGER         NCHAR    ! number of strings returned from PARSCSRC
        INTEGER         YEAR     ! tmp 4-digit year

        REAL            M2FT     ! meters to feet

        LOGICAL, SAVE :: FIRSTIME = .TRUE.
        LOGICAL, SAVE :: LNONRD   = .FALSE.
        LOGICAL ORLCOLS( 7 )
        DATA    ORLCOLS / 6*.TRUE., .FALSE. /

        CHARACTER(4)   CYEAR          !  character 4-digit year
        CHARACTER(128) CHARS( 7 )     !  source fields for output
        CHARACTER(256) POLBUF         !  pollutant list buffer
        CHARACTER(256) MESG           !  message buffer

        CHARACTER(16) :: PROGNAME = 'WRORLOUT' ! program name

C***********************************************************************
C   begin body of subroutine WRORLOUT


C.........  Initializations
        STATUS = 0

C.........  Write warning if activities are in inventory.
        IF( NIACT .GT. 0 ) THEN
            MESG = 'WARNING: Activities in inventory cannot be '//
     &             'written to ORL format because' //CRLF()//
     &             BLANK10 // 'activities are not supported in '//
     &             'ORL format.'
            CALL M3MSG2( MESG )
        END IF

C.........  Write each ORL record
        SELECT CASE( CATEGORY )

C............................................................................
C.........  For area sources...
        CASE ( 'AREA' )

            IF ( FIRSTIME ) THEN

C.................  Retrieve environment variable to optionally write out ORL 
C                   Nonroad format instead of ORL nonpoint.
                MESG = 'Output Nonroad ORL format instead of nonpoint'        
                LNONRD = ENVYN ( 'ORL_NONROAD_OUT', MESG, .FALSE., IOS )

C.................  Write initial header
                IF( LNONRD ) THEN
                    WRITE( RDEV, 93000 ) 
     &                  '#ORL NONROAD',
     &                  '#TYPE     Nonroad Inventory',
     &                  '#DESC        (Output from SMOKE)'
                ELSE
                    WRITE( RDEV, 93000 ) 
     &                  '#ORL NONPOINT',
     &                  '#TYPE     Nonpoint Inventory',
     &                  '#DESC        (Output from SMOKE)'
                END IF

                FIRSTIME = .FALSE.

            END IF

C.............  Write area-source characteristics to output file
            DO C = 1, NREC

                S = SRCID( C )

C.................  Store others in temporary variables
                COID = STR2INT( CIFIP( S ) ) / 100000
                FIP  = STR2INT( CIFIP( S ) ) - COID * 100000
                SIC  = CISIC( S )
                YEAR = INVYR( S )
                SCC  = CSCC ( S )
 
                IF( ASSOCIATED( CMACT ) ) THEN
                    MACT = CMACT( S )
                ELSE
                    MACT = '-9'
                END IF

                IF( ASSOCIATED( CSRCTYP ) ) THEN
                    SRCTYPE = CSRCTYP( S )
                ELSE
                    SRCTYPE = '-9'
                END IF

                IF( ASSOCIATED( CNAICS ) ) THEN
                    NAICS = CNAICS( S )
                ELSE
                    NAICS = '-9'
                END IF

                IF( ASSOCIATED( CEXTORL ) ) THEN
                    CEXT = ADJUSTL( CEXTORL( S ) )
                ELSE
                    CEXT = ''
                END IF

C.................  Account for missing or default codes
                IF ( SCC(1:2) .EQ. '00' ) SCC = SCC(3:SCCLEN3)
                IF ( MACT .EQ. '000000' ) MACT = '-9'
                IF ( MACT(1:2) .EQ. '00' ) MACT = MACT(3:MACLEN3)
                IF ( LEN( TRIM( MACT ) ) .EQ. 0 ) MACT = '-9'
                IF ( LEN( TRIM( SRCTYPE ) ) .EQ. 0 ) SRCTYPE = '-9'
                IF ( NAICS .EQ. '000000' ) NAICS = '-9'
                IF ( LEN( TRIM( NAICS ) ) .EQ. 0 ) NAICS = '-9'
                IF ( LEN( TRIM( CEXT ) ) .EQ. 0 ) CEXT = ''

C.................  Retrieve pollutant code from Inventory Table
C.................  Ensure that only CAS codes where the factor is 1 are used.
C.................  If there are multiple CAS codes where the factor is 1, use the first one
C.................     (accomplished since original INVTABLE list order is maintained already)
                CAS = ' '
                DO I = 1, NINVTBL

                    IF ( DATNAM .EQ. ITNAMA( I ) .AND. ITFACA( I ) .EQ. 1. ) THEN
                        CAS = ITCASA( I )
                        EXIT              ! Exit loop
                    END IF

                END DO           

                IF ( CAS .EQ. ' ' ) THEN
                    MESG = 'ERROR: CAS number for pollutant "'// 
     &                     TRIM( DATNAM ) // '" and factor = 1 in '//
     &                     'INVTABLE is no found.'
                    CALL M3MSG2( MESG )
                    STATUS = 2 
                    CYCLE
                END IF

C.................  Get emissions and emissions-dependent values.
                ANN_EMIS = SRCDAT( C,1 )
                AVD_EMIS = SRCDAT( C,2 )
                CEFF = SRCDAT( C,4 ) * 100.
                REFF = SRCDAT( C,5 ) * 100.
                RPEN = SRCDAT( C,6 ) * 100.

C.................  Write out header
                CALL WRITE_ORL_HEADER( RDEV, IOS )
                IF( IOS .GT. 0 ) CYCLE

C.................  Write out entries for this record
C.................  Write ORL nonroad format
                IF( LNONRD ) THEN
                   WRITE( RDEV,93200 ) FIP, SCC, TRIM(CAS), ANN_EMIS, 
     &                  AVD_EMIS, CEFF, REFF, RPEN, SRCTYPE, TRIM(CEXT)
C.................  Write ORL nonpoint format
                ELSE
                   WRITE( RDEV,93210 ) FIP, SCC, TRIM(SIC), TRIM(MACT), 
     &                  TRIM(SRCTYPE), TRIM(NAICS), TRIM(CAS), ANN_EMIS,
     &                  AVD_EMIS, CEFF, REFF, RPEN, TRIM(CEXT)
                END IF

                LCOID = COID
                LYEAR = YEAR

            END DO  ! loop through sources

C............................................................................
C.........  For mobile sources...
        CASE ( 'MOBILE' )

            IF ( FIRSTIME ) THEN

C.................  Write initial header
                WRITE( RDEV, 93000 ) 
     &              '#ORL',
     &              '#TYPE     Onroad Inventory',
     &              '#DESC        (Output from SMOKE)'
                FIRSTIME = .FALSE.

            END IF

C.............  Write mobile-source characteristics to output file
            DO C = 1, NREC

                S = SRCID( C )

C.................  Store others in temporary variables
                COID = STR2INT( CIFIP( S ) ) / 100000
                FIP  = STR2INT( CIFIP( S ) ) - COID * 100000
                SIC  = CISIC( S )
                YEAR = INVYR( S )

                SCC  = CSCC ( S )
 
                IF( ASSOCIATED( CSRCTYP ) ) THEN
                    SRCTYPE = CSRCTYP( S )
                ELSE
                    SRCTYPE = '-9'
                END IF

                IF( ASSOCIATED( CEXTORL ) ) THEN
                    CEXT = ADJUSTL( CEXTORL( S ) )
                ELSE
                    CEXT = ''
                END IF

C.................  Account for missing or default codes
                IF ( LEN( TRIM( SRCTYPE ) ) .EQ. 0 ) SRCTYPE = '-9'
                IF ( LEN( TRIM( CEXT ) ) .EQ. 0 ) CEXT = ''

C.................  Retrieve pollutant code from Inventory Table
C.................  Ensure that only CAS codes where the factor is 1 are used.
C.................  If there are multiple CAS codes where the factor is 1, use the first one
C.................     (accomplished since original INVTABLE list order is maintained already)
                CAS = ' '
                DO I = 1, NINVTBL

                    IF ( DATNAM .EQ. ITNAMA( I ) .AND. ITFACA( I ) .EQ. 1. ) THEN
                        CAS = ITCASA( I )
                        EXIT              ! Exit loop
                    END IF

                END DO           

                IF ( CAS .EQ. ' ' ) THEN
                    MESG = 'ERROR: CAS number for pollutant "'// 
     &                     TRIM( DATNAM ) // '" and factor = 1 in '//
     &                     'INVTABLE is no found.'
                    CALL M3MSG2( MESG )
                    STATUS = 2 
                    CYCLE
                END IF

C.................  Get emissions and emissions-dependent values.
                ANN_EMIS = SRCDAT( C,1 )
                AVD_EMIS = SRCDAT( C,2 )
                CEFF = SRCDAT( C,4 ) * 100.
                REFF = SRCDAT( C,5 ) * 100.
                RPEN = SRCDAT( C,6 ) * 100.

C.................  Write out header
                CALL WRITE_ORL_HEADER( RDEV, IOS )
                IF( IOS .GT. 0 ) CYCLE

C.................  Write out entries for this record
C.................  Write ORL onroad format
                WRITE( RDEV,93300 ) FIP, SCC, TRIM(CAS), ANN_EMIS, 
     &                 AVD_EMIS, SRCTYPE, TRIM(CEXT)

                LCOID = COID
                LYEAR = YEAR

            END DO  ! loop through sources

C............................................................................
C.........  For point sources...
        CASE ( 'POINT' )

            IF ( FIRSTIME ) THEN

C.................  Write initial header
                WRITE( RDEV, 93000 ) 
     &               '#ORL',
     &               '#TYPE     Point Source Inventory',
     &               '#DESC        (Output from SMOKE)'
                FIRSTIME = .FALSE.

            END IF

C.............  Compute conversion constants
            M2FT  = 1./FT2M

C.............  Write point-source characteristics to output file
            DO C = 1, NREC

                S = SRCID( C )
                CALL PARSCSRC( CSOURC( S ), MXCHRS, SC_BEGP, SC_ENDP,
     &                         ORLCOLS, NCHAR, CHARS )

C.................  Truncate character string variables
                PLANTID  = CHARS( 2 ) 
                POINTID  = CHARS( 3 )
                STACKID  = CHARS( 4 )
                SEGMENT  = CHARS( 5 )
                SCC      = CHARS( 6 )

C.................  Store others in temporary variables
                COID   = STR2INT( CIFIP( S ) ) / 100000
                FIP    = STR2INT( CIFIP( S ) ) - COID * 100000
                SIC    = CISIC( S )(SICLEN3-3:SICLEN3)
                YEAR   = INVYR( S )
                CORS   = CORIS( S )
                CBLR   = CBLRID( S )
                CNEI   = ADJUSTL( CNEIUID( S ) )
                PLNTDESC = CPDESC( S )

                IF( ASSOCIATED( CMACT ) ) THEN
                    MACT = CMACT( S )
                ELSE
                    MACT = '-9'
                END IF

                IF( ALLOCATED( CERPTYP ) ) THEN
                    ERPTYPE = CERPTYP( S )
                ELSE
                    ERPTYPE = '-9'
                END IF

                IF( ASSOCIATED( CSRCTYP ) ) THEN
                    SRCTYPE = CSRCTYP( S )
                ELSE
                    SRCTYPE = '-9'
                END IF

                IF( ASSOCIATED( CNAICS ) ) THEN
                    NAICS = CNAICS( S )
                ELSE
                    NAICS = '-9'
                END IF

                IF( ASSOCIATED( CEXTORL ) ) THEN
                    CEXT = ADJUSTL( CEXTORL( S ) )
                ELSE
                    CEXT = '-9'
                END IF

                XLOC     = XLOCA( S )                
                YLOC     = YLOCA( S )                

C.................  Account for missing or default codes
                IF ( LEN( TRIM( PLANTID ) ) .EQ. 0 ) PLANTID = '0'
                IF ( LEN( TRIM( POINTID ) ) .EQ. 0 ) POINTID = '0'
                IF ( LEN( TRIM( STACKID ) ) .EQ. 0 ) STACKID = '0'
                IF ( LEN( TRIM( SEGMENT ) ) .EQ. 0 ) SEGMENT = '0'
                IF ( SCC(1:2) .EQ. '00' ) SCC = SCC(3:SCCLEN3)
                IF ( MACT .EQ. '000000' ) MACT = '-9'
                IF ( MACT(1:2) .EQ. '00' ) MACT = MACT(3:MACLEN3)
                IF ( LEN( TRIM( MACT ) ) .EQ. 0 ) MACT = '-9'
                IF ( LEN( TRIM( SRCTYPE ) ) .EQ. 0 ) SRCTYPE = '-9'
                IF ( LEN( TRIM( ERPTYPE ) ) .EQ. 0 ) ERPTYPE = '-9'
                IF ( NAICS .EQ. '000000' ) NAICS = '-9'
                IF ( NAICS .EQ. '0000-9' ) NAICS = '-9'
                IF ( LEN( TRIM( NAICS ) ) .EQ. 0 ) NAICS = '-9'
                IF ( LEN( TRIM( CEXT ) ) .EQ. 0 ) CEXT = ''

C.................  Convert units of stack parameters
                STKHGT  = STKHT( S ) * M2FT
                STKDIAM = STKDM( S ) * M2FT
                STKTEMP = ( STKTK( S ) - CTOK ) * CTOF + 32.
                STKVEL  = STKVE( S ) * M2FT
                STKFLOW = STKVEL * 0.25 * PI * STKDIAM * STKDIAM

C.................  Retrieve pollutant code from Inventory Table
                I = INDEX1( DATNAM, NINVTBL, ITNAMA )
                CAS = ITCASA( I )

C.................  Get emissions and emissions-dependent values.
                ANN_EMIS = SRCDAT( C,1 )
                AVD_EMIS = SRCDAT( C,2 )
                CEFF = SRCDAT( C,3 ) * 100.
                REFF = SRCDAT( C,4 ) * 100.

C.................  Set default or placeholder variables
                CTYPE = "L"
                UTMZ = -9
                CPRI = 0
                CSEC = 0

C.................  Write out header
                CALL WRITE_ORL_HEADER( RDEV, IOS )
                IF( IOS .GT. 0 ) CYCLE

C.................  Write out data
                WRITE( RDEV, 93600 ) FIP, TRIM(PLANTID), TRIM(POINTID), 
     &                 TRIM(STACKID), TRIM(SEGMENT), TRIM(PLNTDESC),
     &                 TRIM(SCC), TRIM(ERPTYPE), TRIM(SRCTYPE), STKHGT,
     &                 STKDIAM, STKTEMP, STKFLOW, STKVEL, SIC, 
     &                 TRIM(MACT), TRIM(NAICS), CTYPE, XLOC, YLOC, UTMZ, 
     &                 TRIM(CAS), ANN_EMIS, AVD_EMIS, CEFF, REFF, CPRI, 
     &                 CSEC, TRIM(CNEI), TRIM(CORS), TRIM(CBLR),
     &                 TRIM(CEXT)     ! Extended ORL variables

                LCOID = COID
                LYEAR = YEAR

            END DO  ! loop through sources

        END SELECT  ! for source category

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93200   FORMAT( I5.5, ',', A, ',', A, ',', E15.8, ',', E15.8, 
     &          3( ',', F6.2 ),',', A, A )   ! nonroad

93210   FORMAT( I5.5, ',' A, ',', A, ',' A, ',', A, ',', A, ',',
     &          A, ',', E15.8, ',', E15.8, 3( ',', F6.2 ), A )   ! nonpoint

93300   FORMAT( I5.5, ',', A, ',', A, ',', E15.8, ',', E15.8, 
     &          ',', A, A )   ! onroad

93600   FORMAT( I5.5, 8( ',"',A, '"'), 4( ',', F10.2), ',', F10.4,
     &          ',', A4, 2( ',"',A, '"'),',', A1, 
     &          2( ',', F10.5), ',', I3, ',"', A, '"', 2( ',', E13.6 ),
     &          2( ',', F6.2 ), ',', I2, ',', I2, 3(',"',A,'"'), A )  ! point

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram is for writing the COUNTRY and YEAR
C               parts of the ORL header fields, when the COUNTRY or YEAR are
C               inconsistent with the previous country or year.
            SUBROUTINE WRITE_ORL_HEADER( FDEV, LOCSTAT )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: FDEV    ! file unit no.
            INTEGER, INTENT(OUT) :: LOCSTAT ! exit status

C.............  Subprogram local variables
            INTEGER     K
            CHARACTER( FIPLEN3 ) CTRY

C-------------------------------------------------------------------------

            LOCSTAT = 0

            IF( COID .NE. LCOID .OR. YEAR .NE. LYEAR ) THEN

                WRITE( CYEAR, '(I4)' ) YEAR
                WRITE( CTRY,  '(I6)' ) COID
                CALL PADZERO( CTRY ) 

                K = FINDC( CTRY, NCOUNTRY, CTRYCOD )

                IF( K .GT. 0 ) THEN

                    WRITE( FDEV, 93000 ) 
     &                 '#COUNTRY  ' // CTRYNAM( K ),
     &                 '#YEAR     ' // CYEAR

                ELSE
                    STATUS = 1
                    LOCSTAT = 1
                    WRITE( MESG,94010 ) 'ERROR: Invalid country code', K,
     &                     'found at source', S
                    CALL M3MESG( MESG )

                END IF

            END IF

C---------------------  FORMAT  STATEMENTS   -----------------------------

C...............   Formatted file I/O formats............ 93xxx

93000       FORMAT( A )

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE WRITE_ORL_HEADER

C-------------------------------------------------------------------------
C-------------------------------------------------------------------------

        END SUBROUTINE WRORLOUT

