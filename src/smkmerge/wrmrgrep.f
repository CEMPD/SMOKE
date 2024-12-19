
        SUBROUTINE WRMRGREP( JDATE, JTIME, NIDX )

C***********************************************************************
C  subroutine WRMRGREP body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to write state and county totals
C      with headers and formatted nicely by column width.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 8/99 by M. Houyoux
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: SDATE, STIME, EDATE, ETIME, TSTEP,
     &                      AFLAG, BFLAG, MFLAG, PFLAG, XFLAG, SFLAG, 
     &                      ANMSPC, BNMSPC, MNMSPC, PNMSPC, NMSPC, 
     &                      AEMNAM, BEMNAM, MEMNAM, PEMNAM, 
     &                      ANIPOL, MNIPPA, PNIPOL, NIPPA, 
     &                      AEINAM, MEANAM, PEINAM, 
     &                      AEBCNY, BEBCNY, MEBCNY, PEBCNY, TEBCNY, 
     &                      AEBSTA, BEBSTA, MEBSTA, PEBSTA, TEBSTA, 
     &                      ARDEV,  BRDEV,  MRDEV,  PRDEV,  TRDEV, 
     &                      AUFLAG,         MUFLAG, PUFLAG, TUFLAG, 
     &                      ARFLAG,         MRFLAG, PRFLAG, TRFLAG, 
     &                      AECCNY,         MECCNY, PECCNY, TECCNY, 
     &                      AECSTA,         MECSTA, PECSTA, TECSTA, 
     &                      VGRPCNT, LREPSTA, LREPCTL, LREPCNY, 
     &                      LAVEDAY, EMNAM, EANAM, TOTUNIT, 
     &                      SIINDEX, SPINDEX

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTY, NSTATE, STATNAM, CNTYNAM, CNTYCOD

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)    CRLF
        INTEGER         ENVINT  
        INTEGER         INDEX1  
        CHARACTER(14)   MMDDYY
        INTEGER         WKDAY
        REAL            YR2DAY

        EXTERNAL    CRLF, ENVINT, INDEX1, MMDDYY, WKDAY, YR2DAY

C...........   Subroutine arguments
        INTEGER, INTENT (IN) :: JDATE  ! julian date  (YYYYDDD)
        INTEGER, INTENT (IN) :: JTIME  ! time (HHMMSS)
        INTEGER, INTENT (IN) :: NIDX   ! group index

C...........   Local allocatable arrays
        LOGICAL, ALLOCATABLE :: LUPDATE( : ) ! true: units/names not yet updated

        CHARACTER(30), ALLOCATABLE, SAVE ::  NAMES( : )
        CHARACTER(30), ALLOCATABLE, SAVE :: ANAMES( : )
        CHARACTER(30), ALLOCATABLE, SAVE :: BNAMES( : )
        CHARACTER(30), ALLOCATABLE, SAVE :: MNAMES( : )
        CHARACTER(30), ALLOCATABLE, SAVE :: PNAMES( : )

        CHARACTER(30), ALLOCATABLE, SAVE ::  UNITS( : )
        CHARACTER(30), ALLOCATABLE, SAVE :: AUNITS( : )
        CHARACTER(30), ALLOCATABLE, SAVE :: BUNITS( : )
        CHARACTER(30), ALLOCATABLE, SAVE :: MUNITS( : )
        CHARACTER(30), ALLOCATABLE, SAVE :: PUNITS( : )

C...........   Local group counts
        INTEGER, SAVE :: ACNT = 0       ! area source output vars count
        INTEGER, SAVE :: BCNT = 0       ! biogenic source output vars count
        INTEGER, SAVE :: MCNT = 0       ! mobile source output vars count
        INTEGER, SAVE :: PCNT = 0       ! point source output vars count
        INTEGER, SAVE :: TCNT = 0       ! total source output vars count

C...........   Other local variables

        INTEGER          C, F, I, J, K, L, L2, N, V  ! counters and indices
        INTEGER          IOS            ! i/o status
        INTEGER       :: KA = 0         ! tmp search indices
        INTEGER       :: KB = 0         ! tmp search indices
        INTEGER       :: KM = 0         ! tmp search indices
        INTEGER       :: KP = 0         ! tmp search indices
        INTEGER, SAVE :: MAXCYWID       ! max width of county names
        INTEGER, SAVE :: MAXSTWID       ! max width of state names
        INTEGER, SAVE :: NC             ! tmp no. counties in domain
        INTEGER, SAVE :: NS             ! tmp no. states in domain
        INTEGER, SAVE :: PGRP   = 0     ! group from previous call
        INTEGER, SAVE :: PDATE          ! date for computing PTIME
        INTEGER, SAVE :: PTIME          ! start time of period
        INTEGER, SAVE :: REPTIME        ! report time
        INTEGER, SAVE :: NVPGP          ! no. variables per group

        REAL   , SAVE :: UNITFAC        ! units conversion factor for reports
 
        LOGICAL, SAVE :: FIRSTIME= .TRUE. ! true: first time routine called
        LOGICAL, SAVE :: WARNON  = .TRUE. ! true: report warning

        CHARACTER(3000)    DATFMT     ! format for data
        CHARACTER(3000)    HDRFMT     ! format for header
        CHARACTER(3000)    HEADER     ! header for output files
        CHARACTER(3000)    LINFLD     ! line of dashes
        CHARACTER(300)     MESG       ! message buffer
        CHARACTER(30) SBUF       ! tmp pol or species name
        CHARACTER(30) CBUF       ! tmp units field

        CHARACTER(16) :: PROGNAME = 'WRMRGREP' ! program name

C***********************************************************************
C   begin body of subroutine WRMRGREP

C.........    If first time the routine is called
        IF( FIRSTIME ) THEN

            NC = NCOUNTY
            NS = NSTATE

C.............  Get time for output of reports 
            MESG = 'Report time for country, state, and county totals'
            REPTIME = ENVINT( 'SMK_REPORT_TIME', MESG, 230000, IOS )

C............. NOTE that this would be a good place for an environment variable
C              that permits users to select the write format of the data

C.............  Get the maximum width for the state names
            MAXSTWID = 0 
            DO I = 1, NS
                L = LEN_TRIM( STATNAM( I ) )
                MAXSTWID = MAX( MAXSTWID, L )
            END DO

C.............  Get the maximum width for the county names
            MAXCYWID = 0 
            DO I = 1, NC
                L = LEN_TRIM( CNTYNAM( I ) )
                MAXCYWID = MAX( MAXCYWID, L )
            END DO

            PDATE = SDATE
            PTIME = STIME

C.............  Allocate memory for the units and names based on the
C               largest group size
            I = MAXVAL( VGRPCNT )
            ALLOCATE( NAMES( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NAMES', PROGNAME )            
            ALLOCATE( UNITS( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'UNITS', PROGNAME )            

            ALLOCATE( ANAMES( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ANAMES', PROGNAME )
            ALLOCATE( AUNITS( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AUNITS', PROGNAME )

            ALLOCATE( BNAMES( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BNAMES', PROGNAME )
            ALLOCATE( BUNITS( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BUNITS', PROGNAME )

            ALLOCATE( MNAMES( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MNAMES', PROGNAME )
            ALLOCATE( MUNITS( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MUNITS', PROGNAME )
 
            ALLOCATE( PNAMES( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PNAMES', PROGNAME )
            ALLOCATE( PUNITS( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PUNITS', PROGNAME )

C.............  Allocate memory for the logical update names/units flag
            I = MAX( NIPPA, NMSPC )
            ALLOCATE( LUPDATE( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LUPDATE', PROGNAME )

            FIRSTIME = .FALSE.

        END IF

C.........  If current group is a new group...
        IF( NIDX .NE. PGRP ) THEN

C.............  Initialize units and names
            UNITS  = ' '  ! array
            AUNITS = ' '  ! array
            BUNITS = ' '  ! array
            MUNITS = ' '  ! array
            PUNITS = ' '  ! array
            ANAMES = ' '  ! array
            BNAMES = ' '  ! array
            MNAMES = ' '  ! array
            PNAMES = ' '  ! array

C.............  Initialize logical update states
            LUPDATE = .TRUE.   ! array

C.............  Set no. variables per group. Recall that this number can include
C               multiple species-pollutant combos (e.g., multiple processes for
C               mobile sources).
            NVPGP = VGRPCNT( NIDX )

C.............  Create units labels from variable units for current group
            ACNT = 0
            BCNT = 0
            MCNT = 0
            PCNT = 0
            TCNT = 0
            DO V = 1, NVPGP

C.................  Get indices
                I = SIINDEX( V,NIDX )  ! index to EANAM
                J = I
                SBUF = EANAM( J )

                IF( SFLAG ) THEN
                    J = SPINDEX( V,NIDX )  ! index to EMNAM 
                    SBUF = EMNAM( J )
                END IF

C.................  If this species or pollutant has not already been 
C                   encountered, update the source-category-specific arrays 
C                   so that they will be in the same order as the global lists
                IF( LUPDATE( J ) ) THEN

                    LUPDATE( J ) = .FALSE.

C.....................  Set names and units for totals output
                    IF( SFLAG ) THEN
                        L = LEN_TRIM( TOTUNIT( J ) )
                        CBUF = '[' // TOTUNIT( J )( 1:L ) // ']'
                    ELSE
                        L = LEN_TRIM( TOTUNIT( I ) )
                        CBUF = '[' // TOTUNIT( I )( 1:L ) // ']'
                    END IF

                    TCNT = TCNT + 1
                    NAMES( TCNT ) = SBUF
                    UNITS( TCNT ) = CBUF

C.....................  Get indices for source categories, depending on
C                       speciation or not
                    IF( SFLAG ) THEN
                        IF( AFLAG ) KA = INDEX1( SBUF, ANMSPC, AEMNAM )
                        IF( BFLAG ) KB = INDEX1( SBUF, BNMSPC, BEMNAM )
                        IF( MFLAG ) KM = INDEX1( SBUF, MNMSPC, MEMNAM )
                        IF( PFLAG ) KP = INDEX1( SBUF, PNMSPC, PEMNAM )
                    ELSE
                        IF( AFLAG ) KA = INDEX1( SBUF, ANIPOL, AEINAM )
                        KB = 0
                        IF( MFLAG ) KM = INDEX1( SBUF, MNIPPA, MEANAM )
                        IF( PFLAG ) KP = INDEX1( SBUF, PNIPOL, PEINAM )
                    END IF

C.....................  Set names and units for area sources
                    IF( KA .GT. 0 ) THEN
                        ACNT = ACNT + 1
                        ANAMES( ACNT ) = SBUF
                        AUNITS( ACNT ) = CBUF
                    END IF

C.....................  Set names and units for biogenics
                    IF( KB .GT. 0 ) THEN
                        BCNT = BCNT + 1
                        BNAMES( BCNT ) = SBUF
                        BUNITS( BCNT ) = CBUF
                    END IF

C.....................  Set names and units for mobile sources
                    IF( KM .GT. 0 ) THEN
                        MCNT = MCNT + 1
                        MNAMES( MCNT ) = SBUF
                        MUNITS( MCNT ) = CBUF
                    END IF

C.....................  Set names and units for point sources
                    IF( KP .GT. 0 ) THEN
                        PCNT = PCNT + 1
                        PNAMES( PCNT ) = SBUF
                        PUNITS( PCNT ) = CBUF
                    END IF

                END IF  ! End if global units already defined or not

            END DO      ! End loop on group

            PGRP = NIDX

        END IF   ! End of processing for current group


C.........  Do not report if this time is not appropriate
        IF( JTIME .NE. REPTIME .AND. 
     &    ( JDATE .NE. EDATE .OR. JTIME .NE. ETIME ) ) RETURN

C.............  If required, create and write state totals, either controlled
C               or uncontrolled, depending on sector and which are controlled.
        IF ( LREPSTA ) THEN  

C.............  Controlled area sources
            IF( ( AUFLAG .OR. ARFLAG ) .AND. LREPCTL ) THEN
                CALL CREATE_HEADER( 'Controlled area' )
                CALL CREATE_STATE( NC, NS, ACNT, AECCNY, AECSTA )
                CALL WRITE_STA( ARDEV, NS, ACNT, ANAMES, AUNITS, AECSTA)

C.................  Update state totals
                IF( XFLAG ) THEN
                    CALL TOT_UPDATE( NS, ACNT, ANAMES, AECSTA, TECSTA)
                END IF

C.............  Uncontrolled area sources
            ELSE IF ( AFLAG ) THEN

                CALL CREATE_HEADER( 'Area' )
                CALL CREATE_STATE( NC, NS, ACNT, AEBCNY, AEBSTA )
                CALL WRITE_STA( ARDEV, NS, ACNT, ANAMES, AUNITS, AEBSTA)

C.....................  Update state totals
                IF( XFLAG ) THEN
                    CALL TOT_UPDATE( NS, ACNT, ANAMES, AEBSTA, TEBSTA)
                END IF

            END IF

C.............  Biogenic sources (speciated by definition)
            IF( BFLAG ) THEN

                CALL CREATE_HEADER( 'Biogenic' )
                CALL CREATE_STATE( NC, NS, BCNT, BEBCNY, BEBSTA )
                CALL WRITE_STA( BRDEV, NS, BCNT, BNAMES, BUNITS, BEBSTA)

C.................  Update state totals
                IF( XFLAG ) THEN
                    CALL TOT_UPDATE( NS, BCNT, BNAMES, BEBSTA, TEBSTA )
                END IF

            END IF

C.............  Controlled mobile sources
            IF( ( MUFLAG .OR. MRFLAG ) .AND. LREPCTL ) THEN
                CALL CREATE_HEADER( 'Controlled mobile' )
                CALL CREATE_STATE( NC, NS, MCNT, MECCNY, MECSTA )
                CALL WRITE_STA( MRDEV, NS, MCNT, MNAMES, MUNITS, MECSTA)

C.................  Update state totals
                IF( XFLAG ) THEN
                    CALL TOT_UPDATE( NS, MCNT, MNAMES, MECSTA, TECSTA)
                END IF

C.............  Uncontrolled mobile sources
            ELSE IF( MFLAG ) THEN

                CALL CREATE_HEADER( 'Mobile' )
                CALL CREATE_STATE( NC, NS, MCNT, MEBCNY, MEBSTA )
                CALL WRITE_STA( MRDEV, NS, MCNT, MNAMES, MUNITS, MEBSTA)

C.................  Update state totals
                IF( XFLAG ) THEN
                    CALL TOT_UPDATE( NS, MCNT, MNAMES, MEBSTA, TEBSTA)
                END IF

            END IF

C.............  Controlled point sources
            IF( ( PUFLAG .OR. PRFLAG ) .AND. LREPCTL ) THEN
                CALL CREATE_HEADER( 'Controlled point' )
                CALL CREATE_STATE( NC, NS, PCNT, PECCNY, PECSTA )
                CALL WRITE_STA( PRDEV, NS, PCNT, PNAMES, PUNITS, PECSTA)

C.................  Update state totals
                IF( XFLAG ) THEN
                    CALL TOT_UPDATE( NS, PCNT, PNAMES, PECSTA, TECSTA)
                END IF

C.............  Uncontrolled point sources
            ELSE IF( PFLAG ) THEN

                CALL CREATE_HEADER( 'Point' )
                CALL CREATE_STATE( NC, NS, PCNT, PEBCNY, PEBSTA )
                CALL WRITE_STA( PRDEV, NS, PCNT, PNAMES, PUNITS, PEBSTA)

C.................  Update state totals
                IF( XFLAG ) THEN
                    CALL TOT_UPDATE( NS, PCNT, PNAMES, PEBSTA, TEBSTA )
                END IF
            END IF

C.............  Controlled total sources
            IF( ( TUFLAG .OR. TRFLAG ) .AND. LREPCTL ) THEN
                CALL CREATE_HEADER( 'Controlled total' )
                CALL CREATE_STATE( NC, NS, TCNT, TECCNY, TECSTA )
                CALL WRITE_STA( TRDEV, NS, TCNT, NAMES, UNITS, TECSTA)

C.............  Uncontrolled total sources
            ELSE IF( XFLAG ) THEN

                CALL CREATE_HEADER( 'Total' )
                CALL WRITE_STA( TRDEV, NS, TCNT, NAMES, UNITS, TEBSTA )

            END IF

        END IF

C.........  If required, write county totals
        IF( LREPCNY ) THEN

C.............  Give warning if not already given and if reporting of controlled
C               emissions was selected
            IF ( LREPCTL .AND. WARNON ) THEN
                WARNON = .FALSE.

                MESG = 'WARNING: County reports do not support '//
     &                 'including controlled emissions, but this '//
     &                 'combination was requested.'
                CALL M3MESG( MESG )

            END IF

C.............  Area sources
            IF( AFLAG ) THEN

                CALL CREATE_HEADER( 'Area' )
                CALL WRITE_CNY( ARDEV, NC, ACNT, ANAMES, AUNITS, AEBCNY)

C.................  Update state totals
                IF( XFLAG ) THEN
                    CALL TOT_UPDATE( NC, ACNT, ANAMES, AEBCNY, TEBCNY)
                END IF
            END IF

C.............  Biogenc sources (speciated by definition)
            IF( BFLAG ) THEN

                CALL CREATE_HEADER( 'Biogenic' )
                CALL WRITE_CNY( BRDEV, NC, BCNT, BNAMES, BUNITS, BEBCNY)

C.................  Update state totals
                IF( XFLAG ) THEN
                    CALL TOT_UPDATE( NC, BCNT, BNAMES, BEBCNY, TEBCNY )
                END IF

            END IF

C.............  Mobile sources
            IF( MFLAG ) THEN

                CALL CREATE_HEADER( 'Mobile' )
                CALL WRITE_CNY( MRDEV, NC, MCNT, MNAMES, MUNITS, MEBCNY)

C.................  Update state totals
                IF( XFLAG ) THEN
                    CALL TOT_UPDATE( NC, MCNT, MNAMES, MEBCNY, TEBCNY)
                END IF

            END IF

C.............  Point sources
            IF( PFLAG ) THEN

                CALL CREATE_HEADER( 'Point' )
                CALL WRITE_CNY( PRDEV, NC, PCNT, PNAMES, PUNITS, PEBCNY)

C.................  Update state totals
                IF( XFLAG ) THEN
                    CALL TOT_UPDATE( NC, PCNT, PNAMES, PEBCNY, TEBCNY)
                END IF

            END IF

C.............  Combined sources
            IF( XFLAG ) THEN

                CALL CREATE_HEADER( 'Total' )
                CALL WRITE_CNY( TRDEV, NC, TCNT, NAMES, UNITS, TEBCNY )

            END IF

        END IF

C.........  Intialize summed emissions to zero
        CALL INITSTCY

C........  After output, increment date and time one step to set the start
C          of the next period
        PDATE = JDATE
        PTIME = JTIME
        CALL NEXTIME( PDATE, PTIME, TSTEP )  ! advance one

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C*****************  INTERNAL SUBPROGRAMS   *****************************

        CONTAINS

            SUBROUTINE CREATE_STATE( NC, NS, NDIM, CY_EMIS, ST_EMIS )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: NC
            INTEGER     , INTENT (IN) :: NS
            INTEGER     , INTENT (IN) :: NDIM
            REAL        , INTENT (IN) :: CY_EMIS( NC, NDIM )
            REAL        , INTENT(OUT) :: ST_EMIS( NS, NDIM )

C.............  Local variables
            INTEGER  I, J, N
            
            CHARACTER(FIPLEN3) PSTA, STA

C..............................................................................

            PSTA = ' '
            N = 0
            DO I = 1, NC


                STA = CNTYCOD( I )
                STA( FIPLEN3-2:FIPLEN3 ) = '000'
                IF( STA .NE. PSTA ) THEN
                    N = N + 1
                    PSTA = STA
                END IF

                DO J = 1, NDIM
                    ST_EMIS( N,J ) = ST_EMIS( N,J ) + CY_EMIS( I,J )
                END DO

            END DO

            RETURN

            END SUBROUTINE CREATE_STATE

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

            SUBROUTINE CREATE_HEADER( CATNAME )

C.............  MODULES for public variables
C.............  This module contains the global variables for the 3-d grid
            USE MODGRID, ONLY: GRDNM

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: CATNAME

C.............  Local variables
            INTEGER   K1, K2, L, LD1, LD2

            CHARACTER(10) TYPENAM
            CHARACTER(15) DATANAM

C..............................................................................

            HEADER = '# ' // CATNAME // ' source'
            L  = LEN_TRIM( HEADER )

            DATANAM = ' average'
            IF( LAVEDAY .AND. ( AFLAG .OR. PFLAG ) ) 
     &          DATANAM = ' average day'
            LD1 = LEN_TRIM( DATANAM )
     
            TYPENAM = ' inventory'
            IF( SFLAG ) TYPENAM = ' speciated'

            HEADER = HEADER( 1:L ) // DATANAM( 1:LD1 ) //
     &               TYPENAM // ' emissions'
            L = LEN_TRIM( HEADER )

            IF( JDATE .NE. 0 ) THEN

                K1  = WKDAY( PDATE )
                K2  = WKDAY( JDATE )
                LD1 = LEN_TRIM( DAYS( K1 ) )
                LD2 = LEN_TRIM( DAYS( K2 ) )

                WRITE( HEADER,94010 ) HEADER( 1:L ) // ' from ' //
     &             CRLF() // '#' // BLANK5(1:4) // 
     &             DAYS( K1 )( 1:LD1 ) // ' ' // MMDDYY( PDATE ) //
     &             ' at', PTIME, 'to' // CRLF() // '#' // BLANK5(1:4) // 
     &             DAYS( K2 )( 1:LD2 ) // ' '// MMDDYY( JDATE ) //
     &             ' at', JTIME
                L = LEN_TRIM( HEADER )

            END IF

            HEADER = HEADER( 1:L ) // ' within grid ' // GRDNM

            RETURN

C.......................... FORMAT STATEMENTS ................................

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE CREATE_HEADER

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

            SUBROUTINE CREATE_FORMATS( NDIM, MAXCOL1, INNAMS, INUNIT,
     &                                 WIDTHS, OUTNAMS, OUTUNIT  )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN)     :: NDIM
            INTEGER     , INTENT (IN)     :: MAXCOL1
            CHARACTER(*), INTENT (IN)     :: INNAMS ( NDIM )
            CHARACTER(*), INTENT (IN)     :: INUNIT ( NDIM )
            INTEGER     , INTENT (IN OUT) :: WIDTHS ( 0:NDIM )
            CHARACTER(*), INTENT (IN OUT) :: OUTNAMS( NDIM )
            CHARACTER(*), INTENT (IN OUT) :: OUTUNIT( NDIM )

C.............  Local parameters
            INTEGER, PARAMETER :: EFMTWID = 10  ! minimum width of fields
            INTEGER, PARAMETER :: EFMTDEC = 4   ! number of decimal places

C.............  Local variables
            INTEGER       I1, I2, J, L, L1, L2
            CHARACTER(30)  :: SPACE = ' '
            CHARACTER(3000) :: TMPFMT

C.............................................................................

C.............  Adjust maximum width of numbers in case width of variable
C               names or units is greater than width of numbers
C.............  Also, move the position of the names and units so that output 
C               strings will look right-justified in the file.
            DO J = 1, NDIM
                L1 = LEN_TRIM( INNAMS( J ) )
                WIDTHS ( J ) = MAX( WIDTHS( J ), L1, EFMTWID )
                L2 = LEN_TRIM( INUNIT( J ) )
                WIDTHS ( J ) = MAX( WIDTHS( J ), L2 )

                I1 = WIDTHS( J ) - L1
                IF( I1 .GT. 0 ) THEN
                    WRITE( OUTNAMS( J ), '(A,A)' ) SPACE( 1:I1 ), 
     &                                             INNAMS( J )( 1:L1 )
                ELSE
                    WRITE( OUTNAMS( J ), '(A)' ) INNAMS( J )( 1:L1 )
                END IF

                I2 = WIDTHS( J ) - L2
                IF( I2 .GT. 0 ) THEN
                    WRITE( OUTUNIT( J ), '(A,A)' ) SPACE( 1:I2 ), 
     &                                             INUNIT( J )( 1:L2 )
                ELSE
                    WRITE( OUTUNIT( J ), '(A)' ) INUNIT( J )( 1:L2 )
                END IF

            END DO

C.............  Create format statement for output of header
            WIDTHS( 0 ) = MAXCOL1
            WRITE( HDRFMT, '( "(A",A)' ) ',";"'
            DO J = 1, NDIM
                TMPFMT = HDRFMT
                L = LEN_TRIM( TMPFMT ) 
                WRITE( HDRFMT, '(A, ",1X,A",I2.2,A)' ) 
     &                 TMPFMT(1:L), WIDTHS( J ), ',";"'
            END DO
            TMPFMT = HDRFMT
            L = LEN_TRIM( TMPFMT )
            WRITE( HDRFMT, '(A)' ) TMPFMT( 1:L ) // ')'

C.............  Create format statement for output of emissions
            WRITE( DATFMT, '( "(A",I2.2,A)' ) WIDTHS( 0 ), ',";"'
            DO J = 1, NDIM
                TMPFMT = DATFMT
                L = LEN_TRIM( TMPFMT ) 
                WRITE( DATFMT, '(A, ",1X,E",I2.2,".",I1,A)' ) 
     &                 TMPFMT(1:L), WIDTHS( J ), EFMTDEC, ',";"'
            END DO
            TMPFMT = DATFMT
            L = LEN_TRIM( TMPFMT )
            WRITE( DATFMT, '(A)' ) TMPFMT( 1:L ) // ')'

            RETURN

            END SUBROUTINE CREATE_FORMATS

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

            SUBROUTINE WRITE_STA( FDEV, NS, NDIM, VNAMES, INUNIT,
     &                            ST_EMIS )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: FDEV
            INTEGER     , INTENT (IN) :: NS
            INTEGER     , INTENT (IN) :: NDIM
            CHARACTER(*), INTENT (IN) :: VNAMES ( NDIM )
            CHARACTER(*), INTENT (IN) :: INUNIT ( NDIM )
            REAL        , INTENT (IN) :: ST_EMIS( NS, NDIM )

C.............  Arrays allocated by subprogram argument
            INTEGER            MAXWID ( 0:NDIM )
            CHARACTER(30) OUTNAMS( NDIM )
            CHARACTER(30) OUTUNIT( NDIM )

C.............  Local variables
            INTEGER       I, J, L, L2

            REAL          VAL

            CHARACTER(20) :: HDRBUF  = '#'
            CHARACTER(20) :: STLABEL = '# State'
            CHARACTER(30)    BUFFER

C..............................................................................

C.............  Get maximum width of numbers and 
            DO J = 1, NDIM

                VAL = MAXVAL( ST_EMIS( 1:NS, J ) )
                WRITE( BUFFER, '(F30.1)' ) VAL
                BUFFER = ADJUSTL( BUFFER )
                MAXWID( J ) = LEN_TRIM( BUFFER )

            END DO

C.............  Rearrange labels and units to be in order of master list, which
C               is the order that the emission values themselves will be in

C.............  Get column labels and formats
            CALL CREATE_FORMATS( NDIM, 6+MAXSTWID, VNAMES, INUNIT,
     &                           MAXWID, OUTNAMS, OUTUNIT )

C.............  Create line format
c            L2 = SUM( MAXWID ) + NDIM
c            LINFLD = REPEAT( '-', L2 )

C.............  Write header for state totals
            WRITE( FDEV, '(A)' ) '# '
            WRITE( FDEV, '(A)' ) HEADER( 1:LEN_TRIM( HEADER ) )

C.............  Write line
c            WRITE( FDEV, '(A)' ) LINFLD( 1:L2 )

C.............  Write units for columns
            WRITE( FDEV, HDRFMT ) ADJUSTL( HDRBUF), 
     &                            ( OUTUNIT(J), J=1,NDIM )

C.............  Write column labels
            WRITE( FDEV, HDRFMT ) ADJUSTL( STLABEL ),
     &                          ( OUTNAMS( J ), J=1, NDIM )

C.............  Write state total emissions
            DO I = 1, NSTATE

C.................  Build output format depending on data values
                CALL DYNAMIC_FORMATS( NSTATE, NDIM, I, ST_EMIS,
     &                                MAXWID(1), DATFMT )

C.................  Write out state name and converted emissions
                WRITE( FDEV, DATFMT ) STATNAM( I ), 
     &                                ( ST_EMIS( I,J ), J=1, NDIM )
            END DO

            RETURN

            END SUBROUTINE WRITE_STA

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

            SUBROUTINE WRITE_CNY( FDEV, NC, NDIM, VNAMES, INUNIT,
     &                            CY_EMIS )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: FDEV
            INTEGER     , INTENT (IN) :: NC
            INTEGER     , INTENT (IN) :: NDIM
            CHARACTER(*), INTENT (IN) :: VNAMES ( NDIM )
            CHARACTER(*), INTENT (IN) :: INUNIT ( NDIM )
            REAL        , INTENT (IN) :: CY_EMIS( NC, NDIM )

C.............  Arrays allocated by subprogram argument
            INTEGER            MAXWID ( 0:NDIM )
            CHARACTER(30) OUTNAMS( NDIM )
            CHARACTER(30) OUTUNIT( NDIM )

C.............  Local variables
            INTEGER       I, J, L, L1, L2, N

            REAL          VAL

            CHARACTER(FIPLEN3+8) CDATFIP
            CHARACTER(FIPLEN3) PSTA, STA
            CHARACTER(60) :: HDRBUF  = '#'
            CHARACTER(60) :: STLABEL = '# County'
            CHARACTER(30)    BUFFER

C..............................................................................

C.............  Get maximum width of numbers
            DO J = 1, NDIM

                VAL = MAXVAL( CY_EMIS( 1:NC, J ) )
                WRITE( BUFFER, '(F30.1)' ) VAL
                BUFFER = ADJUSTL( BUFFER )
                MAXWID( J ) = LEN_TRIM( BUFFER )

            END DO

C.............  Get column labels and formats
            CALL CREATE_FORMATS( NDIM, 26+MAXSTWID+MAXCYWID, VNAMES,
     &                           INUNIT, MAXWID, OUTNAMS, OUTUNIT )

C.............  Create line format
c            L2 = SUM( MAXWID ) + NDIM
c            LINFLD = REPEAT( '-', L2 )

C.............  Write header for county totals
            WRITE( FDEV, '(A)' ) '# '
            WRITE( FDEV, '(A)' ) HEADER( 1:LEN_TRIM( HEADER ) )

C.............  Write units for columns
            WRITE( FDEV, HDRFMT ) ADJUSTL( HDRBUF), 
     &                            ( OUTUNIT(J), J=1,NDIM )

C.............  Write column labels
            WRITE( FDEV, HDRFMT ) ADJUSTL( STLABEL ),
     &                            ( OUTNAMS( J ), J=1, NDIM )

C.............  Write line
c            WRITE( FDEV, '(A)' ) LINFLD( 1:L2 )

C.............  Write county total emissions
            PSTA = ' '
            N = 0
            DO I = 1, NC

                STA = CNTYCOD( I )
                STA( FIPLEN3-2:FIPLEN3 ) = '000'
                IF( STA .NE. PSTA ) THEN
                    N = N + 1
                    PSTA = STA
                END IF

C.................  Write out county name and converted emissions
                WRITE( CDATFIP, '(I7.7,1X,A)' ) JDATE, CNTYCOD( I )

C.................  Build output format depending on data values
                CALL DYNAMIC_FORMATS( NC, NDIM, I, CY_EMIS,
     &                                MAXWID(1), DATFMT )

                WRITE( FDEV,DATFMT ) CDATFIP // ' '// STATNAM(N) // 
     &                               CNTYNAM(I), 
     &                               ( CY_EMIS( I,J ), J=1, NDIM )

            END DO

            RETURN

            END SUBROUTINE WRITE_CNY

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

            SUBROUTINE TOT_UPDATE( NU, NDIM, VNAMES, INEMIS, TOTEMIS )

C.............  Subprogram arguments
            INTEGER     , INTENT    (IN) :: NU               ! no. summing units
            INTEGER     , INTENT    (IN) :: NDIM           ! no. pols or species
            CHARACTER(*), INTENT    (IN) :: VNAMES ( NDIM )   ! pol or spc names
            REAL        , INTENT    (IN) :: INEMIS ( NU, NDIM ) ! component emis
            REAL        , INTENT(IN OUT) :: TOTEMIS( NU, * )    ! total emis

C.............  Local variables
            INTEGER       I, K, V

C..............................................................................

C.............  Abort function call if only one source category
            IF( .NOT. XFLAG ) RETURN

C.............  Loop through pollutants or species, and search for
C               local names in master lists to get index to totals array           
            DO V = 1, NDIM

                IF( SFLAG ) THEN
                    K = INDEX1( VNAMES( V ), NMSPC, EMNAM )
                ELSE
                    K = INDEX1( VNAMES( V ), NIPPA, EANAM )
                END IF

C.................  If pollutant or species is found, then add component
C                   emissions to the total
                IF( K .GT. 0 ) THEN

                    DO I = 1, NU

                        TOTEMIS( I,K ) = TOTEMIS( I,K ) + INEMIS( I,V )

                    END DO

                END IF

            END DO

            RETURN

            END SUBROUTINE TOT_UPDATE

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

            SUBROUTINE DYNAMIC_FORMATS( N1, N2, STCNT, EMIS, 
     &                                  WIDTH, DATFMT )

            INTEGER     , INTENT (IN)  :: N1
            INTEGER     , INTENT (IN)  :: N2
            INTEGER     , INTENT (IN)  :: STCNT      ! counter for state/county index
            REAL        , INTENT (IN)  :: EMIS( N1, N2 )  ! state/count emission
            INTEGER     , INTENT (IN)  :: WIDTH( N2 )
            CHARACTER(*), INTENT (OUT) :: DATFMT

C.............  Local variables
            INTEGER  I
            CHARACTER(5) :: FMT

C..............................................................................

C.............  Initialize format array
            DATFMT = '(A'

C............. Determine significant figures that we want to report. 
C............. Use a minimum of 5 significant figures, and a minimum
C              of 1 decimal place.  If value is < 0.1, then use exponential
C              with 4 decimal places
            DO I = 1, N2

                IF ( EMIS( STCNT,I ) .EQ. 0. ) THEN
                   WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.0'

                ELSE IF( EMIS( STCNT,I ) .GE. 1000. ) THEN
                   WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.1'

                ELSE IF ( EMIS( STCNT,I ) .GE. 100. ) THEN
                   WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.2'

                ELSE IF ( EMIS( STCNT,I ) .GE. 10. ) THEN
                   WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.3'

                ELSE IF ( EMIS( STCNT,I ) .GE. 0.1 ) THEN
                   WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.4'

                ELSE IF ( EMIS( STCNT,I ) .LT. 0.1 ) THEN 
                   WRITE( FMT, '(A,I2.2,A)' ) 'E',WIDTH(I),'.4'
                END IF

                DATFMT = TRIM( DATFMT ) // ',"; "' // FMT

            END DO

            DATFMT = TRIM( DATFMT ) // ')'

            RETURN

            END SUBROUTINE DYNAMIC_FORMATS

        END SUBROUTINE WRMRGREP

