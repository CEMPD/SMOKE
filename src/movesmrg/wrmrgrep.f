
        SUBROUTINE WRMRGREP( JDATE, JTIME )

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
     &                      NMSPC, MEBCNY, MEBSTA,
     &                      MRDEV, LREPSTA, LREPCNY, 
     &                      EMNAM, TOTUNIT

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

C...........   Local allocatable arrays
        CHARACTER(IOVLEN3), ALLOCATABLE, SAVE :: MNAMES( : )
        CHARACTER(IOULEN3), ALLOCATABLE, SAVE :: MUNITS( : )

C...........   Other local variables

        INTEGER          C, F, I, J, K, L, L2, N, V  ! counters and indices
        INTEGER          IOS            ! i/o status
        INTEGER          PSTATE         ! previous state
        INTEGER          STATE          ! current state
        INTEGER, SAVE :: MAXCYWID       ! max width of county names
        INTEGER, SAVE :: MAXSTWID       ! max width of state names
        INTEGER, SAVE :: NC             ! tmp no. counties in domain
        INTEGER, SAVE :: NS             ! tmp no. states in domain
        INTEGER, SAVE :: PDATE          ! date for computing PTIME
        INTEGER, SAVE :: PTIME          ! start time of period
        INTEGER, SAVE :: REPTIME        ! report time

        REAL   , SAVE :: UNITFAC        ! units conversion factor for reports
 
        LOGICAL, SAVE :: FIRSTIME= .TRUE. ! true: first time routine called

        CHARACTER(3000)    DATFMT     ! format for data
        CHARACTER(3000)    HDRFMT     ! format for header
        CHARACTER(3000)    HEADER     ! header for output files
        CHARACTER(3000)    LINFLD     ! line of dashes
        CHARACTER(300)     MESG       ! message buffer
        CHARACTER(IOULEN3) CBUF       ! tmp units field

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

C.............  Allocate memory for the units and names
            ALLOCATE( MNAMES( NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MNAMES', PROGNAME )
            ALLOCATE( MUNITS( NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MUNITS', PROGNAME )

C.............  Initialize units and names
            MUNITS = ' '  ! array
            MNAMES = ' '  ! array
    
C.............  Create units labels from species
            DO V = 1, NMSPC
    
C.................  Set names and units for output
                L = LEN_TRIM( TOTUNIT( V ) )
                CBUF = '[' // TOTUNIT( V )( 1:L ) // ']'
    
                MNAMES( V ) = EMNAM( V )
                MUNITS( V ) = CBUF
    
            END DO

            FIRSTIME = .FALSE.
        
        END IF

C.........  Do not report if this time is not appropriate
        IF( JTIME .NE. REPTIME .AND. 
     &    ( JDATE .NE. EDATE .OR. JTIME .NE. ETIME ) ) RETURN

C.............  If required, create and write state totals
        IF ( LREPSTA ) THEN  

            CALL CREATE_HEADER( 'Mobile' )
            
            PSTATE = -9
            N = 0
            DO I = 1, NC
            
                STATE = CNTYCOD( I ) / 1000
                IF( STATE .NE. PSTATE ) THEN
                    N = N + 1
                    PSTATE = STATE
                END IF
                
                DO J = 1, NMSPC
                    MEBSTA( N,J ) = MEBSTA( N,J ) + MEBCNY( I,J )
                END DO
            
            END DO
            
            CALL WRITE_STA( MRDEV, NS, NMSPC, MNAMES, MUNITS, MEBSTA)

        END IF

C.........  If required, write county totals
        IF( LREPCNY ) THEN

            CALL CREATE_HEADER( 'Mobile' )
            CALL WRITE_CNY( MRDEV, NC, NMSPC, MNAMES, MUNITS, MEBCNY)

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
            LD1 = LEN_TRIM( DATANAM )
     
            TYPENAM = ' inventory'
            TYPENAM = ' speciated'

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
            CHARACTER(IOVLEN3) OUTNAMS( NDIM )
            CHARACTER(IOULEN3) OUTUNIT( NDIM )

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
            CHARACTER(IOVLEN3) OUTNAMS( NDIM )
            CHARACTER(IOULEN3) OUTUNIT( NDIM )

C.............  Local variables
            INTEGER       I, J, L, L1, L2, N
            INTEGER       PSTA, STA

            REAL          VAL

            CHARACTER(FIPLEN3+8) CDATFIP
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
            PSTA = -9
            N = 0
            DO I = 1, NC

                STA = CNTYCOD( I ) / 1000
                IF( STA .NE. PSTA ) THEN
                    N = N + 1
                    PSTA = STA
                END IF

C.................  Write out county name and converted emissions
                WRITE( CDATFIP, '(I7.7,1X,I6.6)' ) JDATE, CNTYCOD( I )

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

