
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
     &                      NMSPC, MEBCNY, MEBSTA, MEBSUM, SMATCHK,
     &                      MEBSCC, MEBSTC,
     &                      MRDEV, LREPSTA, LREPCNY, LREPSCC, LREPSRC,
     &                      EMNAM, EANAM, TOTUNIT, GRDUNIT, NMSRC, NIPPA

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTY, NSTATE, STATNAM, STATCOD, 
     &                     CNTYNAM, CNTYCOD, MICNY

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: MISCC

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVSCC, INVSCC

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)    CRLF
        INTEGER         ENVINT, STR2INT
        INTEGER         INDEX1  
        INTEGER         FINDC
        CHARACTER(14)   MMDDYY
        INTEGER         WKDAY
        REAL            YR2DAY
        CHARACTER(16)   MULTUNIT

        EXTERNAL    CRLF, ENVINT, INDEX1, MMDDYY, WKDAY, YR2DAY,
     &              MULTUNIT, FINDC, STR2INT

C...........   Subroutine arguments
        INTEGER, INTENT (IN) :: JDATE  ! julian date  (YYYYDDD)
        INTEGER, INTENT (IN) :: JTIME  ! time (HHMMSS)

C...........   Local allocatable arrays
        CHARACTER(IOVLEN3), ALLOCATABLE, SAVE :: MNAMES( : )
        CHARACTER(IOULEN3), ALLOCATABLE, SAVE :: MUNITS( : )

C...........   Other local variables

        INTEGER          C, F, I, J, K, L, L2, N, V  ! counters and indices
        INTEGER          SRC, STAIDX, CNTYIDX, SCCIDX
        INTEGER          IOS            ! i/o status
        INTEGER, SAVE :: MAXCYWID       ! max width of county names
        INTEGER, SAVE :: MAXSTWID       ! max width of state names
        INTEGER, SAVE :: MAXSCCWID      ! max width of SCCs
        INTEGER, SAVE :: NC             ! tmp no. counties in domain
        INTEGER, SAVE :: NS             ! tmp no. states in domain
        INTEGER, SAVE :: PDATE          ! date for computing PTIME
        INTEGER, SAVE :: PTIME          ! start time of period
        INTEGER, SAVE :: REPTIME        ! report time

        REAL   , SAVE :: UNITFAC        ! units conversion factor for reports
        REAL          :: VAL            ! current emissions value
 
        LOGICAL, SAVE :: FIRSTIME= .TRUE. ! true: first time routine called

        CHARACTER(FIPLEN3) CSTA       ! current country/state
        CHARACTER(FIPLEN3) PSTA       ! prev country/state
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
            
            MAXSCCWID = 10

            PDATE = SDATE
            PTIME = STIME

C.............  Allocate memory for the units and names
            ALLOCATE( MNAMES( NMSPC+NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MNAMES', PROGNAME )
            ALLOCATE( MUNITS( NMSPC+NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MUNITS', PROGNAME )

C.............  Initialize units and names
            MUNITS = ' '  ! array
            MNAMES = ' '  ! array
    
C.............  Create units names and units from species
            DO V = 1, NMSPC
                IF( .NOT. SMATCHK ) THEN
                    GRDUNIT(V) = MULTUNIT( GRDUNIT( V ), 's/day' )
                
                    L = LEN_TRIM( GRDUNIT( V ) )
                    CBUF = '[' // GRDUNIT( V )( 1:L ) // ']'
                ELSE
                    L = LEN_TRIM( TOTUNIT( V ) )
                    CBUF = '[' // TOTUNIT( V )( 1:L ) // ']'
                END IF 
    
                MNAMES( V ) = EMNAM( V )
                MUNITS( V ) = CBUF
    
            END DO

C.............  Create units labels from pollutants
            DO V = 1, NIPPA
    
C.................  Set names and units for output
                L = LEN_TRIM( TOTUNIT( NMSPC+V ) )
                CBUF = '[' // TOTUNIT( NMSPC+V )( 1:L ) // ']'
    
                MNAMES( NMSPC+V ) = EANAM( V )
                MUNITS( NMSPC+V ) = CBUF
    
            END DO

            FIRSTIME = .FALSE.
        
        END IF

C.........  Do not report if this time is not appropriate
        IF( JTIME .NE. REPTIME ) RETURN 

C.............  Create totals as needed
        PSTA = ' '
        STAIDX = 0
        DO SRC = 1, NMSRC
        
            CNTYIDX = MICNY( SRC )
            SCCIDX  = MISCC( SRC )

            CSTA = CNTYCOD( CNTYIDX )( 1:STALEN3 ) // '000'
            IF( CSTA .NE. PSTA ) THEN
                STAIDX = MAX( FINDC( CSTA, NSTATE, STATCOD ), 0 )
                PSTA = CSTA 
            END IF

            DO J = 1, NMSPC+NIPPA
                VAL = MEBSUM( SRC,J )
                
                IF( LREPSCC ) THEN
                    MEBSCC( SCCIDX,J ) = MEBSCC( SCCIDX,J ) + VAL
                END IF

                IF( LREPCNY ) THEN
                    MEBCNY( CNTYIDX,J ) = MEBCNY( CNTYIDX,J ) + VAL
                END IF
                
                IF( LREPSTA ) THEN
                    MEBSTA( STAIDX,J ) = MEBSTA( STAIDX,J ) + VAL
                    IF( LREPSCC ) THEN
                        MEBSTC( STAIDX,SCCIDX,J ) = 
     &                  MEBSTC( STAIDX,SCCIDX,J ) + VAL
                    END IF
                END IF
            END DO
        
        END DO

C.............  If required, write SCC totals
        IF ( LREPSCC ) THEN
        
            CALL CREATE_HEADER( 'Mobile' )
            CALL WRITE_SCC( MRDEV, NINVSCC, NMSPC+NIPPA, MNAMES, MUNITS, MEBSCC )
        
        END IF

C.............  If required, write state totals
        IF ( LREPSTA ) THEN  

            CALL CREATE_HEADER( 'Mobile' )
            CALL WRITE_STA( MRDEV, NS, NMSPC+NIPPA, MNAMES, MUNITS, MEBSTA)

            IF ( LREPSCC ) THEN

                CALL CREATE_HEADER( 'Mobile' )
                CALL WRITE_STASCC( MRDEV, NS, NINVSCC, NMSPC+NIPPA, MNAMES, MUNITS, MEBSTC )

            END IF

        END IF

C.........  If required, write county totals
        IF( LREPCNY ) THEN

            CALL CREATE_HEADER( 'Mobile' )
            CALL WRITE_CNY( MRDEV, NC, NMSPC+NIPPA, MNAMES, MUNITS, MEBCNY)

            IF ( LREPSCC ) THEN
        
                CALL CREATE_HEADER( 'Mobile' )
                CALL WRITE_CNYSCC( MRDEV, NMSRC, NMSPC+NIPPA, MNAMES, MUNITS, MEBSUM )

            END IF

        END IF

        IF ( LREPSRC ) THEN
        
            CALL CREATE_HEADER( 'Mobile' )
            CALL WRITE_CNYSCC( MRDEV, NMSRC, NMSPC+NIPPA, MNAMES, MUNITS, MEBSUM )

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

            DATANAM = ' summary'
            LD1 = LEN_TRIM( DATANAM )
     
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

            SUBROUTINE WRITE_SCC( FDEV, NSCC, NDIM, VNAMES, INUNIT,
     &                            SCC_EMIS )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: FDEV
            INTEGER     , INTENT (IN) :: NSCC
            INTEGER     , INTENT (IN) :: NDIM
            CHARACTER(*), INTENT (IN) :: VNAMES ( NDIM )
            CHARACTER(*), INTENT (IN) :: INUNIT ( NDIM )
            REAL        , INTENT (IN) :: SCC_EMIS( NSCC, NDIM )

C.............  Arrays allocated by subprogram argument
            INTEGER            MAXWID ( 0:NDIM )
            CHARACTER(IOVLEN3) OUTNAMS( NDIM )
            CHARACTER(IOULEN3) OUTUNIT( NDIM )

C.............  Local variables
            INTEGER       I, J, L, L2

            REAL          VAL

            CHARACTER(7)  :: CDATE
            CHARACTER(28) :: HDRBUF  = '#'
            CHARACTER(28) :: STLABEL = '# Date ; SCC                '
            CHARACTER(30)    BUFFER

C..............................................................................

C.............  Get maximum width of numbers and 
            DO J = 1, NDIM

                VAL = MAXVAL( SCC_EMIS( 1:NSCC, J ) )
                WRITE( BUFFER, '(F30.1)' ) VAL
                BUFFER = ADJUSTL( BUFFER )
                MAXWID( J ) = LEN_TRIM( BUFFER )

            END DO

C.............  Rearrange labels and units to be in order of master list, which
C               is the order that the emission values themselves will be in

C.............  Get column labels and formats
            CALL CREATE_FORMATS( NDIM, 6+MAXSCCWID, VNAMES, INUNIT,
     &                           MAXWID, OUTNAMS, OUTUNIT )

C.............  Write header for SCC totals
            WRITE( FDEV, '(A)' ) '# '
            WRITE( FDEV, '(A)' ) HEADER( 1:LEN_TRIM( HEADER ) )

C.............  Write units for columns
            WRITE( FDEV, HDRFMT ) ADJUSTL( HDRBUF), 
     &                            ( OUTUNIT(J), J=1,NDIM )

C.............  Write column labels
            WRITE( FDEV, HDRFMT ) ADJUSTL( STLABEL ),
     &                          ( OUTNAMS( J ), J=1, NDIM )

C.............  Write SCC total emissions
            DO I = 1, NSCC
            
                WRITE( CDATE, '(I7.7)' ) JDATE

C.................  Build output format depending on data values
                CALL DYNAMIC_FORMATS( NSCC, NDIM, I, SCC_EMIS,
     &                                MAXWID(1), DATFMT )

C.................  Write out SCC and converted emissions
                WRITE( FDEV, DATFMT ) CDATE // ';' // INVSCC( I ), 
     &                                ( SCC_EMIS( I,J ), J=1, NDIM )
            END DO

            RETURN

            END SUBROUTINE WRITE_SCC

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

            CHARACTER(7)  :: CDATE
            CHARACTER(28) :: HDRBUF  = '#'
            CHARACTER(28) :: STLABEL = '# Date ; State              '
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
            DO I = 1, NS
            
                WRITE( CDATE, '(I7.7)' ) JDATE

C.................  Build output format depending on data values
                CALL DYNAMIC_FORMATS( NS, NDIM, I, ST_EMIS,
     &                                MAXWID(1), DATFMT )

C.................  Write out state name and converted emissions
                WRITE( FDEV, DATFMT ) CDATE // ';' // STATNAM( I ), 
     &                                ( ST_EMIS( I,J ), J=1, NDIM )
            END DO

            RETURN

            END SUBROUTINE WRITE_STA

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

            SUBROUTINE WRITE_STASCC( FDEV, NSTA, NSCC, NDIM, VNAMES, INUNIT,
     &                               STSCC_EMIS )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: FDEV
            INTEGER     , INTENT (IN) :: NSTA
            INTEGER     , INTENT (IN) :: NSCC
            INTEGER     , INTENT (IN) :: NDIM
            CHARACTER(*), INTENT (IN) :: VNAMES ( NDIM )
            CHARACTER(*), INTENT (IN) :: INUNIT ( NDIM )
            REAL        , INTENT (IN) :: STSCC_EMIS( NSTA, NSCC, NDIM )

C.............  Arrays allocated by subprogram argument
            INTEGER            MAXWID ( 0:NDIM )
            CHARACTER(IOVLEN3) OUTNAMS( NDIM )
            CHARACTER(IOULEN3) OUTUNIT( NDIM )

C.............  Local variables
            INTEGER       I, J, L, L2

            REAL          VAL

            CHARACTER(7)  :: CDATE
            CHARACTER(49) :: HDRBUF  = '#'
            CHARACTER(49) :: STLABEL='# Date ; State              '
     &                           //'; SCC                '
            CHARACTER(30)    BUFFER

C..............................................................................

C.............  Get maximum width of numbers and 
            DO J = 1, NDIM

                VAL = MAXVAL( STSCC_EMIS( 1:NSTA, 1:NSCC, J ) )
                WRITE( BUFFER, '(F30.1)' ) VAL
                BUFFER = ADJUSTL( BUFFER )
                MAXWID( J ) = LEN_TRIM( BUFFER )

            END DO

C.............  Rearrange labels and units to be in order of master list, which
C               is the order that the emission values themselves will be in

C.............  Get column labels and formats
            CALL CREATE_FORMATS( NDIM, 6+MAXSTWID+MAXSCCWID, VNAMES, INUNIT,
     &                           MAXWID, OUTNAMS, OUTUNIT )

C.............  Write header for state/SCC totals
            WRITE( FDEV, '(A)' ) '# '
            WRITE( FDEV, '(A)' ) HEADER( 1:LEN_TRIM( HEADER ) )

C.............  Write units for columns
            WRITE( FDEV, HDRFMT ) ADJUSTL( HDRBUF), 
     &                            ( OUTUNIT(J), J=1,NDIM )

C.............  Write column labels
            WRITE( FDEV, HDRFMT ) ADJUSTL( STLABEL ),
     &                          ( OUTNAMS( J ), J=1, NDIM )

C.............  Write state/SCC total emissions
            DO I = 1, NSTA

                DO J = 1, NSCC
                    WRITE( CDATE, '(I7.7)' ) JDATE

C.....................  Build output format depending on data values
                    CALL DYNAMIC_FORMATS2( NSTA, NSCC, NDIM, I, J, STSCC_EMIS,
     &                                MAXWID(1), DATFMT )

C.....................  Write out state name and converted emissions
                    WRITE( FDEV, DATFMT ) CDATE // ';' // STATNAM( I ) // 
     &                               ';' // INVSCC( J ), 
     &                                ( STSCC_EMIS( I,J,K ), K=1, NDIM )

                END DO

            END DO

            RETURN

            END SUBROUTINE WRITE_STASCC

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

            REAL          VAL

            CHARACTER(FIPLEN3+8) CDATFIP
            CHARACTER(62) :: HDRBUF  = '#'
            CHARACTER(62) :: STLABEL = '# Date ; FIPS       ; State              '
     &                           //'; County              '
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

                CSTA = CNTYCOD( I )( 1:STALEN3 ) // '000'
                IF( CSTA .NE. PSTA ) THEN
                    N = MAX( FINDC( CSTA, NSTATE, STATCOD ), 0 ) 
                    PSTA = CSTA
                END IF

C.................  Write out county name and converted emissions
                WRITE( CDATFIP, '(I7.7,";",A)' ) JDATE, CNTYCOD( I )

C.................  Build output format depending on data values
                CALL DYNAMIC_FORMATS( NC, NDIM, I, CY_EMIS,
     &                                MAXWID(1), DATFMT )

                WRITE( FDEV,DATFMT ) CDATFIP // ';' // STATNAM(N) // 
     &                               ';' // CNTYNAM(I), 
     &                               ( CY_EMIS( I,J ), J=1, NDIM )

            END DO

            RETURN

            END SUBROUTINE WRITE_CNY

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

            SUBROUTINE WRITE_CNYSCC( FDEV, NSRC, NDIM, VNAMES, INUNIT,
     &                            SRC_EMIS )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: FDEV
            INTEGER     , INTENT (IN) :: NSRC
            INTEGER     , INTENT (IN) :: NDIM
            CHARACTER(*), INTENT (IN) :: VNAMES ( NDIM )
            CHARACTER(*), INTENT (IN) :: INUNIT ( NDIM )
            REAL        , INTENT (IN) :: SRC_EMIS( NSRC, NDIM )

C.............  Arrays allocated by subprogram argument
            INTEGER            MAXWID ( 0:NDIM )
            CHARACTER(IOVLEN3) OUTNAMS( NDIM )
            CHARACTER(IOULEN3) OUTUNIT( NDIM )

C.............  Local variables
            INTEGER       I, J, L, L1, L2, N

            REAL          VAL

            CHARACTER(FIPLEN3+8) CDATFIP
            CHARACTER(83) :: HDRBUF  = '#'
            CHARACTER(83) :: STLABEL = '# Date ; FIPS       ; State              '
     &                           //'; County             ; SCC                '
            CHARACTER(30)    BUFFER

C..............................................................................

C.............  Get maximum width of numbers
            DO J = 1, NDIM

                VAL = MAXVAL( SRC_EMIS( 1:NSRC, J ) )
                WRITE( BUFFER, '(F30.1)' ) VAL
                BUFFER = ADJUSTL( BUFFER )
                MAXWID( J ) = LEN_TRIM( BUFFER )

            END DO

C.............  Get column labels and formats
            CALL CREATE_FORMATS( NDIM, 26+MAXSTWID+MAXCYWID, VNAMES,
     &                           INUNIT, MAXWID, OUTNAMS, OUTUNIT )

C.............  Write header for county totals
            WRITE( FDEV, '(A)' ) '# '
            WRITE( FDEV, '(A)' ) HEADER( 1:LEN_TRIM( HEADER ) )

C.............  Write units for columns
            WRITE( FDEV, HDRFMT ) ADJUSTL( HDRBUF), 
     &                            ( OUTUNIT(J), J=1,NDIM )

C.............  Write column labels
            WRITE( FDEV, HDRFMT ) ADJUSTL( STLABEL ),
     &                            ( OUTNAMS( J ), J=1, NDIM )

C.............  Write county total emissions
            PSTA = ' '
            N = 0
            DO I = 1, NSRC

                CSTA = CNTYCOD( MICNY( I ) )( 1:STALEN3 ) // '000' 
                IF( CSTA .NE. PSTA ) THEN
                    N = MAX( FINDC( CSTA, NSTATE, STATCOD ),0 )
                    PSTA = CSTA
                END IF

C.................  Write out county name and converted emissions
                WRITE( CDATFIP, '(I7.7,";",A)' ) JDATE, CNTYCOD( MICNY( I ) )

C.................  Build output format depending on data values
                CALL DYNAMIC_FORMATS( NSRC, NDIM, I, SRC_EMIS,
     &                                MAXWID(1), DATFMT )

                WRITE( FDEV,DATFMT ) CDATFIP // ';' // STATNAM(N) // 
     &                               ';' // CNTYNAM(MICNY(I)) // 
     &                               ';' // INVSCC(MISCC(I)), 
     &                               ( SRC_EMIS( I,J ), J=1, NDIM )

            END DO

            RETURN

            END SUBROUTINE WRITE_CNYSCC

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

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

            SUBROUTINE DYNAMIC_FORMATS2( N1, N2, N3, STCNT, SCCCNT, EMIS, 
     &                                  WIDTH, DATFMT )

            INTEGER     , INTENT (IN)  :: N1
            INTEGER     , INTENT (IN)  :: N2
            INTEGER     , INTENT (IN)  :: N3
            INTEGER     , INTENT (IN)  :: STCNT      ! counter for state/county index
            INTEGER     , INTENT (IN)  :: SCCCNT     ! counter for SCC index
            REAL        , INTENT (IN)  :: EMIS( N1, N2, N3 )  ! state/count emission
            INTEGER     , INTENT (IN)  :: WIDTH( N3 )
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
            DO I = 1, N3

                IF ( EMIS( STCNT,SCCCNT,I ) .EQ. 0. ) THEN
                   WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.0'

                ELSE IF( EMIS( STCNT,SCCCNT,I ) .GE. 1000. ) THEN
                   WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.1'

                ELSE IF ( EMIS( STCNT,SCCCNT,I ) .GE. 100. ) THEN
                   WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.2'

                ELSE IF ( EMIS( STCNT,SCCCNT,I ) .GE. 10. ) THEN
                   WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.3'

                ELSE IF ( EMIS( STCNT,SCCCNT,I ) .GE. 0.1 ) THEN
                   WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.4'

                ELSE IF ( EMIS( STCNT,SCCCNT,I ) .LT. 0.1 ) THEN 
                   WRITE( FMT, '(A,I2.2,A)' ) 'E',WIDTH(I),'.4'
                END IF

                DATFMT = TRIM( DATFMT ) // ',"; "' // FMT

            END DO

            DATFMT = TRIM( DATFMT ) // ')'

            RETURN

            END SUBROUTINE DYNAMIC_FORMATS2

        END SUBROUTINE WRMRGREP

