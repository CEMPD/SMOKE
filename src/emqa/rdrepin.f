
        SUBROUTINE RDREPIN( GDIM, NSLIN, NSSIN, SDEV, GDEV, PDEV, TDEV,
     &                      EDEV, YDEV, NDEV, ENAME, GNAME, LNAME, 
     &                      SLNAME, SSNAME, NX, IX, CX, SSMAT, SLMAT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      The RDREPIN routine reads in the SMOKE intermediate files and other
C      files needed for generating the reports.
C
C  PRECONDITIONS REQUIRED:
C    REPCONFIG file is opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 7/2000 by M Houyoux
C
C***********************************************************************
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
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains Smkreport-specific settings
        USE MODREPRT

C.........  This module contains report arrays for each output bin
        USE MODREPBN

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........  EXTERNAL FUNCTIONS and their descriptions:
        INTEGER    GETFLINE
        INTEGER    INDEX1

        EXTERNAL   GETFLINE, INDEX1

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: GDIM   ! no. mass spec input vars
        INTEGER     , INTENT (IN) :: NSLIN  ! no. mass spec input vars
        INTEGER     , INTENT (IN) :: NSSIN  ! no. mass spec input vars
        INTEGER     , INTENT (IN) :: SDEV   ! unit no.: ASCII inven file
        INTEGER     , INTENT (IN) :: GDEV   ! unit no.: gridding supplemental
        INTEGER     , INTENT (IN) :: PDEV   ! unit no.: speciation supplemental
        INTEGER     , INTENT (IN) :: TDEV   ! unit no.: temporal supplemental
        INTEGER     , INTENT (IN) :: EDEV   ! unit no.: elevated ID file (PELV)
        INTEGER     , INTENT (IN) :: YDEV   ! unit no.: cy/st/co file
        INTEGER     , INTENT (IN) :: NDEV   ! unit no.: SCC descriptions
        CHARACTER(*), INTENT (IN) :: ENAME  ! name for I/O API inven input
        CHARACTER(*), INTENT (IN) :: GNAME  ! gridding matrix name
        CHARACTER(*), INTENT (IN) :: LNAME  ! layer fractions file name
        CHARACTER(*), INTENT (IN) :: SLNAME ! speciation matrix name
        CHARACTER(*), INTENT (IN) :: SSNAME ! speciation matrix name
        INTEGER     , INTENT(OUT) :: NX( NGRID ) ! no. srcs per cell
        INTEGER     , INTENT(OUT) :: IX( NMATX ) ! src IDs
        REAL        , INTENT(OUT) :: CX( NMATX ) ! gridding coefficients
        REAL        , INTENT(OUT) :: SLMAT( NSRC, NSLIN ) ! mole spec coefs
        REAL        , INTENT(OUT) :: SSMAT( NSRC, NSSIN ) ! mass spec coefs
 
C.........  Local allocatable arrays
        INTEGER, ALLOCATABLE :: IBUF   ( : )  ! tmp var for temporal profiles
        REAL   , ALLOCATABLE :: LFRAC1L( : )  ! 1st-layer fraction

C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C...........   Local variables that depend on module variables
        INTEGER    SWIDTH( NCHARS )

C...........   Other local variables
        INTEGER          I, J, K, L, L1, L2, N, V, S, T ! counters and indices

        INTEGER          DIU                ! tmp diurnal profile number
        INTEGER          IOS                ! i/o status
        INTEGER          IREC               ! tmp record number
        INTEGER       :: JDATE = 0          ! Julian date
        INTEGER       :: JTIME = 0          ! time (HHMMSS)
        INTEGER          MON                ! tmp monthly profile number
        INTEGER       :: NINVARR = 0        ! no. actual inventory inputs
        INTEGER          NV                 ! tmp no. variables in temporal suplm
        INTEGER       :: SRGID1             ! tmp primary surrogate IDs
        INTEGER       :: SRGID2             ! tmp fallback surrogate IDs
        INTEGER          WEK                ! tmp weekly profile number

        LOGICAL       :: LRDREGN = .FALSE.  !  true: read region code
        LOGICAL       :: EFLAG   = .FALSE.  !  true: error found
        LOGICAL       :: LTMP    = .FALSE.  !  true: temporary logical

        CHARACTER*1            TTYP         !  temporal profile entry type
        CHARACTER*16  ::       BNAME = ' '  !  name buffer
        CHARACTER*50           BUFFER       !  string buffer
        CHARACTER*300          MESG         !  message buffer
        CHARACTER(LEN=IOVLEN3) CBUF         !  tmp pollutant name
        CHARACTER(LEN=IODLEN3) DBUF         !  tmp variable name
        CHARACTER(LEN=SRCLEN3) CSRC         !  tmp source chars

        CHARACTER*16 :: PROGNAME = 'RDREPIN' ! program name

C***********************************************************************
C   begin body of subroutine RDREPIN

C.........  Set local variables for determining input inventory variables
        LRDREGN = ( ANY_TRUE( NREPORT, ALLRPT%BYCNRY ) .OR.
     &              ANY_TRUE( NREPORT, ALLRPT%BYSTAT ) .OR.
     &              ANY_TRUE( NREPORT, ALLRPT%BYCNTY ) .OR.
     &              ANY_CVAL( NREPORT, ALLRPT%REGNNAM )     )

C.........  Build array of inventory variable names based on report settings
C.........  Region code
        IF( LRDREGN ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'IFIP'
        END IF

C.........  Road class code
        IF( ANY_TRUE( NREPORT, ALLRPT%BYRCL ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'IRCLAS'
        END IF

C.........  SCC code
        IF( ANY_TRUE( NREPORT, ALLRPT%BYSCC ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CSCC'
        END IF

C.........  Source description
        IF( ANY_TRUE( NREPORT, ALLRPT%BYSRC ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CSOURC'
        END IF

C.........  Stack parameters
        IF( ANY_TRUE( NREPORT, ALLRPT%STKPARM ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'STKHT'
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'STKDM'
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'STKTK'
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'STKVE'
        END IF

C.........  Plant name
        IF( ANY_TRUE( NREPORT, ALLRPT%SRCNAM ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CPDESC'
        END IF

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Create unique source characteristic lists
        CALL GENUSLST

C.........  If needed, read in gridding matrix
C.........  Initialize all to 1 for point sources
        IF( GFLAG ) THEN

            CALL RDGMAT( GNAME, NGRID, NMATX, NMATX, NX, IX, CX )

C.............  Initialize part of gridding matrix array for point sources
            IF( CATEGORY .EQ. 'POINT' ) THEN
                CX = 1   ! array            
            END IF

        END IF

C.........  If needed, read in gridding supplementation matrix
        IF( GSFLAG ) THEN

            ALLOCATE( SRGID( NSRC,2 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRGID', PROGNAME )
            SRGID = -9

            MESG = 'Supplemental gridding file'
            N = GETFLINE( GDEV, MESG )

            IREC = 0
            DO I = 1, N

                READ( GDEV, *, END=999, IOSTAT=IOS ) S, SRGID1, SRGID2
                IREC = IREC + 1

                IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                'I/O error', IOS, 
     &                'reading supplemental gridding file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

                SRGID( S,1 ) = SRGID1
                SRGID( S,2 ) = SRGID2
                
            END DO            

        END IF

C.........  If needed, read in speciation matrices
C.........  NOTE that only the variables that are needed are read in 
        IF( SLFLAG .OR. SSFLAG ) THEN

            IF( SLNAME .NE. ' ' ) BNAME = SLNAME
            IF( SSNAME .NE. ' ' ) BNAME = SSNAME

C.............  Get file header for variable names
            IF ( .NOT. DESC3( BNAME ) ) THEN

                MESG = 'Could not get description of file "' //
     &                 BNAME( 1:LEN_TRIM( BNAME ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ENDIF

            K = 0
            DO V = 1, NSVARS

                IF( SPCOUT( V ) ) THEN

                    K = K + 1
                    DBUF = VDESC3D( V )
                    IF( SLFLAG ) THEN
                        CALL RDSMAT( SLNAME, DBUF, SLMAT( 1,K ) )
                    END IF

                    IF( SSFLAG ) THEN
                        CALL RDSMAT( SSNAME, DBUF, SSMAT( 1,K ) )
                    END IF

                END IF

            END DO

        END IF

C.........  If needed, read in speciation supplementation matrix
        IF( PSFLAG ) THEN

            ALLOCATE( SPPROF( NSRC,NSPCPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPPROF', PROGNAME )
            SPPROF = ' '

            MESG = 'supplemental speciation file'
            N = GETFLINE( PDEV, MESG )

            IREC = 0
            
            DO I = 1, N

                READ( PDEV, '(A)', END=999, IOSTAT=IOS ) BUFFER
                IREC = IREC + 1

                IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                'I/O error', IOS, 'reading supplemental ' //
     &                'speciation file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  See if this line is a pollutant name
                L1 = INDEX( BUFFER, '"' )        ! Find start quote

C.................  If pollutant name, figure out which pollutant index and
C                   reset source counter to 0.
                IF ( L1 .GT. 0 ) THEN

                    L2 = LEN_TRIM( BUFFER )     ! Find end quote
                    CBUF = BUFFER( L1+1:L2-1 )  

C.....................  Check if this pollutant is one selected for reporting 
                    V = INDEX1( CBUF, NSPCPOL, SPCPOL )
                    IF ( V .GT. 0 ) THEN
                        LSPCPOL( V ) = .TRUE.
                        S = 0
                    END IF
                        
C.................  If not pollutant name, then continue to read in the 
C                   pollutant codes and store them by source
                ELSE IF ( V .GT. 0 ) THEN
                    S = S + 1

                    BUFFER = ADJUSTL( BUFFER )
                    SPPROF( S,V ) = ADJUSTR( BUFFER( 1:SPNLEN3 ) )

C.....................  Handle case where spaces have been added to file.
                    IF ( S .EQ. NSRC ) V = 0

                END IF
                
            END DO            

        END IF

C.........  If needed, read in temporal supplementation matrix
        IF( TSFLAG ) THEN

            ALLOCATE( IDIU( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IDIU', PROGNAME )
            ALLOCATE( IWEK( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IWEK', PROGNAME )
            ALLOCATE( IMON( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IMON', PROGNAME )
            ALLOCATE( IBUF( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IBUF', PROGNAME )
            IDIU = 0
            IWEK = 0
            IMON = 0
            IBUF = 0

            MESG = 'Supplemental temporal file'
            N = GETFLINE( TDEV, MESG )

C.............  Check file format, assuming that pollutants weren't processed
C               in > 1 groups.  (this routine doesn't handle grouped processing)
            IF ( ( (N-1)/3 ) .NE. NSRC ) THEN
                MESG = 'INTERNAL ERROR: ' // CRL// 'TSUP file has '//
     &                 'inconsistent number of lines with NSRC'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END IF

C.............  Skip file header
            READ( TDEV, * )

            IREC = 1
            DO I = 2, N

                READ( TDEV, *, END=999, IOSTAT=IOS ) 
     &              TTYP, NV, ( IBUF( V ), V=1, NV )
                IREC = IREC + 1

                IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                'I/O error', IOS, 
     &                'reading supplemental temporal file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  NV > 1 is not supported
                IF ( NV .GT. 1 ) THEN
                    MESG = 'Number of variables in ' // CRL //
     &                     'TSUP file is not supported.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                S = INT( ( I-2 ) / 3 ) + 1
                SELECT CASE( TTYP )
                CASE ( 'M' )
                    IF ( NV .EQ. 1 ) THEN
                        IMON( S ) = IBUF( 1 )
                    ELSE
c                   note: this is not supported yet in Temporal
                    ENDIF

                CASE ( 'W' )
                    IF ( NV .EQ. 1 ) THEN
                        IWEK( S ) = IBUF( 1 )
                    ELSE
c                   note: this is not supported yet in Temporal
                    ENDIF

                CASE ( 'H' )
                    IF ( NV .EQ. 1 ) THEN
                        IDIU( S ) = IBUF( 1 )
                    ELSE
c                   note: this is not supported yet in Temporal
                    ENDIF

                END SELECT
              
            END DO            

        END IF

C.........  If needed, read in country, state, county file
        IF( YFLAG ) THEN
            CALL RDSTCY( YDEV, NINVIFIP, INVIFIP )
        END IF

C.........  If needed, read in elevated source indentification file
        IF( VFLAG ) THEN
            LTMP = ( .NOT. LFLAG )
            CALL RDPELV( EDEV, NSRC, LTMP, NMAJOR, NPING )
        END IF

C.........  If needed, read in SCC descriptions file
        IF( NFLAG ) CALL RDSCCDSC( NDEV )

C.........  If needed, read in layer fractions file to identify elevated
C           sources
        IF( LFLAG ) THEN

            IF( .NOT. ALLOCATED( LMAJOR ) ) THEN
                ALLOCATE( LMAJOR( NSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'LMAJOR', PROGNAME )
                LMAJOR = .FALSE.   ! array
            END IF

            IF( .NOT. ALLOCATED( LPING ) ) THEN
                ALLOCATE( LPING( NSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'LPING', PROGNAME )
                LPING  = .FALSE.   ! array
            END IF

            ALLOCATE( LFRAC1L( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LFRAC1L', PROGNAME )

            JDATE = SDATE
            JTIME = STIME
            DO T = 1, NSTEPS

                IF( READ3( LNAME, 'LFRAC', 1,
     &                     JDATE, JTIME, LFRAC1L ) ) THEN
                    DO S = 1, NSRC
                        IF( LFRAC1L( S ) .LT. 1. ) LMAJOR( S ) = .TRUE.
                    END DO

                ELSE  !  Read failed
                    MESG = 'Could not read "LFRAC" from '// LNAME
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                END IF

                CALL NEXTIME( JDATE, JTIME, TSTEP )

            END DO

        END IF

C.........  Reformat source characteristics and set widths.  Do this once
C           for the entire run of the program, so that it doesn't have to be
C           done for each report (it is slow)

        IF( ANY_TRUE( NREPORT, ALLRPT%BYSRC  ) ) THEN

C.............  Determine width of source chararactistic columns over the
C               whole inventory
            SWIDTH = 0   ! initialize array
            DO S = 1, NSRC

                K = 0
                DO J = MINC, NCHARS

                    K  = K + 1
                    L1 = SC_BEGP( J )
                    L2 = SC_ENDP( J )
                    BUFFER = ADJUSTL( CSOURC( S )( L1:L2 ) )
                    SWIDTH( K ) = MAX( SWIDTH( K ), LEN_TRIM( BUFFER ) )

                END DO

            END DO

C.............  Reset CSOURC based on these widths
C.............  Also remove SCC from the source characteristics
            DO S = 1, NSRC

                CSRC = CSOURC( S )

                L = SC_ENDP( MINC-1 ) 
                K  = 0
                DO J = MINC, NCHARS

                    K = K + 1
                    IF( J .NE. JSCC ) THEN
                        L1 = SC_BEGP( J )
                        L2 = SC_ENDP( J )
                        CSOURC( S ) = CSOURC( S )( 1:L ) //
     &                                ADJUSTL( CSRC( L1:L2 ) )
                        L = L + SWIDTH( K )
                    END IF

                END DO

            END DO

C.............  Allocate and initialize arrays for storing new field lengths
            ALLOCATE( LOC_BEGP( NCHARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LOC_BEGP', PROGNAME )
            ALLOCATE( LOC_ENDP( NCHARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LOC_ENDP', PROGNAME )
            LOC_BEGP = SC_BEGP   ! array
            LOC_ENDP = SC_ENDP   ! array

C.............  Set local start and end fields based on new widths
            K = 0
            DO J = MINC, NCHARS
                K = K + 1
                IF( J .EQ. JSCC ) CYCLE
                LOC_BEGP( J ) = SC_ENDP( J-1 ) + 1
                LOC_ENDP( J ) = SC_BEGP( J ) + SWIDTH( K ) - 1
            END DO

            IF( JSCC .GT. 0 ) NCHARS = NCHARS - 1

        END IF

C.........  Deallocate local memory
        IF( ALLOCATED( LFRAC1L ) ) DEALLOCATE( LFRAC1L )

        RETURN

999     MESG = 'Unexpected end of file reached while reading ' //
     &         'supplementary gridding file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS
 
C.............  This internal function scans a logical array for any
C               true values, and if it finds one, returns true.
            LOGICAL FUNCTION ANY_TRUE( NDIM, LOGARR )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: NDIM
            LOGICAL, INTENT (IN) :: LOGARR( NDIM )

C.............  Local variables
            INTEGER   I

C----------------------------------------------------------------------

            ANY_TRUE = .FALSE.

            DO I = 1, NDIM

                IF( LOGARR( I ) ) THEN
                    ANY_TRUE = .TRUE.
                    RETURN
                END IF

            END DO

            RETURN
 
            END FUNCTION ANY_TRUE

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal function scans a character array for any
C               non-blank values, and if it finds one, returns true.
            LOGICAL FUNCTION ANY_CVAL( NDIM, CHARARR )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: NDIM
            CHARACTER(*), INTENT (IN) :: CHARARR( NDIM )

C.............  Local variables
            INTEGER   I

C----------------------------------------------------------------------

            ANY_CVAL = .FALSE.

            DO I = 1, NDIM

                IF( CHARARR( I ) .NE. ' ' ) THEN
                    ANY_CVAL = .TRUE.
                    RETURN
                END IF

            END DO

            RETURN
 
            END FUNCTION ANY_CVAL

        END SUBROUTINE RDREPIN

