
        SUBROUTINE RDREPIN( GDIM, NSLIN, NSSIN, SDEV, EDEV, YDEV, NDEV, 
     &                      ENAME, GNAME, LNAME, SLNAME, SSNAME, 
     &                      NX, IX, CX, SSMAT, SLMAT )

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
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C  
C See file COPYRIGHT for conditions of use.
C  
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C  
C env_progs@mcnc.org
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

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: GDIM   ! no. mass spec input vars
        INTEGER     , INTENT (IN) :: NSLIN  ! no. mass spec input vars
        INTEGER     , INTENT (IN) :: NSSIN  ! no. mass spec input vars
        INTEGER     , INTENT (IN) :: SDEV   ! unit no.: ASCII inven file
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
 
C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C...........   Local variables that depend on module variables
        INTEGER    SWIDTH( NCHARS )

C...........   Other local variables
        INTEGER          I, J, K, L, L1, L2, V, S, T ! counters and indices

        INTEGER       :: JDATE = 0          ! Julian date
        INTEGER       :: JTIME = 0          ! time (HHMMSS)
        INTEGER       :: NINVARR = 0        !  no. actual inventory inputs

        LOGICAL       :: LRDREGN = .FALSE.  !  true: read region code
        LOGICAL       :: EFLAG   = .FALSE.  !  true: error found

         
        CHARACTER*16  ::       BNAME = ' '  !  name buffer
        CHARACTER*50           BUFFER       !  string buffer
        CHARACTER*300          MESG         !  message buffer
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

            IF( CATEGORY .EQ. 'POINT' ) THEN
                CX = 1   ! array
            END IF

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

C.........  If needed, read in country, state, county file
        IF( YFLAG ) THEN
            CALL RDSTCY( YDEV, NINVIFIP, INVIFIP )
        END IF

C.........  If needed, read in elevated source indentification file
        IF( VFLAG ) THEN
            CALL RDPELV( EDEV, NSRC, NMAJOR, NPING )
        END IF

C.........  If needed, read in SCC descriptions file
        IF( NFLAG ) THEN
c           CALL RDSCCDSC( NDEV, NINVSCC, INVSCC )
c note: still must write this
        END IF

C.........  If needed, read in layer fractions file to identify elevated
C           sources
        IF( LFLAG ) THEN

            JDATE = SDATE
            JTIME = STIME
            DO T = 1, NSTEPS

c note: add this later. May need a new array in the MODELEV module?  Maybe can
C    n: use LMAJOR

            END DO

        END IF

C.........  Reformat source characteristics and set widths.  Do this once
C           for the entire run of the program, so that it doesn't have to be
C           done for each report (it is slow)

        IF( ANY_TRUE( NREPORT, ALLRPT%BYSRC  ) ) THEN

C.............  Determine width of source chararactistic columns over the
C               whole inventory
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

C.............  Reset start and end fields based on new widths
            K = 0
            DO J = MINC, NCHARS
                K = K + 1
                IF( J .EQ. JSCC ) CYCLE
                SC_BEGP( J ) = SC_ENDP( J-1 ) + 1
                SC_ENDP( J ) = SC_BEGP( J ) + SWIDTH( K ) - 1
            END DO

            IF( JSCC .GT. 0 ) NCHARS = NCHARS - 1

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

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

