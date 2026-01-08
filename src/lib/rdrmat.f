
        SUBROUTINE RDRMAT( FNAME, NSREAC, SPECNUM, IDX, REPEM,
     &                     PRJFC, MKTPN, RFAC )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads a reactivity matrix for any source category.
C      Certain variable names are expected to be in the files, so if the
C      reactivity matrix writer is changed, then this routine needs to be
C      changed as well.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C****************************************************************************/
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2)    CRLF
C        EXTERNAL        CRLF

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME     ! speciation matrix file name
        INTEGER     , INTENT (IN) :: NSREAC    ! no. reactivity sources
        INTEGER     , INTENT (IN) :: SPECNUM   ! no. source-spec species (all)
        INTEGER     , INTENT(OUT) :: IDX  ( NSREAC ) ! source index
        REAL        , INTENT(OUT) :: REPEM( NSREAC ) ! replacement emissions
        REAL        , INTENT(OUT) :: PRJFC( NSREAC ) ! projection factors
        REAL        , INTENT(OUT) :: MKTPN( NSREAC ) ! market penetration
        REAL        , INTENT(OUT) :: RFAC ( NSREAC, SPECNUM ) ! spc coeffs

C...........   Error message strings
        CHARACTER(23), PARAMETER :: PART1 = 'Error reading variable '
        CHARACTER(23), PARAMETER :: PART3 = ' from REACTIVITY MATRIX'

C.........  Other local variables
        INTEGER            N, V       !  counters and indices

        CHARACTER(300)     MESG    ! message buffer
        CHARACTER(IOVLEN3) INVAR    ! tmp inventory pollutant name

        CHARACTER(16) :: PROGNAME = 'RDRMAT' ! program name

C***********************************************************************
C   begin body of subroutine RDRMAT

C.........  Retrieve file header
        IF ( .NOT. DESC3( FNAME ) ) THEN
            MESG = 'Could not get description of file ' // FNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Loop through variables and read them
        V = 0
        DO N = 1, NVARS3D

            INVAR = VNAME3D( N )

            MESG = PART1 // INVAR( 1:LEN_TRIM( INVAR ) ) // PART3

C.............  Section for reading non-species variables
            SELECT CASE( INVAR )

            CASE( 'SRCID' )

                IF( .NOT. READ3( FNAME,INVAR,ALLAYS3,0,0,IDX ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF
                CYCLE

            CASE( 'REPEMIS' )

                IF( .NOT. READ3( FNAME,INVAR,ALLAYS3,0,0,REPEM ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF
                CYCLE

            CASE( 'PRJFAC' )

                IF( .NOT. READ3( FNAME,INVAR,ALLAYS3,0,0,PRJFC ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF
                CYCLE

            CASE( 'MKTPEN' )

                IF( .NOT. READ3( FNAME,INVAR,ALLAYS3,0,0,MKTPN ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF
                CYCLE

            END SELECT

C.............  Section for reading species variables

            V = V + 1   ! Increment count for species-specific variables

            IF( V .LE. SPECNUM ) THEN

                IF( .NOT. READ3(FNAME,INVAR,ALLAYS3,0,0,RFAC(1,V))) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            END IF
             
        END DO  ! End loop on variables indicated for reading

        IF( V .GT. SPECNUM ) THEN

            WRITE( MESG,94010 ) 
     &             'INTERNAL ERROR: Dimension mismatch. Memory ' // 
     &             'needed for reactivity factors was', V, CRLF() //
     &             BLANK10 // 'but', SPECNUM, 'was allocated.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE RDRMAT
