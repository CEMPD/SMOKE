
        SUBROUTINE FILLATBL( NXREF, ICSIZE, XTYPE, XTCNT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine populates the area-to-point table that contains
C      the table numbers, row numbers, and counts per FIPS code.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 11/02 by M. Houyoux
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
C***************************************************************************

C...........   This module is for cross reference tables
        USE MODXREF, ONLY: ARPT08, ARPT09, IARPTA, INDXTA 

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        EXTERNAL        CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NXREF           ! no. ungrpd x-ref entries
        INTEGER     , INTENT (IN) :: ICSIZE( * )     ! size of x-ref groups
        INTEGER     , INTENT (IN) :: XTYPE ( NXREF ) ! group no. of x-ref entry
        INTEGER     , INTENT (IN) :: XTCNT ( NXREF ) ! pos. in x-ref group

C...........   Other local variables
        INTEGER       I, J, K, T       ! counter and indices
        INTEGER       NTBL             ! tmp table number
        INTEGER       IROW             ! tmp row number
        INTEGER       ICNT             ! tmp counter of FIPS code

        LOGICAL    :: EFLAG = .FALSE.  ! true: error has occurred

        CHARACTER(257)         MESG    ! message buffer

        CHARACTER(16) :: PROGNAME = 'FILLATBL' ! program name

C***********************************************************************
C   begin body of subroutine FILLATBL

C.........  Store the temporal profile codes for each x-ref entry, depending
C           on the group (XTYPE) and the position in that group (XTCNT)

        DO I = 1, NXREF

            J    = INDXTA ( I )
            NTBL = IARPTA( J,1 )
            IROW = IARPTA( J,2 )
            ICNT = IARPTA( J,3 )

            T      = XTYPE ( I )
            K      = XTCNT ( I )

C.................  Populate tables depending on type. Note that tables
C                   are not pollutant-specific
            SELECT CASE ( T )

            CASE( 0 )  ! Skip this x-ref because it is invalid or duplicate

c            CASE( 1 )  
c                ARPT01( 1 ) = NTBL
c                ARPT01( 2 ) = IROW
c                ARPT01( 3 ) = ICNT

c            CASE( 2 )
c                ARPT02( K,1 ) = NTBL
c                ARPT02( K,2 ) = IROW
c                ARPT02( K,3 ) = ICNT

c            CASE( 3 )
c                ARPT03( K,1 ) = NTBL
c                ARPT03( K,2 ) = IROW
c                ARPT03( K,3 ) = ICNT

c            CASE( 4 )
c                ARPT04( K,1 ) = NTBL
c                ARPT04( K,2 ) = IROW
c                ARPT04( K,3 ) = ICNT

c            CASE( 5 )
c                ARPT05( K,1 ) = NTBL
c                ARPT05( K,2 ) = IROW
c                ARPT05( K,3 ) = ICNT

c            CASE( 6 )
c                ARPT06( K,1 ) = NTBL
c                ARPT06( K,2 ) = IROW
c                ARPT06( K,3 ) = ICNT

c            CASE( 7 )
c                ARPT07( K,1 ) = NTBL
c                ARPT07( K,2 ) = IROW
c                ARPT07( K,3 ) = ICNT

            CASE( 8 )
                ARPT08( K,1 ) = NTBL
                ARPT08( K,2 ) = IROW
                ARPT08( K,3 ) = ICNT

            CASE( 9 )
                ARPT09( K,1 ) = NTBL
                ARPT09( K,2 ) = IROW
                ARPT09( K,3 ) = ICNT
                   
            CASE DEFAULT

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'INTERNAL ERROR: Group', T,
     &                 'not valid in subroutine ' // TRIM(PROGNAME)
                CALL M3MSG2( MESG ) 

            END SELECT

        END DO                            ! End loop on sorted x-ref entries

        IF( EFLAG ) THEN
            MESG = 'Problem processing area-to-point records.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE FILLATBL
