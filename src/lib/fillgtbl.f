
        SUBROUTINE FILLGTBL( NXREF, ICSIZE, XTYPE, XTCNT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine populates the gridding surrogate codes part of the 
C      grouped gridding cross-reference tables.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 4/99 by M. Houyoux
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
        USE MODXREF, ONLY: ISRG01, ISRG02, ISRG03, ISRG04, ISRG05,
     &                     ISRG06, ISRG07, ISRG08, ISRG09,
     &                     INDXTA, ISRGCDA

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
        INTEGER       ISRG             ! tmp for gridding surrogate code

        LOGICAL    :: EFLAG = .FALSE.  ! true: error has occurred

        CHARACTER(300)         MESG    ! message buffer

        CHARACTER(16) :: PROGNAME = 'FILLGTBL' ! program name

C***********************************************************************
C   begin body of subroutine FILLGTBL

C.........  Store the temporal profile codes for each x-ref entry, depending
C           on the group (XTYPE) and the position in that group (XTCNT)

        DO I = 1, NXREF

            J    = INDXTA ( I )
            ISRG = ISRGCDA( J )

            T      = XTYPE ( I )
            K      = XTCNT ( I )

C.................  Populate tables depending on type. Note that tables
C                   are not pollutant-specific
            SELECT CASE ( T )

            CASE( 0 )  ! Skip this x-ref because it is invalid or duplicate

            CASE( 1 )  
                ISRG01 = ISRG

            CASE( 2 )
                ISRG02( K ) = ISRG

            CASE( 3 )
                ISRG03( K ) = ISRG

            CASE( 4 )
                ISRG04( K ) = ISRG

            CASE( 5 )
                ISRG05( K ) = ISRG

            CASE( 6 )
                ISRG06( K ) = ISRG

            CASE( 7 )
                ISRG07( K ) = ISRG

            CASE( 8 )
                ISRG08( K ) = ISRG

            CASE( 9 )
                ISRG09( K ) = ISRG
                                                            
            CASE DEFAULT

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'INTERNAL ERROR: Group', T,
     &                 'not valid in subroutine ' // TRIM(PROGNAME)
                CALL M3MESG( MESG ) 

            END SELECT

        END DO                            ! End Loop on sorted x-ref entries

        IF( EFLAG ) THEN
            MESG = 'Problem processing cross-reference records.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE FILLGTBL
