
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
C     Modified 11/01 by Gabe Cano - deterministic mode
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

C...........   This module is for cross reference tables
        USE MODXREF

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
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

        CHARACTER*300          MESG    ! message buffer

        CHARACTER*16 :: PROGNAME = 'FILLGTBL' ! program name

C***********************************************************************
C   begin body of subroutine FILLGTBL

C.........  Store the temporal profile codes for each x-ref entry, depending
C           on the group (XTYPE) and the position in that group (XTCNT)

        DO I = 1, NXREF

            J    = INDXTA ( I )
            ISRG = J
c            ISRG = ISRGCDA( J, 1 )

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
     &                 'not valid in subroutine', PROGNAME
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
