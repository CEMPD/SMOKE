
        SUBROUTINE FILLMTBL( NXREF, ICSIZE, XTYPE, XTCNT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine populates the mobile data part of the grouped VMT Mix
C      
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 2/2000 by M. Houyoux
C
C****************************************************************************/
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
C Last updated: %G 
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
        INTEGER       IMVS             ! tmp for gridding surrogate code

        LOGICAL    :: EFLAG = .FALSE.  ! true: error has occurred

        CHARACTER*300          MESG    ! message buffer

        CHARACTER*16 :: PROGNAME = 'FILLMTBL' ! program name

C***********************************************************************
C   begin body of subroutine FILLMTBL

C.........  Store the position for each unsorted x-ref entry, depending
C           on the group (XTYPE) and the position in that group (XTCNT)

        DO I = 1, NXREF

            J    = INDXTA ( I )

            T      = XTYPE ( I )
            K      = XTCNT ( I )

C.................  Populate tables depending on type. Note that tables
C                   are not pollutant-specific
            SELECT CASE ( T )

            CASE( 0 )  ! Skip this x-ref because it is invalid or duplicate

            CASE( 1 )  
                IMVS01 = J

            CASE( 2 )
                IMVS02( K ) = J

            CASE( 3 )
                IMVS03( K ) = J

            CASE( 4 )
                IMVS04( K ) = J

            CASE( 5 )
                IMVS05( K ) = J

            CASE( 6 )
                IMVS06( K ) = J

            CASE( 7 )
                IMVS07( K ) = J

            CASE( 8 )
                IMVS08( K ) = J

            CASE( 9 )
                IMVS09( K ) = J
                                                            
            CASE( 10 )
                IMVS10( K ) = J
                                                            
            CASE( 11 )
                IMVS11( K ) = J
                                                            
            CASE( 12 )
                IMVS12( K ) = J
                                                            
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

        END SUBROUTINE FILLMTBL
