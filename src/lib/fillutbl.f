
        SUBROUTINE FILLUTBL( NUOVAR, NXREF, ICSIZE, XTYPE, XTCNT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine populates the uncertainty cross-reference table.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 9/2001 by A. Holland
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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

C.........  MODULES for public variables
C.........  This module is for cross reference tables
        USE MODXREF

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: NUOVAR         ! no. pollutants + activities
        INTEGER, INTENT (IN) :: NXREF           ! no. ungrpd x-ref entries
        INTEGER, INTENT (IN) :: ICSIZE( * )     ! size of x-ref groups
        INTEGER, INTENT (IN) :: XTYPE ( NXREF ) ! group no. of x-ref entry
        INTEGER, INTENT (IN) :: XTCNT ( NXREF ) ! pos. in x-ref group

C...........   Other local variables
        INTEGER       I, J, K, T     ! counter and indices
        INTEGER       IDX            ! tmp index to control data table
        INTEGER       IDXPS          ! IDX with ADDPS added for pol-specific
        INTEGER       ISP            ! tmp pollutant/activity position index
        INTEGER       TDIM           ! temporary table dimension

        CHARACTER*16 :: PROGNAME = 'FILLUTBL' ! program name

C***********************************************************************
C   begin body of subroutine FILLUTBL

C.........  Store the temporal profile codes for each x-ref entry, depending
C           on the group (XTYPE) and the position in that group (XTCNT)

        DO I = 1, NXREF

            J      = INDXTA( I )
            ISP    = ISPTA ( J )

            T      = XTYPE ( I )
            K      = XTCNT ( I )

C.............  Skip invalid or duplicate records
            IF( T .LE. 0 ) CYCLE

            TDIM   = ICSIZE( T )

            IDX    = J
            IDXPS  = IDX + ADDPS

C.............  Populate tables depending on type. Note that the pollutant/
C               activity-specific entries are assumed to always come after the
C               non-specific ones (based on the previous sorting).
C.............  The pol-specific entries are stored by adding 90000
C               to the monthly profile number (which has a maximum of 3 
C               digits) so that the pol-specific can be identified later
            SELECT CASE ( T )

            CASE( 1 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL01 )

            CASE( 2 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL02 )

            CASE( 3 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL03 )

            CASE( 4 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL04 )

            CASE( 5 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL05 )

            CASE( 6 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL06 )

            CASE( 7 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL07 )

            CASE( 8 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL08 )

            CASE( 9 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL09 )

            CASE( 10 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL10 )
                 
            CASE( 11 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL11 )
   
            CASE( 12 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL12 )

            CASE( 13 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL13 )

            CASE( 14 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL14 )

            CASE( 15 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL15 )
                      
            CASE( 16 )
                CALL SET_CNTRL_INDEX( TDIM, NUOVAR, ICTL16 )
                      
            CASE DEFAULT

            END SELECT

        END DO                            ! End Loop on sorted x-ref entries

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram stores the appropriate index to the
C               control factors tables.  Note that the subprogram inherits 
C               variables from the main program.
            SUBROUTINE SET_CNTRL_INDEX( N, M, ITABLE )

C.............  Subprogram arguments
            INTEGER   N             ! local size for dimensioning
            INTEGER   M             ! local size for dimensioning
            INTEGER   ITABLE( N,M ) ! index to control data table

C......................................................................

            IF( ISP .EQ. 0 ) THEN
                ITABLE( K,: ) = IDX

            ELSE
                ITABLE( K,ISP ) = IDXPS

            ENDIF

            END SUBROUTINE SET_CNTRL_INDEX

        END SUBROUTINE FILLUTBL
