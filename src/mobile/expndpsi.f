
        SUBROUTINE EXPNDPSI( MXPSI, NP, PSILOC )
   
C***********************************************************************
C  subroutine EXPNDPSI body starts at line < >
C
C  DESCRIPTION:
C      Expands the list of unique PSIs from the orignal EF cross-reference
C      list to an updated one with the component PSIs from the combination
C      emission factors.  Also expands the list to include any root PSIs from
C      a PSI scenario group that are not included in the PSI list from the
C      cross-reference (via the MEFTEMP file).
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
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
C****************************************************************************

C...........   MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C.........  External functions
        CHARACTER*2   CRLF
        INTEGER       FIND1

        EXTERNAL      CRLF, FIND1

C.........  Subroutine arguments
        INTEGER    , INTENT (IN) :: MXPSI
        INTEGER    , INTENT (IN) :: NP               ! no. PSIs for current call
        INTEGER    , INTENT (IN) :: PSILOC( NP )     ! PSIs for current call

C.........  Arrays declared by subroutine arguments

        INTEGER INDEXA ( MXPSI )
        INTEGER PSIALLA( MXPSI )

C.........  Other local variables

        INTEGER         I, J, K, L, M, N

        INTEGER         CHKPSI    ! tmp PSI from combo list
        INTEGER         CINDX     ! tmp pntr to combo PSIs table or pure scnario
        INTEGER         CNTCOMBO  ! tmp number of PSIs in a combo PSI
        INTEGER         PSIPNTR   ! tmp PSI pointer (see below)

        LOGICAL      :: EFLAG = .FALSE.  ! true: error found

        CHARACTER*300   MESG      ! message buffer

        CHARACTER*16 :: PROGNAME = 'EXPNDPSI' ! program name

C***********************************************************************
C   begin body of subroutine EXPNDPSI

C.........  The PDATPNTR points to either the table of contributing PSIs for
C           combination PSIs, or it is the position of the PSI in a multi-
C           scenario group.

C......... Loop through PSIs from emission factor cross-reference and
C          for any of these that are combo factors, add non-combo PSIs that
C          are not already in the list
        K = 0
        DO I = 1, NP

C.............  Find PSI in emission factor data file
            N = FIND1( PSILOC( I ), NPSIDAT, PSIDAT )
            J = PDATINDX( N )
            CNTCOMBO = PDATTYPE( J )
            PSIPNTR  = PDATPNTR( J )  ! Pointer (see above)

C.............  For pure factors that are not the first in a group, make sure
C               that the first PSI in the group is in the list
            IF( CNTCOMBO .EQ. 0 .AND. PSIPNTR .NE. 1 ) THEN

                N = FIND1( PDATROOT( J ), NPSIDAT, PSIDAT )

                IF( N .LE. 0 ) THEN
                    K = K + 1
                    IF( K .LE. MXPSI ) THEN
                        INDEXA ( K ) = K
                        PSIALLA( K ) = PDATROOT( N )
                    END IF
                ELSE
                    K = K + 1
                    IF( K .LE. MXPSI ) THEN
                	INDEXA ( K ) = K
                	PSIALLA( K ) = PSILOC( I )
                    END IF
                END IF

C.............  For pure factors that are the first in a group, simply 
C               store the PSI
            ELSE IF( CNTCOMBO .EQ. 0 ) THEN

                K = K + 1
                IF( K .LE. MXPSI ) THEN
                    INDEXA ( K ) = K
                    PSIALLA( K ) = PSILOC( I )
                END IF

C.............  Otherwise, check the component PSIs in the combo and add these
C               to the unsorted list
            ELSE

C.................  Add the current PSI to the list, and then the component PSIs
                K = K + 1
                IF( K .LE. MXPSI ) THEN
                    INDEXA ( K ) = K
                    PSIALLA( K ) = PSILOC( I )
                END IF

C.................  Retrieve index to combination table
                CINDX = PDATPNTR( J )

C.................  Loop through component PSIs for this combo PSI
                DO L = 1, CNTCOMBO

C..................... Search for component PSI in master list, and if its not
C                      there, then add it to the list
                    CHKPSI = CMBOPSI( L, CINDX )
                    M = FIND1( CHKPSI, NP, PSILOC )

                    IF( M .LE. 0 ) THEN
                        K = K + 1
                        IF( K .LE. MXPSI ) THEN
                            INDEXA ( K ) = K
                            PSIALLA( K ) = CHKPSI
                        END IF
                    END IF

                END DO  ! end loop on component PSIs

            END IF      ! if combination PSI or not
 
        END DO          ! loop through PSI list from EF cross-reference

        NPSIALL = K

C.........  Ensure that memory was properly allocated        
        IF( NPSIALL .GT. MXPSI ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &             'storing expanded list of PSIs was', MXPSI,
     &             CRLF() // BLANK10 // 'but actually needed', NPSIALL
            CALL M3MSG2( MESG )

        END IF

C.........  Abort if an error was encountered
        IF( EFLAG ) THEN

            MESG = 'Problem updating PSI list with component PSIs'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Sort new list of PSIs
        CALL SORTI1( NPSIALL, INDEXA, PSIALLA )

C.........  Store sorted final list of PSIs
        DO I = 1, NPSIALL

            J = INDEXA( I )
            PSIALL( I ) = PSIALLA( J )

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE EXPNDPSI




