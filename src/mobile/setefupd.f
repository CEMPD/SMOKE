
        SUBROUTINE SETEFUPD( NUPDAT, NDREUSE, DIREUSE, NNAME, DNAME, 
     &                       UPDATE, UPDATNDI, UPDATDIU )
   
C***********************************************************************
C  subroutine SETEFUPD body starts at line < >
C
C  DESCRIPTION:
C      This subroutine sets the flags for each PSI for non-diurnal and
C      diurnal as to which PSIs need to be updated/created for a 
C      the circumstances of the run.  The determination is affected by
C      whether there is reuse of the emission factors files, combination
C      emission factors, whether a PSI is the first in a group of 
C      scenarios, and the targeted update inputs.
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
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS
        LOGICAL     CHKEMFAC
        INTEGER     FIND1
        INTEGER     SEC2TIME

        EXTERNAL    CHKEMFAC, FIND1, SEC2TIME

C...........   SUBROUTINE ARGUMENTS
        INTEGER,      INTENT (IN):: NUPDAT              ! no. trgtd update PSIs
        LOGICAL,      INTENT (IN):: NDREUSE             ! true: reusing nondiur
        LOGICAL,      INTENT (IN):: DIREUSE             ! true: reusing diur
        CHARACTER(*), INTENT (IN):: NNAME               ! name of non-diur file
        CHARACTER(*), INTENT (IN):: DNAME               ! name of diur file
        INTEGER,      INTENT (IN):: UPDATE  ( NUPDAT  ) ! targeted update PSIs
        LOGICAL,      INTENT(OUT):: UPDATNDI( NPSIALL ) ! true: update non-diur
        LOGICAL,      INTENT(OUT):: UPDATDIU( NPSIALL ) ! true: update diur

C...........   Local variables
        INTEGER         I, J, K, L, M, N ! counters and indices

        INTEGER         CINDX        ! index to contributing PSIs and factors 
        INTEGER         CNTCOMBO     ! tmp count for combo EFs (or 0 for pure)
        INTEGER         CPSI         ! tmp contributing PSI to a combo PSI
        INTEGER         DDATE        ! julian date for reading diurnal EFs
        INTEGER         IPTIM        ! PSI converted to time for I/O API read
        INTEGER         NDATE        ! julian date for reading non-diur EFs
        INTEGER         NPSISCN      ! no. PSIs in multi-scenario group
        INTEGER         PSI          ! tmp PSI
        INTEGER         PSIROOT      ! root PSI of multi-scenario group

        CHARACTER*300   MESG         ! message buffer

        CHARACTER*16 :: PROGNAME = 'SETEFUPD' ! program name

C***********************************************************************
C   begin body of subroutine SETEFUPD

C.........  Get dates for checking non-diurnal and diurnal emission factors
C           files
        CALL GET_FILE_DATE( NNAME, NDATE )
        CALL GET_FILE_DATE( DNAME, DDATE )

C.........  If new non-diurnal emission factors file, initialize update 
C           array to true for all PSIs
        IF( NDREUSE ) THEN
            UPDATNDI = .FALSE.   ! array
        ELSE
            UPDATNDI = .TRUE.    ! array
        END IF

C.........  If new diurnal emission factors file, initialize update array
C           to true for all PSIs
        IF( DIREUSE ) THEN
            UPDATDIU = .FALSE.   ! array
        ELSE
            UPDATDIU = .TRUE.    ! array
        END IF

C.........  Loop through all PSIs and set update arrays depending on contents
C           of reused emission factor files or the targeted updates arrays
C.........  The PDATPNTR points to either the table of contributing PSIs for
C           combination PSIs, or it is the position of the PSI in a multi-
C           scenario group.
        DO I = 1, NPSIALL

            PSI   = PSIALL  ( I )
            IPTIM = SEC2TIME( PSI )

            J = FIND1( PSI, NPSIDAT, PSIDAT )
            K = PDATINDX( J )

            CNTCOMBO = PDATTYPE( K ) ! pure or combo
            CINDX    = PDATPNTR( K ) ! pointer (see above)
            NPSISCN  = PDATMCNT( K ) ! no. PSIs in multi-scenario group
            PSIROOT  = PDATROOT( K ) ! root PSI for multi-scenario group

C.............  Check non-diurnal factors to see if current PSI is already 
C               there and set flag to update it if it's not
C.............  NOTE - if any non-diurnals are there, they all will be in
C               current version of the code.
            IF( .NOT. CHKEMFAC( 'NON-DIURNAL', NNAME, NDATE, 
     &                          IPTIM, NDREUSE, .FALSE.           ) )
     &              UPDATNDI( I ) = .TRUE.

C.............  Check diurnal factors to see if current PSI is already there
C               and set flag to update it if it's not
            IF( .NOT. CHKEMFAC( 'DIURNAL', DNAME, DDATE, 
     &                          IPTIM, DIREUSE, .FALSE.           ) )
     &              UPDATDIU( I ) = .TRUE.

C.............  Search for PSI in update list  
            N = FIND1( PSI, NUPDAT, UPDATE )

C.............  If PSI is in update list, check various conditions, and set the
C               update arrays accordingly for the affected PSIs.
            IF( N .GT. 0 ) THEN

C.................  If this PSI is not the first in a multi-scenario input of
C                   PSIs, then make sure that the first PSI in the multi-
C                   scenario list will be updated.
                IF( CNTCOMBO .EQ. 0 .AND. PSI .NE. PSIROOT ) THEN

                    M = FIND1( PSIROOT, NPSIALL, PSIALL )
                    UPDATNDI( M ) = .TRUE.
                    UPDATDIU( M ) = .TRUE.

C.................  If this PSI is a combination PSI, make sure that all 
C                   contributing PSIs will be updated as well
                ELSE IF( CNTCOMBO .GT. 0 ) THEN
                
                    DO L = 1, CNTCOMBO

                        CPSI = CMBOPSI( L,CINDX )
                        M   = FIND1( CPSI, NPSIALL, PSIALL )

                	UPDATNDI( M ) = .TRUE.
                	UPDATDIU( M ) = .TRUE.

                    END DO  ! End loop on contributing PSIs

                END IF      ! End if for condition of PSI

C.................  Make sure the PSI itself is updated
                UPDATNDI( I ) = .TRUE.
                UPDATDIU( I ) = .TRUE.

            END IF          ! End if for type of PSI

C.............  No matter what is going on otherwise, cannot process the PSI
C               if there are no temperature indices for the PSI...
            IF( MMTEFIDX( I ) .LE. 0 ) THEN
                UPDATNDI( I ) = .FALSE.
                UPDATDIU( I ) = .FALSE.
            END IF

        END DO              ! End loop on all PSIs

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This subprogram retrieves the date from an i/o api file

            SUBROUTINE GET_FILE_DATE( FNAME, JDATE )

            CHARACTER(*), INTENT (IN) :: FNAME
            INTEGER     , INTENT(OUT) :: JDATE

C----------------------------------------------------------------------------

            IF( .NOT. DESC3( FNAME ) ) THEN

                MESG = 'Could not read header of file "' // FNAME // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF
            JDATE = SDATE3D

            RETURN

            END SUBROUTINE GET_FILE_DATE

        END SUBROUTINE SETEFUPD
