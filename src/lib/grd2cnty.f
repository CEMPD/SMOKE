
        SUBROUTINE GRD2CNTY( IDXINV, IDXSPC, NGRID, NCNTY, CNVFAC,
     &                       GRDDAT, CNYDAT )

C************************************************************************
C  subroutine GRD2CNTY body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to convert gridded values to county-
C      total values for reporting purposes
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 8/99 by M. Houyoux
C
C************************************************************************
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
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the arrays for state and county summaries
        USE MODSTCY

C...........   This module contains the gridding surrogates tables
        USE MODSURG

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        INTEGER         FIND1  

        EXTERNAL   CRLF, ENVINT, FIND1

C...........   Subroutine arguments
        INTEGER, INTENT (IN) :: IDXINV
        INTEGER, INTENT (IN) :: IDXSPC
        INTEGER, INTENT (IN) :: NGRID
        INTEGER, INTENT (IN) :: NCNTY
        REAL   , INTENT (IN) :: CNVFAC
        REAL   , INTENT (IN) :: GRDDAT( NGRID )
        REAL   , INTENT(OUT) :: CNYDAT( NCNTY,* )

C...........   Local allocatable arrays
        REAL, ALLOCATABLE :: SRGSUM( : )  ! dim: ngrid - sum of surrogate fracs

C...........   Other local variables

        INTEGER          C, F, J, K, N     ! counters and indices
        INTEGER          IDX               ! index to 2nd dim of CNYDAT
        INTEGER          IOS               ! i/o status
        INTEGER          SSC               ! 
        INTEGER, SAVE :: SRGID             ! surrogate ID for area surrogate

        REAL             FRAC              ! tmp surrogate fraction

        LOGICAL, SAVE :: FIRSTIME = .TRUE. ! true: first timwe routine called

        CHARACTER*300    MESG              ! message buffer

        CHARACTER*16  :: PROGNAME = 'GRD2CNTY' ! program name

C***********************************************************************
C   begin body of subroutine GRD2CNTY
        
        IF( FIRSTIME ) THEN

C.............  Get surrogate code to use as area surrogate
            MESG = 'Code number for the area surrogate'
            SSC = ENVINT( 'AREA_SURROGATE_NUM', MESG, 60, IOS )

C.............  Find index of surrogate code in list of these
            SRGID = FIND1( SSC, NSRGS, SRGLIST )

            IF( SRGID .LE. 0 ) THEN

                WRITE( MESG,94010 ) 'The AREA_SURROGATE_NUM '//
     &                 'environment variable was set to surrogate '//
     &                 'no.', SSC, CRLF() // BLANK10 //
     &                 'but this does not exist in the surrogates data!'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

C.............  Convert surrogate fractions to grid-to-county fractions...
C.............  Allocatable memory to sum surrogates by cell
            ALLOCATE( SRGSUM( NGRID ), STAT=IOS ) 
            CALL CHECKMEM( IOS, 'SRGSUM', PROGNAME )
            SRGSUM = 0.  ! array

C.............  Sum by grid cell
            DO F = 1, NSRGFIPS
        	DO N = 1, NCELLS( F )

                    C    = FIPCELL( N,F )        ! Retrieve cell number
                    FRAC = SRGFRAC( SRGID,N,F )
                    SRGSUM( C ) = SRGSUM( C ) + FRAC

        	END DO  ! End loop on cells in county
            END DO      ! End loop on counties in domain

C.............  Divide surrogate value by sum on cell
            DO F = 1, NSRGFIPS
        	DO N = 1, NCELLS( F )

                    C = FIPCELL( N,F )
                    IF( SRGFRAC( SRGID,N,F ) .EQ. 0. ) THEN
                        WRITE( MESG,94010 )
     &                         'WARNING: Area surrogate is 0. for ' //
     &                         CRLF() // BLANK10 //
     &                         'Country/state/county', SRGFIPS( F ),
     &                         'and cell', C
                        CALL M3MSG2( MESG )

                    ELSE IF( SRGSUM( C ) .NE. 0. ) THEN
                        SRGFRAC( SRGID,N,F ) = SRGFRAC( SRGID,N,F ) / 
     &                                         SRGSUM( C )
                    END IF

        	END DO  ! End loop on cells in county
            END DO      ! End loop on counties in domain

            FIRSTIME = .FALSE.

        END IF

C.........  Species index
        IF( IDXSPC .GT. 0 ) THEN
            IDX = IDXSPC

C.........  Pollutant index
       ELSE IF( IDXINV .GT. 0 ) THEN
            IDX = IDXINV

C.........  If this pollutant or species is not valid for current call, return
        ELSE
            RETURN

        END IF

C.........  Loop through county codes and compute county total emissions
        DO J = 1, NCOUNTY

C.............  Make sure county is in surrogates file
            F = FIND1( CNTYCOD( J ), NSRGFIPS, SRGFIPS )

C.............  Skip county if its not in surrogates file
            IF( F .LE. 0 ) CYCLE

C.............  Otherwise, loop through cells in county and get total
            DO N = 1, NCELLS( F )

                C    = FIPCELL( N,F )        ! Retrieve cell number
                FRAC = SRGFRAC( SRGID,N,F )

                CNYDAT( J,IDX ) = CNYDAT( J,IDX ) + 
     &                            CNVFAC * GRDDAT( C ) * FRAC

            END DO  ! End loop on cells in county

        END DO      ! End loop on counties in domain

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GRD2CNTY
