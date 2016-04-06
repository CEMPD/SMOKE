
        SUBROUTINE GRD2CNTY( IDXINV, IDXSPC, NCNTY, CNVFAC,
     &                       GRDDAT, CNYDAT, SGFLAG, SGFAC, 
     &                       FIPTOSG, SGDAT )

C************************************************************************
C  subroutine GRD2CNTY body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to convert gridded values to county-
C      total values for reporting purposes. It can also save totals for
C      source apportionment.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 8/99 by M. Houyoux
C       Added source apportionment handling 07/13 by C. Seppanen
C
C************************************************************************
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
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTY, CNTYCOD

C...........   This module contains the gridding surrogates tables
        USE MODSURG, ONLY: NSRGS, SRGLIST, NSRGFIPS, SRGFIPS, 
     &                     SRGFRAC, FIPCELL, NCELLS

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        INTEGER         ENVINT
        INTEGER         FIND1  
        INTEGER         FINDC

        EXTERNAL   CRLF, ENVINT, FIND1, FINDC

C...........   Subroutine arguments
        INTEGER, INTENT (IN) :: IDXINV
        INTEGER, INTENT (IN) :: IDXSPC
        INTEGER, INTENT (IN) :: NCNTY
        REAL   , INTENT (IN) :: CNVFAC
        REAL   , INTENT (IN) :: GRDDAT( NGRID )
        REAL   , INTENT(OUT) :: CNYDAT( NCNTY,* )
        LOGICAL, INTENT (IN) :: SGFLAG           ! true: store totals for source apportionment
        REAL   , INTENT (IN) :: SGFAC            ! conversion factor for gridded emissions
        INTEGER, INTENT (IN) :: FIPTOSG( NCNTY ) ! source group for each FIPS code
        REAL   , INTENT(OUT) :: SGDAT( NGRID,* )

C...........   Local allocatable arrays
        REAL, ALLOCATABLE :: SRGSUM( : )  ! dim: ngrid - sum of surrogate fracs

C...........   Other local variables

        INTEGER          C, F, J, K, N     ! counters and indices
        INTEGER          IDX               ! index to 2nd dim of CNYDAT
        INTEGER          IOS               ! i/o status
        INTEGER          SSC               ! 
        INTEGER, SAVE :: SRGID             ! surrogate ID for area surrogate
        INTEGER          GIDX              ! source group index

        REAL             FRAC              ! tmp surrogate fraction
        REAL             VAL               ! tmp data value

        LOGICAL, SAVE :: FIRSTIME = .TRUE. ! true: first timwe routine called

        CHARACTER(300)   MESG              ! message buffer

        CHARACTER(16) :: PROGNAME = 'GRD2CNTY' ! program name

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
     &                         'Country/state/county' // SRGFIPS( F ) //
     &                         'and cell', C
                        CALL M3MSG2( MESG )

                    ELSE IF( SRGSUM( C ) .NE. 0. ) THEN
                        SRGFRAC( SRGID,N,F ) = SRGFRAC( SRGID,N,F ) / 
     &                                         SRGSUM( C )
                    END IF

                END DO  ! End loop on cells in county
            END DO      ! End loop on counties in domain

            DEALLOCATE( SRGSUM )

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
            F = FINDC( CNTYCOD( J ), NSRGFIPS, SRGFIPS )

C.............  Skip county if its not in surrogates file
            IF( F .LE. 0 ) CYCLE

C.............  Otherwise, loop through cells in county and get total
            DO N = 1, NCELLS( F )

                C    = FIPCELL( N,F )        ! Retrieve cell number
                FRAC = SRGFRAC( SRGID,N,F )

                VAL = GRDDAT( C ) * FRAC

                CNYDAT( J,IDX ) = CNYDAT( J,IDX ) +
     &                            CNVFAC * VAL

C.................  Store source apportionment data
                IF( SGFLAG ) THEN
                    GIDX = FIPTOSG( F )
                    SGDAT( C,GIDX ) = SGDAT( C,GIDX ) + 
     &                                SGFAC * VAL
                END IF

            END DO  ! End loop on cells in county

        END DO      ! End loop on counties in domain

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GRD2CNTY
