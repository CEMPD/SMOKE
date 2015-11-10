
        SUBROUTINE MRGUNITS

C***********************************************************************
C  subroutine MRGUNITS body starts at line 93
C
C  DESCRIPTION:
C      The purpose of this subroutine
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/99 by M. Houyoux
C
C***********************************************************************
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
C***********************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: SMATCHK, GRDFAC, TOTFAC,
     &                      GRDUNIT, TOTUNIT, 
     &                      NMSPC, NIPPA, NUNITS

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: SPCUNIT_L, SPCUNIT_S, GRDENV, TOTENV

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)    CRLF
        CHARACTER(16)   MULTUNIT
        REAL            UNITFAC

        EXTERNAL  CRLF, MULTUNIT, UNITFAC

C...........   Other local variables

        INTEGER         IOS      ! tmp I/O status
        INTEGER         L, V        ! counter

        REAL            EMFAC, FAC1, FAC2

        CHARACTER(300)  BUFFER   ! text buffer
        CHARACTER(300)  MESG     ! message buffer
        
        CHARACTER(IOULEN3) EMUNIT      ! converted emissions unit

        CHARACTER(IOULEN3) GRDUNIT_I   ! initialized gridded outputs units
        CHARACTER(IOULEN3) GDEN_I      ! initialized gridded denominator
        CHARACTER(IOULEN3) GNUM_I      ! initialized gridded numerator
        CHARACTER(IOULEN3) GDEN        ! work gridded denominator
        CHARACTER(IOULEN3) GNUM        ! work gridded numerator
        CHARACTER(IOULEN3) GRDBUF      ! work gridded output units

        CHARACTER(IOULEN3) TOTUNIT_I   ! initialized output totals units
        CHARACTER(IOULEN3) TDEN_I      ! initialized totals denominator
        CHARACTER(IOULEN3) TNUM_I      ! initialized totals numerator
        CHARACTER(IOULEN3) TDEN        ! work totals denominator
        CHARACTER(IOULEN3) TNUM        ! work totals numerator
        CHARACTER(IOULEN3) TOTBUF      ! work output totals  units

        CHARACTER(16) :: PROGNAME = 'MRGUNITS' ! program name

C***********************************************************************
C   begin body of subroutine MRGUNITS

C.........  Allocate memory for units conversion factors and units.

C.........  If using speciation, allocate arrays to number of species;
C           otherwise, use total number of activities and pollutants.
C           This approach assumes that each species has the same units
C           regardless of which pollutant it comes from, i.e.
C           ALD2 from VOC must have the same units as ALD2 from ACETALD
        NUNITS = NMSPC
        
        ALLOCATE( GRDFAC( NUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRDFAC', PROGNAME )
        ALLOCATE( TOTFAC( NUNITS+NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TOTFAC', PROGNAME )
        ALLOCATE( GRDUNIT( NUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRDUNIT', PROGNAME )
        ALLOCATE( TOTUNIT( NUNITS+NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TOTUNIT', PROGNAME )

C.........  Initialize all
        GRDFAC  = 1.   ! array
        TOTFAC  = 1.   ! array
        GRDUNIT = ' '  ! array
        TOTUNIT = ' '  ! array

C.........  Loop through pollutants, and create output units accordingly
C.........  For speciation, use the first unit for the speciation units from
C           a given pollutants speciation factors
        DO V = 1, NUNITS

C.............  Initialize the output units
            CALL UNITMATCH( SPCUNIT_L( V ) )
            CALL UNITMATCH( SPCUNIT_S( V ) )

C.............  Convert emissions units (g/hr) to tons/hr
            EMUNIT = 'tons/hr'
            GRDUNIT_I = MULTUNIT( SPCUNIT_L( V ), EMUNIT )
            TOTUNIT_I = MULTUNIT( MULTUNIT( SPCUNIT_S( V ), EMUNIT ), 'hr/day' )

C.............  Set the trial units
            GRDBUF = GRDUNIT_I
            TOTBUF = TOTUNIT_I
            IF( GRDENV .NE. ' ' ) GRDBUF = GRDENV
            IF( TOTENV .NE. ' ' ) TOTBUF = TOTENV

C.............  Set the numerators and denominators
            L = INDEX( GRDUNIT_I, '/' ) 
            GNUM_I = ADJUSTL( GRDUNIT_I(   1:L-1     ) )
            GDEN_I = ADJUSTL( GRDUNIT_I( L+1:IOULEN3 ) )

            L = INDEX( TOTUNIT_I, '/' ) 
            TNUM_I = ADJUSTL( TOTUNIT_I(   1:L-1     ) )
            TDEN_I = ADJUSTL( TOTUNIT_I( L+1:IOULEN3 ) )

            L = INDEX( GRDBUF, '/' ) 
            GNUM = ADJUSTL( GRDBUF(   1:L-1     ) )
            GDEN = ADJUSTL( GRDBUF( L+1:IOULEN3 ) )

            L = INDEX( TOTBUF, '/' ) 
            TNUM = ADJUSTL( TOTBUF(   1:L-1     ) )
            TDEN = ADJUSTL( TOTBUF( L+1:IOULEN3 ) )

C.............  Get factor for the numerators for the gridded outputs...
            FAC1 = UNITFAC( SPCUNIT_L( V ), GRDBUF, .TRUE. )  ! speciation

C.............  Get factor for the denominators for the gridded outputs
            FAC2 = UNITFAC( EMUNIT, GRDBUF, .FALSE. )

C.............  In case E.V. setting was bogus rebuild output units based
C               on which factors were valid
C.............  Also set negative factors (from unknown conversions) to 1.
            CALL CORRECT_UNITS( GNUM_I, GDEN_I, GNUM, GDEN, GRDBUF )

C.............  When SMAT is used to calcuate model species
            IF( SMATCHK ) THEN
                EMFAC = UNITFAC( 'g/hr', EMUNIT, .TRUE. )
                FAC1 = FAC1 * EMFAC
            END IF

C.............  Set factors for gridded outputs
            GRDFAC( V ) = FAC1 / FAC2

C.............  Get conversion factor for the numerators for totals
            FAC1 = UNITFAC( SPCUNIT_S( V ), TOTBUF, .TRUE. )  ! speciation

C.............  Get factors for the denominators for the totals.  Note that
C               the hourly data are output as daily totals.
            FAC2 = UNITFAC( TOTUNIT_I, TOTBUF, .FALSE. )

C.............  In case E.V. setting was bogus rebuild output units based
C               on which factors were valid
C.............  Also set negative factors (from unknown conversions) to 1.
            CALL CORRECT_UNITS( TNUM_I, TDEN_I, TNUM, TDEN, TOTBUF )

C.............  Set factors for totaled outputs
            TOTFAC( V ) = FAC1 / FAC2

C.............  Set the output units per pollutant/activity
            GRDUNIT( V ) = GRDBUF
            TOTUNIT( V ) = TOTBUF

        END DO

C.........  Set up units for non-speciated emissions report (tons/day for now)
        EMFAC = UNITFAC( 'g/hr', 'tons/day', .TRUE. )
        DO V = 1, NIPPA
            TOTFAC ( NUNITS+V ) = EMFAC
            TOTUNIT( NUNITS+V ) = 'tons/day'
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS
 
C.............  This internal subprogram corrects badly formatted units
            SUBROUTINE CORRECT_UNITS( NUM_I, DEN_I, NUM, DEN, OUTUNIT )

C.............  Subprogram arguments
            CHARACTER(IOULEN3) NUM_I
            CHARACTER(IOULEN3) DEN_I
            CHARACTER(IOULEN3) NUM
            CHARACTER(IOULEN3) DEN
            CHARACTER(IOULEN3) OUTUNIT

C.............  Local variables
            INTEGER L1,L2

C----------------------------------------------------------------------

            IF( FAC1 .LT. 0. ) THEN
                NUM = NUM_I
                FAC1 = 1.
            END IF
            IF( FAC2 .LT. 0. ) THEN
                DEN = DEN_I
                FAC2 = 1.
            END IF

            L1  = LEN_TRIM( NUM )
            L2  = LEN_TRIM( DEN )
            OUTUNIT =  NUM( 1:L1 ) // '/' // DEN( 1:L2 )
 
            END SUBROUTINE CORRECT_UNITS

        END SUBROUTINE MRGUNITS
