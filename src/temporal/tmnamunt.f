
        SUBROUTINE TMNAMUNT

C***********************************************************************
C  subroutine body starts at line 81
C
C  DESCRIPTION:
C       This program creates the temporal emissions output file variable names
C       and associated activities.  It also sets the units and conversion
C       factors for creating the output emission values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 10/99 by M. Houyoux
C
C*************************************************************************
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER   INDEX1
        CHARACTER(LEN=IOULEN3) MULTUNIT
        REAL                   UNITFAC 

        EXTERNAL     INDEX1, MULTUNIT, UNITFAC

C...........   Other local variables
        INTEGER         I, J, K, L, L2, M     !  counters and indices

        INTEGER         IOS               !  i/o status

        REAL            FAC1, FAC2        ! tmp conversion factors

        LOGICAL      :: FIXDESC = .FALSE. ! true: append info to description
        LOGICAL      :: EFLAG = .FALSE.   ! true: error found

        CHARACTER*16           CURRUNIT   !  current unit
        CHARACTER*16           CURRVNAME  !  current variable name
        CHARACTER*300          MESG       !  message buffer
        CHARACTER(LEN=IOVLEN3) CBUF       !  tmp variable name

        CHARACTER*16 :: PROGNAME = 'TMNAMUNT' ! program name

C***********************************************************************
C   begin body of subroutine TMNAMUNT

C.........  Allocate memory for the names of the emission types and associated
C           arrays for all activities in the inventory
        ALLOCATE( EMTUNT( MXETYPE, NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMTUNT', PROGNAME )
        ALLOCATE( EMTDSC( MXETYPE, NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMTDSC', PROGNAME )

C.........  Allocate memory for units conversions for inventory pollutants and
C           activities (stored in MODINFO)
        ALLOCATE( EACNV( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EACNV', PROGNAME )

C.........  Initialize arrays
        EMTUNT = ' '  ! array
        EMTDSC = ' '  ! array
        EACNV  = 1.   ! array

C.........  Loop through the emission types for each activity and determine 
C           their associated emission factors and units for emission factors
C.........  Also, for each pollutant or activity, store the output units and
C           conversion factor to ton/hr (required output units from Temporal)
C.........  The hours adjustment is part of the temporal allocation, which
C           assumes that the input data are annual data. So, if not, add the
C           conversion to annual data here.
        DO I = 1, NIACT

            M = INDEX1( ACTVTY( I ), NIPPA, EANAM )

            DO K = 1, NETYPE( I )

                CURRVNAME = EMTNAM( K,I )

                FIXDESC = .FALSE.
C.................  Check if pollutant is output hydrocarbon
                IF( OUTPUTHC /= ' ' ) THEN
                    J = INDEX( CURRVNAME, TRIM( OUTPUTHC ) )
                    IF( J > 0 ) THEN
                    	FIXDESC = .TRUE.
                        CURRVNAME = CURRVNAME( 1:J-1 ) // INPUTHC
                    END IF
                END IF

C.................  Search for emission type in emission factors
                J = INDEX1( CURRVNAME, NEFS, EFSNAM )

C.................  Store info if this emissions type is found
C.................  For the units, multiply the emission factor units with the
C                   activity units
                IF( J .GT. 0 ) THEN

C.....................  Ensure that emission factor units are consistent
                    CALL UNITMATCH( EFSUNT( J ) )

                    L  = INDEX( EFSDSC( J ), 'for' )
                    L2 = LEN_TRIM( EFSDSC( J ) )

C.....................  Store for emission types
                    EMTUNT( K,I ) = MULTUNIT( EFSUNT( J ), EAUNIT( M ) )
                    EMTDSC( K,I ) = EFSDSC( J )( L+3:L2 )
                    IF( FIXDESC ) THEN
                        EMTDSC( K,I ) = TRIM( EMTDSC( K,I ) ) // 
     &                                  ' (minus HAPS)'
                    END IF
                    EMTDSC( K,I ) = TRIM( EMTDSC( K,I ) ) // 
     &                              ' from ' // ACTVTY( I )

                ELSE

C.....................  Otherwise, build info for user-defined HAPS
                    CURRUNIT = M6UNIT
                    CALL UNITMATCH( CURRUNIT )
                    
                    EMTUNT( K,I ) = MULTUNIT( CURRUNIT, EAUNIT( M ) )
                    
                    EMTDSC( K,I ) = CURRVNAME
                    EMTDSC( K,I ) = TRIM( EMTDSC( K,I ) ) //
     &                              ' from ' // ACTVTY( I )
                END IF

C.................  If emission type has not been associated with an emission
C                   factor, then error
                IF( EMTDSC( K,I ) .EQ. ' ' ) THEN

                    EFLAG = .TRUE.
                    L = LEN_TRIM( EMTNAM( K,I ) )
                    MESG = 'ERROR: Emission type "' // 
     &                     EMTNAM( K,I )( 1:L ) // '" could not be '//
     &                     'associated with an emission factor!'
                    CALL M3MSG2( MESG )

                END IF

            END DO     ! End of loop through emission types

C.............  Store units and convversion factors for activities and output
C.............  NOTE - this assumes that the units of all emission types
C               from one activity are the same.
            CBUF = EMTUNT( 1,I )
            FAC1 = UNITFAC( CBUF, 'tons', .TRUE. )
            FAC2 = UNITFAC( EAUNIT( M ), '1/yr', .FALSE. )

            IF ( FAC1 .LT. 0. ) FAC1 = 1.
            IF ( FAC2 .LT. 0. ) FAC2 = 1.

            EAUNIT( M ) = 'tons/hr'
            EACNV ( M ) = FAC1 / FAC2

        END DO         ! End of loop through activities

C.........  Now loop through pollutants and create units and conversion factors
        DO I = 1, NIPOL

            M = INDEX1( EINAM( I ), NIPPA, EANAM )
            
            CBUF = EAUNIT ( M )
            FAC1 = UNITFAC( CBUF, 'tons', .TRUE. )
            FAC2 = UNITFAC( EAUNIT( M ), '1/yr', .FALSE. )

            IF ( FAC1 .LT. 0. ) FAC1 = 1.
            IF ( FAC2 .LT. 0. ) FAC2 = 1.

            EAUNIT( M ) = 'tons/hr'
            EACNV ( M ) = FAC1 / FAC2

        END DO

C.........  Abort if error was found
        IF ( EFLAG ) THEN
            MESG = 'Problem with emission types or emission factors'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE TMNAMUNT
