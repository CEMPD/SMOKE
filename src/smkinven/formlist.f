
        SUBROUTINE FORMLIST

C**************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine is designed to read SMKINVEN_FORMULA and determine how many formulas and list of variables 
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 5/2012 by B. Baek
C      09/2025 by HT UNC-IE: Use M3UTILIO; Removed MESG format 93000 which is not used anywhere
C
C**************************************************************************
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
C**************************************************************************
        USE M3UTILIO

C...........   Modules for public variables
C...........   This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVTBL, ITNAMA, ITCASA
     
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NMAP,NIACT, NIPPA, NPPOL, EANAM, ACTVTY,
     &                     NCOMP, VAR_FORMULA, CHKPLUS, CHKMINUS, 
     &                     FORMULAS, VIN_A, VIN_B, VNAME
        IMPLICIT NONE

C.........  EXTERNAL FUNCTIONS
C       CHARACTER(2)  CRLF
C       INTEGER       INDEX1
C       EXTERNAL      CRLF, INDEX1

C.........   Subroutine arguments
        
C.........   Local variables
        INTEGER         I, J, F, L, N, V1, V2, VA, VB    ! indices and counters

        INTEGER         IOS                   ! i/o status
        INTEGER         LEQU                  ! position of '=' in formula
        INTEGER         LDIV                  ! position of '-' or '+' in formula
        INTEGER         LMNS                  ! position of '-' in formula
        INTEGER         LPLS                  ! position of '+' in formula

        LOGICAL,SAVE :: FIRSTIME = .TRUE.     ! true: first time 
        LOGICAL      :: TFLAG    = .FALSE.    ! true: current formula has error
        LOGICAL      :: EFLAG    = .FALSE.    ! true: error found

        CHARACTER(300)  MESG                  ! Message buffer

        CHARACTER(16) :: PROGNAME = 'FORMLIST'    !  program name

C***********************************************************************
C   Begin body of subroutine FORMULAS

C.........  Figure out how many variables there are based on the
C           number of commas found in the string.
        NCOMP = 1
        L = LEN_TRIM( VAR_FORMULA )
        DO I = 1, L
            IF( VAR_FORMULA( I:I ) == ',' ) NCOMP = NCOMP + 1
        ENDDO

        NMAP = NMAP + NCOMP  ! NCOMP more variable(s) to map

C.........  Allocate array to store formulas
        IF( FIRSTIME ) THEN
            ALLOCATE( CHKPLUS( NCOMP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHKPLUS', PROGNAME )
            ALLOCATE( CHKMINUS( NCOMP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHKMINUS', PROGNAME )
            ALLOCATE( FORMULAS( NCOMP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FORMULAS', PROGNAME )
            ALLOCATE( VIN_A( NCOMP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VIN_A', PROGNAME )
            ALLOCATE( VIN_B( NCOMP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VIN_B', PROGNAME )
            ALLOCATE( VNAME( NCOMP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VNAME', PROGNAME )
            FIRSTIME = .FALSE.
        END IF
C.........  Split out formulas in string to array
        CALL PARSLINE( VAR_FORMULA, NCOMP, FORMULAS )

C.........  Loop through formulas
        DO F = 1, NCOMP

C.............  Make sure formula makes sense
            LEQU = INDEX( FORMULAS( F ), '=' )
            LPLS = INDEX( FORMULAS( F ), '+' )
            LMNS = INDEX( FORMULAS( F ), '-' )

            CHKPLUS( F )  = ( LPLS .GT. 0 )
            CHKMINUS( F ) = ( LMNS .GT. 0 )

            LDIV = LPLS
            IF( CHKMINUS( F ) ) LDIV = LMNS

            IF( LEQU .LE. 0 .OR.
     &        ( .NOT. CHKPLUS(F) .AND. .NOT. CHKMINUS(F) ) ) THEN
                MESG = 'Could not interpret formula for extra ' //
     &                 'pollutant from ' // 'SMKINVEN_FORMULAS=' //
     &                 TRIM( FORMULAS( F ) )
                CALL M3MSG2( MESG )
                TFLAG = .TRUE.
            END IF

C.............  Extract formula variable names
            L      = LEN_TRIM( FORMULAS( F ) )
            VNAME( F )= ADJUSTL( FORMULAS( F )(      1:LEQU-1 ) )
            VIN_A( F )= ADJUSTL( FORMULAS( F )( LEQU+1:LDIV-1 ) )
            VIN_B( F )= ADJUSTL( FORMULAS( F )( LDIV+1:L      ) )

C.............  Find formula inputs in existing variable list
            J = INDEX1( VIN_A( F ), NINVTBL, ITCASA )
            IF( J < 1 ) THEN
                L = LEN_TRIM( VIN_A( F ) )
                MESG = 'Variable "'// VIN_A( F )( 1:L ) // 
     &             '" from formula was not found in inventory ' //
     &             'pollutant code (CAS nubmer)'
                CALL M3MSG2( MESG )

            ELSE
                VIN_A( F ) = ITNAMA( J )

            END IF

            J = INDEX1( VIN_B( F ), NINVTBL, ITCASA )
            IF( J < 1 ) THEN
                L = LEN_TRIM( VIN_B( F ) )
                MESG = 'Variable "'// VIN_B( F )( 1:L ) // 
     &             '" from formula was not found in inventory ' //
     &             'pollutant code (CAS nubmer)'
                CALL M3MSG2( MESG )

            ELSE
                VIN_B( F ) = ITNAMA( J )

            END IF

            VA = INDEX1( VIN_A( F ), NIPPA, EANAM )
            VB = INDEX1( VIN_B( F ), NIPPA, EANAM )

            IF( VA .LE. 0 ) THEN
                TFLAG = .TRUE.
                L = LEN_TRIM( VIN_A( F ) )
                MESG = 'Variable "'// VIN_A( F )( 1:L ) // 
     &                 '" from formula was not found in inventory.'
                CALL M3MSG2( MESG )
            END IF

            IF( VB .LE. 0 ) THEN
                TFLAG = .TRUE.
                L = LEN_TRIM( VIN_B( F ) )
                MESG = 'Variable "'// VIN_B( F )( 1:L ) // 
     &                 '" from formula was not found in inventory.'
                CALL M3MSG2( MESG )
            END IF

            V1 = INDEX1( VIN_A( F ), NIACT, ACTVTY )
            V2 = INDEX1( VIN_B( F ), NIACT, ACTVTY )

            IF( V1 .GT. 0 ) THEN
                TFLAG = .TRUE.
                L = LEN_TRIM( VIN_A( F ) )
                MESG = 'ERROR: Variable "'//VIN_A(F)(1:L)//'" is an'//
     &                 'activity, which is not allowed in a formula.'
                CALL M3MSG2( MESG )
            END IF

            IF( V2 .GT. 0 ) THEN
                TFLAG = .TRUE.
                L = LEN_TRIM( VIN_B( F ) )
                MESG = 'ERROR: Variable "'//VIN_B(F)(1:L)//'" is an'//
     &                 'activity, which is not allowed in a formula.'
                CALL M3MSG2( MESG )
            END IF

            IF( TFLAG ) THEN
                WRITE( MESG,94010 ) 'ERROR: Problem processing '//
     &                 'formula', F, ': "'//TRIM(FORMULAS(F)) // '"'
                CALL M3MSG2( MESG )
            END IF

            IF ( TFLAG ) EFLAG = .TRUE.

        END DO

        IF ( EFLAG ) THEN
            MESG = 'ERROR: Problem processing formulas. ' //
     &             'See previous error messages.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C******************  FORMAT  STATEMENTS   ******************************

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE FORMLIST
