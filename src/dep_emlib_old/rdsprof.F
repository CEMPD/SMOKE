C ****************************************************************************
C *
C * Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
C * File: @(#)$Id$
C *
C * COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C * All Rights Reserved
C *
C * See file COPYRIGHT for conditions of use.
C *
C * Environmental Programs Group
C * MCNC--North Carolina Supercomputing Center
C * P.O. Box 12889
C * Research Triangle Park, NC  27709-2889
C *
C * env_progs@mcnc.org
C *
C * Pathname: $Source$
C * Last updated: $Date$ 
C *
C ****************************************************************************/
 
        INTEGER FUNCTION RDSPROF( SDEV, NP, PLIST, N, 
     &                            SXREF, SCOEF, IPOL )

C***********************************************************************
C  function body starts at line  107
C
C  FUNCTION:
C       Reads the indicated speciation profile file from DEV, 
C       sorts it, and returns the number of profile entries found.
C       In case of overflow, returns negative of number of profile
C       entries found.
C
C  PRECONDITIONS REQUIRED:
C       speciation profile file already opened on unit SDEV
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       M3WARN, SORTI1
C
C  REVISION  HISTORY:
C       prototype 1/96 by CJC
C       updated for variable species 5/95 by MRH
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'CHDIMS3.EXT'     !  emis chem parms (inventory + model)
      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations

C...........   ARGUMENTS and their descriptions:

        INTEGER     SDEV        !  unit number for speciation profile file
        INTEGER     N           !  dimension for profile codes
        INTEGER     NP          !  number of model species used
        INTEGER     SXREF( N )  !  profile indexes (sorted)
        REAL        SCOEF( N, MXMPOL )!  profile coefficients
        CHARACTER*8 PLIST( NP ) !  actual model species list
        CHARACTER*5 IPOL        !  inventory pollutant for speciation profile


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER     INDEX1
        EXTERNAL    INDEX1


C...........   SAVED LOCAL VARIABLES and their descriptions:
C...........   NOTE:  the ANSI standard requires the use of SAVE statements
C...........   for variables which must retain their values from call to call.

        LOGICAL         FIRSTIME
        DATA            FIRSTIME / .TRUE. /

        SAVE   FIRSTIME

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         I, J, K, L, V
        INTEGER         IREC, IOS
        INTEGER         NSPEC            ! number of model spcs in profile file
        INTEGER         SID, LID
        INTEGER         INDX ( MXMPOL )  ! pointer from SPRO spcs to actual spcs
        INTEGER         INDEXA( MXSPRO ) ! profile subscripts for sorting
        INTEGER         SXREFA( MXSPRO ) ! profile indexes (sorted)
        REAL            SCOEFA( MXSPRO, MXMPOL ) !  profile coefficients
        LOGICAL         EFLAG
        LOGICAL         SORTED
        REAL            FACS ( MXMPOL )

        CHARACTER*5     ISPEC( MXMPOL )  ! temporary for inventory spcs name
        CHARACTER*5     LSPEC            ! temporary for inventory spcs name
        CHARACTER*10    MSPEC( MXMPOL )  ! model spcs names in profile file
        CHARACTER*22    LINE
        CHARACTER*256   MESG


C***********************************************************************
C   begin body of function  RDSPROF

        IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
        END IF          !  if firstime

        EFLAG = .FALSE.

C.........  Read number of species
        READ( SDEV, * ) NSPEC

        IF( NSPEC .GT. MXMPOL ) THEN
            WRITE( MESG,94010 )
     &          'Max model species', MXMPOL,
     &          'exceeded in SPRO file: count=', NSPEC
            CALL M3MESG( MESG )
            RDSPROF = IMISS3
            RETURN    

        ENDIF

C.........  Read list of model species
        DO 44 I = 1, NSPEC

            READ( SDEV, 93010, IOSTAT=IOS ) LINE
            MSPEC( I ) = LINE(  1:10 )
            ISPEC( I ) = LINE( 18:22 )

            J = INDEX1( ISPEC( I ), MXIPOL, INVSPEC )

            IF( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'Error', IOS, 
     &              'reading SPECIATION XREF file at line', I+1
                CALL M3MESG( MESG )

            ELSEIF( J .EQ. 0 ) THEN   ! Invalid pollutant ID

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'ERROR: Pollutant ' // ISPEC( I ) // 
     &              ' in SPRO; no match in valid list'
                CALL M3MESG( MESG )

            ENDIF

44      CONTINUE

C.........  Abort if problem with species in header
        IF( EFLAG ) THEN
            RDSPROF = IMISS3
            RETURN

        ENDIF

C.........  Check to make sure ISPECs are all the same
        DO 50 I = 1, NSPEC
            IF( I .EQ. 1 ) THEN
                LSPEC = ISPEC( I )
            ELSEIF( ISPEC( I ) .NE. LSPEC ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'Inconsistent inventory pollutant ' //
     &              ISPEC( I ) // ' found at line ', I,
     &              'of SPECIATION PROFILE file'
                CALL M3MESG( MESG )
             
            ENDIF
50      CONTINUE

        IPOL = LSPEC
        
C.........  Find and index actual species to profile format
        DO 55 I = 1,NP

            INDX( I ) = INDEX1( PLIST( I ), NSPEC, MSPEC )

            IF( INDX( I ) .EQ. 0 ) THEN

                WRITE( MESG,94010 )
     &              'NOTE: Species ' // PLIST( I ) // 
     &              ' not found in SPRO file'
                CALL M3MESG( MESG )

            ENDIF

55      CONTINUE

C.........  Read SPECIATION PROFILE file
        IREC   =  0
        I      =  0
        EFLAG  = .FALSE.
        SORTED = .TRUE.
        LID    = -1
77      CONTINUE        !  head of the SDEV-read loop

            READ( SDEV, *, END=99, IOSTAT=IOS ) 
     &          SID, ( FACS( V ), V = 1, NSPEC )

            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 'reading SPC PROFs at line', IREC
                CALL M3MESG( MESG )

            ELSE

                I = I + 1
                IF ( I .LE. MXSPRO ) THEN

                    INDEXA( I ) = I
                    SORTED      = ( SORTED .AND. ( LID .LT. SID ) )
                    LID         = SID
                    SXREFA( I ) = SID
                    DO  88  V = 1, NSPEC
                        SCOEFA( I,V ) = FACS( V )
88                  CONTINUE

                ELSE

                    WRITE( MESG, 94010 )
     &                  'Max SPROF table size', MXSPRO, 
     &                  'exceeded at line', IREC
                    CALL M3MESG( MESG )

                END IF

            END IF

            GO TO  77   !  to head of loop

99      CONTINUE        !  end of the SDEV-read loop

        IF ( EFLAG ) THEN

            CALL M3WARN( 'RDSPROF', 0, 0, 
     &                   'Error reading SPECIATION PROFILE file.' )
            RDSPROF = IMISS3
            RETURN

        ELSE IF ( I .GT. MXSPRO ) THEN

            WRITE( MESG,94010 ) 
     &          'Internal Table Overflow.  Dimensioned', MXSPRO,
     &          'Needed', I
            CALL M3WARN( 'RDSPROF', 0, 0, MESG )
            RDSPROF = -I
            RETURN

        ELSE IF ( I .GT. N ) THEN

            WRITE( MESG,94010 ) 
     &          'Caller Table Overflow.  Dimensioned', N,
     &          'Needed', I
            CALL M3WARN( 'RDSPROF', 0, 0, MESG )
            RDSPROF = -I
            RETURN

        END IF

        RDSPROF = I
                

C...........   Now sort the profiles:
        IF ( .NOT. SORTED ) CALL SORTI1( I, INDEXA, SXREFA )
       
C...........   Use sorted index and actual species index to create
C...........   properly ordered and adjusted coefficient table.

        DO  122  K = 1, I

            J = INDEXA( K )
            SXREF( J ) = SXREFA( K )

            DO  111  L = 1, NP

                V = INDX( L )

C.................  Non-VOC 
                IF( V .EQ. 0 ) THEN

                    IF    ( PLIST( L ) .EQ. 'CO' ) THEN
                        SCOEF( J,L ) = COFAC

                    ELSEIF( PLIST( L ) .EQ. 'NO' ) THEN
                        SCOEF( J,L ) = NOFAC

                    ELSEIF( PLIST( L ) .EQ. 'NO2' ) THEN
                        SCOEF( J,L ) = NO2FAC

                    ELSE

                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 
     &                         'Program does not handle species "'
     &                          // PLIST( L ) // '" correctly. ' //
     &                         'See rdsprof library program.'
                        CALL M3WARN( 'RDSPROF', 0, 0, MESG )

                    ENDIF

                ELSE

                    SCOEF( J,L ) = SCOEFA( K,V )

                ENDIF

111         CONTINUE
122     CONTINUE
        
        IF( EFLAG ) THEN
            RDSPROF = IMISS3
            RETURN

        ENDIF
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xx

93010   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )



        END

