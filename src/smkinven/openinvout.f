
        SUBROUTINE OPENINVOUT( MXIDAT, INVDNAM, ENAME, ANAME, SDEV )

C*************************************************************************
C  subroutine body starts at line 118
C
C  DESCRIPTION:
C      This subroutine sets up the header and variables for the I/O API 
C      inventory file, and opens the I/O API and ASCII files for the SMOKE
C      area, mobile, or point source inventory.
C
C  PRECONDITIONS REQUIRED:
C      Correct number of pollutants and names EINAM are set
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines, BLDENAMS
C      Functions: I/O API functions, VERCHAR
C
C  REVISION  HISTORY:
C      Created 4/99 by M. Houyoux
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
C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptionsNRAWIN
        CHARACTER*2     CRLF
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE
        CHARACTER*16    VERCHAR

        EXTERNAL CRLF, INDEX1, PROMPTFFILE, PROMPTMFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT(IN)  :: MXIDAT ! max. inventory data names
        CHARACTER(*), INTENT(IN)  :: INVDNAM( MXIDAT ) ! inv. data names
        CHARACTER(*), INTENT(OUT) :: ENAME  ! emis i/o api inven logical name
        CHARACTER(*), INTENT(OUT) :: ANAME  ! emis ASCII inven logical name
        INTEGER     , INTENT(OUT) :: SDEV   ! ascii output inven file unit no.

C...........   LOCAL PARAMETERS
        CHARACTER*16, PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'
        CHARACTER*50, PARAMETER :: SCCSW  = '@(#)$Id$'  ! SCCS string with vers no.

C...........   Names, Units, types, & descriptions for pollutant-specific 
C              output variables

        CHARACTER(LEN=IOVLEN3) EONAMES( NIPOL+NIACT,NPPOL ) ! Names 
        INTEGER                EOTYPES( NIPOL+NIACT,NPPOL ) ! Types (Real|Int)
        CHARACTER(LEN=IOULEN3) EOUNITS( NIPOL+NIACT,NPPOL ) ! Units  
        CHARACTER(LEN=IODLEN3) EODESCS( NIPOL+NIACT,NPPOL ) ! Dscriptions  

C...........   Other local variables

        INTEGER       I, J, K, L, L1, L2, V     ! counter and indices

        INTEGER       IOS       ! i/o status
        INTEGER       LEQU      ! position of '=' in formula
        INTEGER       LDIV      ! position of '-' or '+' in formula
        INTEGER       LMNS      ! position of '-' in formula
        INTEGER       LPLS      ! position of '+' in formula
        INTEGER       NIOVARS   ! Number of I/O API file non-emis variables
        INTEGER       NDATMAX   ! Max no of pols+activitys, based on I/O API
        INTEGER       NNPVAR    ! No. non-pollutant inventory variables

        LOGICAL    :: CHKPLUS  = .FALSE. ! true: formula uses a + sign
        LOGICAL    :: CHKMINUS = .FALSE. ! true: formula uses a - sign

        CHARACTER*60  VAR_FORMULA
        CHARACTER*300 MESG      ! message buffer 

        CHARACTER(LEN=IOVLEN3 ) VIN_A
        CHARACTER(LEN=IOVLEN3 ) VNAME
        CHARACTER(LEN=IOULEN3)  UNITS  ! tmp units name

        CHARACTER*16 :: PROGNAME = 'OPENINVOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENINVOUT

C.........  Get environment variable settings
        CALL ENVSTR( FORMEVNM, MESG, ' ', VAR_FORMULA, IOS )

C.........  Get output inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Depending on source category, set number of non-pollutant 
C           inventory variables
        SELECT CASE( CATEGORY )
        CASE( 'AREA' )
            NNPVAR = NARVAR3
        CASE( 'MOBILE' )
            NNPVAR = NMBVAR3
        CASE( 'POINT' )
            NNPVAR = NPTVAR3
        END SELECT

C.........  Check number of output variables against I/O API maximum

        NIOVARS = NNPVAR + NPPOL * ( NIPOL + NIACT )
        NDATMAX = INT( ( MXVARS3 - NNPVAR ) / NPPOL )

C.........  If there are too many output variables, reset NIPOL

        IF( NIOVARS .GT. MXVARS3 ) THEN

            WRITE( MESG,94010 ) 
     &             'WARNING: Maximum number of pollutants or ' //
     &             'activities that can be be'// CRLF()// BLANK10//
     &             'written to the I/O API file is', NDATMAX, 
     &             '. This limitation is caused by'// CRLF()// BLANK10//
     &             'the I/O API variable limit of', MXVARS3,'.'
            CALL M3MSG2( MESG )
 
            WRITE( MESG,94010 ) 
     &             'WARNING: Reseting number of output pollutants '//
     &             'and activities to',
     &             NDATMAX
            CALL M3MSG2( MESG )

C.............  If the number of pollutants alone are too much, then reset
C               that number.
            IF( NIPOL .GT. NDATMAX ) THEN
                NIPOL   = NDATMAX
                NIACT   = 0
                NIOVARS = NNPVAR + NPPOL * NIPOL                

C.............  Otherwise, reset the number of activities by subtracting the
C               number of pollutants from the maximum allowed number of 
C               variables
            ELSE
                NIACT = NDATMAX - NIPOL
                NIOVARS = NNPVAR + NPPOL * ( NIPOL + NIACT )

            END IF

        ENDIF

C.........  Determine the base year of the inventory
        CALL GETBASYR( NSRC, BYEAR )

C.........  Set up for opening I/O API output file header

        CALL HDRMISS3  ! Initialize for emissions 

        NVARS3D = NIOVARS
        NROWS3D = NSRC   !  number of rows = # of sources.

        FDESC3D( 1 ) = CATDESC // ' source inventory'
        FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( SCCSW )
        WRITE( FDESC3D( 4 ),94010 ) '/NON POLLUTANT/ ', NNPVAR

        IF( NIPOL .GT. 0 ) THEN
            WRITE( FDESC3D( 5 ),94010 ) '/POLLUTANTS/', NIPOL
            WRITE( FDESC3D( 6 ),94010 ) '/PER POLLUTANT/ ', NPPOL    
        END IF

        IF( NIACT .GT. 0 ) THEN
            WRITE( FDESC3D( 7 ),94010 ) '/ACTIVITIES/', NIACT
            WRITE( FDESC3D( 8 ),94010 ) '/PER ACTIVITY/ ', NPPOL
        END IF

        WRITE( FDESC3D( 9  ),94010 ) '/NUMBER CHARS/ ' , NCHARS 
        WRITE( FDESC3D( 10 ),94010 ) '/SCC POSITION/ ' , JSCC 
        WRITE( FDESC3D( 11 ),94010 ) '/STACK POSITION/ ' , JSTACK 
        WRITE( FDESC3D( 12 ),94010 ) '/BASE YEAR/ '    , BYEAR 

c note: now that FDESC for the inven file goes passed 10 fields, must check
C n:    other programs to avoid conflicts with FDESC

C.........  Define source characteristic variables that are not strings

        J = 1
        VNAME3D( J ) = 'IFIP'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'State and county FIPS code'
        J = J + 1

        VNAME3D( J ) = 'TZONES'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Time zone for site'
        J = J + 1

        VNAME3D( J ) = 'TPFLAG'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'T|2? T|3?'
        VDESC3D( J ) = 'Use week(2), month(3) temporal profiles or not'
        J = J + 1

        VNAME3D( J ) = 'INVYR'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'year AD'
        VDESC3D( J ) = 'Year of inventory for this record'
        J = J + 1

        IF( CATEGORY .EQ. 'MOBILE' ) THEN

            VNAME3D( J ) = 'IRCLAS'
            VTYPE3D( J ) = M3INT
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'Roadway type'
            J = J + 1

            VNAME3D( J ) = 'IVTYPE'
            VTYPE3D( J ) = M3INT
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'Vehicle type code'
            J = J + 1

            VNAME3D( J ) = 'XLOC1'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'degrees'
            VDESC3D( J ) = 'Longitude at beginning of link'
            J = J + 1

            VNAME3D( J ) = 'YLOC1'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'degrees'
            VDESC3D( J ) = 'Latitude at beginning of link'
            J = J + 1

            VNAME3D( J ) = 'XLOC2'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'degrees'
            VDESC3D( J ) = 'Longitude at end of link'
            J = J + 1

            VNAME3D( J ) = 'YLOC2'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'degrees'
            VDESC3D( J ) = 'Latitude at end of link'
            J = J + 1

            VNAME3D( J ) = 'SPEED'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'miles/hr'
            VDESC3D( J ) = 'Average speed of vehicles'
            J = J + 1

        ELSE IF ( CATEGORY .EQ. 'POINT' ) THEN

            VNAME3D( J ) = 'ISIC'
            VTYPE3D( J ) = M3INT
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'Source Industrial Code'
            J = J + 1

            VNAME3D( J ) = 'XLOCA'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'degrees'
            VDESC3D( J ) = 'longitude'
            J = J + 1

            VNAME3D( J ) = 'YLOCA'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'degrees'
            VDESC3D( J ) = 'latitude'
            J = J + 1

            VNAME3D( J ) = 'STKHT'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'm'
            VDESC3D( J ) = 'Stack height'
            J = J + 1

            VNAME3D( J ) = 'STKDM'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'm'
            VDESC3D( J ) = 'Stack diameter'
            J = J + 1

            VNAME3D( J ) = 'STKTK'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'deg K'
            VDESC3D( J ) = 'Stack exhaust temperature'
            J = J + 1

            VNAME3D( J ) = 'STKVE'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'm/s'
            VDESC3D( J ) = 'Stack exhaust velocity'
            J = J + 1

        END IF

C.........  Get names, units, etc. of output pollutant-specific and activity-
C           specific records
        CALL BLDENAMS( CATEGORY, NIPPA, NPPOL, EANAM, 
     &                 EONAMES, EOUNITS, EOTYPES, EODESCS )

        DO V = 1 , NIPPA
            
            DO I = 1, NPPOL ! Loop through number of variables per pollutant

C.................  Reset units from inventory input file, if they are available
                IF( ALLOCATED( INVDUNT ) ) THEN
                    K = INDEX1( EONAMES( V, 1 ), MXIDAT, INVDNAM )
                    UNITS = INVDUNT( K, I )
                ELSE
                    UNITS = EOUNITS( V, I )
                END IF

C.................  Store variable names and information
                VNAME3D( J ) = EONAMES( V, I )
                VTYPE3D( J ) = EOTYPES( V, I )
                UNITS3D( J ) = UNITS
                VDESC3D( J ) = EODESCS( V, I )
                J = J + 1

            END DO    !  end loop on number of variables per pollutant

        END DO        !  end loop on inventory pollutants V

C.........  If there is a computed output variable, add it to the variable list.
C           Base the other variable settings (besides the name) on the first 
C           of the two variables in the formula.
        IF( VAR_FORMULA .NE. ' ' ) THEN

C.............  Make sure not using activities - not supported 
            IF( NIACT .GT. 0 ) THEN
                MESG = 'Formula not support when using activities'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Make sure formula makes sense
            LEQU = INDEX( VAR_FORMULA, '=' )
            LPLS = INDEX( VAR_FORMULA, '+' )
            LMNS = INDEX( VAR_FORMULA, '-' )

            CHKPLUS  = ( LPLS .GT. 0 )
            CHKMINUS = ( LMNS .GT. 0 )

            LDIV = LPLS
            IF( CHKMINUS ) LDIV = LMNS

            IF( LEQU .LE. 0 .OR. 
     &        ( .NOT. CHKPLUS .AND. .NOT. CHKMINUS ) ) THEN

                L = LEN_TRIM( FORMEVNM )
                MESG = 'Could not interpret formula for extra ' //
     &                 'pollutant from environment variable ' //
     &                 CRLF() // BLANK10 // '"' // FORMEVNM( 1:L ) //
     &                 '": ' // VAR_FORMULA
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Extract formula variable names
            VNAME = ADJUSTL ( VAR_FORMULA(      1:LEQU-1 ) )
            VIN_A = ADJUSTL ( VAR_FORMULA( LEQU+1:LDIV-1 ) )

C.............  Find formula inputs in existing variable list
            V = INDEX1( VIN_A, NIPPA, EANAM )

            IF( V .LE. 0 ) THEN
                L = LEN_TRIM( VIN_A )
                MESG = 'Variable "'// VIN_A( 1:L ) // 
     &                 '" from formula was not found in inventory.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            DO I = 1, NPPOL ! Loop through number of variables per pollutant

                L = LEN_TRIM( VIN_A )
                L = INDEX( EONAMES( V, I ), VIN_A(1:L) )
        	VNAME3D( J ) = EONAMES( V, I )( 1:L-1 ) // VNAME
        	VTYPE3D( J ) = EOTYPES( V, I )
        	UNITS3D( J ) = EOUNITS( V, I )
        	VDESC3D( J ) = EODESCS( V, I )
        	J = J + 1

            END DO    !  end loop on number of variables per pollutant

C.............  Update header settings
            NVARS3D = NVARS3D + NPPOL
            WRITE( FDESC3D( 5 ),94010 ) '/POLLUTANTS/', NIPOL + 1

        END IF

C.........  Prompt for and open I/O API output file
        ENAME= PROMPTMFILE( 
     &       'Enter logical name for the I/O API INVENTORY output file',
     &       FSUNKN3, ENAME, PROGNAME )
        
C.........  Prompt for and open ASCII output file
        SDEV= PROMPTFFILE( 
     &      'Enter logical name for the ASCII INVENTORY output file',
     &      .FALSE., .TRUE., ANAME, PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE OPENINVOUT

