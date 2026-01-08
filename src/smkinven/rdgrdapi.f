
        SUBROUTINE RDGRDAPI( FNAME, GRDNM )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads an I/O API NetCDF gridded inventory file.
C
C  PRECONDITIONS REQUIRED:
C      Input file logical name FNAME opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 10/2000 by M. Houyoux
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: POLVAL, TZONES, CIFIP, CELLID, TPFLAG, INVYR,
     &                      NPCNT, CSCC, IPOSCOD, CSOURC

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPOL, NIPPA, NPPOL, CATEGORY, NEM, NDY,
     &                     EIIDX, EINAM, EANAM, EAUNIT, EADESC, NSRC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
C       CHARACTER(2)    CRLF
C       LOGICAL         ENVYN
        INTEGER         GETIFDSC
C       INTEGER         INDEX1

C        EXTERNAL        CRLF, ENVYN, GETIFDSC, INDEX1
        EXTERNAL     GETIFDSC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME     ! logical name of input file
        CHARACTER(*), INTENT(OUT) :: GRDNM     ! grid name for input data

C...........   Local allocatable arrays
        REAL, ALLOCATABLE :: INDATA( : )   ! tmp gridded input data

C...........   Other local variables
        INTEGER         C, ES, J, K, L, S, V     !  counters and indices

        INTEGER         FLEN        !  file name length
        INTEGER         IDX         !  emissions array index (ann or ave day)
        INTEGER      :: INY = 0     !  tmp inventory year
        INTEGER         IOS         !  i/o status
        INTEGER      :: NCELL       !  tmp cell numbers
        INTEGER         NEDIM1      !  1st dimension for sparse emis arrays
        INTEGER         NVARS       !  no. variables in gridded input file
        INTEGER      :: TPF         !  tmp temporal adjustments setting
        INTEGER      :: TZONE = -50 !  tmp time zone
        INTEGER         WKSET       !  setting for wkly profile TPFLAG component

        LOGICAL      :: EFLAG  = .FALSE. ! true: error occured
        LOGICAL      :: DFLAG  = .FALSE. ! true: weekday (not full week) nrmlizr
        LOGICAL      :: TVFLAG = .FALSE. ! true: time zone is a variable

        CHARACTER(300)  MESG        !  message buffer

        CHARACTER(CELLEN3) CCELL    ! tmp cell ID
        CHARACTER(POLLEN3) CCOD     ! character pollutant index
        CHARACTER(FIPLEN3) CFIP     ! character FIP code
        CHARACTER(IOVLEN3) VBUF     ! tmp variable name
        CHARACTER(SCCLEN3) SCCZERO  ! default source category code

        CHARACTER(16) :: PROGNAME =  'RDGRDAPI' ! program name

C***********************************************************************
C   begin body of subroutine RDGRDAPI

        FLEN = LEN_TRIM( FNAME )

C.........  Get setting for interpreting weekly temporal profiles from the
C           environment.
        DFLAG = .FALSE.
        MESG = 'Use weekdays only to normalize weekly profiles'
        DFLAG = ENVYN( 'WKDAY_NORMALIZE', MESG, DFLAG, IOS )

C.........  Set weekly profile interpretation flag...
C.........  Weekday normalized
        IF( DFLAG ) THEN
            WKSET = WDTPFAC
            MESG = 'NOTE: Setting inventory to use weekday '//
     &             'normalizer for weekly profiles'

C.........  Full-week normalized
        ELSE
            WKSET = WTPRFAC
            MESG = 'NOTE: Setting inventory to use full-week '//
     &             'normalizer for weekly profiles'

        END IF

C.........  Write message
        CALL M3MSG2( MESG )

C.........  Set default inventory characteristics (declared in MODINFO)
        CALL INITINFO( IOGFMT )

C.........  Read the header of the gridded I/O API file
        IF( .NOT. DESC3( FNAME ) ) THEN

            MESG = 'Could not get description of file "' //
     &             FNAME( 1:LEN_TRIM( FNAME ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Give error if there are no variables
        IF( NVARS3D .LT. 1 ) THEN

            MESG = 'No variables found in gridded input file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Give warning if time zone is not set in the header and it is not
C           in the variable list
        TZONE = GETIFDSC( FDESC3D, '/TZONE/', .FALSE. )
        K = INDEX1( 'TZONES', NVARS3D, VNAME3D )
        IF( TZONE .LT. 0 .AND. 
     &      K     .LE. 0       ) THEN

            WRITE( MESG,94010 ) 'WARNING: Time zone not set in ' //
     &             'header or variable list.' // CRLF() // BLANK10 //
     &             'Assuming time zone 0 (GMT) for data.' 
            CALL M3MSG2( MESG )
            TZONE = 0

        ELSE IF( K .GT. 0 ) THEN
            TVFLAG = .TRUE.

        END IF

C.........  Give warning if layers are greater than 1
        IF( NLAYS3D .GT. 1 ) THEN

            WRITE( MESG,94010 ) 'WARNING: Only the first layer out ' //
     &             'of', NLAYS3D, 'will be imported.'
            CALL M3MSG2( MESG )

        END IF

C.........  Give warning if time steps are greater than 1
        IF( MXREC3D .GT. 1 ) THEN

            WRITE( MESG,94010 ) 'WARNING: Only the first time step ' //
     &             'out of', MXREC3D, 'will be imported.'
            CALL M3MSG2( MESG )

        END IF

C.........  Give error if year of data is not provided
        IF( SDATE3D .LT. 1900000 ) THEN

            WRITE( MESG,94010 ) 'Cannot determine year from invalid ' //
     &             'start date', SDATE3D, '.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise set the year of the input data
        ELSE
            INY = SDATE3D / 1000

        END IF

C.........  Set temporal flag based on time step in file
        IF( TSTEP3D .EQ. 87600000 ) THEN     ! Annual data

            TPF = MTPRFAC * WKSET
            IDX = NEM

        ELSE IF( TSTEP3D .EQ. 240000 ) THEN  ! Average day data

            TPF = WKSET
            IDX = NDY

C.........  Give error if time step is not one of the recognized values.
        ELSE
            MESG = 'Time step needs to be 87600000 or 240000 HHMMSS ' //
     &             'to determine temporal approach'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Set the number of cells, variables, grid name, etc.
        NCELL  = NROWS3D * NCOLS3D
        NVARS  = NVARS3D
        NIPOL  = NVARS3D
        GRDNM  = GDNAM3D

C.........  Reset variable number if there are any special variables
        IF( TVFLAG ) THEN
            NIPOL = NIPOL - 1
        END IF

C.........  Set dependent values (no. sources, pol/act, and dimension)
        NSRC   = NCELL
        NIPPA  = NIPOL
        NEDIM1 = NSRC * NIPPA

C.........  Allocate memory for variables
C.........  NOTE - Both EINAM and EANAM are created to support WRINVEMIS
        ALLOCATE( EIIDX( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EIIDX', PROGNAME )
        ALLOCATE( EINAM( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EINAM', PROGNAME )
        ALLOCATE( EANAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAM', PROGNAME )
        ALLOCATE( EAUNIT( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EAUNIT', PROGNAME )
        ALLOCATE( EADESC( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EADESC', PROGNAME )

        J = 0
        DO V = 1, NVARS

C.............  Skip special variable names and set output variable info
            IF( VNAME3D( V ) .NE. 'TZONES' ) THEN

                J = J + 1
                EIIDX ( J ) = V
                EINAM ( J ) = VNAME3D( V )
                EANAM ( J ) = VNAME3D( V )
                EAUNIT( J ) = UNITS3D( V )
                EADESC( J ) = VDESC3D( V )

            END IF

        END DO

C.........  Allocate memory for (sorted) output inventory characteristics.
C.........  The sorted arrays can be allocated right away because the only
C           source characteristic for this type of data is the grid cell, and
C           it is sorted already.
        CALL SRCMEM( CATEGORY, 'SORTED', .TRUE., .FALSE., NSRC, 
     &               NEDIM1, NPPOL )

        CALL SRCMEM( CATEGORY, 'SORTED', .TRUE., .TRUE., NSRC, 
     &               NEDIM1, NPPOL )

C.........  Initialize emissions data
        POLVAL = BADVAL3    ! array

C.........  Initialize tmp source characteristics
        CFIP    = REPEAT( '0', FIPLEN3 )
        SCCZERO = REPEAT( '0', SCCLEN3 )
        CCOD    = ' '

C.........  Loop through cells (same as sources) and store data as sources
        ES  = 0
        DO S = 1, NSRC

            CIFIP ( S ) = CFIP
            CELLID( S ) = S
            TPFLAG( S ) = TPF
            INVYR ( S ) = INY
            NPCNT ( S ) = NIPPA
            CSCC  ( S ) = SCCZERO

C.............  If time zone is not a variable in the input file already, then
C               then set it based on the header
            IF( .NOT. TVFLAG ) THEN
                TZONES( S ) = TZONE
            END IF

            WRITE( CCELL, 94130 ) S

            CALL BLDCSRC( CFIP, SCCZERO, CCELL, CHRBLNK3, 
     &                    CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                    CCOD, CSOURC( S ) )

        END DO

C.........  If time zone is a variable, store it
        IF( TVFLAG ) THEN

            IF( .NOT. READ3( FNAME, 'TZONES', 1, SDATE3D, 
     &                       STIME3D, TZONES              ) ) THEN

C.................  Write error if data could not be read
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not read "TZONES" from file ' //
     &                 FNAME( 1:FLEN ) // '".'
                CALL M3MSG2( MESG )

            END IF

        END IF

C.........  Allocate local memory for temporary gridded data file
        ALLOCATE( INDATA( NCELL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDATA', PROGNAME )

C.........  Read data from gridded file and store in appropriate data structure
C           for use by the rest of the program
        DO V = 1, NVARS

            VBUF = EANAM( V )

C.............  Skip special variable "TZONES"
            IF( VBUF .EQ. 'TZONES' ) CYCLE

            IF( .NOT. READ3( FNAME, VBUF, 1, SDATE3D, 
     &                       STIME3D, INDATA          ) ) THEN

C.................  Write error if data could not be read
                EFLAG = .TRUE.
                L = LEN_TRIM( VBUF )
                MESG = 'ERROR: Could not read "' // VBUF( 1:L ) //
     &                 '" from file "' // FNAME( 1:FLEN ) // '".'
                CALL M3MSG2( MESG )

            ELSE

                DO C = 1, NCELL

                    ES = ( C-1 ) * NIPPA + V
                    IPOSCOD( ES )     = V
                    POLVAL ( ES,IDX ) = INDATA( C )

                END DO

            END IF

        END DO

C.........  Abort if there was a reading error
        IF( EFLAG ) THEN
           MESG = 'Problem reading gridded inventory file "' // 
     &            FNAME( 1:FLEN ) // '".'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Deallocate local memory
        DEALLOCATE( INDATA )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94125   FORMAT( I5 )

94130   FORMAT( I8 )

        END SUBROUTINE RDGRDAPI
