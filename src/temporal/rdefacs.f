
        SUBROUTINE RDEFACS( MODELNAM, NNAME, DNAME, ACTNAM )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Subroutine Rdefacs reads the emission factors for the appropriate
C     activities.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 10/99 by M. Houyoux
C
C****************************************************************************/
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

C...........   This module contains the information about the source category
        USE MODINFO

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     INDEX1
        INTEGER     SEC2TIME

        EXTERNAL    INDEX1, SEC2TIME

C.........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: MODELNAM  ! name for EF model
        CHARACTER(*), INTENT (IN) :: NNAME     ! non-diurnal EFs file name
        CHARACTER(*), INTENT (IN) :: DNAME     ! diurnal EFs file name
        CHARACTER(*), INTENT (IN) :: ACTNAM    ! activity name

C.........  Local variables
        INTEGER       I, T           ! counters and indices

        INTEGER       DDATE          ! diurnal file date
        INTEGER       IDX            ! index of ACTNAM in ACTVTY array
        INTEGER       IOS            ! i/o status
        INTEGER       IPTIM          ! PSI in units of seconds for I/O API read
        INTEGER       NDATE          ! non-diurnal file date
        INTEGER       PSI            ! tmp parameter scheme index

        LOGICAL    :: EFLAG = .FALSE.  ! true: error found

        CHARACTER*300 MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDEFACS' ! program name

C***********************************************************************
C   begin body of subroutine RDEFACS

C.........  Get header description of non-diurnal file to get date for reading
        IF( .NOT. DESC3( NNAME ) ) THEN
            MESG = 'Could not get description of file "' // NNAME // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE
            NDATE = SDATE3D
        END IF

C.........  Get header description of diurnal file to get date for reading
        IF( .NOT. DESC3( DNAME ) ) THEN
            MESG = 'Could not get description of file "' // DNAME // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE
            DDATE = SDATE3D
        END IF

C.........  Get emission factor names, units, and descriptions for the non-
C           diurnal emission factors and the model of interest.
C.........  NOTE - although the arguments to this routine are the I/O API
C           arrays for variables, etc., the subroutine sets the public arrays
C           for non-diurnal (or diurnal) names, units, and descriptions.
        CALL EFSETUP( MODELNAM, 'NONDIURNAL', MXVARS3, NNDI,
     &                VNAME3D, UNITS3D, VDESC3D )

C.........  Get emission factor names, units, and descriptions for the diurnal
C           emission factors and the model of interest
        CALL EFSETUP( MODELNAM, 'DIURNAL', MXVARS3, NDIU,
     &                VNAME3D, UNITS3D, VDESC3D )

C.........  Find index of the activity of interest
        IDX = INDEX1( ACTNAM, NIACT, ACTVTY )

C.........  Allocate memory for emission factors
        ALLOCATE( EFDIUALL( NVLDTMM,NVTYPE,NPSI(IDX),NDIU ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFACDIU', PROGNAME )
        ALLOCATE( EFNDIALL( NTMPR,NVTYPE,NPSI(IDX),NNDI ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFACNDI', PROGNAME )

        MESG = 'Reading emission factors...'
        CALL M3MSG2( MESG )

C.........  Read emissions factors for all parameter scheme indices (PSIs) 
C           that were in the emission factor cross reference file.
        DO  T = 1, NPSI( IDX )

            PSI = PSILIST( T, IDX )
            IPTIM = SEC2TIME( PSI )

C.............  Read diurnal emission factors
            DO I = 1, NDIU
                CALL READ_EMFACT( DNAME, DIUNAM(I), DDATE, NVLDTMM, 
     &                            NVTYPE, EFDIUALL( 1, 1, T, I ) )
            END DO

C.............  Read non-diurnal emission factors
            DO I = 1, NNDI
                CALL READ_EMFACT( NNAME, NDINAM(I), NDATE, NTMPR,
     &                            NVTYPE, EFNDIALL( 1, 1, T, I )  )
            END DO

        END DO  ! Loop over PSIs

C.........  Abort if an error was found
        IF( EFLAG ) THEN

            WRITE( MESG,94010 ) 'ERROR: Emission factors ' //
     &             'missing!  Rerun EMISFAC program.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C.........  Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS   ******************************

        CONTAINS

            SUBROUTINE READ_EMFACT( FNAME, VNAME, FDATE, 
     &                              DIM1, DIM2, ARRAY )

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: FNAME      ! logical file name          
            CHARACTER(*), INTENT (IN) :: VNAME      ! EF variable name
            INTEGER     , INTENT (IN) :: FDATE      ! read date
            INTEGER     , INTENT (IN) :: DIM1       ! read date
            INTEGER     , INTENT (IN) :: DIM2       ! read date
            REAL        , INTENT(OUT) :: ARRAY( DIM1, DIM2 ) ! emission factors

C.............  Local variables
            CHARACTER*300  MESG

C--------------------------------------------------------------------------

            IF ( .NOT. READ3( FNAME, VNAME, ALLAYS3, FDATE, IPTIM,
     &                        ARRAY ) ) THEN

                WRITE( MESG,94010 ) 
     &                 'Could not read ' // VNAME // ' for PSI', PSI, 
     &                 'from file "' // FNAME // '".'
                EFLAG = .TRUE.

            END IF

C-----------------  SUBPROGRAM FORMAT STATEMENTS -------------------------

C.........  Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE READ_EMFACT

        END SUBROUTINE RDEFACS
