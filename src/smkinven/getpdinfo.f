
        SUBROUTINE GETPDINFO( FDEV, TZONE, INSTEP, OUTSTEP, TYPNAM, 
     &                        FNAME, SDATE, STIME, NSTEPS, NPDVAR, 
     &                        NPDVSP, MXPDSRC, EAIDX, SPIDX, CFLAG )

C***************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine gets the vital information from the day-specific or 
C      hour-specific input files so that memory can be allocated to read in
C      the data.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 12/99 by M. Houyoux
C      09/2025 by HT UNC-IE:  Use M3UTILIO
C
C***************************************************************************
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
C***************************************************************************
        USE M3UTILIO

C.........  MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CINTGR, INTGRFLAG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, NSPDAT, EANAM, NCOMP, VAR_FORMULA,
     &                     VNAME

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: MXPDPT

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NUNIQCAS, UCASNKEP, UNIQCAS, UCASIDX, ITNAMA,
     &                      SCASIDX

        IMPLICIT NONE

C.........  INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
c       INCLUDE 'PARMS3.EXT'    !  I/O API parameters
c       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
c       INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
c       CHARACTER(2)    CRLF
c       INTEGER         GETFLINE
c       INTEGER         GETFORMT
c       INTEGER         FINDC
c       INTEGER         FIND1
c       INTEGER         INDEX1 
c       INTEGER         INDEXINT1 

c       EXTERNAL        CRLF, GETFLINE, GETFORMT, FIND1, FINDC, INDEX1, INDEXINT1
        INTEGER, EXTERNAL :: GETFLINE
        INTEGER, EXTERNAL :: GETFORMT

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN):: FDEV          ! file unit no.
        INTEGER     , INTENT (IN):: TZONE         ! output time zone
        INTEGER     , INTENT (IN):: INSTEP        ! expected data time step HHMMSS
        INTEGER     , INTENT (IN):: OUTSTEP       ! output time step HHMMSS
        CHARACTER(*), INTENT (IN):: TYPNAM        ! name of processing type
        CHARACTER(*), INTENT (IN):: FNAME         ! logical file name
        INTEGER     , INTENT(OUT):: SDATE         ! Julian start date in TZONE
        INTEGER     , INTENT(OUT):: STIME         ! start time of data in TZONE
        INTEGER     , INTENT(OUT):: NSTEPS        ! no. time steps
        INTEGER     , INTENT(OUT):: NPDVAR        ! no. pol/act variables
        INTEGER     , INTENT(OUT):: NPDVSP        ! no. pol/act/special data
        INTEGER     , INTENT(OUT):: MXPDSRC       ! max. no. srcs over all times
        INTEGER     , INTENT(OUT):: EAIDX( NIPPA )! index to EANAM
        INTEGER     , INTENT(OUT):: SPIDX( MXSPDAT )! index to SPDATNAM
        LOGICAL     , INTENT(OUT):: CFLAG         ! true: CEM data processing

C.........  Local parameters
        CHARACTER(16), PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'

C.........  Local allocatable arrays...
        INTEGER, ALLOCATABLE :: EASTAT( : )   ! true: act/pol present in data

C.........  Local arrats
        INTEGER         SPSTAT( MXSPDAT )     ! true: special data variable used

C.........  Other local variables
        INTEGER         I, J, N, V, NV, IV, L, LL          ! counters and indices
        INTEGER         CIDX, NCIDX        ! tmp data index
        INTEGER         FILFMT           ! format code of files in list
        INTEGER         INVFMT           ! inventory format code
        INTEGER         IOS              ! i/o status
        INTEGER         NLINE            ! number of lines
        INTEGER         NPPCAS           !  no. of pollutants per CAS number
        
        LOGICAL       :: DFLAG    = .FALSE.  ! true: day-specific processing

        CHARACTER(IOVLEN3) INVNAM   ! temporary pollutant name
        CHARACTER(IOVLEN3) POLNAM   ! temporary pollutant name
        CHARACTER(300)  MESG        !  message buffer

        CHARACTER(16) :: PROGNAME = 'GETPDINFO' !  program name

C***********************************************************************
C   begin body of program GETPDINFO

C.........  Allocate memory for logical status array for pol/act
        ALLOCATE( EASTAT( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EASTAT', PROGNAME )

        EASTAT = 0  ! array
        SPSTAT = 0  ! array

        MESG = 'Determining number of time steps for ' // TYPNAM //
     &         '-specific files...'
        CALL M3MSG2( MESG )

C.........  Check no of new calculated pollutants        
        CALL ENVSTR( FORMEVNM, MESG, ' ', VAR_FORMULA, IOS )
        IF( LEN_TRIM( VAR_FORMULA ) > 0 ) CALL FORMLIST

C.........  Perform case-specific settings
        SELECT CASE( TYPNAM )
        CASE( 'day' ) 
            DFLAG = .TRUE.

        CASE( 'hour' )
            DFLAG = .FALSE.

        CASE DEFAULT
            MESG = 'INTERNAL ERROR: Do not know type ' // TYPNAM 
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

C.........  Determine whether combining VOC and HAPs together or not 
        DO V = 1, NIPPA
            IF( INDEX( EANAM(V),'_NOI'  ) > 0 ) INTGRFLAG = .TRUE. 
            IF( INDEX( EANAM(V),'NONHAP') > 0 ) INTGRFLAG = .TRUE.
        END DO

C.........  Ensure that input file is a list-formatted file
        INVFMT = GETFORMT( FDEV, -1 )

        IF( INVFMT .NE. LSTFMT ) THEN
            MESG = TYPNAM// '-specific input file is not provided by '//
     &             'a list of files OR ' // CRLF() // BLANK10 // 
     &             'files in list provided could not be found.'
            CALL M3MSG2( MESG )

            MESG = 'Problem reading inventory file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Get the dates (in the output time zone) from the files, 
C           flag the pollutants of interest, and flag the special variables
C           contained in the file.
        CALL RDLOOPPD( FDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, DFLAG, 
     &                 FNAME, SDATE, STIME, NSTEPS, FILFMT, 
     &                 EASTAT, SPSTAT )

C.........  Allocate memory and initialize for the maximum number of 
C           records per time step
        ALLOCATE( MXPDPT( NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MXPDPT', PROGNAME )
        MXPDPT = 0  ! array

        MESG = 'Determining number of sources for ' // TYPNAM //
     &         '-specific files...'
        CALL M3MSG2( MESG )
    
C.........  Get the maximum number of records per time step - i.e., populate
C           MXPDSRC
        CALL RDLOOPPD( FDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, DFLAG, 
     &                 FNAME, SDATE, STIME, NSTEPS, FILFMT, 
     &                 EASTAT, SPSTAT )

C.........  Check whether processing CEM dataset or not
        NV = 0
        DO V = 1, NIPPA
            NV = NV + EASTAT( V )
        END DO 
        IF( NV == NIPPA ) CFLAG = .TRUE.

        N = 0
        DO V = 1, NIPPA

C.............  Processing CEM data
            IF( CFLAG ) THEN
                N = N + 1
                EAIDX( N ) = V
                CYCLE
            END IF

C.............  Add multiple inventory pollutant(s) with same CAS name
C               Find code corresponding to current pollutant before you add
            IF( EASTAT( V ) > 0 ) THEN
                N = N + 1
                EAIDX( N ) = V
                CIDX   = EASTAT( V )
                NPPCAS = UCASNKEP( CIDX )
                IF( NPPCAS > 1 ) THEN
                  DO J = 2, NPPCAS
                    NCIDX   = UCASIDX( CIDX ) + J - 1
                    POLNAM = ITNAMA( SCASIDX( NCIDX ) )
                    NV = INDEX1( POLNAM, NIPPA, EANAM )
                    IF( INDEXINT1( NV, NIPPA, EAIDX ) < 1 .AND. NV > 0 ) THEN
                        N = N + 1
                        EAIDX( N ) = NV
                    END IF
                  END DO
                END IF
            END IF

        END DO

C............  Add new computed pollutants
        IF( NCOMP > 0 ) THEN
            DO I = 1, NCOMP
                POLNAM = VNAME( I )
                NV = INDEX1( POLNAM, NIPPA, EANAM )
                IV = FIND1( NV, N, EAIDX )
                IF( NV > 0 .AND. IV < 1 ) THEN    ! only add if it doesn't exit in PDAY
                    N = N + 1
                    EAIDX( N ) = NV
                END IF
            END DO
        END IF

        NPDVAR = N

C.........  Create index to special data variable names for current data files
C.........  The idex serves a different purpose from EAIDX and is constructed
C           differently intentionally.
        N = 0
        DO V = 1, MXSPDAT

            IF( SPSTAT( V ) > 0 ) THEN
                N = N + 1
                SPIDX( V ) = N
            END IF

        END DO
        NSPDAT = N

        NPDVSP = NPDVAR + NSPDAT

C.........  Compute the maximum number of sources per time step
C.........  NOTE - MXPDPT is in the MODDAYHR module
        MXPDSRC = MAXVAL( MXPDPT )

C.........  If no sources matched then error
        IF ( MXPDSRC .EQ. 0 ) THEN

            MESG = 'No ' // TYPNAM //'-specific sources matched ' //
     &             'the inventory for the time period processed.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Deallocate local memory
        DEALLOCATE( EASTAT )

C.........  Deallocate global memory that is no longer needed
        DEALLOCATE( MXPDPT )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GETPDINFO


