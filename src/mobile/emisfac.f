
        PROGRAM EMISFAC

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Copied from emisfac.F in mobile module 7/99 by M. Houyoux and modified
C
C  FUTURE NEEDS:
C 
C***********************************************************************
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
C****************************************************************************

C...........   MODULES for public variables
C.........  This module is for mobile-specific data
        USE MODMOBIL

C...........   This module contains emission factor tables and related
        USE MODEMFAC

C...........   This module contains the information about the source category
        USE MODINFO

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

         IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL      CHKEMFAC
        INTEGER      FIND1
        INTEGER      SEC2TIME

        EXTERNAL     CHKEMFAC, FIND1, SEC2TIME

C...........   EXTERNAL BLOCK DATA for MOBILE model:

        EXTERNAL      BD01, BD02, BD03, BD04, BD05, BD06, BD10,
     &                BD11, BD12, BD13, BD14, BD15, BD16, BD17,
     &                BD18, BD19, BD20, BD21, BD22, BD23, BD24,
     &                BD25, BD26, BD27, BD28, BD29, BD30, BD31,
     &                BD32, BD33, BD34, BD35, BD36, BD37, BD38

C...........   Variables for common blocks used to carry over from one MOBILE5
C              call to the next
        INTEGER      IREAD
        INTEGER      NOTUSED

        REAL         CR12HC
        REAL         CR12CO

        CHARACTER*40 TNAME

        COMMON /IMPAR7/ TNAME,IREAD,NOTUSED
        COMMON /IM12HC/ CR12HC(19,20,5,2)
        COMMON /IM12CO/ CR12CO(19,20,5,2)

C...........   Commons linked with MOBILE subroutine
        REAL       TF                     ! temporary temperature
        REAL       TF_MAX                 ! temporary maximum temperature
        REAL       TF_MIN                 ! temporary minimum temperature

        COMMON / TOVRIDE / TF, TF_MIN, TF_MAX

C...........   Local parameters

        CHARACTER*50, PARAMETER :: SCCSW    = '@(#)$Id$'
        INTEGER     , PARAMETER :: PURETYPE = -1
        INTEGER     , PARAMETER :: CMBOTYPE =  1

C...........   Local allocatable arrays....

C...........   Simulation-specific PSIs and related arrays
        LOGICAL, ALLOCATABLE :: NDISTAT ( : )  ! true: non-diu was updated
        LOGICAL, ALLOCATABLE :: DIUSTAT ( : )  ! true: diurnal was updated
        LOGICAL, ALLOCATABLE :: UPDATNDI( : )  ! true: should update non-diu
        LOGICAL, ALLOCATABLE :: UPDATDIU( : )  ! true: should update diurnal

C...........   Variables for targeted update of EFs
        INTEGER                 NUPDAT      ! No. of PSIs targeted for update
        INTEGER, ALLOCATABLE :: UPDATE( : ) ! Sorted targeted PSI list

C...........   Update st

C...........   Temporary emission factor arrays
        REAL, ALLOCATABLE :: EFACT ( :,: )      ! Temporary non-diurnal EFs 
        REAL, ALLOCATABLE :: DFACT ( :,: )      ! Temporary diurnal EFs
        REAL, ALLOCATABLE :: EFSAVND( :,:,:,: ) ! Multi-scenario saved nondi EFs
        REAL, ALLOCATABLE :: EFSAVDI( :,:,:,: ) ! Multi-scenario saved diurn EFs

C...........   Unit numbers and logical names
        INTEGER       EDEV        ! ef-ref/temperature file unit
        INTEGER       LDEV        ! i/o api log unit
        INTEGER       IDEV        ! emission factor cross-reference file unit
        INTEGER       RDEV        ! emission factor data file (MOBILE inputs)
        INTEGER   ::  SDEV  = 0   ! Optional targeted update on EF file unit

        CHARACTER*16  DNAME       ! Output diurnal EF file logical name
        CHARACTER*16  NNAME       ! Output non-diurnal EF file logical name

C...........   Other local variables
        INTEGER    I, J, K, L     ! counters and indices

        INTEGER    CNTCOMBO       ! tmp count of no. PSI in combo PSI
        INTEGER    CPSI           ! tmp contributing PSI
        INTEGER    FDATE          ! date of read/write for EF file(s)
        INTEGER    IERR           ! MOBILE5 error value
        INTEGER    IOS            ! i/o status
        INTEGER    IPM            ! index to PSIALL
        INTEGER    IPTIM          ! PSI converted to time for file operations
        INTEGER    ITEMP          ! indicator of tmpr status for MOBILE5 
        INTEGER    JYEAR          ! year for computing emission factors
        INTEGER    LINECNT        ! position of MOBILE packet in RDEV
c        INTEGER    MXEMIS         ! sum of no. diurnal + no. non-diur EF types
        INTEGER    MXPSI          ! max no. PSIs after expansion
        INTEGER    NNDIP1         ! no. non-diurnal EFs, plus 1
        INTEGER    NPSISCN        ! no. PSIs in scenario for current PSI
        INTEGER    NV             ! local no. vehicle types
        INTEGER    PSI            ! tmp parameter scheme index
        INTEGER    PSIPNTR        ! pure: no. in scenario; combo: pntr to table
        INTEGER    TI             ! index to VLDTMPR (in MODMET)
        INTEGER    TMMI           ! index to VLDTMIN (in MODMET)

        LOGICAL    CFLAG_DIU              ! true: diurnal EFs were calc'd any T
        LOGICAL    CFLAG_NDI              ! true: non-diu EFs were calc'd any T
        LOGICAL    DIREUSE                ! true: reusing diurnal
        LOGICAL :: EFLAG = .FALSE.        ! true: error found
        LOGICAL    ENDLOOP                ! true: end tmpr loop for current PSI 
        LOGICAL    NDREUSE                ! true: reusing non-diurnal
        LOGICAL    NFLAG                  ! null flag (not used)
        LOGICAL    RFLAG                  ! true: more MOBILE5 changed then tmpr
        LOGICAL    TFLAG_DIU              ! true: diurnal EFs were calc'd this T
        LOGICAL    TFLAG_NDI              ! true: non-diu EFs were calc'd this T

        CHARACTER*20           MODELNAM   ! emission factor model name from E.V.
        CHARACTER*300          MESG       ! message buffer
        CHARACTER(LEN=IOVLEN3) ACTNAM     ! tmp activity name from E.V.


        CHARACTER*16 :: PROGNAME = 'EMISFAC'   ! program name

C***********************************************************************
C   begin body of program EMISFAC

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Obtain settings from the environment...
C.........  Get the name of the activity to use for one run
        MESG = 'Emission factor activity'
        CALL ENVSTR( 'SMK_EF_ACTIVITY', MESG, 'VMT', ACTNAM, IOS )

C.........  Get the name of the emission factor model to use for one run
        MESG = 'Emission factor model'
        CALL ENVSTR( 'SMK_EF_MODEL', MESG, 'MOBILE5', MODELNAM, IOS )

C.........  End program if source category is not mobile sources
        IF( CATEGORY .NE. 'MOBILE' ) THEN
            L = LEN_TRIM( PROGNAME )
            MESG = 'Program ' // PROGNAME( 1:L ) // ' does not ' //
     &             'support ' // CATEGORY( 1:CATLEN ) // ' sources.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Set number and name of activity
        ELSE

            NIACT = 1
            ALLOCATE( ACTVTY( NIACT ),STAT=IOS )
            CALL CHECKMEM( IOS, 'ACTVTY', PROGNAME )
            ACTVTY( 1 ) = ACTNAM

        END IF

C.........  Get file names and units; open input and output files
        CALL EMFACIO( MODELNAM, NNAME, DNAME, JYEAR, EDEV, IDEV, RDEV,
     &                SDEV, NDREUSE, DIREUSE )
        NV = NVTYPE

C.........  Pre-process emission factors inputs file.  This involves identifying
C           the model(s) being used, flaging the types, flagging combination
C           emission factors, and allocating and storing various arrays.
        CALL RDEFDAT( RDEV )

C.........  Read the emission-factor index for the purpose of getting a 
C           list of the PSIs that are used in the EF X-ref
C.........  Creates NPSI and PSILIST
        CALL RDEFTMPR( EDEV, .FALSE. )

C.........  Allocate memory for all-PSI list. Make sure that it will be
C           large enough by allocating for all in EF x-ref and in EF data files
        MXPSI = NPSI( 1 ) + NPSIDAT
        ALLOCATE( PSIALL( MXPSI ),STAT=IOS )
        CALL CHECKMEM( IOS, 'PSIALL', PROGNAME )

C.........  Expand unique PSI list to include pure EFs that are not in the 
C           emission factor cross-reference, but that are used in combination
C           factors that are in the cross-reference.
C.........  Also expand list to include the first PSI in a group, if it is not
C           already in the list.
        CALL EXPNDPSI( MXPSI, NPSI( 1 ), PSILIST( 1,1 ) )

C.........  Read the emission-factor index and min/max temperature combinations
C           from the emission factor preprocessing program
        CALL RDEFTMPR( EDEV, .TRUE. )

C.........  Expand the records of the emission-factor index and min/max 
C           temperature combinations so that pure PSIs contain the temperature
C           indices from the combination factors that use them.
        CALL EXPNDEFT

C.........  Process targeted EF update file if it exists... 
        IF( SDEV .GT. 0 ) THEN

C.............  Allocate memory for targeted update array
            ALLOCATE( UPDATE( NPSIALL ),STAT=IOS )
            CALL CHECKMEM( IOS, 'UPDATE', PROGNAME )
            
            CALL RDEFUPD( SDEV, NPSIALL, PSIALL, NPSIALL, 
     &                    NUPDAT, UPDATE )

        END IF

C.........  Allocate memory for emission factors (from MODEMFAC) and emission
C           factor temporary arrays
c        MXEMIS = NNDI + NDIU
        ALLOCATE( EFACNDI( NTMPR, NV, NNDI ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFACNDI', PROGNAME )
        ALLOCATE( EFACDIU( NVLDTMM, NV, NDIU ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFACDIU', PROGNAME )
        ALLOCATE( EFACT( NV, NNDI ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFACT', PROGNAME )
        ALLOCATE( DFACT( NV, NDIU ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DFACT', PROGNAME )
        ALLOCATE( EFSAVND( NV, MXPPGRP, NTMPR, NNDI ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFSAVND', PROGNAME )
        ALLOCATE( EFSAVDI( NV, MXPPGRP, NVLDTMM, NDIU ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFSAVDI', PROGNAME )

C.........  Allocate memory for status arrays by PSI
        ALLOCATE( UPDATNDI( NPSIALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UPDATNDI', PROGNAME )
        ALLOCATE( UPDATDIU( NPSIALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UPDATDIU', PROGNAME )        
        ALLOCATE( NDISTAT( NPSIALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NDISTAT', PROGNAME )
        ALLOCATE( DIUSTAT( NPSIALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DIUSTAT', PROGNAME )        

C.........  Create update status arrays for diurnal and non-diurnal emission
C           factors.  These depend on the type of PSI, the contents of the
C           UPDATE array, and the PSIs for which there are min/max temperatures
C           (i.e., PSIs that are in the domain).
        CALL SETEFUPD( NUPDAT, NDREUSE, DIREUSE, NNAME, DNAME,
     &                 UPDATE, UPDATNDI, UPDATDIU )

C.........  Initialize constants for loops
        FDATE = JYEAR * 1000 + 1
 
C.........  Initialize array for determining whether combo EFs need to be
C           updated based on changes to pure EFs
        NDISTAT = .FALSE.  ! array 
        DIUSTAT = .FALSE.  ! array

C...................................................................
C.........  Create new pure EFs where needed
C...................................................................

C.........  Loop through parameter scheme indices (PSIs) unique list created 
C           from emission factor cross-reference file
        NNDIP1 = NNDI + 1
        DO IPM = 1, NPSIALL

C.............  Initializations for the current PSI
            CFLAG_NDI = .FALSE.  ! Non-diurnal not calculated
            CFLAG_DIU = .FALSE.  ! Diurnal not calculated

C.............  Get PSI and seconds-count equivalent
            PSI   = PSIALL  ( IPM )
            IPTIM = SEC2TIME( PSI )

C.............  Search for PSI in list from the emission factor inputs file
            J = FIND1( PSI, NPSIDAT, PSIDAT )

            IF( J .LE. 0 ) 
     &          CALL REPORT_MISS_PSI( PSI, 'EF x-ref', 'EF settings' )

C.............  Retrieve count of PSIs that current PSI depends on, index to 
C               combo-only PSI list, and no. PSIs per scenario
            K        = PDATINDX( J )
            CNTCOMBO = PDATTYPE( K )
            PSIPNTR  = PDATPNTR( K )
            LINECNT  = PDATLINE( K )
            NPSISCN  = PDATMCNT( K )

C.............  Skip over any PSIs that are for combination emission factors
            IF( CNTCOMBO .GT. 0 ) CYCLE

C.............  Write message to screen if this is the first PSI in a group
            IF( PSIPNTR .EQ. 1 ) THEN

                WRITE( MESG, 94010 ) 'Processing PSI', PSI

                IF( NPSISCN .GT. 1 ) THEN
                    L = LEN_TRIM( MESG ) 
                    WRITE( MESG, 94010 ) MESG( 1:L ) // ' through', 
     &                                   PSI + NPSISCN - 1
                END IF

                L = LEN_TRIM( MESG ) 
                MESG = MESG( 1:L ) // ' ...'

                CALL M3MSG2( MESG )

            END IF

C.............  Check to see any diurnals for this profile indice exist
C.............  If not, initialize diurnal values to missing because they 
C               are not defined for all vehicle types
C.............  Note that NFLAG is a null flag simply used to call routine
            NFLAG = CHKEMFAC( 'DIURNAL', DNAME, FDATE, IPTIM, 
     &                        DIREUSE, .TRUE. )

C.............  Specify that more has changed than just temperature
C.............  This should not be set inside non-diurnal and diurnal
C               sections individually, because ITEMP = 1 will reset TCNT
C               in MOBILE subroutine.
            RFLAG = .TRUE.

C.............  Loop over temperatures for non-diurnal and diurnal emission
C               factors.
C.............  The temperature-handling routine will first provide temperatures
C               for the non-diurnal EFs, with diurnal min/maxs in special
C               cases.  Then the remaining min/max temperatures will be set
C               for computing the remaining diurnal emission factors
            DO

                CALL SETEFTMPR( IPM, ENDLOOP, TI, TMMI, TF, 
     &                          TF_MIN, TF_MAX )

C.................  End the loop if temperatures have been exhausted
                IF( ENDLOOP ) EXIT

C.................  Determine the status for the current temperature. This will
C                   determine whether or not to save the emisison factors
C                   in the output array. When TI = 0, do not store non-diurnal
C                   emission factors.  When TMMI = 0, do not store diurnal ones
                TFLAG_NDI = .FALSE.
                TFLAG_DIU = .FALSE.
                IF( UPDATNDI( IPM ) .AND. TI   .GT. 0 ) TFLAG_NDI=.TRUE.
                IF( UPDATDIU( IPM ) .AND. TMMI .GT. 0 ) TFLAG_DIU=.TRUE.

C.................  Force PSIs to be updated in output file if they have been
C                   computed for any temperatures
                CFLAG_NDI = ( CFLAG_NDI .OR. TFLAG_NDI )
                CFLAG_DIU = ( CFLAG_DIU .OR. TFLAG_DIU )

C.................  Skip iteration if no emission factors need updating
                IF( .NOT. CFLAG_NDI .AND. .NOT. CFLAG_DIU ) CYCLE

C.................  Initialize MOBILE5 error count 
                IERR = 0 

C.................  For MOBILE* input, when the PSI is the first in a
C                   multi-scenario run, skip the appropriate number of lines
C                   in the emission factor data file (RDEV).  Then, create
C                   emission factors and store extras.
C.................  Note that combination types have already been screened out
                IF( PSIPNTR .EQ. 1 ) THEN
                            
                    REWIND( RDEV )                  
                    CALL SKIPL ( RDEV, LINECNT )

                    CALL MOBILE( IERR , TI, TMMI, RDEV, NPSISCN, RFLAG,
     &                           EFACT, DFACT, EFSAVND, EFSAVDI  )
                    RFLAG = .FALSE.

C.................  For MOBILE5* input when the PSI is NOT the first in a
C                   multi-scenario run, retrieve EFs from storage
C.................  TFLAG_* variables screen for TI = 0 or TMMI = 0
                ELSE IF( PSIPNTR .GT. 1 ) THEN

                    IF( TFLAG_NDI ) THEN
                        EFACT( 1:NV,1:NNDI ) = 
     &                         EFSAVND( 1:NV,PSIPNTR,TI,1:NNDI )
                    END IF

                    IF( TFLAG_DIU ) THEN
                        DFACT( 1:NV,1:NDIU ) = 
     &                         EFSAVDI( 1:NV,PSIPNTR,TMMI,1:NDIU )
                    END IF

                END IF

C.................  Put non-diurnal factors from emission factor routine into
C                   non-diurnal arrays for output
C.................  Set EFACNDI to indicate that non-diurnal emission factors
C                   were calculated
                IF( TFLAG_NDI ) THEN

                    EFACNDI( TI, 1:NV, 1:NNDI ) = EFACT( 1:NV, 1:NNDI )

                END IF

C.................  Put diurnal factors from emission factor routine into 
C                   individual diurnal arrays for output
C.................  Set RFLAG to indicate that diurnal emission factors
C                   were calculated
                IF( TFLAG_DIU ) THEN
 
                    EFACDIU( TMMI, 1:NV, 1:NDIU )= DFACT( 1:NV, 1:NDIU )

                END IF

C.................  Update status of whether current PSI has been recalculated
C                   for any temperatures
                NDISTAT( IPM ) = ( NDISTAT( IPM ) .OR. CFLAG_NDI )
                DIUSTAT( IPM ) = ( DIUSTAT( IPM ) .OR. CFLAG_DIU )

            END DO  ! end loop on temperatures for computing emission factors

C.............  Write out emission factors if they were updated, or write 
C               message saying that no update was needed.
            CALL WREMFACS( NNAME, DNAME, FDATE, FDATE, IPTIM, 
     &                     NDISTAT( IPM ), DIUSTAT( IPM ) )

        END DO           ! End loop over PSIs

C...................................................................
C.........  Create new mixed EFs where needed
C...................................................................

        DO IPM = 1, NPSIALL

C.............  Get indice from MPLIST
            PSI   = PSIALL  ( IPM )
            IPTIM = SEC2TIME( PSI )

C.............  Search for PSI in emission factors data PSIs list
            J = FIND1( PSI, NPSIDAT, PSIDAT )

C.............  Get information from emission factors data table
            K        = PDATINDX( J )
            CNTCOMBO = PDATTYPE( K )
            PSIPNTR  = PDATPNTR( K )

C.............  If this PSI is not a combination PSI, go to next iteration
            IF( CNTCOMBO .EQ. 0 ) CYCLE

C.............  Write message to screen
            WRITE( MESG, 94010 ) 'Processing PSI', PSI, '...'
            CALL M3MSG2( MESG )

C.............  Determine if creation/update of combination factor is needed...
C.............  Loop through contributing PSIs for current PSI
            CFLAG_NDI = .FALSE.
            CFLAG_DIU = .FALSE.
            DO I = 1, CNTCOMBO

                CPSI = CMBOPSI( I,PSIPNTR )

C.................  Find contributing PSI in total list of PSIs
                K = FIND1( CPSI, NPSIALL, PSIALL )

C.................  Consider if current PSI is forced update or if contributing
C                   PSI was updating during this run
                IF( UPDATNDI(IPM) .OR. NDISTAT(K) ) CFLAG_NDI = .TRUE.
                IF( UPDATDIU(IPM) .OR. DIUSTAT(K) ) CFLAG_DIU = .TRUE.

C.................  Leave loop if this PSI should be updated
                IF( CFLAG_NDI .OR. CFLAG_DIU ) EXIT

            END DO

C.............  If either non-diurnal or diurnal are required, update both using 
C               EFCOMBO, because running it is fast
            IF( CFLAG_NDI .OR. CFLAG_DIU ) THEN

C.................  Create combo profile through common blocks EFS1 and EFS2
                CALL EFCOMBO( NNAME, DNAME, CNTCOMBO, FDATE, FDATE,
     &                        CMBOPSI(1,PSIPNTR), CMBOFAC(1,PSIPNTR) )

            END IF

C.............  Update status of whether current PSI has been recalculated
            NDISTAT( IPM ) = ( NDISTAT( IPM ) .OR. CFLAG_NDI )
            DIUSTAT( IPM ) = ( DIUSTAT( IPM ) .OR. CFLAG_DIU )

C.............  Store emission factors in master lists
            CALL WREMFACS( NNAME, DNAME, FDATE, FDATE, IPTIM, 
     &                     NDISTAT( IPM ), DIUSTAT( IPM ) )

        END DO           ! End loop on parameter schemes

C.........  Normal program completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )
 
C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C..............  This internal subprogram reports missing PSIs

            SUBROUTINE REPORT_MISS_PSI( PSI, FILDSC1, FILDSC2 )

            INTEGER     , INTENT (IN):: PSI
            CHARACTER(*), INTENT (IN):: FILDSC1
            CHARACTER(*), INTENT (IN):: FILDSC2

C......................................................................

            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 
     &             'Parameter scheme index', PSI, 'from ' // FILDSC1 // 
     &             ' file not found in ' // FILDSC2 // ' file.'
            CALL M3MSG2( MESG )

            RETURN

C--------------------  FORMAT  STATEMENTS   ----------------------------

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE REPORT_MISS_PSI

        END PROGRAM EMISFAC

