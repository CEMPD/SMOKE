
        SUBROUTINE MRGVNAMS

C***********************************************************************
C  subroutine MRGVNAMS body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to merge the pollutant-to-species,
C      pollutant, and species names from all of the open source categories,
C      and populate the global unique lists of these names.
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        INTEGER         GETFLINE
        INTEGER         INDEX1

        EXTERNAL   CRLF, GETFLINE, INDEX1

C...........   Master list of valid inventory pollutants (needed for sorting
C              output variables)
        INTEGER                   :: MXIDAT = 0   ! max no pols/acts
        INTEGER                   :: NDAT   = 0   ! actual no pols/act in master
        INTEGER                   :: NDAT1  = 0   ! actual no pols in master
        INTEGER                   :: NIACT  = 0   ! actial no activites
        INTEGER      , ALLOCATABLE:: INVDCOD( : ) ! codes
        INTEGER      , ALLOCATABLE:: INVSTAT( : ) ! Status (<0 activity; >0 pol)
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: INVDNAM( : ) ! pollutant/act name

C...........   Temporary arrays for building sorted pol-to-species names list
        INTEGER                                MXPTOSPC     ! max for dimension
        INTEGER               , ALLOCATABLE :: INDXA  ( : ) ! sorting index
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE :: TVSORTA( : ) ! basis for sorting
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE :: TVDESCA( : ) ! pol-to-spec concat

C...........   Other local variables
        INTEGER         I, J, J1, J2, L, K, L1, L2, M, N  ! counters and indices
        INTEGER         LJ              ! string length of emis types joiner
        INTEGER         IOS             ! i/o error status
        INTEGER         JA, JM, JP      ! position in rctvty var desc of pol-spc
        INTEGER         NCNT            ! counter

        LOGICAL      :: EFLAG = .FALSE. ! error flag

        CHARACTER(LEN=IOVLEN3)   CPOL     ! tmp pol/act buffer
        CHARACTER*300            MESG     ! message buffer

        CHARACTER*16 :: PROGNAME = 'MRGVNAMS' ! program name

C***********************************************************************
C   begin body of subroutine MRGVNAMS

C.........  Get number of lines of pollutant codes/names file 
        IF( PDEV .GT. 0 ) THEN
            MXIDAT = GETFLINE( PDEV, 'Pollutant codes and names file' )
        END IF
        IF( VDEV .GT. 0 ) THEN
            I = GETFLINE( VDEV, 'Activity names file' )
            MXIDAT = MXIDAT + I
        END IF

C.........  Allocate memory for storing contents of pollutants file
        ALLOCATE( INVDCOD( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDCOD', PROGNAME )
        ALLOCATE( INVDNAM( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDNAM', PROGNAME )
        ALLOCATE( INVSTAT( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVSTAT', PROGNAME )

C.........  Initialize inventory data status
        INVSTAT = 0   ! array

C.........  Read, sort, and store pollutant codes/names file
        IF( PDEV .GT. 0 ) THEN
            CALL RDCODNAM( PDEV, MXIDAT, NDAT, INVDCOD, 
     &                     INVDNAM )
            NDAT1 = NDAT
        END IF

C.........  Read, sort, and store activity codes/names file
        IF( VDEV .GT. 0 ) THEN
            CALL RDCODNAM( VDEV, MXIDAT, NDAT, INVDCOD, 
     &                     INVDNAM )
        END IF

        MXIDAT = NDAT

C.........  Loop through pollutants and/or activities in each source category
C           and update status of pollutants and activities in master list.
        CALL SET_MASTER_POL_STAT( 'area',   ANIPOL, AEINAM, EFLAG )
        CALL SET_MASTER_POL_STAT( 'mobile', MNIPPA, MEANAM, EFLAG )
        CALL SET_MASTER_POL_STAT( 'point',  PNIPOL, PEINAM, EFLAG )

C.........  If any pollutants/activities from matrices not found in master
C           list, exit
        IF( EFLAG ) THEN
            MESG = 'Make sure that master pollutants and/or ' //
     &             'activities files are' // CRLF() // BLANK10 // 
     &             'consistent with those used to process the ' //
     &             'inventory.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Loop through master list count the actual number of total pollutants
C           and activities
        NIPOL = 0
        NIACT = 0
        DO I = 1, MXIDAT
            IF( INVSTAT( I ) .GT. 0 ) NIPOL = NIPOL + 1
            IF( INVSTAT( I ) .LT. 0 ) NIACT = NIACT + 1
        END DO
        NIPPA = NIPOL + NIACT

C.........  Allocate memory for array of sorted pollutants and activities and
C           for pollutants only, and for the input variable names and units
C.........  NOTE - the input variable names could be different from EANAM if
C           user is using emissions from the inventory file, and has selected
C           ozone season emissions instead of annual emissions.
        ALLOCATE( EANAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAM', PROGNAME )
        ALLOCATE( EINAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EINAM', PROGNAME )
        ALLOCATE( TONAMES( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TONAMES', PROGNAME )
        ALLOCATE( TOUNITS( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TOUNITS', PROGNAME )
        ALLOCATE( SPCUNIT( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCUNIT', PROGNAME )

C.........  Initialize all
        EANAM   = ' '   ! array
        EINAM   = ' '   ! array
        TONAMES = ' '   ! array
        TOUNITS = ' '   ! array
        SPCUNIT = ' '   ! array

C.........  Create array of sorted unique pollutants and activities
C.........  Also determine number of pollutants and store pollutants-only array
        J1 = 0
        J2 = 0
        LJ = LEN_TRIM( ETJOIN )
        DO I = 1, MXIDAT

            IF( INVSTAT( I ) .NE. 0 ) THEN
                J1 = J1 + 1
                EANAM( J1 ) = INVDNAM( I )

                K  = INDEX( EANAM( J1 ), ETJOIN )
                L2 = LEN_TRIM( EANAM( J1 ) )
                IF( K .GT. 0 ) THEN
                    CPOL = EANAM( J1 )( K+LJ:L2 )
                ELSE
                    CPOL = EANAM( J1 )
                END IF
            END IF

            IF( INVSTAT( I ) .GT. 0 ) THEN
                J2 = J2 + 1
                EINAM( J2 ) = CPOL
            END IF

        END DO

C.........  Create array of input variable names and units that combines
C           this information from anthropogenic source categories.
        CALL BUILD_NAMES_UNITS( ANIPOL, AEINAM, AONAMES, AOUNITS, 
     &                          TONAMES, TOUNITS )
        CALL BUILD_NAMES_UNITS( MNIPPA, MEANAM, MONAMES, MOUNITS, 
     &                          TONAMES, TOUNITS )
        CALL BUILD_NAMES_UNITS( PNIPOL, PEINAM, PONAMES, POUNITS, 
     &                          TONAMES, TOUNITS )

C.........  Create array of sorted unique pol-to-species, sorted in order of
C           pollutants, and then in alphabetical order by species...
        
C.........  Allocate maximum possible memory needed by summing number of 
C           variables from all source categories

        MXPTOSPC = ANSMATV + ARNMSPC + BNSMATV + MNSMATV + MRNMSPC +
     &             PRNMSPC + PNSMATV

        ALLOCATE( INDXA( MXPTOSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXA', PROGNAME )
        ALLOCATE( TVSORTA( MXPTOSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TVSORTA', PROGNAME )
        ALLOCATE( TVDESCA( MXPTOSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TVDESCA', PROGNAME )

C.........  Loop through source-specific variable descriptions, find position of
C           of pollutant for output, and store concatenated position with
C           species name.  Make sure the same name is not stored twice.
        NCNT = 0
        JA   = ANRMATV - ARNMSPC + 1
        JM   = MNRMATV - MRNMSPC + 1
        JP   = PNRMATV - PRNMSPC + 1
        IF( ARFLAG ) 
     &    CALL BUILD_VDESC_UNSORT( NCNT, ARNMSPC, ARVDESC(JA) ) ! area rctvty
        IF( MRFLAG ) 
     &    CALL BUILD_VDESC_UNSORT( NCNT, MRNMSPC, MRVDESC(JM) ) ! mobile rctvty
        IF( PRFLAG ) 
     &    CALL BUILD_VDESC_UNSORT( NCNT, PRNMSPC, PRVDESC(JP) ) ! point rctvty
        IF( AFLAG ) 
     &    CALL BUILD_VDESC_UNSORT( NCNT, ANSMATV, ASVDESC     ) ! area spectn
        IF( BFLAG ) 
     &    CALL BUILD_VDESC_UNSORT( NCNT, BNSMATV, BSVDESC     ) ! biogenics
        IF( MFLAG ) 
     &    CALL BUILD_VDESC_UNSORT( NCNT, MNSMATV, MSVDESC     ) ! mobile spectn
        IF( PFLAG ) 
     &    CALL BUILD_VDESC_UNSORT( NCNT, PNSMATV, PSVDESC     ) ! point spectn

        NSMATV = NCNT 

C.........  Make sure some species match some of the pollutants
        IF( SFLAG .AND. NSMATV .EQ. 0 ) THEN
            MESG = 'ERROR: No speciation factors match the inventory!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Sort cumulative list
        CALL SORTIC( NSMATV, INDXA, TVSORTA )

C.........  Allocate memory for sorted list
        ALLOCATE( TSVDESC( NSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TSVDESC', PROGNAME )

C.........  Store sorted list
        DO I = 1, NSMATV
            J = INDXA( I )
            TSVDESC( I ) = TVDESCA( J )
        END DO

C.........  Create array of sorted unique species, sorted in order of their
C           associated pollutants, and then in alphabetical order by species...
C
C.........  Allocate memory with the number of variables in the speciation
C           matrices, which will always be >= needed space
C.........  Also allocate memory for the index between the master species names
C           and the master pollutant names
        ALLOCATE( AEMNAM( ANSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AEMNAM', PROGNAME )
        ALLOCATE( BEMNAM( BNSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BEMNAM', PROGNAME )
        ALLOCATE( MEMNAM( MNSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MEMNAM', PROGNAME )
        ALLOCATE( PEMNAM( PNSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PEMNAM', PROGNAME )
        ALLOCATE( EMNAM( NSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMNAM', PROGNAME )
        ALLOCATE( EMIDX( NSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMIDX', PROGNAME )

C.........  Call subprogram to store species names in appropriate order
C           for all source categories and the total.
        IF( AFLAG ) 
     &    CALL BUILD_SPECIES_ARRAY( ANSMATV, ASVDESC, ANMSPC, AEMNAM )
        IF( BFLAG ) 
     &    CALL BUILD_SPECIES_ARRAY( BNSMATV, BSVDESC, BNMSPC, BEMNAM )
        IF( MFLAG ) 
     &    CALL BUILD_SPECIES_ARRAY( MNSMATV, MSVDESC, MNMSPC, MEMNAM )
        IF( PFLAG ) 
     &    CALL BUILD_SPECIES_ARRAY( PNSMATV, PSVDESC, PNMSPC, PEMNAM )
        CALL BUILD_SPECIES_ARRAY(  NSMATV, TSVDESC, NMSPC, EMNAM )

C.........  Create index between master species names and master inventory names
C.........  If there are multiple pollutants per species, the last pollutant
C           will be put in the index.
C.........  Also create the speciation units per pollutant
        DO I = 1, NSMATV

C.............  Get positions of pollutants and model species
            L1 = INDEX( TSVDESC( I ), SPJOIN )
            L2 = LEN_TRIM( TSVDESC( I ) )
            J  = INDEX1( TSVDESC( I )(    1:L1-1 ), NIPPA, EANAM )
            K  = INDEX1( TSVDESC( I )( L1+1:L2   ), NMSPC, EMNAM )

C.............  Store position of pollutant for each species
            EMIDX( K ) = J

C.............  Find speciation name in one of the speciation matrices and
C               set units accordingly.  Set it based on the first one found.
            IF( SPCUNIT( J ) .EQ. ' ' ) THEN

                IF( AFLAG ) THEN
                    M = INDEX1( TSVDESC( I ), ANSMATV, ASVDESC )
                    IF( M .GT. 0 ) SPCUNIT( J ) = ASVUNIT( M )
                END IF

                IF( MFLAG ) THEN
                    M = INDEX1( TSVDESC( I ), MNSMATV, MSVDESC )
                    IF( M .GT. 0 ) SPCUNIT( J ) = MSVUNIT( M )
                END IF

                IF( PFLAG ) THEN
                    M = INDEX1( TSVDESC( I ), PNSMATV, PSVDESC )
                    IF( M .GT. 0 ) SPCUNIT( J ) = PSVUNIT( M )
                END IF

            END IF

        END DO

C.........  Deallocate temporary arrays
        DEALLOCATE( INDXA, TVSORTA, TVDESCA, INVDCOD, INVDNAM, INVSTAT )

C.........  Abort if error(s) found
        IF( EFLAG ) THEN
            MESG = 'Problem processing variable names.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS
 
C.............  This internal subprogram updates the status for the master
C               list of inventory pollutants
            SUBROUTINE SET_MASTER_POL_STAT( CATDESC, N, PNAMES, EFLAG )

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: CATDESC     ! src category desc
            INTEGER     , INTENT (IN) :: N           ! number of src-spec pols
            CHARACTER(*), INTENT (IN) :: PNAMES( N ) ! names of src-spec pols
            LOGICAL     , INTENT(OUT) :: EFLAG       ! error flag

C.............  Local subprogram variables
            INTEGER      I, J, K, L      ! counters and indices        

            CHARACTER(LEN=IOVLEN3) VBUF

C----------------------------------------------------------------------

            DO I = 1, N

                VBUF = PNAMES( I )
                J = INDEX1( VBUF, NDAT1, INVDNAM )
                K = 0
                IF( MXIDAT .NE. NDAT1)
     &              K = INDEX1( VBUF, MXIDAT-NDAT1, INVDNAM( NDAT1+1 ) )

C................. Evaluate for pollutant
                IF( J .GT. 0 ) THEN
                    INVSTAT( J ) = 1

C................. Evaluate for activity
                ELSE IF( K .GT. 0 ) THEN
                    INVSTAT( K ) = -1 

                ELSE
                    EFLAG = .TRUE.
                    L = LEN_TRIM( VBUF )
                    MESG = 'ERROR: Variable "' // VBUF( 1:L ) // 
     &                     '" from ' // CATDESC // 
     &                     ' inventory file is not in ' // CRLF() //
     &                     BLANK10 // 'master pollutants or ' //
     &                     'activites files.'
                    CALL M3MSG2( MESG )

                END IF

            END DO
 
            END SUBROUTINE SET_MASTER_POL_STAT

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram builds the sorted list of input
C               emissions/activity variable names and their units
            SUBROUTINE BUILD_NAMES_UNITS( NINVARS, REFNAMS, INNAMS, 
     &                                    INUNITS, OUTNAMS, OUTUNITS )

C.............  Subprogram arguments

            INTEGER     , INTENT (IN) :: NINVARS            ! no. input vars
            CHARACTER(*), INTENT (IN) :: REFNAMS( NINVARS ) ! reference names
            CHARACTER(*), INTENT (IN) :: INNAMS ( NINVARS ) ! input var names
            CHARACTER(*), INTENT (IN) :: INUNITS( NINVARS ) ! input units
            CHARACTER(*), INTENT(OUT) :: OUTNAMS( NIPPA   ) ! out variable names
            CHARACTER(*), INTENT(OUT) :: OUTUNITS( NIPPA  ) ! out units

C.............  Local subprogram variables
            INTEGER        K, L, V     ! counters and indices
            CHARACTER*300  MESG        ! message buffer

C----------------------------------------------------------------------

C.............  Loop through input variable names
            DO V = 1, NINVARS            

C.................  Look for reference variable name in master list
            
                K = INDEX1( REFNAMS( V ), NIPPA, EANAM )

C.................  When output names and units have already been set, check
C                   for consistency. 
                IF( OUTNAMS( K ) .NE. ' ' ) THEN

                    IF( OUTUNITS( K ) .NE. INUNITS( V ) ) THEN
                        EFLAG = .TRUE.
                        L = LEN_TRIM( INNAMS  ( V ) )
                        MESG = 'ERROR: Units for input variable "' //
     &                         INNAMS( V )( 1:L ) // '" are not ' //
     &                         'consistent among the input files.'
                        CALL M3MSG2( MESG )

                    END IF

C.................  Otherwise, set the output names
                ELSE
                    OUTNAMS ( K ) = INNAMS ( V )
                    OUTUNITS( K ) = INUNITS( V )
                END IF

            END DO

            END SUBROUTINE BUILD_NAMES_UNITS

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram builds the unsorted list of unique
C               pollutant-to-species speciation variable descriptions.
            SUBROUTINE BUILD_VDESC_UNSORT( NCNT, NVARS, VDESCS )

C.............  Subprogram arguments

            INTEGER     , INTENT(IN OUT) :: NCNT            ! running count
            INTEGER     , INTENT    (IN) :: NVARS           ! no of var descs
            CHARACTER(*), INTENT    (IN) :: VDESCS( NVARS ) ! spec var descs

C.............  Local subprogram variables
            INTEGER      I, K1, K2, L, L2     ! counters and indices
            INTEGER      LJ, LS               ! string lengths of separators

            CHARACTER*5                IBUF   ! tmp variable position to
            CHARACTER*5                PBUF   ! tmp pollutant position
            CHARACTER(LEN=IOVLEN3)     POLNAM ! tmp pollutant name
            CHARACTER(LEN=IOVLEN3)     SPCNAM ! tmp species name
            CHARACTER(LEN=IODLEN3)     TDESC  ! tmp combined pollutant & species

C----------------------------------------------------------------------

            LJ = LEN_TRIM( ETJOIN )
            LS = LEN_TRIM( SPJOIN )
            DO I = 1, NVARS

                TDESC  = VDESCS( I )
                L      = INDEX( TDESC, SPJOIN )  ! find species separator
                IF( L .LE. 0 ) THEN
                    MESG = 'ERROR: Variable descriptions do not ' //
     &                     'contain proper separator in spec matrix.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                L2     = INDEX( TDESC, ETJOIN )  ! find emission type separator

C.................  Extract pollutant name (without emissions type, if applies)
                IF( L2 .GT. 0 ) THEN
                    POLNAM = TDESC( L2+LJ:L-1 )
                ELSE
                    POLNAM = TDESC( 1:L-1 )
                END IF
                
                L2     = LEN_TRIM( TDESC )
                SPCNAM = TDESC( L+LS:L2  )     ! extract species name

                K1 = INDEX1( POLNAM, NIPOL, EINAM )  ! find pol in sorted list
                K2 = INDEX1( TDESC, NCNT, TVDESCA )  ! find combo

C.................  When pollutant is found, store sorting variable
C.................  Sorting variable is designed to sort by pollutants (e.g.CO),
C                   then species, and then emission type (e.g. EXH_CO)
                IF( K1 .GT. 0 .AND. K2 .LE. 0 ) THEN
                    NCNT = NCNT + 1
 
                    WRITE( PBUF, '(I5.5)' ) K1
                    WRITE( IBUF, '(I5.5)' ) I
                    INDXA  ( NCNT ) = NCNT
                    TVSORTA( NCNT ) = PBUF // SPCNAM // IBUF
                    TVDESCA( NCNT ) = TDESC
                ELSE IF( K1 .LE. 0 ) THEN
                    L1 = LEN_TRIM( TDESC( 1:L-1 ) )
                    L2 = LEN_TRIM( SPCNAM )
                    MESG = 'NOTE: Speciation factor for "' // 
     &                      TDESC( 1:L1 ) // '-to-' // SPCNAM( 1:L2 )//
     &                     '" was not associated with ' // 
     &                     CRLF()// BLANK10// 'an inventory pollutant.'
                    CALL M3MSG2( MESG )

                END IF
   
            END DO

            END SUBROUTINE BUILD_VDESC_UNSORT

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.........  Build list of species
C.........  To do this, must only condense the pol-to-species list, in case
C           multiple pollutants are creating the same species.  Condense by
C           removing later-appearing duplicates
            SUBROUTINE BUILD_SPECIES_ARRAY( NSMATV,TSVDESC,NMSPC,EMNAM )

C.............  Subprogram arguments

            INTEGER     , INTENT (IN) :: NSMATV             ! input count
            CHARACTER(*), INTENT (IN) :: TSVDESC( NSMATV )  ! spec var descs
            INTEGER     , INTENT(OUT) :: NMSPC              ! output count
            CHARACTER(*), INTENT(OUT) :: EMNAM  ( NSMATV )  ! species list

C.............  Local subprogram variables
            INTEGER      I, J, K, N, NCNT, L1, L2      ! counters and indices
            INTEGER      IPTR                          ! tmp position in array

C----------------------------------------------------------------------

C.............  Populate array by searching remaining list of pollutants-to-
C               species for current iteration's value.  
            NCNT = 0
            DO I = 1, NSMATV 
                N = NSMATV - I

C.................  Search remaining list of pollutants-to-species for current
C                   iteration's value.
                IPTR = MIN( I+1,NSMATV )  
                J = INDEX1( TSVDESC( I ), N, TSVDESC( IPTR ) )

C.................  See if species is not already in the list of species
                L1 = INDEX( TSVDESC( I ), SPJOIN )
                L2 = LEN_TRIM( TSVDESC( I ) )
                K  = INDEX1( TSVDESC( I )( L1+1:L2 ), NCNT, EMNAM )

C.................  If pollutant-to- species is not found, make sure that 
C                   species only is not 
                IF( J .LE. 0 .AND. K .LE. 0 ) THEN
                    NCNT = NCNT + 1
                    EMNAM( NCNT ) = TSVDESC( I )( L1+1:L2 )
                ENDIF
            END DO
            NMSPC = NCNT

            END SUBROUTINE BUILD_SPECIES_ARRAY

        END SUBROUTINE MRGVNAMS
