
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE M3UTILIO

        USE MODMERGE, ONLY: PDEV,
     &                      MNIPPA, NIPPA,
     &                      MEANAM, EINAM,
     &                      EMNAM,
     &                      EANAM, EMIDX, 
     &                      NSMATV, 
     &                      TSVDESC, 
     &                      NMSPC

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: MXIDAT, INVDNAM, INVSTAT, INVDCOD, INVDVTS

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: NHAP, HAPNAM, MNSMATV_L, MSVDESC_L,
     &                       MSVUNIT_L, MSVUNIT_S, SPCUNIT_L, SPCUNIT_S,
     &                       GRDENV, TOTENV

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
C       CHARACTER(2)    CRLF
        INTEGER         GETFLINE
C       INTEGER         INDEX1

C        EXTERNAL   CRLF, GETFLINE, INDEX1
        EXTERNAL     GETFLINE

C...........   Temporary arrays for building sorted pol-to-species names list
        INTEGER           , ALLOCATABLE :: INDXA  ( : ) ! sorting index
        CHARACTER(PLSLEN3), ALLOCATABLE :: TVSORTA( : ) ! basis for sorting
        CHARACTER(PLSLEN3), ALLOCATABLE :: TVDESCA( : ) ! pol-to-spec concat

C...........   Other local variables
        INTEGER         I, J, J1, J2, L, K, L1, L2, M, N, V  ! counters and indices
        INTEGER         LJ              ! string length of emis types joiner
        INTEGER         IOS             ! i/o error status
        INTEGER         NCNT            ! counter

        LOGICAL      :: EFLAG = .FALSE. ! error flag
        LOGICAL      :: FOUND = .FALSE. ! true: found HAP

        CHARACTER(IOVLEN3)   CPOL     ! tmp pol/act buffer
        CHARACTER(IOULEN3)   BNUM     ! tmp biogenic units numerator
        CHARACTER(300)       MESG     ! message buffer

        CHARACTER(16) :: PROGNAME = 'MRGVNAMS' ! program name

C***********************************************************************
C   begin body of subroutine MRGVNAMS

C.........  Read, sort, and store pollutant codes/names file
        CALL RDCODNAM( PDEV )

C.........  Check if list of pollutants from MEPROC file contains
C           any HAPs
        DO I = 1, MNIPPA
            J = INDEX1( MEANAM( I ), MXIDAT, INVDNAM )
            
            IF( J .GT. 0 ) THEN
                IF( INVDVTS( J ) == 'V' .OR.
     &              INVDVTS( J ) == 'T' ) THEN
                    FOUND = .TRUE.
                    EXIT
                END IF
            END IF

        END DO

C.........  If any HAP was found, update TOG to NONHAPTOG
        IF( FOUND ) THEN
            DO I = 1, MNIPPA
                IF( MEANAM( I ) == 'TOG' ) THEN
                    MEANAM( I ) = 'NONHAP'//TRIM( MEANAM(I) )
                END IF
            END DO
        END IF

C.........  Loop through emission process/pollutant combinations 
C           and update status of entry in master list.
C           Also count number of pollutants to subtract when calculating
C           NONHAPTOG
        NHAP = 0
        DO I = 1, MNIPPA
        
            CPOL = MEANAM( I )
            J = INDEX1( CPOL, MXIDAT, INVDNAM )
            IF( J .GT. 0 ) THEN
                INVSTAT( J ) = INVSTAT( J ) * 2
                
                IF( INVDVTS( J ) == 'V' .OR.
     &              INVDVTS( J ) == 'T' ) THEN
                    NHAP = NHAP + 1
                END IF
            ELSE
                EFLAG = .TRUE.
                L = LEN_TRIM( CPOL )
                MESG = 'ERROR: Variable "' // CPOL( 1:L ) //
     &                 '" from MEPROC file is not in ' // CRLF() //
     &                 BLANK10 // 'inventory table.'
                CALL M3MSG2( MESG )
            END IF

        END DO

C.........  If any process/pollutant names are not found in master
C           list, exit
        IF( EFLAG ) THEN
            MESG = 'Make sure that master pollutants and/or ' //
     &             'activities files are' // CRLF() // BLANK10 // 
     &             'consistent with those used to process the ' //
     &             'inventory.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, convert INVSTAT back from 2/-2 to 1/-1, and zero for
C           not used
        ELSE

            INVSTAT = INVSTAT / 2  ! integer math for array

        END IF

C.........  Allocate memory for list of HAPS
        ALLOCATE( HAPNAM( NHAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HAPNAM', PROGNAME )

C.........  Loop through master list count the actual number of 
C           process/pollutant combinations
        NIPPA = 0
        NHAP = 0
        DO I = 1, MXIDAT
            IF( INVSTAT( I ) .GT. 0 ) THEN
                NIPPA = NIPPA + 1
                
                IF( INVDVTS( I ) == 'V' .OR.
     &              INVDVTS( I ) == 'T' ) THEN
                    NHAP = NHAP + 1
                    HAPNAM( NHAP ) = INVDNAM( I )
                END IF
            END IF
        END DO

C.........  Allocate memory for array of sorted process/pollutant names and
C           for pollutants only
        ALLOCATE( EANAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAM', PROGNAME )
        ALLOCATE( EINAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EINAM', PROGNAME )

C.........  Initialize all
        EANAM   = ' '   ! array
        EINAM   = ' '   ! array

C.........  Create array of process/pollutant names matching order of INVTABLE
C.........  Also store pollutants-only array
        J1 = 0
        DO I = 1, MXIDAT

            IF( INVSTAT( I ) .GT. 0 ) THEN
                J1 = J1 + 1
                EANAM( J1 ) = INVDNAM( I )
                EINAM( J1 ) = EANAM( J1 )
            END IF

        END DO

C.........  Create array of sorted unique pol-to-species, sorted in order of
C           pollutants, and then in alphabetical order by species...

        ALLOCATE( INDXA( MNSMATV_L ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXA', PROGNAME )
        ALLOCATE( TVSORTA( MNSMATV_L ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TVSORTA', PROGNAME )
        ALLOCATE( TVDESCA( MNSMATV_L ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TVDESCA', PROGNAME )

C.........  Loop through variable descriptions, find position of
C           of pollutant for output, and store concatenated position with
C           species name.  Make sure the same name is not stored twice.
        NCNT = 0
        CALL BUILD_VDESC_UNSORT( NCNT, MNSMATV_L, MSVDESC_L )

        NSMATV = NCNT 

C.........  Make sure some species match some of the pollutants
        IF( NSMATV .EQ. 0 ) THEN
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
        ALLOCATE( EMNAM( NSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMNAM', PROGNAME )
        ALLOCATE( EMIDX( NSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMIDX', PROGNAME )
        ALLOCATE( SPCUNIT_L( NSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCUNIT_L', PROGNAME )
        ALLOCATE( SPCUNIT_S( NSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCUNIT_S', PROGNAME )
        EMNAM     = ' '  ! array
        SPCUNIT_L = ' '  ! array
        SPCUNIT_S = ' '  ! array

C.........  Call subprogram to store species names in appropriate order
        CALL BUILD_SPECIES_ARRAY(  NSMATV, TSVDESC,  NMSPC,  EMNAM )

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
            IF( SPCUNIT_L( K ) .EQ. ' ' ) THEN

                M = INDEX1( TSVDESC( I ), MNSMATV_L, MSVDESC_L )
                IF( M .GT. 0 ) THEN
                    SPCUNIT_L( K ) = MSVUNIT_L( M )
                    SPCUNIT_S( K ) = MSVUNIT_S( M )
                END IF

            END IF

        END DO

C.........  Switch speciation matrix between moles and mass based on
C           the setting of  MRG_GRDOUT_UNIT and MRG_TOTOUT_UNIT 
        IF( INDEX( GRDENV, 'mole' ) < 1 ) THEN
            SPCUNIT_L = SPCUNIT_S
        END IF

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

C.............  This internal subprogram builds the unsorted list of unique
C               pollutant-to-species speciation variable descriptions.
            SUBROUTINE BUILD_VDESC_UNSORT( NCNT, NVARS, VDESCS )

C            INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C.............  Subprogram arguments

            INTEGER     , INTENT(IN OUT) :: NCNT            ! running count
            INTEGER     , INTENT    (IN) :: NVARS           ! no of var descs
            CHARACTER(*), INTENT    (IN) :: VDESCS( NVARS ) ! spec var descs

C.............  Local subprogram arrays
            INTEGER,            ALLOCATABLE, SAVE :: EMBIN ( : )
            CHARACTER(IOVLEN3), ALLOCATABLE, SAVE :: EMLIST( : )

C.............  Local subprogram variables
            INTEGER, SAVE :: EMCNT = 0            ! count of species
            INTEGER          I, K1, K2, K3, L, L2 ! counters and indices
            INTEGER          LJ, LS               ! string lengths of separators

            CHARACTER(5)           CBIN   ! tmp pollutant bin
            CHARACTER(5)           CPOL   ! tmp pollutant index
            CHARACTER(5)           IBUF   ! tmp variable position to INVDNAM
            CHARACTER(IOVLEN3)     POLNAM ! tmp pollutant name
            CHARACTER(IOVLEN3)     SPCNAM ! tmp species name
            CHARACTER(IODLEN3)     TDESC  ! tmp combined pollutant & species

C----------------------------------------------------------------------

C.............  Allocate local memory
            IF ( .NOT. ALLOCATED( EMBIN ) ) THEN
                ALLOCATE( EMBIN( NVARS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'EMBIN', PROGNAME )
                ALLOCATE( EMLIST( NVARS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'EMLIST', PROGNAME )
                EMBIN  = 0
                EMLIST = ' '
            END IF

            LS = LEN_TRIM( SPJOIN )
            DO I = 1, NVARS

                TDESC  = VDESCS( I )
                L      = INDEX( TDESC, SPJOIN )  ! find species separator
                IF( L .LE. 0 ) THEN
                    MESG = 'ERROR: Variable descriptions do not ' //
     &                     'contain proper separator in spec matrix.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                POLNAM = TDESC( 1:L-1 )
                
                L2     = LEN_TRIM( TDESC )
                SPCNAM = TDESC( L+LS:L2  )     ! extract species name

                K1 = INDEX1( POLNAM, NIPPA, EINAM )  ! find pol in sorted list
                K2 = INDEX1( TDESC, NCNT, TVDESCA )  ! find combo

C.................  Look for species name in temporary list. If found, assign
C                   pollutant bin number, if not, initialize pollutant bin 
C                   number.                
                K3 = INDEX1( SPCNAM, EMCNT, EMLIST )

                IF( K3 .GT. 0 ) THEN
                    K3 = EMBIN( K3 )
                    
                ELSE IF ( EMCNT + 1 <= NVARS ) THEN
                    EMCNT = EMCNT + 1
                    EMLIST( EMCNT ) = SPCNAM
                    EMBIN ( EMCNT ) = K1
                    K3 = K1

                ELSE
                    MESG = 'INTERNAL ERROR: Array overflow for EMBIN.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                                        
                END IF

C.................  When pollutant is found, store sorting variable
C.................  Sorting variable is designed to sort by pollutant (e.g.CO),
C                   then species, and then emission type (e.g. EXH_CO).
C.................  If multiple pollutants contribute to the same species, the
C                   sorting must ensure the species are together. An example
C                   of this is ROG and VOC creating the same species. The
C                   CBIN variable ensures that this is the case.
                IF( K1 .GT. 0 .AND. K2 .LE. 0 ) THEN
                    NCNT = NCNT + 1
 
                    WRITE( CBIN, '(I5.5)' ) K3
                    WRITE( CPOL, '(I5.5)' ) K1
                    WRITE( IBUF, '(I5.5)' ) I
                    INDXA  ( NCNT ) = NCNT
                    TVSORTA( NCNT ) = CBIN // SPCNAM // CPOL // IBUF
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
            SUBROUTINE BUILD_SPECIES_ARRAY( LSMATV,LOCDESC,LMSPC,LOCNAM)

C.............  Subprogram arguments

            INTEGER     , INTENT (IN) :: LSMATV             ! input count
            CHARACTER(*), INTENT (IN) :: LOCDESC( LSMATV )  ! spec var descs
            INTEGER     , INTENT(OUT) :: LMSPC              ! output count
            CHARACTER(*), INTENT(OUT) :: LOCNAM ( LSMATV )  ! species list

C.............  Local subprogram variables
            INTEGER      I, J, K, N, NCNT, L1, L2      ! counters and indices
            INTEGER      IPTR                          ! tmp position in array

C----------------------------------------------------------------------

C.............  Populate array by looping through master list of pol-to-species in
C               output order, and if that entry is there for local subroutine
C               arguments, and if the species hasn't been added yet, then add it.  
            NCNT = 0
            DO I = 1, NSMATV 

C.................  Check if master pol-to-species is in local list. If not cycle
                J = INDEX1( TSVDESC( I ), LSMATV, LOCDESC )
                IF( J .LE. 0 ) CYCLE
                N = NSMATV - I

C.................  Search remaining list of pollutants-to-species for current
C                   iteration's value.
                IPTR = MIN( I+1,NSMATV )  
                J = INDEX1( TSVDESC( I ), N, TSVDESC( IPTR ) )

C.................  See if species is not already in the list of species
                L1 = INDEX( TSVDESC( I ), SPJOIN )
                L2 = LEN_TRIM( TSVDESC( I ) )
                K  = INDEX1( TSVDESC( I )( L1+1:L2 ), NCNT, LOCNAM )

C.................  If pollutant-to- species is not found, make sure that 
C                   species only is not 
                IF( J .LE. 0 .AND. K .LE. 0 ) THEN
                    NCNT = NCNT + 1
                    LOCNAM( NCNT ) = TSVDESC( I )( L1+1:L2 )
                ENDIF
            END DO
            LMSPC = NCNT

            END SUBROUTINE BUILD_SPECIES_ARRAY

        END SUBROUTINE MRGVNAMS
