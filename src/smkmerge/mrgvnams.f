
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
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
        INTEGER                                MXIPOL       ! max no pollutants
        INTEGER                             :: NPOL   = 0   ! actual no pols
        INTEGER               , ALLOCATABLE :: INVPCOD( : ) ! codes
        LOGICAL               , ALLOCATABLE :: INVSTAT( : ) ! true: used in run
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: INVPNAM( : ) ! pollutant name

C...........   Temporary arrays for building sorted pol-to-species names list
        INTEGER                                MXPTOSPC     ! max for dimension
        INTEGER               , ALLOCATABLE :: INDXA  ( : ) ! sorting index
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE :: TVSORTA( : ) ! basis for sorting
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE :: TVDESCA( : ) ! pol-to-spec concat

C...........   Other local variables
        INTEGER         I, J, L, L2, N  ! counters and indices
        INTEGER         IOS             ! i/o error status
        INTEGER         JA, JM, JP      ! position in rctvty var desc of pol-spc
        INTEGER         NCNT            ! counter

        LOGICAL      :: EFLAG = .FALSE. ! error flag

        CHARACTER*300   MESG     ! message buffer

        CHARACTER*16 :: PROGNAME = 'MRGVNAMS' ! program name

C***********************************************************************
C   begin body of subroutine MRGVNAMS

C.........  Get number of lines of pollutant codes/names file 
        MXIPOL = GETFLINE( PDEV, 'Pollutant codes and names file')

C.........  Allocate memory for storing contents of pollutants file
        ALLOCATE( INVPCOD( MXIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVPCOD', PROGNAME )
        ALLOCATE( INVPNAM( MXIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVPNAM', PROGNAME )
        ALLOCATE( INVSTAT( MXIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVSTAT', PROGNAME )

        INVSTAT = .FALSE.  ! array

C.........  Read and sort pollutant codes/names file
        CALL RDCODNAM( PDEV, MXIPOL, NPOL, INVPCOD, INVPNAM )

C.........  Loop through pollutants in each source category and update status
C           of pollutants in master list.
        CALL SET_MASTER_POL_STAT( 'area',   ANIPOL, AEINAM, EFLAG )
        CALL SET_MASTER_POL_STAT( 'mobile', MNIPOL, MEINAM, EFLAG )
        CALL SET_MASTER_POL_STAT( 'point',  PNIPOL, PEINAM, EFLAG )

C.........  If any pollutants from matrices not found in master list, exit
        IF( EFLAG ) THEN
            MESG = 'Make sure that master pollutants file is ' //
     &             'consistent with the one used to ' // CRLF() //
     &             BLANK10 // 'process the inventory.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Loop through master list count the actual number of total pollutants
        NIPOL = 0
        DO I = 1, MXIPOL
            IF( INVSTAT( I ) ) NIPOL = NIPOL + 1
        END DO

C.........  Allocate memory for array of sorted polllutants
        ALLOCATE( EINAM( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EINAM', PROGNAME )
        
C.........  Create array of sorted unique pollutants
        J = 0
        DO I = 1, MXIPOL
            IF( INVSTAT( I ) ) THEN
                J = J + 1
                EINAM( J ) = INVPNAM( I )
            END IF
        END DO
        
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
        CALL BUILD_VDESC_UNSORT( NCNT, ANSMATV, ASVDESC       ) ! area spectn
        CALL BUILD_VDESC_UNSORT( NCNT, ARNMSPC, ARVDESC( JA ) ) ! area rctvty
        CALL BUILD_VDESC_UNSORT( NCNT, BNSMATV, BSVDESC       ) ! biogenics
        CALL BUILD_VDESC_UNSORT( NCNT, MNSMATV, MSVDESC       ) ! mobile spectn
        CALL BUILD_VDESC_UNSORT( NCNT, MRNMSPC, MRVDESC( JM ) ) ! mobile rctvty
        CALL BUILD_VDESC_UNSORT( NCNT, PNSMATV, PSVDESC       ) ! point spectn
        CALL BUILD_VDESC_UNSORT( NCNT, PRNMSPC, PRVDESC( JP ) ) ! point rctvty
        NSMATV = NCNT

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

C.........  Call subprogram to store species names in appropriate order
C           for all source categories and the total.
        CALL BUILD_SPECIES_ARRAY( ANSMATV, ASVDESC, ANMSPC, AEMNAM )
        CALL BUILD_SPECIES_ARRAY( BNSMATV, BSVDESC, BNMSPC, BEMNAM )
        CALL BUILD_SPECIES_ARRAY( MNSMATV, MSVDESC, MNMSPC, MEMNAM )
        CALL BUILD_SPECIES_ARRAY( PNSMATV, PSVDESC, PNMSPC, PEMNAM )
        CALL BUILD_SPECIES_ARRAY(  NSMATV, TSVDESC, NMSPC, EMNAM )

C.........  Deallocate temporary arrays
        DEALLOCATE( INDXA, TVSORTA, TVDESCA, INVPCOD, INVPNAM, INVSTAT )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS
 
C.............  This internal subprogram updates the status for the master
C               list of inventory pollutants
            SUBROUTINE SET_MASTER_POL_STAT( CATDESC, N, PNAMES, EFLAG )

C.............  Subprogram arguments
            CHARACTER(*) CATDESC       ! source category description
            INTEGER      N             ! number of source-specific pollutants
            CHARACTER(*) PNAMES( N )   ! names of source-specific pollutants
            LOGICAL      EFLAG         ! error flag

C.............  Local subprogram variables
            INTEGER      I, J, L      ! counters and indices        

            CHARACTER(LEN=IOVLEN3) POLNAM

C----------------------------------------------------------------------

            DO I = 1, N

                POLNAM = PNAMES( I )
                J = INDEX1( POLNAM, MXIPOL, INVPNAM )
                IF( J .GT. 0 ) THEN
                    INVSTAT( J ) = .TRUE.

                ELSE
                    EFLAG = .TRUE.
                    L = LEN_TRIM( POLNAM )
                    MESG = 'ERROR: Pollutant "' // POLNAM( 1:L ) // 
     &                     '" from ' // CATDESC // 
     &                     ' source speciation matrix is ' // CRLF() //
     &                     BLANK10 // 'not in pollutants file.'
                    CALL M3MSG2( MESG )

                END IF

            END DO
 
            END SUBROUTINE SET_MASTER_POL_STAT
C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram builds the unsorted list of unique
C               pollutant-to-species speciation variable descriptions.
            SUBROUTINE BUILD_VDESC_UNSORT( NCNT, NVARS, VDESCS )

C.............  Subprogram arguments

            INTEGER      NCNT             ! running count
            INTEGER      NVARS            ! number of variable descriptions
            CHARACTER(*) VDESCS( NVARS )  ! speciation variable descriptions

C.............  Local subprogram variables
            INTEGER      I, K1, K2, L, L2      ! counters and indices
     
            CHARACTER*5                PBUF   ! tmp pollutant position
            CHARACTER(LEN=IOVLEN3)     POLNAM ! tmp pollutant name
            CHARACTER(LEN=IOVLEN3)     SPCNAM ! tmp species name
            CHARACTER(LEN=IOVLEN3*2+1) TDESC  ! tmp combined pollutant & species

C----------------------------------------------------------------------

            DO I = 1, NVARS

                L      = INDEX( VDESCS( I ), '_' )  ! find parsing '_'
                POLNAM = VDESCS( I )(   1:L-1 )     ! extract pollutant name

                L2     = LEN_TRIM( VDESCS( I ) )
                SPCNAM = VDESCS( I )( L+1:L2  )     ! extract species name

                TDESC  = POLNAM( 1:IOVLEN3 ) // '_' // SPCNAM

                K1 = INDEX1( POLNAM, NIPOL, EINAM )  ! find pol in sorted list
                K2 = INDEX1( TDESC, NCNT, TVDESCA )  ! find combo

C.................  When pollutant is found, store sorting variable
                IF( K1 .GT. 0 .AND. K2 .LE. 0 ) THEN
                    NCNT = NCNT + 1
 
                    WRITE( PBUF, '(I5.5)' ) K1
                    INDXA  ( NCNT ) = NCNT
                    TVSORTA( NCNT ) = PBUF // SPCNAM
                    TVDESCA( NCNT ) = TDESC
                END IF
   
            END DO

            END SUBROUTINE BUILD_VDESC_UNSORT

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.........  To do this, must only condense the pol-to-species list, in case
C           multiple pollutants are creating the same species.  Condense by
C           removing later-appearing duplicates
            SUBROUTINE BUILD_SPECIES_ARRAY( NSMATV,TSVDESC,NMSPC,EMNAM )

C.............  Subprogram arguments

            INTEGER      NSMATV             ! input count
            CHARACTER(*) TSVDESC( NSMATV )  ! speciation variable descriptions
            INTEGER      NMSPC              ! output count
            CHARACTER(*) EMNAM( NSMATV  )   ! species list

C.............  Local subprogram variables
            INTEGER      I, J, N, NCNT, L, L2      ! counters and indices
     
C----------------------------------------------------------------------

C.............  Populate array by searching remaining list for current 
C               iteration's value
            NCNT = 0
            DO I = 1, NSMATV
                N = NSMATV - I
                J = INDEX1( TSVDESC( I ), N, TSVDESC( I+1 ) )

                IF( J .LE. 0 ) THEN
                    NCNT = NCNT + 1
                    L  = INDEX( TSVDESC( I ), '_' )
                    L2 = LEN_TRIM( TSVDESC( I ) )
                    EMNAM( NCNT ) = TSVDESC( I )( L+1:L2 )
                ENDIF
            END DO
            NMSPC = NCNT

            END SUBROUTINE BUILD_SPECIES_ARRAY

        END SUBROUTINE MRGVNAMS
