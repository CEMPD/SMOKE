
        SUBROUTINE BLDMRGIDX( MXGRP, MXPOLPGP, MXSPPOL, NGRP )

C***********************************************************************
C  subroutine BLDMRGIDX body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to allocate and populate indicator
C      arrays that say which pollutants and species are present for each
C      type of input file (inventory, speciation matrix, multiplicative
C      control matrix, reactivity matrix, etc.) for each source category
C      (area, biogenic, mobile, point).  The indicator arrays store the
C      position in the list of i/o api file variables that the given pollutant
C      or species matches.
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

C...........  EXTERNAL FUNCTIONS
        INTEGER         INDEX1
        EXTERNAL        INDEX1

C...........  SUBROUTINE ARGUMENTS

        INTEGER, INTENT  (IN) :: MXGRP    ! Maximum no. groups
        INTEGER, INTENT  (IN) :: MXPOLPGP ! Max no. pollutants per group
        INTEGER, INTENT  (IN) :: MXSPPOL  ! Max no. species per pollutant
        INTEGER, INTENT (OUT) :: NGRP     ! Actual number of groups

C...........   Other local variables
        INTEGER         J, K, L, N, V    !  counters and indices

        INTEGER         ACNT     ! area src var counter
        INTEGER         BCNT     ! biogenic src var counter
        INTEGER         GCNT     ! group counter
        INTEGER         IOS      ! i/o error status
        INTEGER         ISDIM    ! tmp dimension for speciation
        INTEGER         JA, JM, JP ! position in rctvty var desc of pol-spc
        INTEGER         MCNT     ! mobile src var counter
        INTEGER         NLOOP    ! tmp loop total
        INTEGER         PCNT     ! pollutant counter OR point src var counter
        INTEGER         SCNT     ! species counter
        INTEGER         TGRP     ! tmp group number

        LOGICAL         NEXTGRP  ! true: time to increment group number cntr

        CHARACTER(LEN=IOVLEN3) :: CPOL ! tmp pollutant buffer 
        CHARACTER(LEN=IOVLEN3) :: CSPC ! tmp species buffer 
        CHARACTER(LEN=IOVLEN3) :: PPOL ! tmp previous pollutant  
        CHARACTER(LEN=IOVLEN3) :: VBUF ! tmp variable name buffer 
        CHARACTER(LEN=PLSLEN3) :: SVBUF ! tmp speciation name buffer 

        CHARACTER*16 :: PROGNAME = 'BLDMRGIDX' ! program name

C***********************************************************************
C   begin body of subroutine BLDMRGIDX

C.........  Allocate memory for all non-speciation indices and initialize all 
C           to zero...
C.........  Global pollutants
        ALLOCATE( PLCNT( MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PLCNT', PROGNAME )
        ALLOCATE( IDPGP( MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDPGP', PROGNAME )
        ALLOCATE( PLNAMES( MXPOLPGP, MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PLNAMES', PROGNAME )

C.........  Compute dimension for allocating EXIST arrays for speciation
        ISDIM = MAX( 1, MXSPPOL )

C.........  Allocate memory for all *EXIST arrays for a source category if
C           that source category is present.  This is necessary to simplify
C           using these arrays, and it will be a small memory penalty.
C.........  Inventory emissions: e.g., A_EXIST
C.........  Multiplicative control matrices: e.g., AU_EXIST
C.........  Additive control matrices: e.g., AA_EXIST
C.........  Reactivity control matrices: e.g., AS_EXIST

C.........  Allocate and initialize to zero...
        IF( AFLAG ) THEN   ! area

            ALLOCATE( A_EXIST( MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'A_EXIST', PROGNAME )
            ALLOCATE( AU_EXIST( MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AU_EXIST', PROGNAME )
            ALLOCATE( AA_EXIST( MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AA_EXIST', PROGNAME )
            ALLOCATE( AR_EXIST( ISDIM, MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AR_EXIST', PROGNAME )
            ALLOCATE( AS_EXIST( ISDIM, MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AS_EXIST', PROGNAME )
            A_EXIST  = 0  ! array
            AU_EXIST = 0  ! array
            AA_EXIST = 0  ! array
            AR_EXIST = 0  ! array
            AS_EXIST = 0  ! array
        ENDIF

        IF( BFLAG ) THEN   ! biogenics
            ALLOCATE( BS_EXIST( MXSPPOL,MXPOLPGP,MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BS_EXIST', PROGNAME )
            BS_EXIST = 0  ! array
        ENDIF

        IF( MFLAG ) THEN   ! mobile
            ALLOCATE( M_EXIST( MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'M_EXIST', PROGNAME )
            ALLOCATE( MU_EXIST( MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MU_EXIST', PROGNAME )
            ALLOCATE( MA_EXIST( MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MA_EXIST', PROGNAME )
            ALLOCATE( MR_EXIST( ISDIM, MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MR_EXIST', PROGNAME )
            ALLOCATE( MS_EXIST( ISDIM, MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MS_EXIST', PROGNAME )
            M_EXIST  = 0  ! array
            MU_EXIST = 0  ! array
            MA_EXIST = 0  ! array
            MR_EXIST = 0  ! array
            MS_EXIST = 0  ! array
        ENDIF

        IF( PFLAG ) THEN   ! point
            ALLOCATE( P_EXIST( MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'P_EXIST', PROGNAME )
            ALLOCATE( PU_EXIST( MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PU_EXIST', PROGNAME )
            ALLOCATE( PA_EXIST( MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PA_EXIST', PROGNAME )
            ALLOCATE( PR_EXIST( ISDIM, MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PR_EXIST', PROGNAME )
            ALLOCATE( PS_EXIST( ISDIM, MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PS_EXIST', PROGNAME )
            P_EXIST  = 0  ! array
            PU_EXIST = 0  ! array
            PA_EXIST = 0  ! array
            PR_EXIST = 0  ! array
            PS_EXIST = 0  ! array
        ENDIF
 
C.........  Allocate memory for speciation-based indices for global 
C           pol-to-species and species
C.........  Allocate the NSMPPG even if there is no speciation so that this
C           array will be valid in SMKMERGE no matter what
        ALLOCATE( NSMPPG( MXPOLPGP, MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NSMPPG', PROGNAME )
        NSMPPG  = 1  ! array

        IF( SFLAG ) THEN

            ALLOCATE( SMINDEX( MXSPPOL, MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SMINDEX', PROGNAME )
            ALLOCATE( SPINDEX( MXSPPOL, MXPOLPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPINDEX', PROGNAME )
            SMINDEX = 0  ! array
            SPINDEX = 0  ! array

        END IF

C.........  Create pollutant names array in groups structure, and store counts
C.........  For speciation, also build group indices for pol-to-species...

C.........  For speciation or not, set number of iterations for upcoming loop
        IF ( SFLAG ) THEN
            NLOOP = NSMATV
        ELSE
            NLOOP = NIPOL
            SCNT = 1      ! (will be > MXSPPOL=0)
        END IF

C.........  Loop through pollutants (if no speciation) or pol-to-species names
        PPOL = EMCMISS3  ! initialize previous pollutant
        PCNT = 0         ! initialize pollutant count
        GCNT = 1         ! initialize group count
        DO V = 1, NLOOP

C.............  Extract pollutant name and species name (for speciation only)
            IF( SFLAG ) THEN 
                L = INDEX( TSVDESC( V ), '_' )
                CPOL = TSVDESC( V )(   1:L-1 )
                CSPC = TSVDESC( V )( L+1:LEN_TRIM( TSVDESC( V ) ) )
            ELSE
                CPOL = EINAM( V )
            ENDIF

C.............  If same pollutant as previous, and species county is not
C               up to the maximum yet, increment species count
            IF( CPOL .EQ. PPOL .AND. SCNT .LT. MXSPPOL ) THEN

                SCNT = SCNT + 1

C.............  If different pollutant as previous, and pollutant count
C               is not up the the maximum for the group, increment
C               pollutant count and initialize species count
            ELSEIF( CPOL .NE. PPOL .AND. PCNT .LT. MXPOLPGP ) THEN

                PCNT = PCNT + 1        ! increment pol count
                SCNT = 1               ! initialize spc count

C.............  If different pollutant as previous, and pollutant count
C               is already at the maximum, increment group, and initialize 
C               pollutant and species counts
            ELSE
                GCNT = GCNT + 1        ! increment group
                PCNT = 1               ! re-initialize pol count
                SCNT = 1               ! re-initialize spc count

            END IF

C.............  Store pollutant names and number in group structure
            PLCNT  (       GCNT ) = PCNT ! store no. pols in group
            PLNAMES( PCNT, GCNT ) = CPOL ! store rearranged pol names

C.............  For speciation, store indices and counts
            IF( SFLAG ) THEN
                NSMPPG (       PCNT, GCNT ) = SCNT ! store no. spc in pol/group
                SMINDEX( SCNT, PCNT, GCNT ) = V    ! store pol-to-spc index

                J = INDEX1( CSPC, NMSPC, EMNAM )   ! get position of species
                SPINDEX( SCNT, PCNT, GCNT ) = J    ! store species index
            ENDIF

            PPOL = CPOL

        END DO

C.........  Store actual number of groups
        NGRP = GCNT

C.........  Create pollutant-group array by looping through groups and checking
C           if the pollutants in the current group are the same as those in the
C           previous group...

        IDPGP( 1 ) = 1       ! initialize group number for the first group
        DO N = 2, NGRP

            NEXTGRP = .FALSE.    ! reset indicator for incrementing group
            TGRP = IDPGP( N-1 )  ! set temporary group number

C.............  If any of the pollutants in this group are different from the
C               previous group, turn on indicator for incrementing group
            DO V = 1, PLCNT( N )
                IF( PLNAMES( V,N ) .NE. PLNAMES( V,N-1 ) ) 
     &              NEXTGRP = .TRUE.
            END DO

C.............  If indicator is true, increment group
            IF( NEXTGRP ) THEN
                TGRP = TGRP + 1
            ENDIF

C.............  Set group number for current group
            IDPGP( N ) = TGRP

        END DO

C.........  Precompute the position of the species variables in the reactivity
C           matrices
        JA   = ANRMATV - ARNMSPC + 1
        JM   = MNRMATV - MRNMSPC + 1
        JP   = PNRMATV - PRNMSPC + 1
 
C.........  Create source-categeory-specific indicies for speciation, inventory
C           emissions, and control matrices...

C.........  If speciation, loop through the number of pol-to-spec per
C           pollutant/group and store the count for each source category.
C.........  The counts are kept for each source category because each
C           source category has only it's own speciation factors stored
C           for only the pollutants in that category.
C.........  The counts are used because of the pollutant groups.  If the
C           position in the variable list was stored, this would not be
C           the correct index for a pollutant within a group.  The counts
C           work because all of the pollutant/species are sorted in the same
C           order for all source categories and files.
C.........  If there is reactivity, then also store the position in the
C           reactivity matrix

        IF ( SFLAG ) THEN

C.............  Loop through all groups, then number of pollutants per group
            DO N = 1, NGRP
                    
                ACNT = 0 
                BCNT = 0 
                MCNT = 0 
                PCNT = 0 

                DO V = 1, PLCNT( N )

C.....................  Loop through the number of pol-to-species
                    DO J = 1, NSMPPG ( V,N )

                        SVBUF = TSVDESC( SMINDEX( J,V,N ) )

                        IF ( AFLAG ) THEN    ! Area sources
                            K = INDEX1( SVBUF, ANSMATV, ASVDESC )
                            IF( K .GT. 0 ) THEN
                                ACNT = ACNT + 1
                                AS_EXIST( J,V,N ) = ACNT
                            END IF
                        END IF

                        IF ( ARFLAG ) THEN    ! Area reactivity
                            K = INDEX1( SVBUF, ARNMSPC, ARVDESC( JA ) )
                            AR_EXIST( J,V,N ) = K
                        END IF

                        IF ( BFLAG ) THEN    ! Biogenics sources
                            K = INDEX1( SVBUF, BNSMATV, BSVDESC )
                            IF( K .GT. 0 ) THEN
                                BCNT = BCNT + 1
                                BS_EXIST( J,V,N ) = BCNT
                            END IF
                        END IF

                        IF ( MFLAG ) THEN    ! Mobile sources
                            K = INDEX1( SVBUF, MNSMATV, MSVDESC )
                            IF( K .GT. 0 ) THEN
                                MCNT = MCNT + 1
                                MS_EXIST( J,V,N ) = MCNT
                            END IF
                        END IF

                        IF ( MRFLAG ) THEN    ! Mobile reactivity
                            K = INDEX1( SVBUF, MRNMSPC, MRVDESC( JM ) )
                            MR_EXIST( J,V,N ) = K
                        END IF

                        IF ( PFLAG ) THEN    ! Point sources
                            K = INDEX1( SVBUF, PNSMATV, PSVDESC )
                            IF( K .GT. 0 ) THEN
                                PCNT = PCNT + 1
                                PS_EXIST( J,V,N ) = PCNT
                            END IF
                        END IF

                        IF ( PRFLAG ) THEN    ! Point reactivity
                            K = INDEX1( SVBUF, PRNMSPC, PRVDESC( JP ) )
                            PR_EXIST( J,V,N ) = K
                        END IF

                    END DO  ! end loop on species per pollutant and group
                END DO      ! end loop on pollutants per group
            END DO          ! end loop on groups
        END IF              ! For speciation

C.........  Now build indices for inventory emissions...
C.........  For each source category, store per-group position for 
C           inventory pollutants...

C.........  Loop through all groups, then number of pollutants per group
        DO N = 1, NGRP

            ACNT = 0 
            BCNT = 0 
            MCNT = 0 
            PCNT = 0 

            DO V = 1, PLCNT( N )

                VBUF = PLNAMES( V,N )

                IF ( AFLAG ) THEN    ! Area sources
                    K = INDEX1( VBUF, ANIPOL, AEINAM )
                    IF( K .GT. 0 ) THEN
                        ACNT = ACNT + 1
                        A_EXIST( V,N ) = ACNT
                    END IF
                END IF

                IF ( MFLAG ) THEN    ! Mobile sources
                    K = INDEX1( VBUF, MNIPOL, MEINAM )
                    IF( K .GT. 0 ) THEN
                        MCNT = MCNT + 1
                        M_EXIST( V,N ) = MCNT
                    END IF
                END IF

                IF ( PFLAG ) THEN    ! Point sources
                    K = INDEX1( VBUF, PNIPOL, PEINAM )
                    IF( K .GT. 0 ) THEN
                        PCNT = PCNT + 1
                        P_EXIST( V,N ) = PCNT
                    END IF
                END IF
            END DO      ! end loop on pollutants per group
        END DO          ! end loop on groups

C.........  Now build indices for mulitplicative controls...
C.........  For each source category, store per-group position for 
C           inventory pollutants...

C.........  Loop through all groups, then number of pollutants per group
        DO N = 1, NGRP

            ACNT = 0 
            BCNT = 0 
            MCNT = 0 
            PCNT = 0 

            DO V = 1, PLCNT( N )

                VBUF = PLNAMES( V,N )

                IF ( AUFLAG ) THEN    ! Area sources
                    K = INDEX1( VBUF, ANIPOL, AEINAM )
                    IF( K .GT. 0 ) THEN
                        ACNT = ACNT + 1
                        AU_EXIST( V,N ) = ACNT
                    END IF
                END IF

                IF ( MUFLAG ) THEN    ! Mobile sources
                    K = INDEX1( VBUF, MNIPOL, MEINAM )
                    IF( K .GT. 0 ) THEN
                        MCNT = MCNT + 1
                        MU_EXIST( V,N ) = MCNT
                    END IF
                END IF

                IF ( PUFLAG ) THEN    ! Point sources
                    K = INDEX1( VBUF, PNIPOL, PEINAM )
                    IF( K .GT. 0 ) THEN
                        PCNT = PCNT + 1
                        PU_EXIST( V,N ) = PCNT
                    END IF
                END IF

            END DO      ! end loop on pollutants per group
        END DO          ! end loop on groups

C.........  Now build indices for additive controls...
C.........  For each source category, store per-group position for 
C           inventory pollutants...

C.........  Loop through all groups, then number of pollutants per group
        DO N = 1, NGRP

            ACNT = 0 
            BCNT = 0 
            MCNT = 0 
            PCNT = 0 

            DO V = 1, PLCNT( N )

                VBUF = PLNAMES( V,N )

                IF ( AAFLAG ) THEN    ! Area sources
                    K = INDEX1( VBUF, ANIPOL, AEINAM )
                    IF( K .GT. 0 ) THEN
                        ACNT = ACNT + 1
                        AA_EXIST( V,N ) = ACNT
                    END IF
                END IF

                IF ( MAFLAG ) THEN    ! Mobile sources
                    K = INDEX1( VBUF, MNIPOL, MEINAM )
                    IF( K .GT. 0 ) THEN
                        MCNT = MCNT + 1
                        MA_EXIST( V,N ) = MCNT
                    END IF
                END IF

                IF ( PAFLAG ) THEN    ! Point sources
                    K = INDEX1( VBUF, PNIPOL, PEINAM )
                    IF( K .GT. 0 ) THEN
                        PCNT = PCNT + 1
                        PA_EXIST( V,N ) = PCNT
                    END IF
                END IF

            END DO      ! end loop on pollutants per group
        END DO          ! end loop on groups

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE BLDMRGIDX
