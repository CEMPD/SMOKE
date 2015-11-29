
        SUBROUTINE BLDMRGIDX( MXGRP, MXVARPGP, NGRP )

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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: AFLAG, BFLAG, MFLAG, PFLAG,
     &                      AUFLAG, MUFLAG, PUFLAG,
     &                      ARFLAG, MRFLAG, PRFLAG, SFLAG,
     &                      ANIPOL, PNIPOL, MNIPPA, NIPPA,
     &                      ANSMATV, BNSMATV, MNSMATV, PNSMATV, NSMATV, 
     &                      ANUMATV, MNUMATV, PNUMATV,
     &                      ANRMATV, MNRMATV, PNRMATV,
     &                      ARNMSPC, MRNMSPC, PRNMSPC, NMSPC,
     &                      AEINAM, MEANAM, PEINAM,
     &                      ASVDESC, BSVDESC, MSVDESC, PSVDESC, TSVDESC,
     &                      ARVDESC, MRVDESC, PRVDESC,
     &                      AUVNAMS, MUVNAMS, PUVNAMS,
     &                      SIINDEX, SPINDEX,
     &                      A_EXIST, M_EXIST, P_EXIST,
     &                      AU_EXIST, MU_EXIST, PU_EXIST,
     &                      AR_EXIST, MR_EXIST, PR_EXIST,
     &                      AS_EXIST, BS_EXIST, MS_EXIST, PS_EXIST,
     &                      EMNAM, EANAM, VGRPCNT,
     &                      IDVGP, GVNAMES, GVLOUT

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  EXTERNAL FUNCTIONS
        CHARACTER(2)    CRLF   
        INTEGER         INDEX1
        EXTERNAL        CRLF, INDEX1

C...........  SUBROUTINE ARGUMENTS

        INTEGER, INTENT  (IN) :: MXGRP    ! Maximum no. groups
        INTEGER, INTENT  (IN) :: MXVARPGP ! Max no. pollutants per group
        INTEGER, INTENT (OUT) :: NGRP     ! Actual number of groups

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: PPSCNT( : )  ! count of pols per spc or emtype
        INTEGER, ALLOCATABLE :: KAU   ( : )  ! helper array for area mult matrix indices
        INTEGER, ALLOCATABLE :: KMU   ( : )  ! helper array for mobile mult matrix indices
        INTEGER, ALLOCATABLE :: KPU   ( : )  ! helper array for point mult matrix indices

C...........   Call allocated arrays
C...........   Group index counter for each source-category-specific list of 
C              pollutants and activities.
        INTEGER  KA( NIPPA )    !  area
        INTEGER  KM( NIPPA )    !  mobile
        INTEGER  KP( NIPPA )    !  point

C...........   Other local variables
        INTEGER         J, K, L1, L2, N, V    !  counters and indices

        INTEGER         ACNT     ! area src var counter
        INTEGER         BCNT     ! biogenic src var counter
        INTEGER         GCNT     ! group counter
        INTEGER         IOS      ! i/o error status
        INTEGER         JA, JM, JP ! position in rctvty var desc of pol-spc
        INTEGER         MAXPPS   ! max no. pols per species or per emis type
        INTEGER         MCNT     ! mobile src var counter
        INTEGER         NLOOP    ! tmp loop total
        INTEGER         NSIZE    ! tmp size of PPSCNT
        INTEGER         PCNT     ! point src var counter
        INTEGER         PGCNT    ! previous iteration group count
        INTEGER         PJ       ! previous iteration J index
        INTEGER         PVCNT    ! previous iteration variable count
        INTEGER         VCNT     ! variable counter
        INTEGER         TGRP     ! tmp group number

        LOGICAL      :: EFLAG = .FALSE.  ! true: error found
        LOGICAL         NEXTGRP  ! true: time to increment group number cntr

        CHARACTER(300)     :: MESG ! message buffer
        CHARACTER(IODLEN3) :: CBUF ! tmp pol-to-species buffer
        CHARACTER(IOVLEN3) :: CPOL ! tmp pollutant buffer 
        CHARACTER(IOVLEN3) :: CSPC ! tmp species buffer 
        CHARACTER(IOVLEN3) :: PSPC ! tmp previous species  
        CHARACTER(IOVLEN3) :: PPOL ! tmp previous pollutant  
        CHARACTER(IOVLEN3) :: VBUF ! tmp variable name buffer 
        CHARACTER(PLSLEN3) :: SVBUF ! tmp speciation name buffer 

        CHARACTER(16) :: PROGNAME = 'BLDMRGIDX' ! program name

C***********************************************************************
C   begin body of subroutine BLDMRGIDX

C.........  Ensure that the max no. of variables per group is large 
C           enough so that all pol-to-species combos for the same species
C           are in one group...

C.........  For speciation or not, set number of iterations for upcoming loop
C           and the dimension for the counter
        IF ( SFLAG ) THEN
            NLOOP = NSMATV
            NSIZE = NMSPC
        ELSE
            NLOOP = NIPPA
            NSIZE = NIPPA
        END IF

C.........  Allocate tmp pol per species count (or dummy ol-per-pol)
        ALLOCATE( PPSCNT( NSIZE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PPSCNT', PROGNAME )
        PPSCNT = 0.  ! array

C.........  Set the number of pollutants per species or emission type
        DO V = 1, NLOOP

C.............  For speciation, get position of species
            IF( SFLAG ) THEN
                CBUF = TSVDESC( V )
                L1 = INDEX   ( CBUF, SPJOIN )
                L2 = LEN_TRIM( CBUF )
                CSPC = CBUF( L1+1:L2   )
                J = INDEX1( CSPC, NMSPC, EMNAM )   ! get position of species

C.............  Otherwise, determine position in the list of pols/activities
            ELSE
                CPOL = EANAM( V )
                CBUF = CPOL
                L2 = LEN_TRIM( CBUF )
                J = INDEX1( CPOL, NIPPA, EANAM )   ! get position of pol/act

            END IF

C.............  Make sure variable was found in its list
            IF( J .LE. 0 ) THEN

                EFLAG = .TRUE.
                MESG = 'ERROR: variable "'// CBUF( 1:L2 ) // 
     &                 '" not found in master list.'
                CALL M3MSG2( MESG )

C.............  Count the number per species or pollutant
            ELSE
                PPSCNT( J ) = PPSCNT( J ) + 1
 
            END IF

        END DO

        IF( EFLAG ) THEN
             MESG = 'Problem building indices'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Get the maximum PPSCNT value
        MAXPPS = MAXVAL( PPSCNT )

C.........  Make sure MXVARPGP is large enough for MAXPPS
        IF( MAXPPS .GT. MXVARPGP ) THEN

            WRITE( MESG,94010 )  'ERROR: the maximum variables '//  
     &             'per group of', MXVARPGP, CRLF() // BLANK10 //
     &             'could not support the maximum pollutants ' //
     &             'per species of', MAXPPS, CRLF() // BLANK10 //
     &             'Not enough memory available to run for ' //
     &             'selected settings.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF
 
C.........  Allocate memory the group variables...
        ALLOCATE( VGRPCNT( MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VGRPCNT', PROGNAME )
        ALLOCATE( IDVGP( MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDVGP', PROGNAME )
        ALLOCATE( GVNAMES( MXVARPGP, MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GVNAMES', PROGNAME )
        ALLOCATE( GVLOUT ( MXVARPGP, MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GVLOUT', PROGNAME )

C.........  Initialize group variables
        VGRPCNT = 0        ! array
        IDVGP   = 0        ! array
        GVNAMES = ' '      ! array
        GVLOUT  = .FALSE.  ! array

C.........  Allocate memory for all *EXIST arrays for a source category if
C           that source category is present.  This is necessary to simplify
C           using these arrays, and it will be a small memory penalty.
C.........  NOTE - all arrays are dimensioned by MXVARGRP, but not all
C           need that many (esp. the arrays that contain pollutant info when
C           program is run with speciation)
C.........  Inventory emissions: e.g., A_EXIST
C.........  Multiplicative control matrices: e.g., AU_EXIST
C.........  Reactivity control matrices: e.g., AR_EXIST
C.........  Speciation matrices: e.g., AS_EXIST

C.........  Allocate and initialize to zero...
        IF( AFLAG ) THEN   ! area

            ALLOCATE( A_EXIST ( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'A_EXIST', PROGNAME )
            ALLOCATE( AU_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AU_EXIST', PROGNAME )
            ALLOCATE( AR_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AR_EXIST', PROGNAME )
            ALLOCATE( AS_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AS_EXIST', PROGNAME )
            A_EXIST  = 0  ! array
            AU_EXIST = 0  ! array
            AR_EXIST = 0  ! array
            AS_EXIST = 0  ! array
        ENDIF

        IF( BFLAG ) THEN   ! biogenics
            ALLOCATE( BS_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BS_EXIST', PROGNAME )
            BS_EXIST = 0  ! array
        END IF

        IF( MFLAG ) THEN   ! mobile
            ALLOCATE( M_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'M_EXIST', PROGNAME )
            ALLOCATE( MU_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MU_EXIST', PROGNAME )
            ALLOCATE( MR_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MR_EXIST', PROGNAME )
            ALLOCATE( MS_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MS_EXIST', PROGNAME )
            M_EXIST  = 0  ! array
            MU_EXIST = 0  ! array
            MR_EXIST = 0  ! array
            MS_EXIST = 0  ! array
        ENDIF

        IF( PFLAG ) THEN   ! point
            ALLOCATE( P_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'P_EXIST', PROGNAME )
            ALLOCATE( PU_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PU_EXIST', PROGNAME )
            ALLOCATE( PR_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PR_EXIST', PROGNAME )
            ALLOCATE( PS_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PS_EXIST', PROGNAME )
            P_EXIST  = 0  ! array
            PU_EXIST = 0  ! array
            PR_EXIST = 0  ! array
            PS_EXIST = 0  ! array
        ENDIF

C.........  Allocate memory for variable to pollutant and variable to species
C           indexes. Variable can be either pollutant or pol-to-species
C.........  NOTE - SPINDEX will not be used if there is no speciation
        ALLOCATE( SIINDEX( MXVARPGP, MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SIINDEX', PROGNAME )
        ALLOCATE( SPINDEX( MXVARPGP, MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPINDEX', PROGNAME )
        SIINDEX = 0  ! array
        SPINDEX = 0  ! array

C.........  Create pollutant names array in groups structure, and store counts
C.........  For speciation, also build group indices for pol-to-species...

C.........  Loop through pollutants (if no speciation) or pol-to-species names
        VCNT = 0         ! initialize variable count
        GCNT = 1         ! initialize group count
        PJ   = 0
        PSPC = ' '
        DO V = 1, NLOOP

C.............  Extract pollutant name and species name (for speciation only)
            IF( SFLAG ) THEN 
                CBUF = TSVDESC( V )
                L1 = INDEX   ( CBUF, SPJOIN )
                L2 = LEN_TRIM( CBUF )
                CPOL = CBUF(    1:L1-1 )
                CSPC = CBUF( L1+1:L2   )
            ELSE
                CPOL = EANAM( V )
                CBUF = CPOL
            END IF

C.............  Determine location in the list of pollutants/activities
            J = INDEX1( CPOL, NIPPA, EANAM )
            K = J
            IF( SFLAG ) J = INDEX1( CSPC, NMSPC, EMNAM )

C.............  Ensure that all pollutants for the same species or emission
C               type will be in the same group
            IF( J .NE. PJ .AND. VCNT .GT. MXVARPGP ) THEN
                GCNT = GCNT + 1
                VCNT = 1

C.............  Group count is not up to the maximum yet, increment variable
C               count
            ELSE IF( VCNT .LT. MXVARPGP ) THEN

                VCNT = VCNT + 1

C.............  Otherwise increment group count and initialize count
            ELSE
                GCNT = GCNT + 1 
                VCNT = 1       

            END IF

C.............  Store variable names and number in group structure
            VGRPCNT(       GCNT ) = VCNT ! store no. variables in group
            GVNAMES( VCNT, GCNT ) = CBUF ! store rearranged variable names

C.............  Store pol/act and species indices
            SIINDEX( VCNT, GCNT ) = K             ! store pol/act index
            IF( SFLAG ) SPINDEX( VCNT, GCNT ) = J ! store species index

C.............  Determine the iterations on which output should occur
            IF( V     .NE. 1    .AND. 
     &          SFLAG           .AND. 
     &          CSPC  .NE. PSPC       ) THEN

                GVLOUT( PVCNT, PGCNT ) = .TRUE.

            ELSE IF( V    .NE.  1     .AND. 
     &                    .NOT. SFLAG .AND. 
     &               CPOL .NE.  PPOL        ) THEN

                GVLOUT( PVCNT, PGCNT ) = .TRUE.

            END IF

            PSPC  = CSPC
            PPOL  = CPOL
            PVCNT = VCNT
            PGCNT = GCNT
            PJ    = J

        END DO

C.........  Store actual number of groups
        NGRP = GCNT

C.........  Ensure output on last variable in all groups
        GVLOUT( VGRPCNT( NGRP ), NGRP ) = .TRUE.

C.........  Create pollutant-group array by looping through groups and checking
C           if the pollutants in the current group are the same as those in the
C           previous group...
C.........  This array will be used in helping to decide whether pollutant
C           data need to be read for a given group
        IDVGP( 1 ) = 1       ! initialize group number for the first group
        DO N = 2, NGRP

            NEXTGRP = .FALSE.    ! reset indicator for incrementing group
            TGRP = IDVGP( N-1 )  ! set temporary group number

C.............  If any of the pollutants in this group are different from the
C               previous group, turn on indicator for incrementing group
            DO V = 1, VGRPCNT( N )
                J = SIINDEX( V,N )
                K = SIINDEX( V,N-1 )
                IF( K .LE. 0 ) THEN
                    NEXTGRP = .TRUE.
                ELSE IF( EANAM( J ) .NE. EANAM( K ) ) THEN
                    NEXTGRP = .TRUE.
                END IF
            END DO

C.............  If indicator is true, increment group
            IF( NEXTGRP ) THEN
                TGRP = TGRP + 1
            ENDIF

C.............  Set group number for current group
            IDVGP( N ) = TGRP

        END DO

C.........  Precompute the position of the species variables in the reactivity
C           matrices
        JA   = ANRMATV - ARNMSPC + 1
        JM   = MNRMATV - MRNMSPC + 1
        JP   = PNRMATV - PRNMSPC + 1
 
C.........  Create source-categeory-specific indicies for speciation, inventory
C           emissions, and control matrices...

C.........  If speciation, loop through the number of variables per
C           group and store the count for each source category.
C.........  The counts are kept for each source category because each
C           source category has only its own speciation factors stored
C           for only the pollutants in that source category.
C.........  The counts are used because of the variable groups.  If the
C           position in the variable list were stored, this would not be
C           the correct index for a pollutant within a group.  The counts
C           work because all of the pollutant/species are sorted in the same
C           order for all source categories and files.
C.........  If there is reactivity, then also store the position in the
C           reactivity matrix

        IF ( SFLAG ) THEN

C.............  Loop through all groups, then number of variables per group
            DO N = 1, NGRP
                    
                ACNT = 0 
                BCNT = 0 
                MCNT = 0 
                PCNT = 0 

                DO V = 1, VGRPCNT( N )

                    SVBUF = GVNAMES( V,N )

                    IF ( AFLAG ) THEN    ! Area sources
                        K = INDEX1( SVBUF, ANSMATV, ASVDESC )
                        IF( K .GT. 0 ) THEN
                            ACNT = ACNT + 1
                            AS_EXIST( V,N ) = ACNT
                        END IF
                    END IF

                    IF ( ARFLAG ) THEN    ! Area reactivity
                        K = INDEX1( SVBUF, ARNMSPC, ARVDESC( JA ) )
                        AR_EXIST( V,N ) = K
                    END IF

                    IF ( BFLAG ) THEN    ! Biogenic sources
                        K = INDEX1( SVBUF, BNSMATV, BSVDESC )
                        IF( K .GT. 0 ) THEN
                            BCNT = BCNT + 1
                            BS_EXIST( V,N ) = BCNT
                        END IF
                    END IF

                    IF ( MFLAG ) THEN    ! Mobile sources
                        K = INDEX1( SVBUF, MNSMATV, MSVDESC )
                        IF( K .GT. 0 ) THEN
                            MCNT = MCNT + 1
                            MS_EXIST( V,N ) = MCNT
                        END IF
                    END IF

                    IF ( MRFLAG ) THEN    ! Mobile reactivity
                        K = INDEX1( SVBUF, MRNMSPC, MRVDESC( JM ) )
                        MR_EXIST( V,N ) = K
                    END IF

                    IF ( PFLAG ) THEN    ! Point sources
                        K = INDEX1( SVBUF, PNSMATV, PSVDESC )
                        IF( K .GT. 0 ) THEN
                            PCNT = PCNT + 1
                            PS_EXIST( V,N ) = PCNT
                        END IF
                    END IF

                    IF ( PRFLAG ) THEN    ! Point reactivity
                        K = INDEX1( SVBUF, PRNMSPC, PRVDESC( JP ) )
                        PR_EXIST( V,N ) = K
                    END IF

                END DO      ! end loop on variables per group
            END DO          ! end loop on groups
        END IF              ! For speciation

C.........  Now build indices for inventory emissions...
C.........  For each source category, store per-group position for
C           inventory pollutants...

C.........  Loop through all groups, then number of variables per group
        DO N = 1, NGRP

            ACNT = 0
            MCNT = 0
            PCNT = 0
            KA   = 0  ! array
            KM   = 0  ! array
            KP   = 0  ! array

            DO V = 1, VGRPCNT( N )

                VBUF = EANAM( SIINDEX( V,N ) )

                IF ( AFLAG ) THEN    ! Area sources
                    K = INDEX1( VBUF, ANIPOL, AEINAM )
                    IF( K .GT. 0 ) THEN

C.........................  If index has already been set for the pollutant,
C                           then reuse it, otherwise update the counter and
C                           store
                        IF( KA( K ) .NE. 0 ) THEN
                            A_EXIST( V,N ) = KA( K )
                        ELSE
                            ACNT = ACNT + 1
                            A_EXIST( V,N ) = ACNT
                            KA( K ) = ACNT
                        END IF
                    END IF
                END IF

                IF ( MFLAG ) THEN    ! Mobile sources
                    K = INDEX1( VBUF, MNIPPA, MEANAM )
                    IF( K .GT. 0 ) THEN
                        IF( KM( K ) .NE. 0 ) THEN
                            M_EXIST( V,N ) = KM( K )
                        ELSE
                            MCNT = MCNT + 1
                            M_EXIST( V,N ) = MCNT
                            KM( K ) = MCNT
                        END IF
                    END IF
                END IF

                IF ( PFLAG ) THEN    ! Point sources
                    K = INDEX1( VBUF, PNIPOL, PEINAM )
                    IF( K .GT. 0 ) THEN
                        IF( KP( K ) .NE. 0 ) THEN
                            P_EXIST( V,N ) = KP( K )
                        ELSE
                            PCNT = PCNT + 1
                            P_EXIST( V,N ) = PCNT
                            KP( K ) = PCNT
                        END IF
                    END IF
                END IF
            END DO      ! end loop on pollutants per group
        END DO          ! end loop on groups

C.........  Now build indices for multiplicative controls...
C.........  For each source category, store per-group position for 
C           inventory pollutants...

C.........  Allocate helper arrays for setup of mulitplicative control matrix
        IF( AUFLAG ) THEN
            ALLOCATE( KAU( ANUMATV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'KAU', PROGNAME )
        END IF

        IF( MUFLAG ) THEN
            ALLOCATE( KMU( MNUMATV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'KMU', PROGNAME )
        END IF

        IF( PUFLAG ) THEN
            ALLOCATE( KPU( PNUMATV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'KPU', PROGNAME )
        END IF

C.........  Loop through all groups, then number of pollutants per group
        DO N = 1, NGRP

            ACNT = 0 
            MCNT = 0 
            PCNT = 0 
            IF( AUFLAG ) KAU  = 0   ! array
            IF( MUFLAG ) KMU  = 0   ! array
            IF( PUFLAG ) KPU  = 0   ! array

            DO V = 1, VGRPCNT( N )

                VBUF = EANAM( SIINDEX( V,N ) )

                IF ( AUFLAG ) THEN    ! Area sources
                    K = INDEX1( VBUF, ANUMATV, AUVNAMS )
                    IF( K .GT. 0 ) THEN
                        IF( KAU( K ) .NE. 0 ) THEN
                            AU_EXIST( V,N ) = KAU( K )
                        ELSE
                            ACNT = ACNT + 1
                            AU_EXIST( V,N ) = ACNT
                            KAU( K ) = ACNT
                        END IF
                    END IF
                END IF

                IF ( MUFLAG ) THEN    ! Mobile sources
                    K = INDEX1( VBUF, MNUMATV, MUVNAMS )
                    IF( K .GT. 0 ) THEN
                        IF( KMU( K ) .NE. 0 ) THEN
                            MU_EXIST( V,N ) = KMU( K )
                        ELSE
                            MCNT = MCNT + 1
                            MU_EXIST( V,N ) = MCNT
                            KMU( K ) = MCNT
                        END IF
                    END IF
                END IF

                IF ( PUFLAG ) THEN    ! Point sources
                    K = INDEX1( VBUF, PNUMATV, PUVNAMS )
                    IF( K .GT. 0 ) THEN
                        IF( KPU( K ) .NE. 0 ) THEN
                            PU_EXIST( V,N ) = KPU( K )
                        ELSE
                            PCNT = PCNT + 1
                            PU_EXIST( V,N ) = PCNT
                            KPU( K ) = PCNT
                        END IF
                    END IF
                END IF

            END DO      ! end loop on pollutants per group
        END DO          ! end loop on groups

C......... Deallocate local memory
    	DEALLOCATE( PPSCNT )
        IF( AUFLAG ) DEALLOCATE( KAU )
        IF( MUFLAG ) DEALLOCATE( KMU )
        IF( PUFLAG ) DEALLOCATE( KPU )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE BLDMRGIDX
