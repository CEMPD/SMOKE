
        SUBROUTINE ALLOCMRG( MXGRP, MXVARPGP )

C***********************************************************************
C  subroutine ALLOCMRG body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to allocate the memory for all of the 
C      major arrays in the merge program.
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

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C.........  SUBROUTINE ARGUMENTS and their descriptions:

        INTEGER, INTENT(OUT) :: MXGRP    ! max possible no. of processing groups
        INTEGER, INTENT(OUT) :: MXVARPGP ! maximum number of variables per group

C...........   EXTERNAL FUCNTIONS:

        CHARACTER*2    CRLF
        EXTERNAL       CRLF
   
C...........   Temporary matrix array sizes

        INTEGER         APOLSIZ ! work area inventory emissions array size
        INTEGER         MPOLSIZ ! work mobile inventory emissions array size
        INTEGER         PPOLSIZ ! work point inventory emissions array size

        INTEGER         ASPCSIZ ! work area speciation matrix array size
        INTEGER         MSPCSIZ ! work mobile speciation matrix array size
        INTEGER         PSPCSIZ ! work point speciation matrix array size

        INTEGER         AMULSIZ ! work area multipl control matrix array size
        INTEGER         MMULSIZ ! work mobile multipl control matrix array size
        INTEGER         PMULSIZ ! work point multipl control matrix array size

        INTEGER         AADDSIZ ! work area additive control matrix array size
        INTEGER         MADDSIZ ! work mobile additive control matrix array size
        INTEGER         PADDSIZ ! work point additive control matrix array size

C...........   Array of allocation statuses
        INTEGER         IOSA( 100 )

C...........   Other local variables

        INTEGER         I, J         ! counters and indices
        INTEGER         IOS          ! i/o status
        INTEGER         MXVARPGP_SAV ! saved MXVARPGP value
        INTEGER         NCNY         ! tmp no. counties
        INTEGER         NDIM         ! tmp dimension        
        INTEGER         NSTA         ! tmp no. states

        LOGICAL      :: RESET = .FALSE. ! true: mem alloc fail, try new config

        CHARACTER*300   MESG     ! message buffer

        CHARACTER*16 :: PROGNAME = 'ALLOCMRG' ! program name

C***********************************************************************
C   begin body of subroutine ALLOCMRG

       IF( LREPANY ) THEN
           NCNY = NCOUNTY  ! from modstcy
       ELSE
           NCNY = 0
       END IF

       NSTA = NSTATE   ! from modstcy

C....................................................................
C........  Allocate memory for fixed-sized arrays ...................
C....................................................................

C........  Allocate memory for fixed-size area source arrays        
       IF( AFLAG ) THEN
            J = NGRID + 2 * ANGMAT
            ALLOCATE( AGMATX( J ), STAT=IOS )      ! contiguous gridding matrix
            CALL CHECKMEM( IOS, 'AGMATX', PROGNAME )

            ALLOCATE( ARINFO( NASRC,2 ), STAT=IOS )        ! tmp ar react. data
            CALL CHECKMEM( IOS, 'ARINFO', PROGNAME )

            ALLOCATE( AEMGRD( NGRID ), STAT=IOS )  ! gridded area emissions
            CALL CHECKMEM( IOS, 'AEMGRD', PROGNAME )

            NDIM = 0
            IF( LREPINV ) NDIM = NDIM + ANIPOL
            IF( LREPSPC ) NDIM = NDIM + ANMSPC

            IF( LREPSTA ) THEN
                ALLOCATE( AEBSTA( NSTA,NDIM ), STAT=IOS )         ! state total
                CALL CHECKMEM( IOS, 'AEBSTA', PROGNAME )

                IF( LREPCTL .AND. AUFLAG ) THEN
                    ALLOCATE( AEUSTA( NSTA,NDIM ), STAT=IOS )  ! state mult tot
                    CALL CHECKMEM( IOS, 'AEUSTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. AAFLAG ) THEN
                    ALLOCATE( AEASTA( NSTA,NDIM ), STAT=IOS )   ! state add tot
                    CALL CHECKMEM( IOS, 'AEASTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. ARFLAG ) THEN
                    ALLOCATE( AERSTA( NSTA,NDIM ), STAT=IOS )  ! state reac tot
                    CALL CHECKMEM( IOS, 'AERSTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. 
     &            ( AUFLAG .OR. AAFLAG .OR. ARFLAG ) ) THEN
                    ALLOCATE( AECSTA( NSTA,NDIM ), STAT=IOS )  ! state ctrl tot
                    CALL CHECKMEM( IOS, 'AECSTA', PROGNAME )
                ENDIF

            END IF

            ALLOCATE( AEBCNY( NCNY,NDIM ), STAT=IOS )        ! county total 
            CALL CHECKMEM( IOS, 'AEBCNY', PROGNAME )

            IF( AUFLAG ) THEN
                ALLOCATE( AEUCNY( NCNY,NDIM ), STAT=IOS ) ! county mult tot 
            ELSE
                ALLOCATE( AEUCNY( 0,0 ), STAT=IOS )
            END IF
            CALL CHECKMEM( IOS, 'AEUCNY', PROGNAME )

            IF( AAFLAG ) THEN
                ALLOCATE( AEACNY( NCNY,NDIM ), STAT=IOS )  ! county add tot 
            ELSE
                ALLOCATE( AEACNY( 0,0 ), STAT=IOS )
            END IF
            CALL CHECKMEM( IOS, 'AEACNY', PROGNAME )

C.............  The reactivity controls sum is used when there is speciation,
C               even if no reactivity controls are applied to prevent IF
C               statements in the merge loops in MRGMULT.
            IF( SFLAG ) THEN 
                ALLOCATE( AERCNY( NCNY,NDIM ), STAT=IOS ) ! county reac tot 
            ELSE
                ALLOCATE( AERCNY( 0,0 ), STAT=IOS )
            END IF
            CALL CHECKMEM( IOS, 'AERCNY', PROGNAME )

            IF( AUFLAG .OR. AAFLAG .OR. SFLAG ) THEN
                ALLOCATE( AECCNY( NCNY,NDIM ), STAT=IOS ) ! county ctrl tot 
            ELSE
                ALLOCATE( AECCNY( 0,0 ), STAT=IOS )
            END IF
            CALL CHECKMEM( IOS, 'AECCNY', PROGNAME )


            IF( ARFLAG ) THEN
                ALLOCATE( ACRIDX( ANSREAC ), STAT=IOS )      ! reactivity index
                CALL CHECKMEM( IOS, 'ACRIDX', PROGNAME )
                ALLOCATE( ACRREPEM( ANSREAC ), STAT=IOS ) ! react. replace emis
                CALL CHECKMEM( IOS, 'ACRREPEM', PROGNAME )
                ALLOCATE( ACRPRJFC( ANSREAC ), STAT=IOS )  ! react. projctn fac
                CALL CHECKMEM( IOS, 'ACRPRJFC', PROGNAME )
                ALLOCATE( ACRMKTPN( ANSREAC ), STAT=IOS )   ! react. mkt pentrn
                CALL CHECKMEM( IOS, 'ACRMKTPN', PROGNAME )
                ALLOCATE( ACRFAC( ANSREAC,ANSMATV ), STAT=IOS )       ! factors
                CALL CHECKMEM( IOS, 'ACRFAC', PROGNAME )
            ENDIF

        END IF

C.........  Biogenic source fixed-size arrays
        IF( BFLAG ) THEN
            ALLOCATE( BEMGRD( NGRID ), STAT=IOS )  ! gridded biogenic emissions
            CALL CHECKMEM( IOS, 'BEMGRD', PROGNAME )

            IF( LREPSTA ) THEN
                ALLOCATE( BEBSTA( NSTA,BNMSPC ), STAT=IOS )    ! state total
                CALL CHECKMEM( IOS, 'BEBSTA', PROGNAME )
            ENDIF

            ALLOCATE( BEBCNY( NCNY,BNMSPC ), STAT=IOS )        ! county total 
            CALL CHECKMEM( IOS, 'BEBCNY', PROGNAME )

        END IF

C.........  Mobile source fixed-size arrays
        IF( MFLAG ) THEN
            J = NGRID + 2 * MNGMAT
            ALLOCATE( MGMATX( J ), STAT=IOS )      ! contiguous gridding matrix
            CALL CHECKMEM( IOS, 'MGMATX', PROGNAME )

            ALLOCATE( MRINFO( NMSRC,2 ), STAT=IOS )        ! tmp mb react. data
            CALL CHECKMEM( IOS, 'MRINFO', PROGNAME )

            ALLOCATE( MEMGRD( NGRID ), STAT=IOS )! gridded mobile emissions
            CALL CHECKMEM( IOS, 'MEMGRD', PROGNAME )

            NDIM = 0
            IF( LREPINV ) NDIM = NDIM + MNIPPA
            IF( LREPSPC ) NDIM = NDIM + MNMSPC

            IF( LREPSTA ) THEN
                ALLOCATE( MEBSTA( NSTA,NDIM ), STAT=IOS )         ! state total
                CALL CHECKMEM( IOS, 'MEBSTA', PROGNAME )

                IF( LREPCTL .AND. MUFLAG ) THEN
                    ALLOCATE( MEUSTA( NSTA,NDIM ), STAT=IOS )  ! state mult tot
                    CALL CHECKMEM( IOS, 'MEUSTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. MAFLAG ) THEN
                    ALLOCATE( MEASTA( NSTA,NDIM ), STAT=IOS )   ! state add tot
                    CALL CHECKMEM( IOS, 'MEASTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. MRFLAG ) THEN
                    ALLOCATE( MERSTA( NSTA,NDIM ), STAT=IOS )  ! state reac tot
                    CALL CHECKMEM( IOS, 'MERSTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. 
     &            ( MUFLAG .OR. MAFLAG .OR. MRFLAG ) ) THEN
                    ALLOCATE( MECSTA( NSTA,NDIM ), STAT=IOS )  ! state ctrl tot
                    CALL CHECKMEM( IOS, 'MECSTA', PROGNAME )
                ENDIF

            ENDIF

            ALLOCATE( MEBCNY( NCNY,NDIM ), STAT=IOS )        ! county total 
            CALL CHECKMEM( IOS, 'MEBCNY', PROGNAME )

            IF( MUFLAG ) THEN
                ALLOCATE( MEUCNY( NCNY,NDIM ), STAT=IOS ) ! county mult tot 
            ELSE
                ALLOCATE( MEUCNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'MEUCNY', PROGNAME )

            IF( MAFLAG ) THEN
                ALLOCATE( MEACNY( NCNY,NDIM ), STAT=IOS )  ! county add tot 
            ELSE
                ALLOCATE( MEACNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'MEACNY', PROGNAME )

            IF( SFLAG ) THEN   ! See note on AERCNY definition
                ALLOCATE( MERCNY( NCNY,NDIM ), STAT=IOS ) ! county reac tot 
            ELSE
                ALLOCATE( MERCNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'MERCNY', PROGNAME )

            IF( MUFLAG .OR. MAFLAG .OR. SFLAG ) THEN
                ALLOCATE( MECCNY( NCNY,NDIM ), STAT=IOS ) ! county ctrl tot 
            ELSE
                ALLOCATE( MECCNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'MECCNY', PROGNAME )


            IF( MRFLAG ) THEN
                ALLOCATE( MCRIDX( MNSREAC ), STAT=IOS )      ! reactivity index
                CALL CHECKMEM( IOS, 'MCRIDX', PROGNAME )
                ALLOCATE( MCRREPEM( MNSREAC ), STAT=IOS ) ! react. replace emis
                CALL CHECKMEM( IOS, 'MCRREPEM', PROGNAME )
                ALLOCATE( MCRPRJFC( MNSREAC ), STAT=IOS )  ! react. projctn fac
                CALL CHECKMEM( IOS, 'MCRPRJFC', PROGNAME )
                ALLOCATE( MCRMKTPN( MNSREAC ), STAT=IOS )   ! react. mkt pentrn
                CALL CHECKMEM( IOS, 'MCRMKTPN', PROGNAME )
                ALLOCATE( MCRFAC( MNSREAC,MNSMATV ), STAT=IOS )       ! factors
                CALL CHECKMEM( IOS, 'MCRFAC', PROGNAME )
            END IF

        END IF

C.........  Point source fixed-size arrays
        IF( PFLAG ) THEN
            J = NGRID + 2 * NPSRC
            ALLOCATE( PGMATX( J ), STAT=IOS )      ! contiguous gridding matrix
            CALL CHECKMEM( IOS, 'PGMATX', PROGNAME )

            ALLOCATE( PRINFO( NPSRC,2 ), STAT=IOS )        ! tmp pt react. data
            CALL CHECKMEM( IOS, 'PRINFO', PROGNAME )

            ALLOCATE( PEMGRD( NGRID,EMLAYS ), STAT=IOS ) ! gridded point emissions
            CALL CHECKMEM( IOS, 'PEMGRD', PROGNAME )

            NDIM = 0
            IF( LREPINV ) NDIM = NDIM + PNIPOL
            IF( LREPSPC ) NDIM = NDIM + PNMSPC

            IF( LREPSTA ) THEN
                ALLOCATE( PEBSTA( NSTA,NDIM ), STAT=IOS )         ! state total
                CALL CHECKMEM( IOS, 'PEBSTA', PROGNAME )

                IF( LREPCTL .AND. PUFLAG ) THEN
                    ALLOCATE( PEUSTA( NSTA,NDIM ), STAT=IOS )  ! state mult tot
                    CALL CHECKMEM( IOS, 'PEUSTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. PAFLAG ) THEN
                    ALLOCATE( PEASTA( NSTA,NDIM ), STAT=IOS )   ! state add tot
                    CALL CHECKMEM( IOS, 'PEASTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. PRFLAG ) THEN
                    ALLOCATE( PERSTA( NSTA,NDIM ), STAT=IOS )  ! state reac tot
                    CALL CHECKMEM( IOS, 'PERSTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. 
     &            ( PUFLAG .OR. PAFLAG .OR. PRFLAG ) ) THEN
                    ALLOCATE( PECSTA( NSTA,NDIM ), STAT=IOS )  ! state ctrl tot
                    CALL CHECKMEM( IOS, 'PECSTA', PROGNAME )
                ENDIF

            ENDIF

            ALLOCATE( PEBCNY( NCNY,NDIM ), STAT=IOS )        ! county total 
            CALL CHECKMEM( IOS, 'PEBCNY', PROGNAME )

            IF( PUFLAG ) THEN
                ALLOCATE( PEUCNY( NCNY,NDIM ), STAT=IOS ) ! county mult tot 
            ELSE
                ALLOCATE( PEUCNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'PEUCNY', PROGNAME )

            IF( PAFLAG ) THEN
                ALLOCATE( PEACNY( NCNY,NDIM ), STAT=IOS )  ! county add tot 
            ELSE
                ALLOCATE( PEACNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'PEACNY', PROGNAME )

            IF( SFLAG ) THEN   ! See note on AERCNY definition
                ALLOCATE( PERCNY( NCNY,NDIM ), STAT=IOS ) ! county reac tot 
            ELSE
                ALLOCATE( PERCNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'PERCNY', PROGNAME )

            IF( PUFLAG .OR. PAFLAG .OR. SFLAG ) THEN
                ALLOCATE( PECCNY( NCNY,NDIM ), STAT=IOS ) ! county ctrl tot 
            ELSE
                ALLOCATE( PECCNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'PECCNY', PROGNAME )

            IF( PRFLAG ) THEN
                ALLOCATE( PCRIDX( PNSREAC ), STAT=IOS )      ! reactivity index
                CALL CHECKMEM( IOS, 'PCRIDX', PROGNAME )
                ALLOCATE( PCRREPEM( PNSREAC ), STAT=IOS ) ! react. replace emis
                CALL CHECKMEM( IOS, 'PCRREPEM', PROGNAME )
                ALLOCATE( PCRPRJFC( PNSREAC ), STAT=IOS )  ! react. projctn fac
                CALL CHECKMEM( IOS, 'PCRPRJFC', PROGNAME )
                ALLOCATE( PCRMKTPN( PNSREAC ), STAT=IOS )   ! react. mkt pentrn
                CALL CHECKMEM( IOS, 'PCRMKTPN', PROGNAME )
                ALLOCATE( PCRFAC( PNSREAC,PNSMATV ), STAT=IOS )       ! factors
                CALL CHECKMEM( IOS, 'PCRFAC', PROGNAME )
            END IF

            IF( LFLAG ) THEN
                ALLOCATE( LFRAC( NPSRC,EMLAYS ), STAT=IOS )   ! layer fractions
                CALL CHECKMEM( IOS, 'LFRAC', PROGNAME )
            END IF

        END IF

C.........  Total gridded emissions.  Always allocate this, even if there
C           us only one source category because it will simplify the merging as
C           we won't have to check if it is allocated or not.

        ALLOCATE( TEMGRD( NGRID,EMLAYS ), STAT=IOS ) ! gridded out emis
        CALL CHECKMEM( IOS, 'TEMGRD', PROGNAME )
        TEMGRD = 0.

C.........  Total emissions, fixed-size arrays. 
        IF( XFLAG ) THEN

            NDIM = 0
            IF( LREPINV ) NDIM = NDIM + NIPPA
            IF( LREPSPC ) NDIM = NDIM + NMSPC

            IF( LREPSTA ) THEN
                ALLOCATE( TEBSTA( NSTA,NDIM ), STAT=IOS )         ! state total
                CALL CHECKMEM( IOS, 'TEBSTA', PROGNAME )

                IF( LREPCTL .AND. TUFLAG ) THEN
                    ALLOCATE( TEUSTA( NSTA,NDIM ), STAT=IOS )  ! state mult tot
                    CALL CHECKMEM( IOS, 'TEUSTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. TAFLAG ) THEN
                    ALLOCATE( TEASTA( NSTA,NDIM ), STAT=IOS )   ! state add tot
                    CALL CHECKMEM( IOS, 'TEASTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. TRFLAG ) THEN
                    ALLOCATE( TERSTA( NSTA,NDIM ), STAT=IOS )  ! state reac tot
                    CALL CHECKMEM( IOS, 'TERSTA', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. 
     &            ( TUFLAG .OR. TAFLAG .OR. TRFLAG ) ) THEN
                    ALLOCATE( TECSTA( NSTA,NDIM ), STAT=IOS )  ! state ctrl tot
                    CALL CHECKMEM( IOS, 'TECSTA', PROGNAME )
                ENDIF

            ENDIF

            ALLOCATE( TEBCNY( NCNY,NDIM ), STAT=IOS )        ! county total 
            CALL CHECKMEM( IOS, 'TEBCNY', PROGNAME )

            IF( TUFLAG ) THEN
                ALLOCATE( TEUCNY( NCNY,NDIM ), STAT=IOS ) ! county mult tot 
            ELSE
                ALLOCATE( TEUCNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'TEUCNY', PROGNAME )

            IF( TAFLAG ) THEN
                ALLOCATE( TEACNY( NCNY,NDIM ), STAT=IOS )  ! county add tot 
            ELSE
                ALLOCATE( TEACNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'TEACNY', PROGNAME )

            IF( SFLAG ) THEN   ! See note on AERCNY definition
                ALLOCATE( TERCNY( NCNY,NDIM ), STAT=IOS ) ! county reac tot 
            ELSE
                ALLOCATE( TERCNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'TERCNY', PROGNAME )

            IF( TUFLAG .OR. TAFLAG .OR. SFLAG ) THEN
                ALLOCATE( TECCNY( NCNY,NDIM ), STAT=IOS ) ! county ctrl tot 
            ELSE
                ALLOCATE( TECCNY( 0,0 ), STAT=IOS )
            ENDIF
            CALL CHECKMEM( IOS, 'TECCNY', PROGNAME )

        END IF

C....................................................................
C........  Allocate memory for variable-sized arrays ...................
C....................................................................

C.........  Initialize size for multiplicative control pollutants as the actual
C           number in in matrix for each source category
        AMULSIZ = MIN( ANIPOL, ANUMATV )
        MMULSIZ = MIN( MNIPOL, MNUMATV )
        PMULSIZ = MIN( PNIPOL, PNUMATV )

C.........  Initialize size for additive control pollutants as the actual
C           number in in matrix for each source category
        AADDSIZ = MIN( ANIPOL, ANAMATV )
        MADDSIZ = MIN( MNIPOL, MNAMATV )
        PADDSIZ = MIN( PNIPOL, PNAMATV )

C.........  Initialize size for species as all pol-to-species combos for each 
C           source category.
        ASPCSIZ = ANSMATV
        MSPCSIZ = MNSMATV
        PSPCSIZ = PNSMATV

C.........  Initialize size for inventory pollutants, activities, and/or
C           emission types for each source category
        APOLSIZ = ANIPOL
        MPOLSIZ = MNIPPA
        PPOLSIZ = PNIPOL

C.........  Initialize maximum number of variables per group using all
C           pollutants-to-species combinations in the input data, or all 
C           pollutants in the inventory.
        IF( SFLAG ) THEN 
            MXVARPGP = NSMATV
        ELSE
            MXVARPGP = NIPPA
        END IF
        MXVARPGP_SAV = MXVARPGP

C.........  Head of loop for allocating memory.  When speciation or certain
C           controls are not being used, the value of the variable used for
C           dimensioning will be 0, so no need to use IFs in many cases.
C.........  Note that the matricies below are allocated whether or not they are
C           used. If they are not used, the dimensions will be zero, and they
C           will not actually use memory.  The allocations are necessary
C           because the arrays are always referenced in the main program, even
C           if their data values are not accessed.
        DO

C.............  Allocate speciation matrices...
C.............  Area
            J = 1
            ASPCSIZ = MIN( ASPCSIZ, MXVARPGP )
            ALLOCATE( ASMATX( NASRC,ASPCSIZ ), STAT=IOSA( J ) )
	    CALL CHECKMEM( IOSA( J ), 'ASMATX', PROGNAME )

C.............  Mobile
            J = J + 1
            MSPCSIZ = MIN( MSPCSIZ, MXVARPGP )
            ALLOCATE( MSMATX( NMSRC,MSPCSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'MSMATX', PROGNAME )

C.............  Point
            J = J + 1
            PSPCSIZ = MIN( PSPCSIZ, MXVARPGP )
            ALLOCATE( PSMATX( NPSRC,PSPCSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'PSMATX', PROGNAME )

C.............  Allocate multiplicative control matrices
C.............  Area
            J = J + 1
            AMULSIZ = MIN( AMULSIZ, MXVARPGP )
            ALLOCATE( ACUMATX( NASRC,AMULSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'ACUMATX', PROGNAME )

C.............  Mobile
            J = J + 1
            MMULSIZ = MIN( MMULSIZ, MXVARPGP )
            ALLOCATE( MCUMATX( NMSRC,MMULSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'MCUMATX', PROGNAME )

C.............  Point
            J = J + 1
            PMULSIZ = MIN( PMULSIZ, MXVARPGP )
            ALLOCATE( PCUMATX( NPSRC,PMULSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'PCUMATX', PROGNAME )

C.............  Allocate additive control matrices
C.............  Area
            J = J + 1
            AADDSIZ = MIN( AADDSIZ, MXVARPGP )
            ALLOCATE( ACAMATX( NASRC,AADDSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'ACAMATX', PROGNAME )

C.............  Mobile
            J = J + 1
            MADDSIZ = MIN( MADDSIZ, MXVARPGP )
            ALLOCATE( MCAMATX( NMSRC,MADDSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'MCAMATX', PROGNAME )

C.............  Point
            J = J + 1
            PADDSIZ = MIN( PADDSIZ, MXVARPGP )
            ALLOCATE( PCAMATX( NPSRC,PADDSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'PCAMATX', PROGNAME )

C.............  Allocate emissions arrays
C.............  Area
            J = J + 1
            APOLSIZ = MIN( APOLSIZ, MXVARPGP )
            ALLOCATE( AEMSRC( NASRC,APOLSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'AEMSRC', PROGNAME )

C.............  Mobile
            J = J + 1
            MPOLSIZ = MIN( MPOLSIZ, MXVARPGP )
            ALLOCATE( MEMSRC( NMSRC,MPOLSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'MEMSRC', PROGNAME )

C.............  Point
            J = J + 1
            PPOLSIZ = MIN( PPOLSIZ, MXVARPGP )
            ALLOCATE( PEMSRC( NPSRC,PPOLSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'PEMSRC', PROGNAME )

C.............  Allocate emissions arrays needed for reactivity
C.............  Area
            IF( ARFLAG ) THEN
                J = J + 1
                APOLSIZ = MIN( APOLSIZ, MXVARPGP )
                ALLOCATE( AEISRC( NASRC,APOLSIZ ), STAT=IOSA( J ) )
                CALL CHECKMEM( IOSA( J ), 'AEMSRC', PROGNAME )
            END IF

C.............  Mobile
            IF( MRFLAG ) THEN
                J = J + 1
                MPOLSIZ = MIN( MPOLSIZ, MXVARPGP )
                ALLOCATE( MEISRC( NMSRC,MPOLSIZ ), STAT=IOSA( J ) )
                CALL CHECKMEM( IOSA( J ), 'MEMSRC', PROGNAME )
            END IF

C.............  Point
            IF( PRFLAG ) THEN
                J = J + 1
                PPOLSIZ = MIN( PPOLSIZ, MXVARPGP )
                ALLOCATE( PEISRC( NPSRC,PPOLSIZ ), STAT=IOSA( J ) )
                CALL CHECKMEM( IOSA( J ), 'PEMSRC', PROGNAME )
            END IF

C.............  Check IOSA values to see if any allocations failed
            DO I = 1, J
                IF( IOSA( I ) .GT. 0 ) RESET = .TRUE.
            END DO

C.............  If any memory allocations failed, reset a new size, starting
C               with speciation matrices, and then moving to control and 
C               inventory emission matrices sizes
            IF ( RESET ) THEN
                
C.................  If there is still room to make it smaller, reduce the 
C                   maximum number of variables per group
                IF( MXVARPGP .GT. 1 ) THEN

                    MXVARPGP = MXVARPGP / 2 + MOD( MXVARPGP, 2 )

C.................  Not enough memory available
                ELSE
                    MESG = 'Could not allocate enough memory to run ' //
     &                     'the program for the options selected. The'//
     &                     CRLF() // BLANK10 // 'program might run ' //
     &                     'using fewer source categories or ' //
     &                     'matrices. ' // CRLF() // BLANK10 // 
     &                     'Otherwise the computer must be '//
     &                     'reconfigured for more memory or swap space.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF 

C.................  Deallocate existing allocations to prepare for next
C                   iteration 
                DEALLOCATE( ASMATX , MSMATX , PSMATX  )
                DEALLOCATE( ACUMATX, MCUMATX, PCUMATX )
                DEALLOCATE( ACAMATX, MCAMATX, PCAMATX )
                DEALLOCATE( AEMSRC , MEMSRC , PEMSRC  )

C.............  Memory allocation was successfull
            ELSE
                EXIT  ! Exit loop

            END IF    ! End check on memory allocation

        END DO  ! End of memory allocation loop

C.........  Set the maximum number of processing groups
C.........  Consider worse case in which each pollutant-species combo is a 
C           new pollutant.  MXVARPGP has been set such that this case will be
C           handled, and the maximum number of groups can simply be set as
C           the original number of pol-to-spec or pollutants divided by the
C           current number (plus 1 if there is a remainder)

        MXGRP = MXVARPGP_SAV / MXVARPGP
        IF( MOD( MXVARPGP_SAV, MXVARPGP ) .GT. 0 ) MXGRP = MXGRP + 1

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALLOCMRG
