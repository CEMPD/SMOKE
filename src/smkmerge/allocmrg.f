
        SUBROUTINE ALLOCMRG( MXGRP, MXPOLPGP, MXSPPOL )

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

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C.........  SUBROUTINE ARGUMENTS and their descriptions:

        INTEGER, INTENT(OUT) :: MXGRP    ! max possible no. of processing groups
        INTEGER, INTENT(OUT) :: MXPOLPGP ! maximum number of pols per group
        INTEGER, INTENT(OUT) :: MXSPPOL  ! maximum number of species per pol

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
        INTEGER         MXSPPOL_SAV  ! saved MXSPPOL value
        INTEGER         MXPOLPGP_SAV ! saved MXPOLPGP value
        INTEGER         NDIM         ! tmp dimension

        LOGICAL      :: RESET = .FALSE. ! true: mem alloc fail, try new config

        CHARACTER*300   MESG     ! message buffer

        CHARACTER*16 :: PROGNAME = 'ALLOCMRG' ! program name

C***********************************************************************
C   begin body of subroutine ALLOCMRG

C....................................................................
C........  Allocate memory for fixed-sized arrays ...................
C....................................................................

C........  Allocate memory for fixed-size area source arrays        
       IF( AFLAG ) THEN
            J = NGRID + 2 * ANGMAT
            ALLOCATE( AGMATX( J ), STAT=IOS )      ! contiguous gridding matrix
            CALL CHECKMEM( IOS, 'AGMATX', PROGNAME )
            ALLOCATE( AIFIP( NASRC ), STAT=IOS )   ! country/state/county codes
            CALL CHECKMEM( IOS, 'AIFIP', PROGNAME )

            ALLOCATE( ARINFO( NASRC,2 ), STAT=IOS )        ! tmp ar react. data
            CALL CHECKMEM( IOS, 'ARINFO', PROGNAME )

            IF( LGRDOUT ) THEN
                ALLOCATE( AEMGRD( NGRID ), STAT=IOS )  ! gridded area emissions
                CALL CHECKMEM( IOS, 'AEMGRD', PROGNAME )
            ENDIF

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

            ENDIF

            IF( LREPCNY ) THEN
                ALLOCATE( AEBCNY( NCNY,NDIM ), STAT=IOS )        ! county total 
                CALL CHECKMEM( IOS, 'AEBCNY', PROGNAME )

                IF( LREPCTL .AND. AUFLAG ) THEN
                    ALLOCATE( AEUCNY( NCNY,NDIM ), STAT=IOS ) ! county mult tot 
                    CALL CHECKMEM( IOS, 'AEUCNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. AAFLAG ) THEN
                    ALLOCATE( AEACNY( NCNY,NDIM ), STAT=IOS )  ! county add tot 
                    CALL CHECKMEM( IOS, 'AEACNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. ARFLAG ) THEN
                    ALLOCATE( AERCNY( NCNY,NDIM ), STAT=IOS ) ! county reac tot 
                    CALL CHECKMEM( IOS, 'AERCNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. 
     &            ( AUFLAG .OR. AAFLAG .OR. ARFLAG ) ) THEN
                    ALLOCATE( AECCNY( NCNY,NDIM ), STAT=IOS ) ! county ctrl tot 
                    CALL CHECKMEM( IOS, 'AECCNY', PROGNAME )
                ENDIF

            ENDIF

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
                ALLOCATE( BEBSTA( NSTA,NDIM ), STAT=IOS )         ! state total
                CALL CHECKMEM( IOS, 'BEBSTA', PROGNAME )
            ENDIF

            IF( LREPCNY ) THEN
                ALLOCATE( BEBCNY( NCNY,NDIM ), STAT=IOS )        ! county total 
                CALL CHECKMEM( IOS, 'BEBCNY', PROGNAME )
            ENDIF

        END IF

C.........  Mobile source fixed-size arrays
        IF( MFLAG ) THEN
            J = NGRID + 2 * MNGMAT
            ALLOCATE( MGMATX( J ), STAT=IOS )      ! contiguous gridding matrix
            CALL CHECKMEM( IOS, 'MGMATX', PROGNAME )
            ALLOCATE( MIFIP( NMSRC ), STAT=IOS )   ! country/state/county codes
            CALL CHECKMEM( IOS, 'MIFIP', PROGNAME )

            ALLOCATE( MRINFO( NMSRC,2 ), STAT=IOS )        ! tmp mb react. data
            CALL CHECKMEM( IOS, 'MRINFO', PROGNAME )

            IF( LGRDOUT ) THEN
                ALLOCATE( MEMGRD( NGRID ), STAT=IOS )! gridded mobile emissions
                CALL CHECKMEM( IOS, 'MEMGRD', PROGNAME )
            ENDIF

            NDIM = 0
            IF( LREPINV ) NDIM = NDIM + MNIPOL
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

            IF( LREPCNY ) THEN
                ALLOCATE( MEBCNY( NCNY,NDIM ), STAT=IOS )        ! county total 
                CALL CHECKMEM( IOS, 'MEBCNY', PROGNAME )

                IF( LREPCTL .AND. MUFLAG ) THEN
                    ALLOCATE( MEUCNY( NCNY,NDIM ), STAT=IOS ) ! county mult tot 
                    CALL CHECKMEM( IOS, 'MEUCNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. MAFLAG ) THEN
                    ALLOCATE( MEACNY( NCNY,NDIM ), STAT=IOS )  ! county add tot 
                    CALL CHECKMEM( IOS, 'MEACNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. MRFLAG ) THEN
                    ALLOCATE( MERCNY( NCNY,NDIM ), STAT=IOS ) ! county reac tot 
                    CALL CHECKMEM( IOS, 'MERCNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. 
     &            ( MUFLAG .OR. MAFLAG .OR. MRFLAG ) ) THEN
                    ALLOCATE( MECCNY( NCNY,NDIM ), STAT=IOS ) ! county ctrl tot 
                    CALL CHECKMEM( IOS, 'MECCNY', PROGNAME )
                ENDIF

            ENDIF

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
            ALLOCATE( PIFIP( NPSRC ), STAT=IOS )   ! country/state/county codes
            CALL CHECKMEM( IOS, 'PIFIP', PROGNAME )

            ALLOCATE( PRINFO( NPSRC,2 ), STAT=IOS )        ! tmp pt react. data
            CALL CHECKMEM( IOS, 'PRINFO', PROGNAME )

            IF( LGRDOUT ) THEN
                ALLOCATE( PEMGRD( NGRID ), STAT=IOS ) ! gridded point emissions
                CALL CHECKMEM( IOS, 'PEMGRD', PROGNAME )
            ENDIF

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

            IF( LREPCNY ) THEN
                ALLOCATE( PEBCNY( NCNY,NDIM ), STAT=IOS )        ! county total 
                CALL CHECKMEM( IOS, 'PEBCNY', PROGNAME )

                IF( LREPCTL .AND. PUFLAG ) THEN
                    ALLOCATE( PEUCNY( NCNY,NDIM ), STAT=IOS ) ! county mult tot 
                    CALL CHECKMEM( IOS, 'PEUCNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. PAFLAG ) THEN
                    ALLOCATE( PEACNY( NCNY,NDIM ), STAT=IOS )  ! county add tot 
                    CALL CHECKMEM( IOS, 'PEACNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. PRFLAG ) THEN
                    ALLOCATE( PERCNY( NCNY,NDIM ), STAT=IOS ) ! county reac tot 
                    CALL CHECKMEM( IOS, 'PERCNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. 
     &            ( PUFLAG .OR. PAFLAG .OR. PRFLAG ) ) THEN
                    ALLOCATE( PECCNY( NCNY,NDIM ), STAT=IOS ) ! county ctrl tot 
                    CALL CHECKMEM( IOS, 'PECCNY', PROGNAME )
                ENDIF

            ENDIF

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

        IF( LGRDOUT ) THEN
            ALLOCATE( TEMGRD( NGRID,EMLAYS ), STAT=IOS ) ! gridded out emis
            CALL CHECKMEM( IOS, 'TEMGRD', PROGNAME )
        ENDIF

C.........  Total emissions, fixed-size arrays. 
        IF( XFLAG ) THEN

            NDIM = 0
            IF( LREPINV ) NDIM = NDIM + NIPOL
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

            IF( LREPCNY ) THEN
                ALLOCATE( TEBCNY( NCNY,NDIM ), STAT=IOS )        ! county total 
                CALL CHECKMEM( IOS, 'TEBCNY', PROGNAME )

                IF( LREPCTL .AND. TUFLAG ) THEN
                    ALLOCATE( TEUCNY( NCNY,NDIM ), STAT=IOS ) ! county mult tot 
                    CALL CHECKMEM( IOS, 'TEUCNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. TAFLAG ) THEN
                    ALLOCATE( TEACNY( NCNY,NDIM ), STAT=IOS )  ! county add tot 
                    CALL CHECKMEM( IOS, 'TEACNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. TRFLAG ) THEN
                    ALLOCATE( TERCNY( NCNY,NDIM ), STAT=IOS ) ! county reac tot 
                    CALL CHECKMEM( IOS, 'TERCNY', PROGNAME )
                ENDIF

                IF( LREPCTL .AND. 
     &            ( TUFLAG .OR. TAFLAG .OR. TRFLAG ) ) THEN
                    ALLOCATE( TECCNY( NCNY,NDIM ), STAT=IOS ) ! county ctrl tot 
                    CALL CHECKMEM( IOS, 'TECCNY', PROGNAME )
                ENDIF

            ENDIF

        END IF

C....................................................................
C........  Allocate memory for variable-sized arrays ...................
C....................................................................

C.........  Initialize size for species as all pol-to-species combos for each 
C           source category.
        ASPCSIZ = ANSMATV
        MSPCSIZ = MNSMATV
        PSPCSIZ = PNSMATV

C.........  Initialize maximum number of species per pollutant using all
C           pollutant-species combinations in the inventory because when this
C           value is ultimately used, it will be to dimension and index
C           across all species in the inventory.
        MXSPPOL = NSMATV
        MXSPPOL_SAV = MXSPPOL

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

C.........  Initialize size for inventory pollutants as all pollutants for each 
C           source category
        APOLSIZ = ANIPOL
        MPOLSIZ = MNIPOL
        PPOLSIZ = PNIPOL

C.........  Initialize maximum number of pollutants per group using all
C           pollutants in the inventory because when this value is ultimately 
C           used, it will be to dimension and indexacross all pollutants in the
C           inventory.
        MXPOLPGP = NIPOL
        MXPOLPGP_SAV = MXPOLPGP

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
            ASPCSIZ = MIN( ASPCSIZ, MXSPPOL )
            ALLOCATE( ASMATX( NASRC,ASPCSIZ ), STAT=IOSA( J ) )
	    CALL CHECKMEM( IOSA( J ), 'ASMATX', PROGNAME )

C.............  Mobile
            J = J + 1
            MSPCSIZ = MIN( MSPCSIZ, MXSPPOL )
            ALLOCATE( MSMATX( NMSRC,MSPCSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'MSMATX', PROGNAME )

C.............  Point
            J = J + 1
            PSPCSIZ = MIN( PSPCSIZ, MXSPPOL )
            ALLOCATE( PSMATX( NPSRC,PSPCSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'PSMATX', PROGNAME )

C.............  Allocate multiplicative control matrices
C.............  Area
            J = J + 1
            AMULSIZ = MIN( AMULSIZ, MXPOLPGP )
            ALLOCATE( ACUMATX( NASRC,AMULSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'ACUMATX', PROGNAME )

C.............  Mobile
            J = J + 1
            MMULSIZ = MIN( MMULSIZ, MXPOLPGP )
            ALLOCATE( MCUMATX( NMSRC,MMULSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'MCUMATX', PROGNAME )

C.............  Point
            J = J + 1
            PMULSIZ = MIN( PMULSIZ, MXPOLPGP )
            ALLOCATE( PCUMATX( NPSRC,PMULSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'PCUMATX', PROGNAME )

C.............  Allocate additive control matrices
C.............  Area
            J = J + 1
            AADDSIZ = MIN( AADDSIZ, MXPOLPGP )
            ALLOCATE( ACAMATX( NASRC,AADDSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'ACAMATX', PROGNAME )

C.............  Mobile
            J = J + 1
            MADDSIZ = MIN( MADDSIZ, MXPOLPGP )
            ALLOCATE( MCAMATX( NMSRC,MADDSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'MCAMATX', PROGNAME )

C.............  Point
            J = J + 1
            PADDSIZ = MIN( PADDSIZ, MXPOLPGP )
            ALLOCATE( PCAMATX( NPSRC,PADDSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'PCAMATX', PROGNAME )

C.............  Allocate emissions arrays
C.............  Area
            J = J + 1
            APOLSIZ = MIN( APOLSIZ, MXPOLPGP )
            ALLOCATE( AEMSRC( NASRC,APOLSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'AEMSRC', PROGNAME )

C.............  Mobile
            J = J + 1
            MPOLSIZ = MIN( MPOLSIZ, MXPOLPGP )
            ALLOCATE( MEMSRC( NMSRC,MPOLSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'MEMSRC', PROGNAME )

C.............  Point
            J = J + 1
            PPOLSIZ = MIN( PPOLSIZ, MXPOLPGP )
            ALLOCATE( PEMSRC( NPSRC,PPOLSIZ ), STAT=IOSA( J ) )
            CALL CHECKMEM( IOSA( J ), 'PEMSRC', PROGNAME )

C.............  Allocate emissions arrays needed for reactivity
C.............  Area
            IF( ARFLAG ) THEN
                J = J + 1
                APOLSIZ = MIN( APOLSIZ, MXPOLPGP )
                ALLOCATE( AEISRC( NASRC,APOLSIZ ), STAT=IOSA( J ) )
                CALL CHECKMEM( IOSA( J ), 'AEMSRC', PROGNAME )
            END IF

C.............  Mobile
            IF( MRFLAG ) THEN
                J = J + 1
                MPOLSIZ = MIN( MPOLSIZ, MXPOLPGP )
                ALLOCATE( MEISRC( NMSRC,MPOLSIZ ), STAT=IOSA( J ) )
                CALL CHECKMEM( IOSA( J ), 'MEMSRC', PROGNAME )
            END IF

C.............  Point
            IF( PRFLAG ) THEN
                J = J + 1
                PPOLSIZ = MIN( PPOLSIZ, MXPOLPGP )
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
C                   maximum number of pol-to-spec combinations, and make sure
C                   that both odd and even initial numbers will go to 1
                IF( MXSPPOL .GT. 1 ) THEN

                    MXSPPOL  = MXSPPOL / 2 + MOD( MXSPPOL, 2 )
                    MXPOLPGP = MIN( MXPOLPGP, MXSPPOL )

C.................  If there is still room to make it smaller, reduce the 
C                   maximum number of pollutants
                ELSEIF( MXPOLPGP .GT. 1 ) THEN

                    MXPOLPGP = MXPOLPGP / 2 + MOD( MXPOLPGP, 2 )

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
C           new pollutant.  MXPOLPGP has been set such that this case will be
C           handled, and the maximum number of groups can simply be set as
C           the original number of pol-to-spec or pollutants divided by the
C           current number (plus 1 if there is a remainder)
        IF( SFLAG ) THEN
            MXGRP = MXSPPOL_SAV / MXSPPOL
            IF( MOD( MXSPPOL_SAV, MXSPPOL ) .GT. 0 ) MXGRP = MXGRP + 1

        ELSE
            MXGRP = MXPOLPGP_SAV / MXPOLPGP
            IF( MOD( MXPOLPGP_SAV, MXPOLPGP ) .GT. 0 ) MXGRP=MXGRP + 1

        END IF


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ALLOCMRG
