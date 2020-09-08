
        SUBROUTINE ASGNBINS( RCNT )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       The ASGNBINS routine is responsible for assigning a bin number to each
C       output record.  A bin is a group of records for which the emissions,
C       activity, or emission-type data will be summed to provide a single
C       record in a report.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 7/2000 by M Houyoux
C
C     Version 9/2014 by C Coats:  promote MXOUTREC to INTEGER*8 for CARB;
C     OpenMP parallel; incremental construction of SORTBUF; use SORTINC8();
C     major cleanup of post-sort data reorganization; construction of
C     binning matrices MODREPBN:<NBINS,ISRCB,GFACB>
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
C***********************************************************************

C.......   MODULES for public variables:
C...........   MODSOURC contains the inventory arrays
C............  MODLISTS contains the lists of unique source characteristics
C............  MODREPRT contains Smkreport-specific settings
C............  MODREPBN contains report arrays for each output bin
C............  MODGRID contains the global variables for the 3-d grid
C............  MODELEV contains arrays for plume-in-grid and major sources
C............  MODSTCY contains the arrays for state and county summaries
C............  MODINFO contains the information about the source category

        USE MODSOURC, ONLY: CSOURC, CIFIP, CSCC, IRCLAS, SRGID, CMON,
     &                      CWEK, CDOM, CMND, CTUE, CWED, CTHU, CFRI,
     &                      CSAT, CSUN, CMET, SPPROF, CISIC, CMACT,
     &                      CNAICS, CSRCTYP, CORIS, CINTGR, CERPTYP,
     &                      XLOCA, YLOCA, STKHT, STKDM, STKTK, STKVE,
     &                      FUGHGT, FUGWID, FUGLEN, FUGANG,
     &                      SPPNLO, SPPNHI, NSPFRC, SPPROF, NGSPRO, GSPROID

        USE MODLISTS, ONLY: NINVSCC, INVSCC, NINVSIC, INVSIC, NINVMACT,
     &                      INVMACT, NINVNAICS, INVNAICS

        USE MODREPRT, ONLY: RPT_, LREGION, AFLAG, ALLRPT, NSPCPOL,
     &                      SPCPOL, STKX, STKY, LOC_BEGP, LOC_ENDP

        USE MODREPBN, ONLY: NOUTREC, NOUTBINS, NBINS, ISRCB, ISPRO, GFACB,
     &                      OUTGFAC, OUTSPRO, OUTSFAC, BINBAD, BINCOIDX,
     &                      BINSTIDX, BINCYIDX, BINREGN, BINSMKID,
     &                      BINSCC, BINSRGID1, BINSRGID2, BINSNMIDX,
     &                      BINRCL, BINMONID, BINWEKID, BINDOMID,
     &                      BINMNDID, BINTUEID, BINWEDID, BINTHUID,
     &                      BINFRIID, BINSATID, BINSUNID, BINMETID,
     &                      BINSPCID, BINPLANT, BINX, BINY, BINELEV,
     &                      BINPOPDIV, OUTBIN, OUTCELL, OUTSRC,
     &                      BINSIC, BINSICIDX, BINMACT, BINMACIDX,
     &                      BINNAICS, BINNAIIDX, BINSRCTYP, BINORIS,
     &                      BINORSIDX, BINSTKGRP, BININTGR, BINGEO1IDX,
     &                      BINERPTYP, BINSPCIDX

        USE MODGRID, ONLY: NCOLS

        USE MODELEV, ONLY: LPING, LMAJOR, GROUPID

        USE MODSTCY, ONLY: NCOUNTRY, CTRYCOD, NSTATE, STATCOD, NCOUNTY,
     &                     CNTYCOD, CTRYPOPL, STATPOPL, CNTYPOPL,
     &                     NORIS, ORISLST, NGEOLEV1, GEOLEV1COD

        USE MODINFO, ONLY: CATEGORY

        USE MODSPRO, ONLY: NSPROF, SPROFN

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS

        INTEGER, INTENT (IN) :: RCNT    ! current report number

C...........   INCLUDE Files:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  EXTERNAL FUNCTIONS and their descriptions:

        INTEGER, EXTERNAL :: INDEX1
        INTEGER, EXTERNAL :: FIND1
        INTEGER, EXTERNAL :: FINDC
        LOGICAL, EXTERNAL :: USEEXPGEO

C...........   Local parameters:

        INTEGER, PARAMETER :: BUFLEN =  85 + SCCLEN3 + SICLEN3 + SPNLEN3
     &                                     + MACLEN3 + NAILEN3 + STPLEN3
     &                                     + ORSLEN3 + TMPLEN3 + TMPLEN3
     &                                     + TMPLEN3 + TMPLEN3 + TMPLEN3
     &                                     + TMPLEN3 + TMPLEN3 + TMPLEN3
     &                                     + TMPLEN3 + TMPLEN3 + TMPLEN3
     &                                     + ERPLEN3 + 26 + 65 + 64
        INTEGER, PARAMETER :: PTSCCLEV( NSCCLV3 ) = (/ 1, 3, 6,  8, 9 /)
        INTEGER, PARAMETER :: ARSCCLEV( NSCCLV3 ) = (/ 2, 4, 7, 10, 9 /)

        CHARACTER(1),  PARAMETER :: BLANK = ' '
        CHARACTER(16), PARAMETER :: PROGNAME = 'ASGNBINS' ! program name

C...........   Sorting arrays

        INTEGER          , ALLOCATABLE :: BOUTIDX( : )
        INTEGER(8)       , ALLOCATABLE :: SORTIDX( : )
        CHARACTER(BUFLEN), ALLOCATABLE :: SORTBUF( : )

C...........   Local variables

        INTEGER         B, C, F, I, II, IJ, IS, IV, J, K, L, LB, S

        INTEGER         COL             ! tmp column number
        INTEGER         IOS             ! i/o status
        INTEGER         NDATA           ! no. output data columns for current
        INTEGER         PREVSRCID       ! previous source ID
        INTEGER         RCL             ! tmp road class code
        INTEGER         ROW             ! tmp row number
        INTEGER         SRCID           ! tmp source ID
        INTEGER         SRGID1          ! tmp primary surrogate ID
        INTEGER         SRGID2          ! tmp fallback surrogate ID
        INTEGER         STKGRP          ! tmp stack group ID

        INTEGER(8)      M               ! format-length as INTEGER*8 for SORTINC8()
        INTEGER(8)      N               ! NOUTREC as INTEGER*8 for SORTINC8()
        INTEGER(8)      MXOUTREC        ! max output rec (NOUTREC*BUFLEN)

        LOGICAL         EFLAG

        CHARACTER       ESTAT           ! tmp elevated status
        CHARACTER(300)  MESG            ! message buffer

        CHARACTER(5)       SCCTYPE      ! tmp determination of SCC type
        CHARACTER(BUFLEN)  BUFFER       ! sorting info buffer
        CHARACTER(BUFLEN)  LBUF         ! previous sorting info buffer
        CHARACTER(SCCLEN3) SCC          ! tmp SCC
        CHARACTER(SICLEN3) SIC          ! tmp SIC
        CHARACTER(ERPLEN3) ERPTYP       ! tmp ERPTYP
        CHARACTER(INTLEN3) INTGR        ! tmp INTEGRATE
        CHARACTER(MACLEN3) MACT         ! tmp MACT
        CHARACTER(NAILEN3) NAICS        ! tmp NAICS
        CHARACTER(ORSLEN3) ORIS         ! tmp ORIS
        CHARACTER(STPLEN3) SRCTYP       ! tmp SRCTYP
        CHARACTER(PLTLEN3) PLANT        ! tmp plant ID
        CHARACTER(PLTLEN3) PREVPLT      ! previous plant ID
        CHARACTER(FIPLEN3) CFIP         ! tmp country/state/county
        CHARACTER(FIPLEN3) CCNTRY       ! tmp country
        CHARACTER(FIPLEN3) CSTA         ! tmp country/state
        CHARACTER(FIPLEN3) PREVFIP      ! previous FIPs code

C***********************************************************************
C   begin body of subroutine ASGNBINS

C.........  Set report-specific local settings
        NDATA   = ALLRPT( RCNT )%NUMDATA
        RPT_    = ALLRPT( RCNT )
        LREGION = ( RPT_%BYGEO1 .OR. RPT_%BYCNTY .OR. RPT_%BYSTAT .OR. RPT_%BYCNRY )
        N        = NOUTREC          !!  force INT*8 arithmetic
        MXOUTREC = N * BUFLEN
        B        = 0

C.........  Memory check to check exceeding integer4 maxval=2,147,483,647
        IF( MXOUTREC < 1 ) THEN
            MESG = 'ERROR: Problem processing the size of inventory'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Consistency checking:  inventory vs report

        EFLAG = .FALSE.

        IF( RPT_%BYSPC ) THEN
            IV = INDEX1( RPT_%SPCPOL, NSPCPOL, SPCPOL )
            IF ( IV .LE. 0 ) THEN
                MESG = 'INTERNAL ERROR: Pollutant "'// RPT_%SPCPOL//
     &                     'not found in list created from REPCONFIG.'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF
        END IF

        IF ( RPT_%BYORIS .AND. .NOT. ALLOCATED( CORIS ) ) THEN
            MESG = 'ERROR: BY ORIS is requested, but ' //
     &             'ORIS is not present in ASCII inventory file'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
        END IF

        IF( RPT_%BYMACT .AND. .NOT. ASSOCIATED( CMACT ) ) THEN
            MESG = 'ERROR: BY MACT is requested, but ' //
     &             'MACT is not present in ASCII inventory file'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
        END IF

        IF( RPT_%BYNAICS .AND. .NOT. ASSOCIATED( CNAICS ) ) THEN
            MESG = 'ERROR: BY NAICS is requested, but ' //
     &             'NAICS is not present in ASCII inventory file'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
        END IF

        IF( RPT_%BYSRCTYP .AND. .NOT. ASSOCIATED( CSRCTYP ) ) THEN
            MESG = 'ERROR: BY SRCTYP is requested, but ' //
     &             'SRCTYP code is not present in ASCII  inventory file'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
        END IF

        IF( RPT_%BYINTGR  .AND. .NOT. ASSOCIATED( CINTGR ) ) THEN
            MESG = 'ERROR: BY INTEGRATE is requested, but ' //
     &             'Integrate flag is not present in ASCII '//
     &             'inventory file'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
        END IF

        IF( RPT_%BYERPTYP .AND. .NOT. ALLOCATED( CERPTYP ) ) THEN
            MESG = 'ERROR: BY ERPTYP is requested, but ' //
     &             'ERPTYP code is not present in ASCII  inventory file'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
        END IF

        IF ( EFLAG ) THEN
            MESG = 'Inconsistent INVENTORY-vs-REPORT specification(s)'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate (and deallocate) memory for sorting arrays

        IF( ALLOCATED( NBINS ) ) DEALLOCATE( NBINS )
        IF( ALLOCATED( ISRCB ) ) DEALLOCATE( ISRCB )
        IF( ALLOCATED( ISPRO ) ) DEALLOCATE( ISPRO )
        IF( ALLOCATED( GFACB ) ) DEALLOCATE( GFACB )

        ALLOCATE( BOUTIDX( NOUTREC ),
     &            SORTIDX( NOUTREC ),
     &            SORTBUF( NOUTREC ),
     &            NBINS( 0:NOUTREC ),       !!  zero-based cumulative counts
     &              ISRCB( NOUTREC ),
     &              ISPRO( NOUTREC ),
     &              GFACB( NOUTREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BOUTIDX...GFACB', PROGNAME )



C.........  Create sorting arrays SORTIDX & SORTBUF for all output records
C.........  Add fields to SORTBUF in sorting-order
C.........  [Really wish this could be one parallel loop 1...NOUTREC
C.........  but the previous-plant dependency prevents this. 
C.........  parallel-loop; plant-loop; parallel-loop instead;-( ]

!$OMP   PARALLEL DO DEFAULT( SHARED ),
!$OMP&              PRIVATE( I, L, ROW, COL, SRCID, CFIP, SCC, SIC, 
!$OMP&                       SRGID1, SRGID2 ), 
!$OMP&          LASTPRIVATE( II, IJ, IS )

        DO I = 1, NOUTREC

            SORTIDX( I ) = I

            IF( RPT_%BYCELL ) THEN

                IF( .NOT. AFLAG ) THEN
                    WRITE( SORTBUF( I )( 1:8 ) , '( I8 )' ) OUTCELL( I )
                    II = 9          !!  start of next field
                    IS = 9          !!  start of SRCID field
                ELSE
                    ROW = STKY( I )
                    COL = STKX( I )
                    WRITE( SORTBUF( I )( 1:16 ) , '( 2 I8 )' ) COL, ROW
                    II = 17         !!  start of next field
                    IS = 17         !!  start of SRCID field
                END IF

            ELSE

                II = 1          !!  start of next field
                IS = 1          !!  start of SRCID field

            END IF              !!  if report-by-cell or not


            IF( RPT_%BYSRC ) THEN
                IF ( AFLAG ) THEN
                    IJ = II + 7 + FIPLEN3
                    SRCID = OUTSRC( I )
                    CFIP  = CIFIP( SRCID )
                    WRITE( SORTBUF( I )( II:IJ ), '( I8.8, A )' ) SRCID, CFIP
                    II = IJ + 1
                ELSE IF ( .NOT. RPT_%BYSIC ) THEN
                    IJ = II + 7 + FIPLEN3 + SCCLEN3
                    SRCID = OUTSRC( I )
                    CFIP  = CIFIP( SRCID )
                    SCC   = CSCC( SRCID )
                    WRITE( SORTBUF( I )( II:IJ ), '( I8.8, 2 A )' ) SRCID, CFIP, SCC
                    II = IJ + 1
                ELSE
                    IJ = II + 7 + FIPLEN3 + SCCLEN3 + SICLEN3
                    SRCID = OUTSRC( I )
                    CFIP  = CIFIP( SRCID )
                    SCC   =  CSCC( SRCID )
                    SIC   = CISIC( SRCID )
                    WRITE( SORTBUF( I )( II:IJ ), '( I8.8, 3 A )' ) SRCID, CFIP, SCC, SIC
                    II = IJ + 1
                END IF

            ELSE

C.................  Reporting by plant doesn't assume that the plant ID is unique across counties.
C.................  When building the sorting array, the SRCID for the first source of a plant/county
C.................  combination will be used for all sources at that plant. This happens later in the
C.................  code, so for now save space for the SRCID.
                IF( RPT_%BYPLANT ) THEN
                    IJ = II + 7
                    SRCID = OUTSRC( I )
                    WRITE( SORTBUF( I )( II:IJ ), '( I8.8 )' ) SRCID
                    II = IJ + 1
                END IF

                IF( RPT_%BYCNTY ) THEN
                    IJ = II + FIPLEN3 - 1
                    CFIP  = CIFIP( OUTSRC( I ) )
                    SORTBUF( I )( II:IJ ) = CFIP
                    II = IJ + 1
                ELSE IF( RPT_%BYSTAT ) THEN
                    IJ = II + FIPLEN3 - 1
                    CFIP  = CIFIP( OUTSRC( I ) ) ( 1:STALEN3 ) // '000'
                    SORTBUF( I )( II:IJ ) = CFIP
                    II = IJ + 1
                ELSE IF( RPT_%BYCNRY ) THEN
                    IJ = II + FIPLEN3 - 1
                    IF( USEEXPGEO() ) THEN
                        CFIP  = CIFIP( OUTSRC( I ) )( 1:FIPEXPLEN3 ) // '000000'
                        SORTBUF( I )( II:IJ ) = CFIP
                    ELSE
                        CFIP  = CIFIP( OUTSRC( I ) )( 1:FIPEXPLEN3+1 ) // '00000'
                        SORTBUF( I )( II:IJ ) = CFIP
                    END IF
                    II = IJ + 1
                ELSE IF( RPT_%BYGEO1 ) THEN
                    IJ = II + FIPLEN3 - 1
                    CFIP  = CIFIP( OUTSRC( I ) )( 1:3 ) // '000000000'
                    SORTBUF( I )( II:IJ ) = CFIP
                    II = IJ + 1
                END IF  ! End by county, state, or country

                IF( RPT_%BYSCC .AND. RPT_%SCCRES .NE. 4 ) THEN

                    IJ = II + SCCLEN3 - 1
                    SCC = CSCC( OUTSRC( I ) )
                    IF( SCC( 1:2 ) == '00' ) SCC = SCC(3:SCCLEN3) // '  '
                    L   = LEN_TRIM( SCC )
                    IF ( L .EQ. 10 ) THEN           !!  area
                        SCC = SCC( 1:ARSCCLEV( RPT_%SCCRES ) )
                    ELSE IF ( L .EQ. 8 ) THEN       !!  point
                        SCC = SCC( 1:PTSCCLEV( RPT_%SCCRES ) )
                    ELSE IF ( CATEGORY .EQ. 'POINT' ) THEN
                        SCC = SCC( 1:PTSCCLEV( RPT_%SCCRES ) )
                    ELSE IF ( CATEGORY .EQ. 'AREA' ) THEN
                        SCC = SCC( 1:ARSCCLEV( RPT_%SCCRES ) )
                    ELSE IF ( CATEGORY .EQ. 'MOBILE' ) THEN
                        SCC = SCC( 1:ARSCCLEV( RPT_%SCCRES ) )
                    END IF
                    SORTBUF( I )( II:IJ ) = SCC
                    II = IJ + 1

                ELSE  IF( RPT_%BYSCC ) THEN

                    IJ = II + SCCLEN3 - 1
                    SORTBUF( I )( II:IJ ) = CSCC( OUTSRC( I ) )
                    II = IJ + 1

                END IF

            END IF          !!  if report-by-source, or not


            IF( RPT_%BYSRG .AND. RPT_%SRGRES .EQ. 1 ) THEN
                IJ = II + 15
                SRGID1 = SRGID( OUTSRC( I ), 1 )
                SRGID2 = SRGID( OUTSRC( I ), 2 )
                WRITE( SORTBUF( I )( II:IJ ), '( 2 I8 )' ) SRGID1, SRGID2
                II = IJ + 1
            ELSE IF( RPT_%BYSRG ) THEN
                IJ = II + 7
                SRGID2 = SRGID( OUTSRC( I ), 2 )
                WRITE( SORTBUF( I )( II:IJ ), '( I8 )' ) SRGID2
                II = IJ + 1
            END IF          !!  if report-by-surrogate

            IF( RPT_%BYMON ) THEN
                IJ = II + TMPLEN3 - 1
                SORTBUF( I )( II:IJ ) = CMON( OUTSRC( I ) )
                II = IJ + 1
            END IF          !! if report-by-month

            IF( RPT_%BYWEK ) THEN
                IJ = II + TMPLEN3 - 1
                SORTBUF( I )( II:IJ ) = CWEK( OUTSRC( I ) )
                II = IJ + 1
            END IF          !! if report-by-week

            IF( RPT_%BYDOM ) THEN
                IJ = II + TMPLEN3 - 1
                SORTBUF( I )( II:IJ ) = CDOM( OUTSRC( I ) )
                II = IJ + 1
            END IF          !! if report-by-day-of-month

            IF( RPT_%BYMND ) THEN
                IJ = II + TMPLEN3 - 1
                SORTBUF( I )( II:IJ ) = CMND( OUTSRC( I ) )
                II = IJ + 1
            END IF          !! if report-by-Monday

            IF( RPT_%BYTUE ) THEN
                IJ = II + TMPLEN3 - 1
                SORTBUF( I )( II:IJ ) = CTUE( OUTSRC( I ) )
                II = IJ + 1
            END IF          !! if report-by-Tuesday

            IF( RPT_%BYWED ) THEN
                IJ = II + TMPLEN3 - 1
                SORTBUF( I )( II:IJ ) = CWED( OUTSRC( I ) )
                II = IJ + 1
            END IF          !! if report-by-Wednesday

            IF( RPT_%BYTHU ) THEN
                IJ = II + TMPLEN3 - 1
                SORTBUF( I )( II:IJ ) = CTHU( OUTSRC( I ) )
                II = IJ + 1
            END IF          !! if report-by-Thursday

            IF( RPT_%BYFRI ) THEN
                IJ = II + TMPLEN3 - 1
                SORTBUF( I )( II:IJ ) = CFRI( OUTSRC( I ) )
                II = IJ + 1
            END IF          !! if report-by-Friday

            IF( RPT_%BYSAT ) THEN
                IJ = II + TMPLEN3 - 1
                SORTBUF( I )( II:IJ ) = CSAT( OUTSRC( I ) )
                II = IJ + 1
            END IF          !! if report-by-Saturday

            IF( RPT_%BYSUN ) THEN
                IJ = II + TMPLEN3 - 1
                SORTBUF( I )( II:IJ ) = CSUN( OUTSRC( I ) )
                II = IJ + 1
            END IF          !! if report-by-Sunday

            IF( RPT_%BYMET ) THEN
                IJ = II + TMPLEN3 - 1
                SORTBUF( I )( II:IJ ) = CMET( OUTSRC( I ) )
                II = IJ + 1
            END IF          !! if report-by-met-based

            IF( RPT_%BYSPC ) THEN
                IJ = II + SPNLEN3 - 1
                SORTBUF( I )( II:IJ ) = OUTSPRO( I )
                II = IJ + 1
            END IF          !!  if report-by-species

        END DO      !!  end first parallel loop constructing SORTIDX and SORTBUF

        IF( RPT_%BYPLANT ) THEN

            PREVPLT   = '????????'
            PREVFIP   = '????????'
            PREVSRCID = -9999
            IJ = II + PLTLEN3 - 1

            DO I = 1, NOUTREC
                S     = OUTSRC( I )
                PLANT = CSOURC( S ) (LOC_BEGP(2):LOC_ENDP(2))
                IF ( CIFIP( S ) .EQ. PREVFIP .AND.
     &               PLANT      .EQ. PREVPLT       ) THEN
                    SRCID = PREVSRCID
                    WRITE( SORTBUF( I )( IS:IS+7 ), '( I8.8 )' ) PREVSRCID
                ELSE
                    SRCID     = S
                    PREVFIP   = CIFIP( S )
                    PREVPLT   = PLANT
                    PREVSRCID = S
                END IF
                SORTBUF( I )( II:IJ ) = PLANT
            END DO
            II = IJ + 1

        END IF          !!  if report-by-plant

        IS = II

!$OMP   PARALLEL DO DEFAULT( SHARED ),
!$OMP&              PRIVATE( I, ESTAT ),
!$OMP&         FIRSTPRIVATE( II, IJ, IS )

        DO I = 1, NOUTREC       !!  second parallel loop constructing SORTBUF
        
            II = IS

            IF ( RPT_%BYORIS ) THEN
                IJ = II + ORSLEN3 - 1
                SORTBUF( I )( II:IJ ) = CORIS( OUTSRC( I ) )
                II = IJ + 1
            END IF          !!  if report-by-oris


            IF( RPT_%BYRCL ) THEN
                IJ = II + 7
                WRITE( SORTBUF( I )( II:IJ ), '( I8 )' ) IRCLAS( OUTSRC( I ) )
                II = IJ + 1
            END IF          !!  if report-by-roadclass


            IF ( RPT_%BYELEV ) THEN
                IF( LPING( OUTSRC( I ) ) ) THEN         !!  PinG
                    ESTAT = 'P'
                ELSE IF( LMAJOR( OUTSRC( I ) ) ) THEN   !!  Elevated
                    ESTAT = 'E'
                ELSE                                    !!  Low-level
                    ESTAT = 'L'
                END IF
                SORTBUF( I )( II:II ) = ESTAT
                II = II + 1
            END IF          !!  if report-by-elevstat


            IF ( RPT_%BYMACT ) THEN
                IJ = II + MACLEN3 - 1
                SORTBUF( I )( II:IJ ) = CMACT( OUTSRC( I ) )
                II = IJ + 1
            END IF          !!  if report-by-mact


            IF( RPT_%BYNAICS ) THEN
                IJ = II + NAILEN3 - 1
                SORTBUF( I )( II:IJ ) =  CNAICS( OUTSRC( I ) )
                II = IJ + 1
            END IF          !!  if report-by-naics


            IF ( RPT_%BYSRCTYP ) THEN
                IJ = II + STPLEN3 - 1
                SORTBUF( I )( II:IJ ) =  CSRCTYP( OUTSRC( I ) )
                II = IJ + 1
            END IF          !!  if report-by-sourcetype


            IF ( RPT_%BYELEV .AND. RPT_%ELVSTKGRP ) THEN
                IJ = II + 7
                WRITE( SORTBUF( I )( II:IJ ), '( I8 )' ) GROUPID( OUTSRC( I ) )
                II = IJ + 1
            END IF


            IF ( RPT_%BYINTGR ) THEN
                IJ = II + INTLEN3 - 1
                SORTBUF( I )( II:IJ ) = CINTGR( OUTSRC( I ) )
                II = IJ + 1
            END IF          !!  if report-by-integrate

            IF ( RPT_%BYERPTYP ) THEN
                IJ = II + ERPLEN3 - 1
                SORTBUF( I )( II:IJ ) =  CERPTYP( OUTSRC( I ) )
                II = IJ + 1
            END IF          !!  if report-by-emissions-release-point-type

            IF ( RPT_%BYLATLON ) THEN
                IJ = II + 26 - 1
                WRITE( SORTBUF( I )( II:IJ ), '( F13.8, F13.8 )' ) XLOCA( OUTSRC( I ) ), YLOCA( OUTSRC( I ) )
                II = IJ + 1
            END IF          !!  if report-by-latlon


            IF ( RPT_%BYSTKPARM ) THEN
                IJ = II + 80 - 1
                WRITE( SORTBUF( I )( II:IJ ), '( F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5 )' )
     &              STKHT( OUTSRC( I ) ), STKDM( OUTSRC( I ) ), STKTK( OUTSRC( I ) ),
     &              STKVE( OUTSRC( I ) ), FUGHGT( OUTSRC( I ) ), FUGWID( OUTSRC( I ) ),
     &              FUGLEN( OUTSRC( I ) ), FUGANG( OUTSRC( I ) )
                II = IJ + 1
            END IF          !!  if report-by-stack-params

            SORTBUF( I )( II: ) = ' '

        END DO                  !!  end second parallel loop constructing SORTBUF
 

C.........  Sort sorting array
        N = NOUTREC     !!  INTEGER*8 arguments...
        M = II
        CALL SORTINC8( N, M, SORTIDX, SORTBUF )

C.........  Assign bins to output records based on sorting array
C.........  NOTE:  sequential dependency for B

        NBINS( 0 ) = 0
        LBUF = BLANK
        B    = 0

        IF( RPT_%USEGMAT ) THEN     !!  construct cumulative-count incidence matrices:

            DO I = 1, NOUTREC

                J = SORTIDX( I )
                IF( SORTBUF( J ) .NE. LBUF ) THEN
                    B = B + 1
                    LBUF = SORTBUF( J )
                    BOUTIDX( B ) = J
                END IF
                NBINS( B )  = I
                ISPRO( I )  = J      ! added for GSPRO split factors
                ISRCB( I )  =  OUTSRC( J )
                GFACB( I )  = OUTGFAC( J )
                OUTBIN( J ) = B

            END DO

        ELSE        !!  not usegmat

            DO I = 1, NOUTREC

                J = SORTIDX( I )
                IF( SORTBUF( J ) .NE. LBUF ) THEN
                    B = B + 1
                    LBUF = SORTBUF( J )
                    BOUTIDX( B ) = J
                END IF

                NBINS( B )  = I
                ISPRO( I )  = J      ! added for GSPRO split factors
                ISRCB( I )  = OUTSRC( J )
                OUTBIN( J ) = B
            END DO

        END IF      !!  if usegmat, or not

        NOUTBINS = B

        DEALLOCATE( SORTIDX, SORTBUF )

C.........  If memory is allocated for bin arrays, then deallocate

        IF( ALLOCATED( BINBAD    ) ) DEALLOCATE( BINBAD )
        IF( ALLOCATED( BINGEO1IDX) ) DEALLOCATE( BINGEO1IDX )
        IF( ALLOCATED( BINCOIDX  ) ) DEALLOCATE( BINCOIDX )
        IF( ALLOCATED( BINSTIDX  ) ) DEALLOCATE( BINSTIDX )
        IF( ALLOCATED( BINCYIDX  ) ) DEALLOCATE( BINCYIDX )
        IF( ALLOCATED( BINREGN   ) ) DEALLOCATE( BINREGN )
        IF( ALLOCATED( BINSMKID  ) ) DEALLOCATE( BINSMKID )
        IF( ALLOCATED( BINSCC    ) ) DEALLOCATE( BINSCC )
        IF( ALLOCATED( BINSIC    ) ) DEALLOCATE( BINSIC )
        IF( ALLOCATED( BININTGR  ) ) DEALLOCATE( BININTGR )
        IF( ALLOCATED( BINMACT   ) ) DEALLOCATE( BINMACT )
        IF( ALLOCATED( BINMACIDX ) ) DEALLOCATE( BINMACIDX )
        IF( ALLOCATED( BINNAICS  ) ) DEALLOCATE( BINNAICS )
        IF( ALLOCATED( BINNAIIDX ) ) DEALLOCATE( BINNAIIDX )
        IF( ALLOCATED( BINORIS   ) ) DEALLOCATE( BINORIS )
        IF( ALLOCATED( BINORSIDX ) ) DEALLOCATE( BINORSIDX )
        IF( ALLOCATED( BINSRCTYP ) ) DEALLOCATE( BINSRCTYP )
        IF( ALLOCATED( BINSRGID1 ) ) DEALLOCATE( BINSRGID1 )
        IF( ALLOCATED( BINSRGID2 ) ) DEALLOCATE( BINSRGID2 )
        IF( ALLOCATED( BINSNMIDX ) ) DEALLOCATE( BINSNMIDX )
        IF( ALLOCATED( BINRCL    ) ) DEALLOCATE( BINRCL )
        IF( ALLOCATED( BINMONID  ) ) DEALLOCATE( BINMONID )
        IF( ALLOCATED( BINWEKID  ) ) DEALLOCATE( BINWEKID )
        IF( ALLOCATED( BINDOMID  ) ) DEALLOCATE( BINDOMID )
        IF( ALLOCATED( BINMNDID  ) ) DEALLOCATE( BINMNDID )
        IF( ALLOCATED( BINTUEID  ) ) DEALLOCATE( BINTUEID )
        IF( ALLOCATED( BINWEDID  ) ) DEALLOCATE( BINWEDID )
        IF( ALLOCATED( BINTHUID  ) ) DEALLOCATE( BINTHUID )
        IF( ALLOCATED( BINFRIID  ) ) DEALLOCATE( BINFRIID )
        IF( ALLOCATED( BINSATID  ) ) DEALLOCATE( BINSATID )
        IF( ALLOCATED( BINSUNID  ) ) DEALLOCATE( BINSUNID )
        IF( ALLOCATED( BINMETID  ) ) DEALLOCATE( BINMETID )
        IF( ALLOCATED( BINSPCID  ) ) DEALLOCATE( BINSPCID )
        IF( ALLOCATED( BINSPCIDX ) ) DEALLOCATE( BINSPCIDX )
        IF( ALLOCATED( BINPLANT  ) ) DEALLOCATE( BINPLANT )
        IF( ALLOCATED( BINX      ) ) DEALLOCATE( BINX )
        IF( ALLOCATED( BINY      ) ) DEALLOCATE( BINY )
        IF( ALLOCATED( BINELEV   ) ) DEALLOCATE( BINELEV )
        IF( ALLOCATED( BINSTKGRP ) ) DEALLOCATE( BINSTKGRP )
        IF( ALLOCATED( BINPOPDIV ) ) DEALLOCATE( BINPOPDIV )
        IF( ALLOCATED( BINERPTYP ) ) DEALLOCATE( BINERPTYP )

C.........  Allocate memory for bins

        ALLOCATE( BINBAD( NOUTBINS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BINBAD', PROGNAME )
        BINBAD = 0    ! array

        IF( RPT_%BYGEO1NAM ) THEN
            ALLOCATE( BINGEO1IDX ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINGEO1IDX', PROGNAME )
        ENDIF
        IF( RPT_%BYCONAM ) THEN
            ALLOCATE( BINCOIDX ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINCOIDX', PROGNAME )
        ENDIF
        IF( RPT_%BYSTNAM ) THEN
            ALLOCATE( BINSTIDX ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSTIDX', PROGNAME )
        ENDIF
        IF( RPT_%BYCYNAM ) THEN
            ALLOCATE( BINCYIDX ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINCYIDX', PROGNAME )
        ENDIF
        IF( LREGION      ) THEN
            ALLOCATE( BINREGN  ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINREGN', PROGNAME )
        ENDIF
        IF( RPT_%BYSRC .OR. RPT_%BYPLANT ) THEN
            ALLOCATE( BINSMKID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSMKID', PROGNAME )
        ENDIF
        IF( RPT_%BYSCC   ) THEN
            ALLOCATE( BINSCC   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSCC', PROGNAME )
        ENDIF
        IF( RPT_%SCCNAM  ) THEN
            ALLOCATE( BINSNMIDX( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSNMIDX', PROGNAME )
        ENDIF
        IF( RPT_%BYSIC   ) THEN
            ALLOCATE( BINSIC   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSIC', PROGNAME )
        ENDIF
        IF( RPT_%SICNAM   ) THEN
            ALLOCATE( BINSICIDX( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSICIDX', PROGNAME )
        ENDIF

        IF( RPT_%BYINTGR   ) THEN
            ALLOCATE( BININTGR ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BININTGR', PROGNAME )
        ENDIF

        IF( RPT_%BYMACT   ) THEN
            ALLOCATE( BINMACT   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINMACT', PROGNAME )
        ENDIF

        IF( RPT_%MACTNAM   ) THEN
            ALLOCATE( BINMACIDX( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINMACIDX', PROGNAME )
        ENDIF

        IF( RPT_%BYNAICS   ) THEN
            ALLOCATE( BINNAICS   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINNAICS', PROGNAME )
        ENDIF

        IF( RPT_%NAICSNAM   ) THEN
            ALLOCATE( BINNAIIDX( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINNAIIDX', PROGNAME )
        ENDIF

        IF( RPT_%BYORIS   ) THEN
            ALLOCATE( BINORIS   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINORIS', PROGNAME )
        ENDIF

        IF( RPT_%ORISNAM   ) THEN
            ALLOCATE( BINORSIDX( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINORSIDX', PROGNAME )
        ENDIF

        IF( RPT_%BYSRCTYP   ) THEN
            ALLOCATE( BINSRCTYP   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSRCTYP', PROGNAME )
        ENDIF

        IF( RPT_%SRGRES .EQ. 1 ) THEN
            ALLOCATE( BINSRGID1( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSRGID1', PROGNAME )
        ENDIF
        IF( RPT_%SRGRES .GE. 1 ) THEN
            ALLOCATE( BINSRGID2( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSRGID2', PROGNAME )
        ENDIF
        IF( RPT_%BYMON   ) THEN
            ALLOCATE( BINMONID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINMONID', PROGNAME )
        ENDIF
        IF( RPT_%BYWEK   ) THEN
            ALLOCATE( BINWEKID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINWEKID', PROGNAME )
        ENDIF
        IF( RPT_%BYDOM   ) THEN
            ALLOCATE( BINDOMID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINDOMID', PROGNAME )
        ENDIF
        IF( RPT_%BYMND   ) THEN
            ALLOCATE( BINMNDID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINMNDID', PROGNAME )
        ENDIF
        IF( RPT_%BYTUE   ) THEN
            ALLOCATE( BINTUEID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINTUEID', PROGNAME )
        ENDIF
        IF( RPT_%BYWED   ) THEN
            ALLOCATE( BINWEDID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINWEDID', PROGNAME )
        ENDIF
        IF( RPT_%BYTHU   ) THEN
            ALLOCATE( BINTHUID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINTHUID', PROGNAME )
        ENDIF
        IF( RPT_%BYFRI   ) THEN
            ALLOCATE( BINFRIID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINFRIID', PROGNAME )
        ENDIF
        IF( RPT_%BYSAT   ) THEN
            ALLOCATE( BINSATID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSUNID', PROGNAME )
        ENDIF
        IF( RPT_%BYSUN   ) THEN
            ALLOCATE( BINSUNID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSUNID', PROGNAME )
        ENDIF
        IF( RPT_%BYMET   ) THEN
            ALLOCATE( BINMETID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINMETID', PROGNAME )
        ENDIF
        IF( RPT_%BYSPC   ) THEN
            ALLOCATE( BINSPCID ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSPCID', PROGNAME )
        ENDIF
        IF( RPT_%GSPRONAM   ) THEN
            ALLOCATE( BINSPCIDX( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSPCIDX', PROGNAME )
        ENDIF

        IF( RPT_%BYPLANT ) THEN
            ALLOCATE( BINPLANT ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINPLANT', PROGNAME )
        ENDIF
        IF( RPT_%BYRCL   ) THEN
            ALLOCATE( BINRCL   ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINRCL', PROGNAME )
        ENDIF
        IF( RPT_%BYCELL  ) THEN
            ALLOCATE( BINX     ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINX', PROGNAME )
        ENDIF
        IF( RPT_%BYCELL  ) THEN
            ALLOCATE( BINY     ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINY', PROGNAME )
        ENDIF
        IF( RPT_%BYELEV  ) THEN
            ALLOCATE( BINELEV  ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINELEV', PROGNAME )
        ENDIF
        IF( RPT_%ELVSTKGRP ) THEN
            ALLOCATE( BINSTKGRP  ( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINSTKGRP', PROGNAME )
        ENDIF
        IF( RPT_%BYERPTYP ) THEN
            ALLOCATE( BINERPTYP( NOUTBINS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BINERPTYP', PROGNAME )
        ENDIF

        ALLOCATE( BINPOPDIV( NOUTBINS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BINPOPDIV', PROGNAME )

C.........  Populate the bin characteristic arrays (not the data array)

!$OMP   PARALLEL DO DEFAULT( SHARED ),
!$OMP&              PRIVATE( B, J, C, K, ROW, COL, S, CFIP, SCC,
!$OMP&                       CCNTRY, CSTA )

        DO B = 1, NOUTBINS

            BINPOPDIV( B )  = 1.0

            J = BOUTIDX( B )

            IF( RPT_%BYCELL .AND. .NOT.AFLAG ) THEN
                C = OUTCELL( J ) - 1
                ROW = 1 + INT( C / NCOLS )
                COL = 1 + MOD( C , NCOLS )
            ELSE IF( RPT_%BYCELL ) THEN
                ROW = STKY( J )
                COL = STKX( J )
            END IF
            S    = OUTSRC( J )
            CFIP =  CIFIP( S )
            SCC  =   CSCC( S )

            IF( RPT_%BYSRC )     BINSMKID( B )  = S
            IF( RPT_%BYSCC )       BINSCC( B )  = SCC
            IF( RPT_%BYSIC )       BINSIC( B )  =   CISIC( S )
            IF( RPT_%BYINTGR )   BININTGR( B )  =  CINTGR( S )
            IF( RPT_%BYMACT  )    BINMACT( B )  =   CMACT( S )
            IF( RPT_%BYNAICS )   BINNAICS( B )  =  CNAICS( S )
            IF( RPT_%BYORIS  )    BINORIS( B )  =   CORIS( S )
            IF( RPT_%BYSRCTYP ) BINSRCTYP( B )  = CSRCTYP( S )
            IF( RPT_%BYMON )     BINMONID( B )  =    CMON( S )
            IF( RPT_%BYWEK )     BINWEKID( B )  =    CWEK( S )
            IF( RPT_%BYDOM )     BINDOMID( B )  =    CDOM( S )
            IF( RPT_%BYMND )     BINMNDID( B )  =    CMND( S )
            IF( RPT_%BYTUE )     BINTUEID( B )  =    CTUE( S )
            IF( RPT_%BYWED )     BINWEDID( B )  =    CWED( S )
            IF( RPT_%BYTHU )     BINTHUID( B )  =    CTHU( S )
            IF( RPT_%BYFRI )     BINFRIID( B )  =    CFRI( S )
            IF( RPT_%BYSAT )     BINSATID( B )  =    CSAT( S )
            IF( RPT_%BYSUN )     BINSUNID( B )  =    CSUN( S )
            IF( RPT_%BYMET )     BINMETID( B )  =    CMET( S )
            IF( RPT_%BYSPC )     BINSPCID( B )  = OUTSPRO( J )
            IF( RPT_%BYERPTYP )  BINERPTYP( B ) = CERPTYP( S )

            IF( LREGION ) THEN
                IF( RPT_%BYCNTY ) THEN
                    ! no change, report uses full FIPS code
                ELSE IF( RPT_%BYSTAT ) THEN
                    CFIP = CFIP( 1:STALEN3 ) // '000'
                ELSE IF( RPT_%BYCNRY ) THEN
                    IF( USEEXPGEO() ) THEN
                        CFIP = CFIP( 1:FIPEXPLEN3 ) // '000000'
                    ELSE
                        CFIP = CFIP( 1:FIPEXPLEN3+1 ) // '00000'
                    END IF
                ELSE IF( RPT_%BYGEO1 ) THEN
                    CFIP = CFIP( 1:3 ) // '000000000'
                END IF
                BINREGN( B ) = CFIP
            END IF

            IF( RPT_%BYSRG .AND. RPT_%SRGRES .EQ. 1 ) THEN
                BINSRGID1( B ) = SRGID( S, 1 )
                BINSRGID2( B ) = SRGID( S, 2 )
            ELSE IF( RPT_%BYSRG ) THEN
                BINSRGID2( B ) = SRGID( S, 2 )
            END IF

C.................  Store geocode level 1 index.

            IF( RPT_%BYGEO1NAM ) THEN
                K = FINDC( CFIP( 1:3 ) // '000000000', NGEOLEV1, GEOLEV1COD )
                BINGEO1IDX( B ) = K
            ENDIF

C.................  Store country name index. Note for population that some
C                   form of "by region is required"

            IF( RPT_%BYCONAM ) THEN
                IF( USEEXPGEO() ) THEN
                    CCNTRY = CFIP( 1:FIPEXPLEN3 ) // '000000'
                ELSE
                    CCNTRY = CFIP( 1:FIPEXPLEN3+1 ) // '00000'
                END IF
                K = FINDC( CCNTRY, NCOUNTRY, CTRYCOD )
                BINCOIDX( B ) = K

C.....................  If using population normalization, initialize with
C                       country population

                IF ( RPT_%NORMPOP ) THEN
                    IF( CTRYPOPL( K ) . GT. 0. ) THEN
                        BINPOPDIV( B ) = 1. / CTRYPOPL(K)

C.........................  If population data unavailable, then flags bins
C                           that will not be able to have normalization by pop.
                    ELSE
                        BINBAD   ( B ) = 100  ! cpde for bad pop
                        BINPOPDIV( B ) = 1.
                    END IF
                END IF

            END IF

C.................  Store state name index. Note for population that some
C                   form of "by region is required"

            IF( RPT_%BYSTNAM ) THEN
                CSTA = CFIP( 1:STALEN3 ) // '000'
                K    = FINDC( CSTA, NSTATE, STATCOD )
                BINSTIDX( B ) = K

C.....................  If using population normalization, reset with state population

                IF ( RPT_%NORMPOP ) THEN
                    IF( STATPOPL( K ) . GT. 0. ) THEN
                        BINPOPDIV( B )= 1. / STATPOPL(K)
C.........................  If population data unavailable, then flags bins
C                           that will not be able to have normalization by pop.
                    ELSE
                        BINBAD   ( B ) = 100  ! code for bad pop
                        BINPOPDIV( B ) = 1.
                    END IF
                END IF

            END IF

C.................  Store county name index. Note for population that some
C                   form of "by region is required"

            IF( RPT_%BYCYNAM ) THEN
                K = FINDC( CFIP, NCOUNTY, CNTYCOD )
                BINCYIDX( B ) = K

C.....................  If using population normalization, reset with
C                       county population

                IF ( RPT_%NORMPOP ) THEN
                    IF( CNTYPOPL( K ) . GT. 0. ) THEN
                        BINPOPDIV( B )= 1. / CNTYPOPL(K)
C.........................  If population data unavailable, then flags bins
C                           that will not be able to have normalization by pop.
                    ELSE
                        BINBAD   ( B ) = 100  ! code for bad pop
                        BINPOPDIV( B ) = 1.
                    END IF
                END IF

            END IF

C.................  Store SCC, SIC, MACT, NAICS, ORIS name index
C.................  (for full name, regardless of SCC truncation.
C.................  Note: have confirmed that using the OUTSRC(J)
C.................  index with CSCC maps to SCC from BUFFER properly)

            IF( RPT_%SCCNAM ) THEN
                K = FINDC( CSCC( S ), NINVSCC, INVSCC )
                BINSNMIDX( B ) = K
            END IF

            IF( RPT_%GSPRONAM ) THEN
                K = INDEX1( BINSPCID( B ), NGSPRO, GSPROID )
                BINSPCIDX( B ) = K
            END IF

            IF( RPT_%SICNAM ) THEN
                K = FINDC( BINSIC( B ), NINVSIC, INVSIC )
                BINSICIDX( B ) = K
            END IF

            IF( RPT_%MACTNAM ) THEN
                K = FINDC( BINMACT( B ), NINVMACT, INVMACT )
                BINMACIDX( B ) = K
            END IF

            IF( RPT_%NAICSNAM ) THEN
                K = FINDC( BINNAICS( B ), NINVNAICS, INVNAICS )
                BINNAIIDX( B ) = K
            END IF

            IF( RPT_%ORISNAM ) THEN
                K = FINDC( BINORIS( B ), NORIS, ORISLST )
                BINORSIDX( B ) = K
            END IF

C.................  Store plant ID code
            IF( RPT_%BYPLANT ) THEN
                BINPLANT( B ) = CSOURC( S )( LOC_BEGP(2) : LOC_ENDP(2) )
                BINSMKID( B ) = S           !! Needed for plant names
            END IF

C.................  Store x-cell and y-cell
            IF( RPT_%BYCELL ) BINX( B ) = COL
            IF( RPT_%BYCELL ) BINY( B ) = ROW

C.................  Store Elevated status
            IF( RPT_%BYELEV ) THEN
                IF( LPING( S ) ) THEN       ! PinG
                    BINELEV( B ) = 'P'
                ELSE IF( LMAJOR( S ) ) THEN ! Elevated
                    BINELEV( B ) = 'E'
                ELSE                                  ! Low-level
                    BINELEV( B ) = 'L'
                END IF
                IF( RPT_%ELVSTKGRP ) BINSTKGRP( B ) = GROUPID( S )
            END IF

        END DO   ! End loop I over output records

        DEALLOCATE( BOUTIDX )

        RETURN

        END SUBROUTINE ASGNBINS

