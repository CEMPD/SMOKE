
        SUBROUTINE RDSSUP ( PDEV )

C***********************************************************************
C  subroutine body starts at line  84
C
C  DESCRIPTION:
C       Read new-format SSUP files and places content into MODSOURC
C       sparse matrix profiles-and-fractins data structures
C       <NSPFRC,SPPOLS(:),SPPNS(:,:),SPPROF(:),SPFRAC(:)>
C
C  PRECONDITIONS REQUIRED:
C       setenv [AMP]SSUP   <path>
C
C  REVISION  HISTORY:
C     Created 34/2020 by Carlie J. Coats, Jr.
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
C***********************************************************************

        USE M3UTILIO

C...........   MODULES for public variables
C.........  This module contains the speciation-profiles matrix, among other things.
        USE MODSOURC, ONLY : NSPFRC, SPPNLO, SPPNHI, SPPROF, SPFRAC

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CRL, NSRC, NIPPA

C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: RPT_, NSPCPOL, SPCPOL, LSPCPOL

        USE MODSPRO, ONLY: MXSPEC

        IMPLICIT NONE

C.........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C..........  External Function
        INTEGER     GETFLINE
        EXTERNAL    GETFLINE

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: PDEV   ! unit no.: speciation supplemental

C.........  Local variables
        INTEGER             IOS, IS1, IS2
        CHARACTER(IOVLEN3)  ::   CBUF

        CHARACTER(512) ::   LINE
        CHARACTER(512) ::   MESG

        INTEGER     IREC, PCNT
        INTEGER     K, L, L1, L2, N, N1, N2, M, I, S, PV, V

        CHARACTER(16), PARAMETER :: PNAME = 'RDSSUP' ! subroutine name

C***********************************************************************
C   begin body of subroutine RDSSUP

C.........  Open SSUP speciation supplemental file
        MESG = 'Supplemental speciation file'
        N = GETFLINE( PDEV, MESG )

        IREC   = 0
        NSPFRC = 0

        DO  I = 1, N

            READ( PDEV, '(A)', IOSTAT=IOS ) LINE 
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 )
     &               'I/O error', IOS, 'reading supplemental ' //
     &                'speciation file at line', IREC
                CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            END IF

C.............  See if this line is a pollutant name
            L1 = INDEX( LINE, '"' )        ! Find start quote

C.............  If pollutant name, figure out which pollutant index and
C               reset source counter to 0.
            IF ( L1 .GT. 0 ) THEN

                L2 = LEN_TRIM( LINE )     ! Find end quote
                CBUF = LINE( L1+1:L2-1 )

C.................  Check if this pollutant is one selected for reporting
                V = INDEX1( CBUF, NSPCPOL, SPCPOL )
                IF ( V .GT. 0 ) LSPCPOL( V ) = .TRUE.

C.............  If not pollutant name, then continue to read in the
C                   pollutant codes and store them by source
            ELSE IF ( V .GT. 0 ) THEN

                N1 = INDEX( LINE, 'NFRAC=' )
                N2 = N1 + LEN( 'NFRAC=' )
                READ( LINE(N2: ), *, IOSTAT=IS1 ) PCNT
                IF( N1 > 0 ) THEN
                    S = S + 1
                    NSPFRC = NSPFRC + PCNT
                END IF

                IF( S .EQ. NSRC ) V = 0

            END IF

        END DO

        !!........  Allocate output arrays

        ALLOCATE( SPPNLO( NSRC ),
     &            SPPNHI( NSRC ),
     &            SPPROF( NSPFRC ),
     &            SPFRAC( NSPFRC ), STAT = IOS )
        CALL CHECKMEM( IOS, 'SPPNLO...SPFRAC', PNAME )
        SPFRAC = 1.0

        !!........  Read/Build output arrays for fractions, fractional profiles
        REWIND( PDEV )

        V    = 0
        M    = 0
        S    = 0
        IREC = 0

        DO I = 1, N

            READ( PDEV, '(A)', IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 )
     &               'I/O error', IOS, 'reading supplemental ' //
     &                'speciation file at line', IREC
                CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            END IF

C.............  See if this line is a pollutant name
            L1 = INDEX( LINE, '"' )        ! Find start quote

C.............  If pollutant name, figure out which pollutant index and
C               reset source counter to 0.
            IF ( L1 .GT. 0 ) THEN

                L2 = LEN_TRIM( LINE )     ! Find end quote
                CBUF = LINE( L1+1:L2-1 )

C.................  Check if this pollutant is one selected for reporting
                V = INDEX1( CBUF, NSPCPOL, SPCPOL )
                IF ( V .GT. 0 ) THEN
                    LSPCPOL( V ) = .TRUE.
                    S = 0
                END IF

C.............  If not pollutant name, then continue to read in the
C                   pollutant codes and store them by source
            ELSE IF ( V .GT. 0 ) THEN

                N1 = INDEX( LINE, 'NFRAC=' )
                N2 = N1 + LEN( 'NFRAC=' )
                READ( LINE(N2: ), *, IOSTAT=IS1 ) PCNT
                IF( N1 > 0 ) THEN
                    S = S + 1
                    SPPNLO( S ) = M + 1
                    SPPNHI( S ) = M + PCNT
                    M = M + PCNT

                ELSE
                    READ( LINE, * ) ( SPPROF( K ), SPFRAC( K ), K = SPPNLO(S), SPPNHI(S) )

                END IF

                IF( S == NSRC+1 ) V = 0

            END IF

        END DO

        REWIND( PDEV )

        RETURN

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE RDSSUP
