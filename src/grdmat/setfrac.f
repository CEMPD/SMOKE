
        SUBROUTINE SETFRAC( SRCID, SRGIDX, TGTSRG, CELIDX, FIPIDX, NC, 
     &                      REPORT, CSRC, DEFSRG, SRGFLAG, OUTID1,
     &                      OUTID2, FRAC, SFLAG )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine sets the surrogate fraction from the surrogates
C      module. It tests to make sure that the surrogate for the county is not
C      zero, and if it is, it applies the default surrogate. It also evaluates
C      the default surrogate to use from the environment.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 5/99
C
C**************************************************************************
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
C***************************************************************************

C...........   MODULES for public variables
C...........   This module contains the gridding surrogates tables
        USE MODSURG, ONLY: NSRGS, SRGLIST, SRGCSUM, SRGFRAC

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: ASRGID 

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        INTEGER         ENVINT
        INTEGER         ENVYN
        INTEGER         FIND1

        EXTERNAL        CRLF, ENVINT, FIND1

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: SRCID         ! source ID
        INTEGER     , INTENT (IN) :: SRGIDX        ! surrogate index
        INTEGER     , INTENT (IN) :: TGTSRG        ! target surrogate
        INTEGER     , INTENT (IN) :: CELIDX        ! cell index
        INTEGER     , INTENT (IN) :: FIPIDX        ! FIPS code index
        INTEGER     , INTENT (IN) :: NC            ! no. src chars for msg
        LOGICAL     , INTENT (IN) :: REPORT        ! true: okay to report
        CHARACTER(*), INTENT (IN) :: CSRC          ! source chars
        INTEGER     , INTENT (IN) :: DEFSRG        ! defaul surrogate code
        LOGICAL     , INTENT (IN) :: SRGFLAG       ! true: using fallback surrogate
        INTEGER     , INTENT(OUT) :: OUTID1        ! primary srg ID
        INTEGER     , INTENT(OUT) :: OUTID2        ! secondary srg ID
        REAL        , INTENT(OUT) :: FRAC          ! surrogate fraction
        LOGICAL     , INTENT (IN) :: SFLAG         ! true: called by sizgmat, false: by gengmat

C...........   Local allocatable arrays...

        INTEGER          L2       !  indices and counters.
        INTEGER          IOS      !  i/o status

        CHARACTER(300)  BUFFER    !  source fields buffer
        CHARACTER(300)  MESG      !  message buffer 

        CHARACTER(SRCLEN3), SAVE :: CSRC2  = ' ' ! abridged CSRC
        CHARACTER(SRCLEN3), SAVE :: LCSRC  = ' ' ! prev call src chars string
        CHARACTER(SRCLEN3), SAVE :: LCSRC2 = ' ' ! prev call src chars string

        CHARACTER(16) :: PROGNAME = 'SETFRAC' ! program name

C***********************************************************************
C   begin body of subroutine SETFRAC

C.........  Create abridged name for warning messages
        IF( CSRC /= ' ' ) THEN
            CSRC2 = CSRC( 1:VIDPOS3-1 )
        ELSE
            CSRC2 = ' '
        END IF

C.........  Check if surrogate selected by cross-reference for this
C           source is non-zero in the country/state/county code of interest

        IF( SRGFLAG ) THEN 

            IF( SRGCSUM( SRGIDX,FIPIDX ) .EQ. 0. )THEN
           
C.................  Write note about changing surrogate used for current
C                   source if it has not yet been written
                IF( SFLAG ) THEN
                IF( TGTSRG .NE. DEFSRG ) THEN

                    IF( REPORT .AND. CSRC2 .NE. LCSRC2 ) THEN
                        CALL FMTCSRC( CSRC2, NC, BUFFER, L2 )
                        WRITE( MESG,94010 ) 
     &                     'WARNING: Using fallback surrogate ',DEFSRG,
     &                     ' to prevent zeroing by orig surrogate for:'
     &                     // CRLF()// BLANK10// BUFFER( 1:L2 ) //
     &                     ' Surrogate ID: ', TGTSRG
                        CALL M3MESG( MESG )
                    END IF

                END IF
                END IF

C.....................  Write warning for default fraction of zero
                IF( .NOT. SFLAG ) THEN
                IF( TGTSRG .EQ. DEFSRG ) THEN

                    CALL FMTCSRC( CSRC2, NC, BUFFER, L2 )
                    MESG ='WARNING: Zero fallback surrogate data will'//
     &                    ' cause zero emissions inside the grid for:'//
     &                    CRLF() // BLANK10 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )
           
                END IF
                END IF
           
C.................  Set surrogate fraction using default surrogate
                FRAC   = 1.0E-36
                OUTID1 = TGTSRG
                OUTID2 = DEFSRG

C.............  Set surrogate fraction with cross-reference-selected
C               surrogate
            ELSE

                FRAC   = SRGFRAC( SRGIDX, CELIDX, FIPIDX )
                OUTID1 = TGTSRG
                OUTID2 = 0

            END IF
            
        ELSE
C.........  Check if surrogate selected by cross-reference for this
C           source is non-zero in the country/state/county code of interest

            IF( SRGCSUM( SRGIDX,FIPIDX ) .LE. 0. )THEN
           
C.................  Write warning for default fraction of zero
                IF( .NOT. SFLAG ) THEN

                    CALL FMTCSRC( CSRC2, NC, BUFFER, L2 )
                    MESG ='WARNING: Zero fallback surrogate data will'//
     &                    ' cause zero emissions inside the grid for:'//
     &                    CRLF() // BLANK10 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )

                END IF

C.................  Set surrogate fraction using default surrogate
                FRAC   = 1.0E-36
                OUTID1 = TGTSRG
                OUTID2 = 0

C.............  Set surrogate fraction with cross-reference-selected
C               surrogate
            ELSE

                FRAC   = SRGFRAC( SRGIDX, CELIDX, FIPIDX )
                OUTID1 = TGTSRG
                OUTID2 = 0

            END IF

        END IF            

        LCSRC  = CSRC   ! Store source info for next iteration
        LCSRC2 = CSRC2  ! Store abridged source info

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE SETFRAC

