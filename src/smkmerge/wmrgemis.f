
        SUBROUTINE WRMRGGRD( VNAME, JDATE, JTIME )

C***********************************************************************
C  subroutine WRMRGGRD body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to write out all NetCDF files
C      coming from the merge program.  It is expected that this routine
C      will be called for each time step.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 3/99 by M. Houyoux
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

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        INTEGER         INDEX1

        EXTERNAL        CRLF , INDEX1

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: VNAME   ! variable name to output
        INTEGER     , INTENT (IN) :: JDATE   ! Julian date to output (YYYYDDD)
        INTEGER     , INTENT (IN) :: JTIME   ! time to output (HHMMSS)

C.........  Other local variables
        INTEGER         J, L, L1, L2

        LOGICAL      :: AOUTFLAG = .FALSE.  ! true: output area sources
        LOGICAL      :: MOUTFLAG = .FALSE.  ! true: output mobile sources
        LOGICAL      :: POUTFLAG = .FALSE.  ! true: output point sources

        CHARACTER(LEN=IOVLEN3) FILNAM       ! tmp logical file name
        CHARACTER*300          MESG         ! message buffer

        CHARACTER*16 :: PROGNAME = 'WRMRGGRD' ! program name

C***********************************************************************
C   begin body of subroutine WRMRGGRD

! NOTE: For now, this routine will only write a single output file. In the
!       future, it will write out multiple files.

C.........  Exit subroutine if gridded outputs has been turned off
        IF( .NOT. LGRDOUT ) RETURN

C.........  Determine the source categories that are valid for output
C.........  If the subroutine call is for speciated output, use different
C           indicator arrays for determining output or not.
        IF( SFLAG ) THEN

            IF( AFLAG ) THEN
                J = INDEX1( VNAME, ANMSPC, AEMNAM )
                AOUTFLAG = ( J .GT. 0 )
            END IF

            IF( MFLAG ) THEN
                J = INDEX1( VNAME, MNMSPC, MEMNAM )
                MOUTFLAG = ( J .GT. 0 )
            END IF

            IF( PFLAG ) THEN
                J = INDEX1( VNAME, PNMSPC, PEMNAM )
                POUTFLAG = ( J .GT. 0 )
            END IF

        ELSE

            IF( AFLAG ) THEN
                J = INDEX1( VNAME, ANIPOL, AEINAM )
                AOUTFLAG = ( J .GT. 0 )
            END IF

            IF( MFLAG ) THEN
                J = INDEX1( VNAME, MNIPOL, MEINAM )
                MOUTFLAG = ( J .GT. 0 )
            END IF

            IF( PFLAG ) THEN
                J = INDEX1( VNAME, PNIPOL, PEINAM )
                POUTFLAG = ( J .GT. 0 )
            END IF

        END IF

C.........  For area sources, output file...
        IF( AOUTFLAG ) THEN

            FILNAM = AONAME
            IF( .NOT. WRITE3( FILNAM, VNAME,
     &                        JDATE, JTIME, AEMGRD ) ) THEN

                L  = LEN_TRIM( VNAME )
                L2 = LEN_TRIM( FILNAM )
                MESG = 'Could not write "' // VNAME( 1:L ) //
     &                 '" to file "'// FILNAM( 1:L2 ) //'"'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            ENDIF

        END IF

C.........  For mobile sources, output file...
        IF( MOUTFLAG ) THEN

            FILNAM = MONAME
            IF( .NOT. WRITE3( FILNAM, VNAME,
     &                        JDATE, JTIME, MEMGRD ) ) THEN

                L  = LEN_TRIM( VNAME )
                L2 = LEN_TRIM( FILNAM )
                MESG = 'Could not write "' // VNAME( 1:L ) //
     &                 '" to file "'// FILNAM( 1:L2 ) //'"'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            ENDIF

        END IF

C.........  For point sources, output file...
        IF( POUTFLAG ) THEN

            FILNAM = PONAME
            IF( .NOT. WRITE3( FILNAM, VNAME,
     &                        JDATE, JTIME, PEMGRD ) ) THEN

                L  = LEN_TRIM( VNAME )
                L2 = LEN_TRIM( FILNAM )
                MESG = 'Could not write "' // VNAME( 1:L ) //
     &                 '" to file "'// FILNAM( 1:L2 ) //'"'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            ENDIF

        END IF

C.........  For multiple source categories, output totals file...
        IF( XFLAG ) THEN

            FILNAM = TONAME
            IF( .NOT. WRITE3( FILNAM, VNAME,
     &                        JDATE, JTIME, TEMGRD ) ) THEN

                L  = LEN_TRIM( VNAME )
                L2 = LEN_TRIM( FILNAM )
                MESG = 'Could not write "' // VNAME( 1:L ) //
     &                 '" to file "'// FILNAM( 1:L2 ) //'"'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            ENDIF

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE WRMRGGRD
