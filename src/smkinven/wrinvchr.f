
        SUBROUTINE WPNTSCHR( ENAME, SDEV )

C***********************************************************************
C  subroutine body starts at line 123
C
C  DESCRIPTION:
C      This subroutine writes the point source characteristics to the
C      I/O API and ASCII point source inventory files
C
C  PRECONDITIONS REQUIRED:
C      Logical file ENAME opened
C      File unit number SDEV opened
C      Correct number of sources NSRC
C      Output arrays populated
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines, PARSCSRC
C      Functions: I/O API functions
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 11/98
C
C****************************************************************************/
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     LBLANK
        EXTERNAL    LBLANK

C.........   LOCAL PARAMETERS and their descriptions:
        INTEGER, PARAMETER :: NASCII = 10 !  no. of columns in the ASCII file

C.........  SUBROIUTINE ARGUMENTS and their descriptions:
        CHARACTER(*), INTENT (IN) :: ENAME  !  I/O API file name
        INTEGER     , INTENT (IN) :: SDEV   !  ASCII file unit

C.........  Arrays for column formatting
        INTEGER       COLWID( NASCII ) ! width of source info columns

        LOGICAL       LF    ( NASCII ) !  true if column should be output

        CHARACTER*300 CHARS ( NASCII ) !  source fields for output

        CHARACTER*20 :: HDRFLDS( NASCII+1 ) = ( 
     &                                      / 'SMOKE Source ID     ',  
     &                                      'Cntry/St/Co FIPS    ',
     &                                      'Plant code          ', 
     &                                      'Source char 1       ', 
     &                                      'Source char 2       ',
     &                                      'Source char 3       ',  
     &                                      'Source char 4       ', 
     &                                      'Source char 5       ',
     &                                      'SCC                 ',  
     &                                      'Boiler code         ',  
     &                                      'Facility description' / )

C.........  Other local variables
        INTEGER       COLWID0          !  width of SMOKE source ID column
        INTEGER       COLMAX           !  no. of the most specific plant char 
        INTEGER       I, J, K, S       !  counters and indices
        INTEGER       L1, L2, NL       !  counters and indices
        INTEGER       NC               !  number of output fields

        CHARACTER*100 BUFFER              !  test buffer
        CHARACTER*300 OUTFMT              !  output format
        CHARACTER*300 MESG                !  message buffer

        CHARACTER*16 :: PROGNAME = 'WPNTSCHR' !  program name

C***********************************************************************
C   begin body of subroutine WPNTSCHR

C.........  Get the maximum column width for each of the columns in ASCII file
        COLWID = 0    ! array
        DO S = 1, NSRC

            DO K = 1, 7   ! Loop through source characteristics
                L1 = PTBEGL3( K )
                L2 = PTENDL3( K )
                BUFFER = ADJUSTL( CSOURC( S )( L1:L2 ) )
                J = LEN_TRIM( BUFFER )                      ! could be blank
                IF( BUFFER .NE. ' '    .AND.
     &              J .GT. COLWID( K )      ) COLWID( K ) = J 
            ENDDO

            J = LEN_TRIM( CSCC( S ) )                       ! could be blank
            IF( CSCC( S ) .NE. ' ' .AND.
     &          J .GT. COLWID( 8 ) ) COLWID( 8 ) = J

            J = LEN_TRIM( CBLRID( S ) )                     ! could be blank
            IF( CBLRID( S ) .NE. ' ' .AND.
     &          J .GT. COLWID( 9 ) ) COLWID( 9 ) = J

            J = LEN_TRIM( CPDESC( S ) )                     ! could be blank
            IF( CPDESC( S ) .NE. ' ' .AND.
     &          J .GT. COLWID( 10 ) ) COLWID( 10 ) = J

        ENDDO   ! End loop on sources to get maximum column widths

        WRITE( MESG, * ) NSRC   ! Column width for source IDs
        L1 = LBLANK( MESG ) + 1
        L2 = LEN_TRIM( MESG )
        COLWID0 = L2 - L1 + 1

C.........  It is possible that a source characteristic (such as segment) may
C           never be used, while SCC is defined and also considered a plant 
C           characteristic.  So, make sure that if there are no gaps in the
C           list of source characteristics, even though this means outputting
C           blank fields.
C.........  Find the most specific defined plant characteristic
        DO COLMAX = 7, 3, -1
            IF( COLWID( COLMAX ) .GT. 0 ) EXIT
        END DO

        DO I = 3, COLMAX
            COLWID( I ) = MAX( COLWID( I ), 1 )
        END DO

C.........  Consider that every source might not have the same thing defined!
C.........  If a column is defined for _any_ source, then it must be
C           written for all sources.  Set flags for which are defined.
        DO K = 1, NASCII
            LF( K ) = ( COLWID( K ) .GT. 0 )  
        ENDDO

C.........  Set number of columns for ASCII file (initialize +1 b/c NASCII does
C           not include the SMOKE source ID column)
        NL = NASCII + 1
        DO K = 1, NASCII
            IF( .NOT. LF( K ) ) NL = NL - 1
        ENDDO

C.........  Store the format statement in OUTFMT for well-formatted output
        WRITE( OUTFMT, 94100 ) '(I', COLWID0

        DO K = 1, NASCII
            L1 = LEN_TRIM( OUTFMT )
            IF( LF( K ) ) THEN
               WRITE( MESG, 94100 ) ',1X,A', COLWID( K )
               L2 = LEN_TRIM( MESG )
               OUTFMT = OUTFMT( 1:L1 ) // MESG( 1:L2 )
            ENDIF
        ENDDO
        L1 = LEN_TRIM( OUTFMT )
        OUTFMT = OUTFMT( 1:L1 ) // ')'
              
C.........  Write the ASCII file header
        WRITE( SDEV,93100 ) NL, OUTFMT( 1:LEN_TRIM( OUTFMT ) )
        WRITE( SDEV,93000 ) HDRFLDS( 1 ) 
        DO K = 1, NASCII
            IF( LF( K ) ) WRITE( SDEV,93000 ) HDRFLDS( K+1 )
        ENDDO

C.........  Write the ASCII file data
        DO S = 1, NSRC

            CALL PARSCSRC( CSOURC( S ), LF, CHARS, NC )

            IF( LF( 8 ) ) THEN
                NC = NC + 1
                CHARS( NC ) = CSCC( S )
            ENDIF

            IF( LF( 9 ) ) THEN
                NC = NC + 1
                CHARS( NC ) = CBLRID( S )
            ENDIF

            IF( LF( 10 ) ) THEN
                NC = NC + 1
                CHARS( NC ) = CPDESC( S )
            ENDIF
                      
            WRITE( SDEV, OUTFMT ) S, ( CHARS( I ), I = 1, NC )
 
        ENDDO   ! End loop on sources for writing ASCII file

C.........  Write the I/O API file, one variable at a time

        MESG = 'Error writing output file "' // ENAME // '"'

        IF ( .NOT. WRITE3( ENAME, 'IFIP', 0, 0, IFIP ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'ISIC', 0, 0, ISIC ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'IORIS', 0, 0, IORIS ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'TZONES', 0, 0, TZONES ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'TPFLAG', 0, 0, TPFLAG ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'INVYR', 0, 0, INVYR ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'XLOCA', 0, 0, XLOCA ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'YLOCA', 0, 0, YLOCA ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'STKHT', 0, 0, STKHT ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'STKDM', 0, 0, STKDM ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'STKTK', 0, 0, STKTK ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'STKVE', 0, 0, STKVE ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93100   FORMAT( I2, ', "', A, '"' )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94100   FORMAT( 9( A, I2.2 ) )

        END SUBROUTINE WPNTSCHR
