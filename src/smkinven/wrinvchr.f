
        SUBROUTINE WPNTSCHR( ENAME, SDEV, NPSRC, IFIP, ISCC, ISIC, 
     &                       IORIS, TZONES, TPFLAG, INVYR, XLOCA, 
     &                       YLOCA, STKHT, STKDM, STKTK, STKVE, 
     &                       CBLRID, CPDESC, CSOURC )

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
C      Correct number of sources NPSRC
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

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         LBLANK
        INTEGER         TRIMLEN

        EXTERNAL LBLANK, TRIMLEN

C...........   LOCAL PARAMETERS and their descriptions:
        INTEGER, PARAMETER :: NASCII = 10 !  no. of columns in the ASCII file

C.........  Sorted list of point sources for SMOKE inventory file
        CHARACTER(LEN=NAMLEN3) ENAME  !  I/O API file name
        INTEGER      SDEV             !  ASCII file unit
        INTEGER      NPSRC            !  actual source count
        INTEGER      IFIP  ( NPSRC )  !  source FIPS (county) ID
        INTEGER      ISCC  ( NPSRC )  !  source SCC
        INTEGER      ISIC  ( NPSRC )  !  source SIC
        INTEGER      IORIS ( NPSRC )  !  source ORIS ID code
        INTEGER      TZONES( NPSRC )  !  time zones
        INTEGER      TPFLAG( NPSRC )  !  temporal profile types
        INTEGER      INVYR ( NPSRC )  !  inventory year for this record
        INTEGER      IDXSCC( NPSRC )  !  subscript table for SORTI1() for SCC
        REAL         XLOCA ( NPSRC )  !  UTM X-location (m)
        REAL         YLOCA ( NPSRC )  !  UTM Y-location (m)
        REAL         STKHT ( NPSRC )  !  stack height   (m)
        REAL         STKDM ( NPSRC )  !  stack diameter (m)
        REAL         STKTK ( NPSRC )  !  exhaust temperature (deg K)
        REAL         STKVE ( NPSRC )  !  exhaust velocity    (m/s)

        CHARACTER(LEN=BLRLEN3) CBLRID( NPSRC )  !  boiler ID
        CHARACTER(LEN=DSCLEN3) CPDESC( NPSRC )  !  plant description
        CHARACTER(LEN=SRCLEN3) CSOURC( NPSRC )  !  concatonated source chars

C...........   Other local variables
        INTEGER       COLWID0
        INTEGER       COLWID( 9 )
        INTEGER       I, J, K
        INTEGER       L1, L2, NL
        INTEGER       NC               !  number of output fields
        INTEGER       S

        LOGICAL       LF( 9 )          !  true if column should be output

        CHARACTER*100 BUFFER           !  test buffer
        CHARACTER*300 CHARS( 9 )       !  source fields for output
        CHARACTER*300 OUTFMT           !  output format
        CHARACTER*300 MESG             !  message buffer

        CHARACTER*20 :: HDRFLDS( 10 ) = ( / 'SMOKE Source ID     ',  
     &                                      'Cntry/St/Co FIPS    ',
     &                                      'Plant code          ', 
     &                                      'Source char 1       ', 
     &                                      'Source char 2       ',
     &                                      'Source char 3       ',  
     &                                      'Source char 4       ', 
     &                                      'Source char 5       ',  
     &                                      'Boiler code         ',  
     &                                      'Facility description' / )

        CHARACTER*16 :: PROGNAME = 'WPNTSCHR' !  program name

C***********************************************************************
C   begin body of subroutine WPNTSCHR

C.........  Get the maximum column width for each of the columns in ASCII file
        COLWID = 0    ! array
        DO S = 1, NPSRC

            DO K = 1, 7
                L1 = LENARR3( K )
                L2 = LENARR3( K+1 ) - 1
                BUFFER = ADJUSTL( CSOURC( S )( L1:L2 ) )
                J = TRIMLEN( BUFFER )                      ! could be blank
                IF( BUFFER .NE. ' '    .AND.
     &              J .GT. COLWID( K )      ) COLWID( K ) = J 
            ENDDO

            J = TRIMLEN( CBLRID( S ) )                     ! could be blank
            IF( CBLRID( S ) .NE. ' ' .AND.
     &          J .GT. COLWID( 8 ) ) COLWID( 8 ) = J

            J = TRIMLEN( CPDESC( S ) )                     ! could be blank
            IF( CPDESC( S ) .NE. ' ' .AND.
     &          J .GT. COLWID( 9 ) ) COLWID( 9 ) = J

        ENDDO   ! End loop on sources to get maximum column widths

        WRITE( MESG, * ) NPSRC   ! Column width for source IDs
        L1 = LBLANK( MESG ) + 1
        L2 = TRIMLEN( MESG )
        COLWID0 = L2 - L1 + 1

C.........  Consider that every source might not have the same thing defined!
C.........  If a column is defined for _any_ source, then it must be
C           written for all sources.  Set flags for which are defined.
        DO K = 1, 9
            LF( K ) = ( COLWID( K ) .GT. 0 )  
        ENDDO

C.........  Set number of columns for ASCII file
        NL = NASCII
        DO K = 1, 9
            IF( .NOT. LF( K ) ) NL = NL - 1
        ENDDO

C.........  Write the ASCII file header
        WRITE( SDEV,93100 ) NL
        WRITE( SDEV,93000 ) HDRFLDS( 1 ) 
        DO K = 1, 9
            IF( LF( K ) ) WRITE( SDEV,93000 ) HDRFLDS( K+1 )
        ENDDO

C.........  Store the format statement in OUTFMT for well-formatted output
        WRITE( OUTFMT, 94100 ) '(I', COLWID0

        DO K = 1, 9
            L1 = TRIMLEN( OUTFMT )
            IF( LF( K ) ) THEN
               WRITE( MESG, 94100 ) ',1X,A', COLWID( K )
               L2 = TRIMLEN( MESG )
               OUTFMT = OUTFMT( 1:L1 ) // MESG( 1:L2 )
            ENDIF
        ENDDO
        L1 = TRIMLEN( OUTFMT )
        OUTFMT = OUTFMT( 1:L1 ) // ')'
              
C.........  Write the ASCII file data
        DO S = 1, NPSRC

            CALL PARSCSRC( CSOURC( S ), LF, CHARS, NC )

            IF( LF( 8 ) ) THEN
                NC = NC + 1
                CHARS( NC ) = CBLRID( S )
            ENDIF

            IF( LF( 9 ) ) THEN
                NC = NC + 1
                CHARS( NC ) = CPDESC( S )
            ENDIF
                      
            WRITE( SDEV, OUTFMT ) S, ( CHARS( I ), I=1,NC )
 
        ENDDO   ! End loop on sources for writing ASCII file

C.........  Write the I/O API file, one variable at a time

        MESG = 'Error writing output file "' //
     &         ENAME( 1:TRIMLEN( ENAME ) ) // '"'

        IF ( .NOT. WRITE3( ENAME, 'IFIP', 0, 0, IFIP ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'ISCC', 0, 0, ISCC ) ) THEN
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

93100   FORMAT( I2 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94100   FORMAT( 9( A, I2.2 ) )

        END
