
        SUBROUTINE WRINVCHR( ENAME, SDEV, A2PFLAG )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine writes the source characteristics to the
C      I/O API and ASCII area, mobile, or point source inventory files
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 4/99
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
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

C.........  This module is required by the FileSetAPI
        USE MODFILESET
        
        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

C.........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     LBLANK
        EXTERNAL    LBLANK

C.........  SUBROUTINE ARGUMENTS and their descriptions:
        CHARACTER(*), INTENT (IN) :: ENAME   ! I/O API file name
        INTEGER     , INTENT (IN) :: SDEV    ! ASCII file unit
        LOGICAL     , INTENT (IN) :: A2PFLAG ! true: using area-to-point processing

C.........  Arrays for column formatting
        INTEGER       COLWID( MXCHRS+4 ) !  width of source info columns

        LOGICAL       LF    ( MXCHRS+4 ) !  true if column should be output

        CHARACTER*300 CHARS ( MXCHRS+4 ) !  source fields for output

C.........  Source-specific header arrays
        CHARACTER*20 :: ARHEADRS( MXARCHR3+1 ) = 
     &                                    ( / 'SMOKE Source ID     ',  
     &                                        'Cntry/St/Co FIPS    ',
     &                                        'SCC                 ',
     &                                        'Cell                ' / )

        CHARACTER*20 :: MBHEADRS( MXMBCHR3+2 ) = 
     &                                    ( / 'SMOKE Source ID     ',  
     &                                        'Cntry/St/Co FIPS    ',
     &                                        'Roadway type code   ',
     &                                        'Link ID             ',
     &                                        'Vehicle type code   ',
     &                                        'SCC                 ', 
     &                                        'Vehicle Type Name   ' / ) 

        CHARACTER*20 :: PTHEADRS( MXPTCHR3+5 ) = 
     &                                    ( / 'SMOKE Source ID     ',  
     &                                        'Cntry/St/Co FIPS    ',
     &                                        'Plant code          ', 
     &                                        'Source char 1       ', 
     &                                        'Source char 2       ',
     &                                        'Source char 3       ',  
     &                                        'Source char 4       ', 
     &                                        'Source char 5       ',
     &                                        'SCC                 ',  
     &                                        'DOE plant ID        ',  
     &                                        'Boiler code         ',  
     &                                        'Facility description' / )

C.........  Allocatable header arrays
        CHARACTER*20, ALLOCATABLE :: HDRFLDS( : )

C.........  Other local variables
        INTEGER       COLWID0          !  width of SMOKE source ID column
        INTEGER       COLMAX           !  no. of the most specific plant char 
        INTEGER       I, J, K, S       !  counters and indices
        INTEGER       IOS              !  memory allocation status
        INTEGER       L1, L2, NL       !  counters and indices
        INTEGER       M1, M2, M3, M4   !  positions after src characteristics
        INTEGER       NASCII           !  number of posible output fields
        INTEGER       NC               !  number of output fields

        CHARACTER*100 BUFFER              !  test buffer
        CHARACTER*300 OUTFMT              !  output format
        CHARACTER*300 MESG                !  message buffer

        CHARACTER*16 :: PROGNAME = 'WRINVCHR' !  program name

C***********************************************************************
C   begin body of subroutine WRINVCHR

C.........  Write the I/O API file, one variable at a time

        L1 = LEN_TRIM( ENAME )
        MESG = 'Error writing output file "' // ENAME(1:L1) // '"'

        IF ( .NOT. WRITESET( ENAME, 'IFIP', ALLFILES, 
     &                       0, 0, IFIP ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITESET( ENAME, 'TZONES', ALLFILES, 
     &                       0, 0, TZONES ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITESET( ENAME, 'TPFLAG', ALLFILES, 
     &                       0, 0, TPFLAG ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITESET( ENAME, 'INVYR', ALLFILES, 
     &                       0, 0, INVYR ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        SELECT CASE( CATEGORY )
        CASE( 'AREA' )
        
            IF( A2PFLAG ) THEN
                IF ( .NOT. WRITESET( ENAME, 'XLOCA', ALLFILES, 
     &                               0, 0, XLOCA ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
                IF ( .NOT. WRITESET( ENAME, 'YLOCA', ALLFILES, 
     &                               0, 0, YLOCA ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            END IF
            
            IF ( .NOT. WRITESET( ENAME, 'CELLID', ALLFILES, 
     &                           0, 0, CELLID ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        CASE( 'MOBILE' )

            IF ( .NOT. WRITESET( ENAME, 'IRCLAS', ALLFILES, 
     &                           0, 0, IRCLAS ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. WRITESET( ENAME, 'IVTYPE', ALLFILES, 
     &                           0, 0, IVTYPE ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. WRITESET( ENAME, 'XLOC1', ALLFILES, 
     &                           0, 0, XLOC1 ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. WRITESET( ENAME, 'YLOC1', ALLFILES, 
     &                           0, 0, YLOC1 ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. WRITESET( ENAME, 'XLOC2', ALLFILES, 
     &                           0, 0, XLOC2 ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. WRITESET( ENAME, 'YLOC2', ALLFILES, 
     &                           0, 0, YLOC2 ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        CASE( 'POINT' )

            IF ( .NOT. WRITESET( ENAME, 'ISIC', ALLFILES, 
     &                           0, 0, ISIC ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. WRITESET( ENAME, 'XLOCA', ALLFILES, 
     &                           0, 0, XLOCA ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. WRITESET( ENAME, 'YLOCA', ALLFILES, 
     &                           0, 0, YLOCA ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. WRITESET( ENAME, 'STKHT', ALLFILES, 
     &                           0, 0, STKHT ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. WRITESET( ENAME, 'STKDM', ALLFILES, 
     &                           0, 0, STKDM ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. WRITESET( ENAME, 'STKTK', ALLFILES, 
     &                           0, 0, STKTK ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. WRITESET( ENAME, 'STKVE', ALLFILES, 
     &                           0, 0, STKVE ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END SELECT

C.........  End subroutine if ASCII file is not to be written
        IF( SDEV .LE. 0 ) RETURN

C.........  Set the number of potential ASCII columns in SDEV output file
        SELECT CASE( CATEGORY )
        CASE( 'AREA' )
            NASCII = MXARCHR3
        CASE( 'MOBILE' )
            NASCII = MXMBCHR3+1
        CASE( 'POINT' )
            NASCII = MXPTCHR3+4
        END SELECT

C.........  Allocate memory for and populate the output header fields from
C           the source-specific fields
        ALLOCATE( HDRFLDS( NASCII+1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HDRFLDS', PROGNAME )
        SELECT CASE( CATEGORY )
        CASE( 'AREA' )
            HDRFLDS = ARHEADRS  ! array
        CASE( 'MOBILE' )
            HDRFLDS = MBHEADRS  ! array
        CASE( 'POINT' )
            HDRFLDS = PTHEADRS  ! array
        END SELECT
 
C.........  Set positions for outputs to ASCII file that are after the source
C           characteristics fields
        M1 = MXCHRS + 1
        M2 = MXCHRS + 2
        M3 = MXCHRS + 3
        M4 = MXCHRS + 4

C.........  Get the maximum column width for each of the columns in ASCII file
        COLWID = 0    ! array
        DO S = 1, NSRC

            DO K = 1, NCHARS   ! Loop through source characteristics
                L1 = SC_BEGP( K )
                L2 = SC_ENDP( K )
                BUFFER = ADJUSTL( CSOURC( S )( L1:L2 ) )
                J = LEN_TRIM( BUFFER )                      ! could be blank
                IF( BUFFER .NE. ' '    .AND.
     &              J .GT. COLWID( K )      ) COLWID( K ) = J 
            END DO

            SELECT CASE ( CATEGORY )
            CASE( 'MOBILE' ) 

                J = LEN_TRIM( CVTYPE( S ) )                   ! could be blank
                IF( CVTYPE( S ) .NE. ' ' .AND.
     &              J .GT. COLWID( M1 ) ) COLWID( M1 ) = J

            CASE( 'POINT' )

                J = LEN_TRIM( CSCC( S ) )                   ! could be blank
                IF( CSCC( S ) .NE. ' ' .AND.
     &              J .GT. COLWID( M1 ) ) COLWID( M1 ) = J

                J = LEN_TRIM( CORIS( S ) )                  ! could be blank
                IF( CORIS( S ) .NE. ' ' .AND.
     &              J .GT. COLWID( M2 ) ) COLWID( M2 ) = J

                J = LEN_TRIM( CBLRID( S ) )                 ! could be blank
                IF( CBLRID( S ) .NE. ' ' .AND.
     &              J .GT. COLWID( M3 ) ) COLWID( M3 ) = J

                J = LEN_TRIM( CPDESC( S ) )                 ! could be blank
                IF( CPDESC( S ) .NE. ' ' .AND.
     &              J .GT. COLWID( M4 ) ) COLWID( M4 ) = J

            END SELECT

        END DO   ! End loop on sources to get maximum column widths

        WRITE( MESG, * ) NSRC   ! Column width for source IDs
        L1 = LBLANK( MESG ) + 1
        L2 = LEN_TRIM( MESG )
        COLWID0 = L2 - L1 + 1

C.........  It is possible that a source characteristic (such as segment or
C           VTYPE) may never be used, while SCC or link is defined and also
C           considered a source characteristic.  So, make sure that if there 
C           are no gaps in the list of source characteristics, even though 
C           this means outputting blank fields.
C.........  Find the most specific defined source characteristic
        DO COLMAX = NCHARS, 2, -1
            IF( COLWID( COLMAX ) .GT. 0 ) EXIT
        END DO

        DO I = 2, COLMAX
            COLWID( I ) = MAX( COLWID( I ), 1 )
        END DO

C.........  Consider that every source might not have the same thing defined!
C.........  If a column is defined for _any_ source, then it must be
C           written for all sources.  Set flags for which are defined.
        DO K = 1, NASCII
            LF( K ) = ( COLWID( K ) .GT. 0 )  
        END DO

C.........  Set number of columns for ASCII file (initialize +1 b/c MXARCHR3 
C           does not include the SMOKE source ID column)
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

C.............  Split source characteristics into separate fields (CHARS)
            CALL PARSCSRC( CSOURC( S ), MXCHRS, SC_BEGP, SC_ENDP, LF, 
     &                     NC, CHARS )

C.............  Store remaining source attributes in separate fields (CHARS)
            SELECT CASE ( CATEGORY )          
            CASE( 'MOBILE' )
                IF( LF( M1 ) ) THEN
                    NC = NC + 1
                    CHARS( NC ) = CVTYPE( S )
                END IF

            CASE( 'POINT' ) 
                IF( LF( M1 ) ) THEN
                    NC = NC + 1
                    CHARS( NC ) = CSCC( S )
                END IF

                IF( LF( M2 ) ) THEN
                    NC = NC + 1
                    CHARS( NC ) = CORIS( S )
                END IF

                IF( LF( M3 ) ) THEN
                    NC = NC + 1
                    CHARS( NC ) = CBLRID( S )
                END IF

                IF( LF( M4 ) ) THEN
                    NC = NC + 1
                    CHARS( NC ) = CPDESC( S )
                END IF

            END SELECT

C.............  Write source characteristics and attributes
            WRITE( SDEV, OUTFMT ) S, ( CHARS( I ), I = 1, NC )
 
        END DO   ! End loop on sources for writing ASCII file

C.........  Deallocate locally allocated memory
        IF( ALLOCATED( HDRFLDS ) ) DEALLOCATE( HDRFLDS )

C.........  Close ASCII output file
        CLOSE( SDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93100   FORMAT( I2, ', "', A, '"' )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94100   FORMAT( 9( A, I2.2 ) )

        END SUBROUTINE WRINVCHR
