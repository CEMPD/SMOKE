
        INTEGER FUNCTION GETIDASZ( FDEV, CATEGORY, OUTTYPE )

C***********************************************************************
C  function body starts at line 
C
C  DESCRIPTION:
C      This function returns an exact number of records or records times
C      pollutants (depending on the value of OUTTYPE) for a raw IDA input
C      file opened on unit FDEV.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 1/99
C
C****************************************************************************/
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         GETNLIST

        EXTERNAL        GETNLIST

C...........   SUBROUTINE ARGUMENTS
        INTEGER       FDEV        !  unit number of input file
        CHARACTER*(*) CATEGORY    !  description of source category
        INTEGER       OUTTYPE     !  type of output of the function:
                                  !    1 = input records; 2 = records times pols

C...........   Other local variables
        INTEGER         I, L   !  counters and indices

        INTEGER         IOS         !  i/o status
        INTEGER      :: IREC    = 0 !  input line counter
        INTEGER         FILFMT      !  file format code
        INTEGER, SAVE:: NLINES  = 0 !  number of lines
        INTEGER, SAVE:: NLINEBP = 0 !  number of lines times pollutants
        INTEGER      :: NPOL    = 0 !  number of pollutants at line in file

        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first time routine called

        CHARACTER*300   BUFFER      !  temporary buffer
        CHARACTER*300   LINE        !  input file line buffer
        CHARACTER*300   MESG        !  message buffer

        CHARACTER*16 :: PROGNAME = 'GETIDASZ' ! program name

C***********************************************************************
C   begin body of function GETIDASZ

        IF( FIRSTIME ) THEN

            FIRSTIME = .FALSE.

            DO   ! Head of file read loop

C.................  Read in part of line of file (enough for the header)
                READ( FDEV,93000, END=111, IOSTAT=IOS ) LINE
                IREC = IREC + 1
 
C.................  Check I/O error status
                IF( IOS .GT. 0 ) THEN
                    WRITE( MESG, 94010 )
     &                     'Error', IOS,  'reading IDA input file ' // 
     &                     'as character strings at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

C.................  Skip blank lines
                IF( LINE .EQ. ' ' ) CYCLE

C.................  Scan for header lines
                IF( LINE( 1:1 ) .EQ. '#' ) THEN

C.....................  Scan for pollutant header field
                    L = LEN_TRIM( LINE )
                    I = INDEX( LINE, 'POLID' )

                    IF( I .GT. 0 ) THEN
                        I = I + 5
                        BUFFER = LINE( I:L )
                        L = L - I - 1
  
                        CALL UPCASE( BUFFER )
                        NPOL = GETNLIST( L, BUFFER )

                    END IF

                    CYCLE  ! to end of loop

C.................  Otherwise, count the lines and lines times pollutants
                ELSE

C.....................  First, check to ensure header was there for area and 
C                       point sources
C.....................  For mobile sources, header will not be there, so set
C                       NPOL to 1 for the VMT
                    IF( NPOL .EQ. 0 ) THEN

                        IF( CATEGORY .EQ. 'MOBILE' ) THEN
                            NPOL = 1
                        ELSE
                            WRITE( MESG,94010 ) 
     &                             'No #POLID header in IDA file', FDEV
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        END IF

                    END IF

                    NLINES  = NLINES  + 1
                    NLINEBP = NLINEBP + NPOL

                END IF

            END DO

111         CONTINUE  ! Exit from read loop

            IF( NLINES .EQ. 0 ) THEN
                MESG = 'IDA-formatted inventory file has no valid ' //
     &                 'lines of inventory data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END IF  ! end of first time

        REWIND( FDEV )

        IF( OUTTYPE .EQ. 1 ) THEN
            GETIDASZ = NLINES

        ELSE IF( OUTTYPE .EQ. 2 ) THEN
            GETIDASZ = NLINEBP

        ELSE
            MESG = 'INTERNAL ERROR: Bad output type in call to ' // 
     &             'subroutine ' // PROGNAME     
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
