
        SUBROUTINE GETIDASZ( FDEV, CATEGORY, NLINEBP )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine returns an exact number of records times
C      pollutants (depending on the value of OUTTYPE) for a raw inventory input
C      file opened on unit FDEV.  The file must use the #POLID or #DATA header
C      to indicate how many data variables are contained on each line of the
C      file.  This routine will work with EMS-95 mobile files as well as
C      IDA files.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 1/99
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

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         GETNLIST

        EXTERNAL        GETNLIST

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN ) :: FDEV     !  unit number of input file
        CHARACTER(*), INTENT (IN ) :: CATEGORY !  description of source category
        INTEGER     , INTENT (OUT) :: NLINEBP   !  no. input recs (srcs x poll)

C...........   Other local variables
        INTEGER         I1, I2, L   !  counters and indices

        INTEGER         IOS         !  i/o status
        INTEGER      :: IREC    = 0 !  input line counter
        INTEGER         FILFMT      !  file format code
        INTEGER      :: NVAR    = 0 !  number of data vars at line in file

        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first time routine called

        CHARACTER(300)  BUFFER      !  temporary buffer
        CHARACTER(300)  LINE        !  input file line buffer
        CHARACTER(300)  MESG        !  message buffer

        CHARACTER(16) :: PROGNAME = 'GETIDASZ' ! program name

C***********************************************************************
C   begin body of subroutine GETIDASZ

C.........  Initialize counters
        NLINEBP  = 0

        DO   ! Head of file read loop

C.............  Read in part of line of file (enough for the header)
            READ( FDEV,93000, END=111, IOSTAT=IOS ) LINE
            IREC = IREC + 1
 
C.............  Check I/O error status
            IF( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010 )
     &                     'Error', IOS,  'reading inventory file ' // 
     &                     'as character strings at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) CYCLE

C.............  Scan for header lines
            IF( LINE( 1:1 ) .EQ. CINVHDR ) THEN

C.................  Scan for pollutant header field
                L = LEN_TRIM( LINE )
                I1 = INDEX( LINE, 'POLID' )
                I2 = INDEX( LINE, 'DATA' )

                IF( I1 .GT. 0 ) THEN
                    I1 = I1 + 5
                    BUFFER = ADJUSTL( LINE( I1:L ) )
                    L = LEN_TRIM( BUFFER )
  
                    CALL UPCASE( BUFFER )
                    NVAR = GETNLIST( L, BUFFER )

                ELSE IF( I2 .GT. 0 ) THEN
                    I2 = I2 + 4
                    BUFFER = ADJUSTL( LINE( I2:L ) )
                    L = LEN_TRIM( BUFFER )
  
                    CALL UPCASE( BUFFER )
                    NVAR = GETNLIST( L, BUFFER )

                END IF

                CYCLE  ! to end of loop

C.............  Otherwise, count the lines and lines times pollutants
            ELSE

C.................  First, check to ensure header was there for area and 
C                   point sources
C.................  For mobile sources, header might not be there, so the
C                   default value is 1.
                IF( NVAR .EQ. 0 ) THEN

                    IF( CATEGORY .EQ. 'MOBILE' ) THEN
                        NVAR = 1
                    ELSE
                        WRITE( MESG,94010 ) 
     &                             'No #POLID or #DATA header in ' //
     &                             'inventory file unit number', FDEV
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                END IF

                NLINEBP = NLINEBP + NVAR

            END IF

        END DO

111     CONTINUE  ! Exit from read loop

        IF( NLINEBP .EQ. 0 ) THEN
            MESG = 'Inventory file has no valid ' //
     &             'lines of inventory data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        REWIND( FDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
