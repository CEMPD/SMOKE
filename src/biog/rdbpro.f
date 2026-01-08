
        SUBROUTINE RDBPRO( FDEV, SPPRO, NIPOL, EINAM, MSPCS, EMSPC, 
     &                     MLFAC, MSFAC )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C
C  Find speciation profile to use for speciation of biogenic emissions.
C  Calculate mole and mass factors for each species based on this profile.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C        11/99 by Jeff Vukovich
C       09/2025 by HT UNC-IE: Use M3UTILIO
C
C****************************************************************************/
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

        IMPLICIT NONE

C...........   Include files

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
c       INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
c       CHARACTER(2)    CRLF
c       INTEGER         GETFLINE
c       INTEGER         INDEX1
c       REAL            STR2REAL
c       LOGICAL         BLKORCMT

c       EXTERNAL        CRLF, GETFLINE, INDEX1, STR2REAL, BLKORCMT
        INTEGER, EXTERNAL :: GETFLINE
        LOGICAL, EXTERNAL :: BLKORCMT

C...........   Subroutine arguments (note- outputs MXSPFUL, MXSPEC, and SPCNAMES
C              passed via module MODSPRO)

        INTEGER     , INTENT  (IN) :: FDEV            ! file unit number
        INTEGER     , INTENT  (IN) :: NIPOL           ! number of pollutants
        INTEGER     , INTENT  (IN) :: MSPCS           ! max number of species
        CHARACTER(IOVLEN3), INTENT  (IN) :: EINAM( NIPOL )  ! pollutant names
        CHARACTER(*), INTENT  (IN) :: SPPRO           ! biogenic profile to find 
        CHARACTER(*), INTENT  (IN) :: EMSPC( MSPCS ) 

        REAL, INTENT (OUT) :: MLFAC( MSPCS, NIPOL ) 
        REAL, INTENT (OUT) :: MSFAC( MSPCS, NIPOL ) 

        INTEGER, PARAMETER :: MXSEG = 6        ! # of potential line segments
        LOGICAL      :: EFLAG = .FALSE.   ! error flag

C...........   Other arrays

        CHARACTER(32) SEGMENT( MXSEG )             ! Segments of parsed lines
                
C...........   Local variables

        INTEGER        I, J, K    ! counters and indices
        INTEGER        INPRFTP    ! tmp. profile number
        INTEGER        IOS        ! i/o status
        INTEGER        IREC       ! record counter
        INTEGER        NLINES     ! number of lines in data file
        INTEGER        PPOS       ! tmp position (from INDEX1) of pol in POLNAMA
        INTEGER        SPOS       ! tmp position (from INDEX1) of pol in SPECNMA

        REAL           SPLTFAC, SDIV, SMFAC ! tmp speciation profile factors
        LOGICAL      :: FOUND_FLAG = .FALSE.   ! flag for requested profile
        CHARACTER(SPNLEN3)  TMPPRF     ! tmp profile number
        CHARACTER(IOVLEN3)  POLNAM     ! pollutant name
        CHARACTER(IOVLEN3)  SPECNM     ! tmp species name
        CHARACTER(300) LINE       ! buffer for profile data
        CHARACTER(256) MESG       ! message buffer
        
        CHARACTER(16) :: PROGNAME = 'RDBPRO' ! program name

C***********************************************************************
C   Begin body of subroutine RDBPRO

C...........  Make sure routine arguments are valid
        IF( FDEV .LE. 0 .OR. NIPOL .LE. 0 ) THEN
            MESG = 'INTERNAL ERROR: Invalid subroutine arguments'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
      
        MLFAC = 0.0   ! array
        MSFAC = 0.0   ! array
 
C...........   Determine length of input file 

        NLINES = GETFLINE( FDEV, 'Speciation profile file' )

C...........   Initialize species count per pollutant and flag for indicating
C              true molar conversions (NOTE - for some pollutants like PM10,
C              there is no mole-based factor and outputu should be in units
C              of gm/mole into the mole-base speciation matrix)

C...........   Read through input file to determine the total number
C              of pollutants in the input file, to determine the
C              number of profiles per pollutant, to store the unique 
C              species names, and to store the units for mass-based and
C              mole-based conversions
        IREC   = 0
        DO I = 1, NLINES

            READ( FDEV,93100,END=999,IOSTAT=IOS ) LINE 
      
            IREC = IREC + 1
             
            IF( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 'reading speciation profile '//
     &              'file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Separate the line of data into each part
            CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.............  Left-justify character strings and convert factors to reals
            TMPPRF = ADJUSTL ( SEGMENT( 1 ) )

            IF( TMPPRF .NE. SPPRO ) CYCLE

            FOUND_FLAG = .TRUE.

            POLNAM = ADJUSTL ( SEGMENT( 2 ) )
            SPECNM = ADJUSTL ( SEGMENT( 3 ) )
            SPLTFAC = STR2REAL( SEGMENT( 4 ) )
            SDIV    = STR2REAL( SEGMENT( 5 ) )
            SMFAC   = STR2REAL( SEGMENT( 6 ) )

C.............  Search for pollutant in list of valid names, and go to the end
C               of the loop if not found (skip entry)

            J = INDEX1( POLNAM, NIPOL, EINAM )
    
            IF( J .EQ. 0 ) CYCLE
            
            SPOS = INDEX1( SPECNM, MSPCS, EMSPC )

            IF ( SPOS .GT. 0 ) THEN

                IF ( SDIV .EQ. 0 ) THEN
                    EFLAG = .TRUE.
                    CYCLE
                ENDIF

                MLFAC( SPOS, J ) = SPLTFAC / SDIV
                MSFAC( SPOS, J ) = SMFAC * GM2TON
             
            END IF
       

        END DO

        IF( .NOT. FOUND_FLAG ) THEN
            MESG =  SPPRO // ' profile not found in GSPRO'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( EFLAG ) THEN
            MESG = 'At least one of the divisors was zero in GSPRO.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
 
       
       RETURN
       
C......... Error message for reaching the end of file too soon
999    MESG = 'End of file reached unexpectedly. ' //
     &        'Check format of speciation' // CRLF() // BLANK5 //
     &        'profile file.'
       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

       
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93100   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )
        
       END SUBROUTINE RDBPRO                                                                            
