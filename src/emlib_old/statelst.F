
        SUBROUTINE  STATELST( UNIT , HEADER, MAXROW, NROW, NCOL, 
     &                        FMT, NAMES, STATES, RVAL           )

C***********************************************************************
C  subroutine body starts at line  109
C
C  FUNCTION:  
C     Writes out a report of values to opened file UNIT, from array 
C     RVAL for each state in the STATES list using the format in 
C     FMT and the HEADER
C             
C  PRECONDITIONS REQUIRED:  
C       UNIT must be a valid and opened unit number
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       prototype 4/97 by M Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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

      IMPLICIT NONE

      INCLUDE 'PARMS3.EXT'
      INCLUDE 'IODECL3.EXT'
      INCLUDE 'FDESC3.EXT'

C...........   ARGUMENTS and their descriptions:
        
        INTEGER         UNIT           ! length of input vector
        CHARACTER*256   HEADER         ! header for top of output file
        INTEGER         MAXROW         ! max length of output vector
        INTEGER         NROW           ! length of output vector
        INTEGER         NCOL           ! length of input vector
        CHARACTER*(*)   FMT            ! format to use for out values
        CHARACTER*(*)   NAMES ( NCOL ) ! names for columns
        INTEGER         STATES( NROW ) ! length of output vector
        REAL            RVAL  ( MAXROW, NCOL )   !  reported values

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         STR2INT
        INTEGER         TRIMLEN

        EXTERNAL        STR2INT, TRIMLEN
        
C...........   LOCAL VARIABLES and their descriptions:
        
        REAL            TOTAL( 100 )           ! columns totals

        INTEGER         C, R                   ! counters
        INTEGER         L2, LE, LF, LG, LI, LP
        INTEGER         LENLAB                 ! length of column label fmt
        INTEGER         TOT1, TOT2             ! 1st and 2nd part of totals fmt
        INTEGER         WIDTH

        CHARACTER*1     TOTR                   ! temporary totals format type
        CHARACTER*50    LABELFMT               ! output columns lamel format
        CHARACTER*50    OUTFMT                 ! output table format
        CHARACTER*50    TOTFMT                 ! output totals format

        LOGICAL         INTFLAG                ! true if integer FMT on input

        CHARACTER*20    STATNAM( 99 )   ! Names of states sorted by FIPs code
        CHARACTER*20    TMPNAME

        DATA STATNAM /  'Alabama', 'Alaska', ' ', 'Arizona', 'Arkansas',
     &                  'California', ' ', 'Colorado', 'Connecticut',
     &                  'Delaware', 'District of Columbia', 'Florida',
     &                  'Georgia', ' ', 'Hawaii', 'Idaho', 'Illinois',
     &                  'Indiana', 'Iowa', 'Kansas', 'Kentucky', 
     &                  'Louisiana', 'Maine', 'Maryland', 
     &                  'Massachusetts', 'Michigan', 'Minnesota',
     &                  'Mississippi', 'Missouri', 'Montana', 
     &                  'Nebraska', 'Nevada', 'New Hampshire', 
     &                  'New Jersey', 'New Mexico', 'New York',
     &                  'North Carolina', 'North Dakota', 'Ohio',
     &                  'Oklahoma', 'Oregon', 'Pennsylvania',
     &                  'Puerto Rico', 'Rhode Island', 'South Carolina',
     &                  'South Dakota', 'Tennessee', 'Texas', 'Utah',
     &                  'Vermont', 'Virginia', ' ', 'Washington',
     &                  'West Virginia', 'Wisconsin',
     &                  19*' ', 'Canada', 23*' ', 'Off Shore'  /
 
C***********************************************************************
C   begin body of subroutine  STATELST

C.........  Perform input-format-dependent changes

        INTFLAG = .FALSE.
        LENLAB  = 0

        L2 = TRIMLEN( FMT )     ! Length of format statement
        LI = INDEX( FMT, 'I' )  ! Location of I (for integer format)
        LF = INDEX( FMT, 'F' )  ! Location of F (for float   format)
        LE = INDEX( FMT, 'E' )  ! Location of E (for expontial format)
        LG = INDEX( FMT, 'G' )  ! Location of G (for general format)
        LP = INDEX( FMT, '.' )  ! Location of . ( for F, E, G format)

        IF( LP .EQ. 0 ) LP = L2 ! If no '.', then use end of format

        IF( LI .GT. 0 ) THEN

            INTFLAG = .TRUE.

            LENLAB = STR2INT( FMT( LI+1:L2 ) )

        ELSEIF( LF .GT. 0 ) THEN

            LENLAB = STR2INT( FMT( LF+1:LP-1 ) )
            TOTR   = FMT( LF:LF )
            TOT1   = LENLAB
            TOT2   = MAX( 0, STR2INT( FMT( LP+1:L2 ) ) - 1 )
           
        ELSEIF( LE .GT. 0 ) THEN

            LENLAB = STR2INT( FMT( LE+1:LP-1 ) )
            TOTR   = FMT( LE:LE )
            TOT1   = LENLAB
            TOT2   = STR2INT( FMT( LP+1:L2 ) )

        ELSEIF( LG .GT. 0 ) THEN

            LENLAB = STR2INT( FMT( LG+1:LP-1 ) )
            TOTR   = FMT( LG:LG )
            TOT1   = LENLAB
            TOT2   = MAX( 0, STR2INT( FMT( LP+1:L2 ) ) - 1 )

        ENDIF
       
        IF( LENLAB .EQ. 0 ) THEN
            CALL M3MSG2( 'Bad format in STATELST, using F5.0' )
            LENLAB = 5
            FMT    = 'F5.0'
            TOT1   = LENLAB
            TOT2   = MAX( 0, STR2INT( FMT( LP+1:L2 ) ) - 1 )
        ENDIF

C.........  Determine width
        WIDTH = 16 + LENLAB * NCOL

C.........  Create format statement for column labels

        WRITE( LABELFMT,94010 ) '( A' , LENLAB,  ', 1X, $ )'

C.........  Create format statement for output data rows

        WRITE( OUTFMT,94010 ) '( I2.2, 2X, ', NCOL, '(' // 
     &                        FMT( 1:L2 ) // ', 1X ), 2X, A )'

C.........  Create format statement for output totals

        WRITE( TOTFMT,94010 ) '( I2.2, 2X, ', NCOL, '(' //
     &                        TOTR, TOT1, '.', TOT2, ', 1X ), 2X, A )'

C.........  Write header

        WRITE( UNIT, 93030 ) HEADER( 1:TRIMLEN( HEADER ) )

C.........  Write column labels
        WRITE( UNIT, 93000 ) 'STATE'
        WRITE( UNIT, 93040 ) 'CODE   '

        DO 66 C = 1, NCOL
            WRITE( UNIT, LABELFMT ) NAMES ( C )
66      CONTINUE

        WRITE( UNIT, 93000 ) 'STATE NAME'

C.........  Initialize totals
        DO 88 C = 1, NCOL
            TOTAL( C ) =  0.
88      CONTINUE 

C.........  Write data
        IF( INTFLAG ) THEN

            DO 100 R = 1, NROW

                WRITE( UNIT,OUTFMT ) STATES( R ),
     &                              ( INT( RVAL( R,C ) ), C = 1, NCOL ),
     &                               STATNAM( STATES( R ) )

                DO 99 C = 1, NCOL
                    TOTAL( C ) = TOTAL( C ) + RVAL( R,C )
99              CONTINUE 

100         CONTINUE

        ELSE 

            DO 111 R = 1, NROW

                TMPNAME = STATNAM( STATES( R ) )
                WRITE( UNIT,OUTFMT ) STATES( R ),
     &                              ( RVAL( R,C ), C = 1, NCOL ),
     &                               TMPNAME( 1:TRIMLEN( TMPNAME ) )

                DO 110 C = 1, NCOL
                    TOTAL( C ) = TOTAL( C ) + RVAL( R,C )
110             CONTINUE 

111         CONTINUE

        ENDIF

        WRITE( UNIT, 93050 )

        IF( INTFLAG ) THEN

            WRITE( UNIT,TOTFMT ) 0, ( INT( TOTAL( C ) ), 
     &                                     C = 1, NCOL   ), 'Total'
 
        ELSE 

            WRITE( UNIT,TOTFMT ) 0, ( TOTAL( C ), C = 1, NCOL ), 'Total'

        ENDIF

        WRITE( UNIT,93000 ) ' '
        WRITE( UNIT,93000 ) ' '

        RETURN

C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Formatted file I/O formats............ 93xxx
 
93000   FORMAT( A )
 
93030   FORMAT( A, //, 22('-'), ' By State ', 22('-') )
 
93040   FORMAT( A, $ )

93050   FORMAT( 1X, <WIDTH>('-') )

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10( A, :, I3, :, 2X ) )
 
C******************  FORMAT  STATEMENTS   ******************************

        END

