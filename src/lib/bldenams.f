
        SUBROUTINE BLDENAMS( CATEGORY, NIPOL, NPPOL, EINAM, 
     &                       OUTNAMES, OUTUNITS, OUTTYPES, OUTDESCS )

C***********************************************************************
C  subroutine body starts at line 141
C
C  DESCRIPTION:
C      This subroutine builds names for pollutant-specific inventory variables
C      such as rule effectiveness and control efficiency, while ensuring that 
C      no names are duplicated.  Since inventory pollutant names are allowed
C      to be 16 characters, there needs to be a truncation and check for
C      duplicates, and then an insertion of the prefixes.
C
C  PRECONDITIONS REQUIRED:
C      Source category specified correctly
C      Pollutant list and number provided
C      Memory allocated for output arrays
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     Subroutines: Models-3 subroutines
C     Functions: Models-3 functions
C
C  REVISION  HISTORY:
C     Created 12/98 by M. Houyoux
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

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         TRIMLEN

        EXTERNAL        CRLF, TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=*)       CATEGORY                !  source category
        INTEGER                NIPOL                   !  number of pollutants
        INTEGER                NPPOL                   !  number vars per pol
        CHARACTER(LEN=*)       EINAM   ( NIPOL )       !  pollutant names
        CHARACTER(LEN=*)       OUTNAMES( NIPOL,NPPOL ) !  output names
        CHARACTER(LEN=*)       OUTUNITS( NIPOL,NPPOL ) !  output units
        INTEGER                OUTTYPES( NIPOL,NPPOL ) !  output types int/real
        CHARACTER(LEN=*)       OUTDESCS( NIPOL,NPPOL ) !  output descriptions

C.........  Data arrays for variable names...

C.........  Point source variable name parameters
        CHARACTER(LEN=CPRTLEN3) PTPREFIX( 2:NPTPPOL3 )
        DATA PTPREFIX / OZNSEART, CTLEFFRT, RULEFFRT, 
     &                  EMISFCRT, CECOD1RT, CECOD2RT /

        CHARACTER(LEN=IOULEN3) PTUNITS( NPTPPOL3 )
        DATA PTUNITS / 'tons/year', 'tons/day', '%', '%',
     &                 'SCC units', 'n/a', 'n/a' /

        INTEGER PTTYPES( NPTPPOL3 )
        DATA PTTYPES / M3REAL, M3REAL, M3REAL, M3REAL, M3REAL, 
     &                 M3INT, M3INT /

        CHARACTER(LEN=IODLEN3) PTDESCS( NPTPPOL3 )
        DATA PTDESCS / 'Annual Emissions'
     &               , 'Ozone Season Emissions'
     &      , 'Control efficiency (in [0,100], or "MISSING": < -9.0E36)'
     &      , 'Rule effectiveness (in [0,100], or "MISSING": < -9.0E36)'
     &               , 'Emission Factors'
     &               , 'Primary Control Equipment Code'
     &               , 'Secondary Control Equipment Code'
     &               /

C.........  Area source variable name parameters
        CHARACTER(LEN=CPRTLEN3) ARPREFIX( 2:NARPPOL3 )
        DATA ARPREFIX / OZNSEART, EMISFCRT, CTLEFFRT, 
     &                  RULEFFRT, RULPENRT /

        CHARACTER(LEN=IOULEN3) ARUNITS( NARPPOL3 )
        DATA ARUNITS / 'tons/year', 'tons/day', 
     &                 'SCC units', '%', '%', '%' /

        INTEGER ARTYPES( NARPPOL3 )
        DATA ARTYPES / M3REAL, M3REAL, M3REAL, M3REAL, M3REAL, M3REAL /

        CHARACTER(LEN=IODLEN3) ARDESCS( NARPPOL3 )
        DATA ARDESCS / 'Annual Emissions'
     &               , 'Ozone Season Emissions'
     &               , 'Emission Factors'
     &      , 'Control efficiency (in [0,100], or "MISSING": < -9.0E36)'
     &      , 'Rule effectiveness (in [0,100], or "MISSING": < -9.0E36)'
     &      , 'Rule penetration (in [0,100], or "MISSING": < -9.0E36)'
     &               /

C...........   Unsorted pollutant records
        INTEGER       INDEXA ( NIPOL )
        CHARACTER(LEN=IOVLEN3) ABRNAMA( NIPOL )    !  pollutant names

C...........   Other local variables
        INTEGER         COD     !  tmp for pollutant code
        INTEGER         I, J    !  counters and indices
        INTEGER         IDIF    !  abridged name length
        INTEGER         L, L2   !  length indices
        INTEGER         LCNT    !  same-abrigded-name counter

        CHARACTER(LEN=IOVLEN3)  LNAM  !  previous pollutant name
        CHARACTER(LEN=IOVLEN3)  NAM   !  current pollutant name
        CHARACTER*10    BUFFER        !  tmp string buffer
        CHARACTER*300   MESG          !  message buffer

        CHARACTER*16 :: PROGNAME = 'BLDENAMS' ! program name

C***********************************************************************
C   begin body of subroutine BLDENAMS

C.........  Check to see if someone has modified include file in calling prog
        IF( NPPOL .NE. NPTPPOL3 ) THEN
             WRITE( MESG,94010 ) 'INTERNAL ERROR: ' //
     &         'Subroutine ' // PROGNAME // ' needs to be modified ' //
     &         CRLF() // BLANK16 // 'to be consistent with ' //
     &         'changed value of NPTPPOL3 in EMCNST3.EXT' //
     &         CRLF() // BLANK16 // 'in include file EMDIMS3.EXT'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Truncate variable names based on original length and length of
C.........  fields to be inserted

        IDIF = IOVLEN3 - CPRTLEN3
        DO I = 1, NIPOL
           INDEXA ( I ) = I        
           ABRNAMA( I ) = EINAM( I )( 1:IDIF )
        ENDDO

C.........  Sort to set up for duplicates search
        CALL SORTIC( NIPOL, INDEXA, ABRNAMA )

C.........  Search for duplicates and rename them
        LNAM = EMCMISS3
        LCNT = 0
        DO I = 1, NIPOL

            J = INDEXA( I )
            NAM = ABRNAMA( J )

            IF( NAM .EQ. LNAM ) THEN
                LCNT = LCNT + 1
                WRITE( BUFFER, '(I10)' ) LCNT  ! turn in to string
                BUFFER = ADJUSTL( BUFFER )     ! left justify

                L  = TRIMLEN( BUFFER )         ! get length of trailer
                L2 = MIN( IDIF -L -1, TRIMLEN( NAM ) )  ! remaining name
                ABRNAMA( J ) = NAM( 1:L2 ) // '_' // BUFFER( 1:L ) ! store name
            ELSE
               LCNT = 0
               LNAM = NAM
            ENDIF

        ENDDO

C.........  Build output names
        DO I = 1, NIPOL

            OUTNAMES( I,1 ) = EINAM( I )  ! Annual emis name is same as pol

            L  = TRIMLEN( ABRNAMA( I ) )
            L2 = TRIMLEN( CATEGORY )

            IF( CATEGORY(1:L2) .EQ. 'POINT' ) THEN

                OUTUNITS( I,1 ) = PTUNITS( 1 )          
                OUTTYPES( I,1 ) = PTTYPES( 1 )           
                OUTDESCS( I,1 ) = PTDESCS( 1 )
         
                DO J = 2, NPTPPOL3
                    OUTNAMES( I,J ) = PTPREFIX( J )//ABRNAMA( I )( 1:L )
                    OUTUNITS( I,J ) = PTUNITS( J )          
                    OUTTYPES( I,J ) = PTTYPES( J )           
                    OUTDESCS( I,J ) = PTDESCS( J )          
                ENDDO

            ELSEIF( CATEGORY(1:L2) .EQ. 'AREA' ) THEN

                OUTUNITS( I,1 ) = ARUNITS( 1 )          
                OUTTYPES( I,1 ) = ARTYPES( 1 )           
                OUTDESCS( I,1 ) = ARDESCS( 1 )
         
                DO J = 2, NARPPOL3
                    OUTNAMES( I,J ) = ARPREFIX( J )//ABRNAMA( I )( 1:L )
                    OUTUNITS( I,J ) = ARUNITS( J )          
                    OUTTYPES( I,J ) = ARTYPES( J )           
                    OUTDESCS( I,J ) = ARDESCS( J )          
                ENDDO

            ELSE
                MESG = 'INTERNAL ERROR: Do not know how to build ' //
     &                 'names for category ' // CATEGORY( 1:L2 )
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            ENDIF

        ENDDO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )   

        END
