
        SUBROUTINE BLDENAMS( CATEGORY_L, NIPPA_L, NPPOA, EANAM_L, 
     &                       OUTNAMES, OUTUNITS, OUTTYPES, OUTDESCS )

C***********************************************************************
C  subroutine body starts at line 154
C
C  DESCRIPTION:
C      This subroutine builds names for pol/act-specific inventory variables
C      such as rule effectiveness and control efficiency, while ensuring that 
C      no names are duplicated.  Since inventory pol/act names are allowed
C      to be 16 characters, there needs to be a truncation and check for
C      duplicates, and then an insertion of the prefixes.  If the UNITS
C      header option has been used in an inventory input file, it assigns 
C      variable units based on these settings.
C
C  PRECONDITIONS REQUIRED:
C      Source CATEGORY_L specified correctly
C      Pollutant/activity list and number provided
C      Memory allocated for output arrays
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     Subroutines: Models-3 subroutines
C     Functions: Models-3 functions
C
C  REVISION  HISTORY:
C     Created 12/98 by M. Houyoux
C
C*************************************************************************
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

C.........  This module contains the information about the source CATEGORY
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        EXTERNAL        CRLF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY_L               ! source CATEGORY
        INTEGER     , INTENT (IN) :: NIPPA_L                  ! no. of pol/act
        INTEGER     , INTENT (IN) :: NPPOA                    ! no. vars per
                                                              ! pol or activity
        CHARACTER(*), INTENT (IN) :: EANAM_L ( NIPPA_L )      ! pol/act names
        CHARACTER(*), INTENT(OUT) :: OUTNAMES( NIPPA_L,NPPOA )! var names
        CHARACTER(*), INTENT(OUT) :: OUTUNITS( NIPPA_L,NPPOA )! var units
        INTEGER     , INTENT(OUT) :: OUTTYPES( NIPPA_L,NPPOA )! var typ:int/real
        CHARACTER(*), INTENT(OUT) :: OUTDESCS( NIPPA_L,NPPOA )! var descriptions

C.........  Data arrays for variable names...

C.........  Area source variable name parameters
        CHARACTER(LEN=CPRTLEN3) ARPREFIX( 2:NARPPOL3 )
        DATA ARPREFIX / OZNSEART, EMISFCRT, CTLEFFRT, 
     &                  RULEFFRT, RULPENRT /

        CHARACTER(LEN=IOULEN3) ARUNITS( NARPPOL3 )
        DATA ARUNITS / 'tons/yr', 'tons/day', 
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

C.........  Mobile source variable name parameters
        CHARACTER(LEN=CPRTLEN3) MBPREFIX( 2:NMBPPOL3 )
        DATA MBPREFIX / OZNSEART /

        CHARACTER(LEN=IOULEN3) MBUNITS( NMBPPOL3 )
        DATA MBUNITS / 'tons/yr', 'tons/day' /

        INTEGER MBTYPES( NMBPPOL3 )
        DATA MBTYPES / M3REAL, M3REAL /

        CHARACTER(LEN=IODLEN3) MBDESCS( NMBPPOL3 )
        DATA MBDESCS / 'Annual Data', 'Ozone Season Data' /

C.........  Point source variable name parameters
        CHARACTER(LEN=CPRTLEN3) PTPREFIX( 2:NPTPPOL3 )
        DATA PTPREFIX / OZNSEART, CTLEFFRT, RULEFFRT, 
     &                  EMISFCRT, CECOD1RT, CECOD2RT /

        CHARACTER(LEN=IOULEN3) PTUNITS( NPTPPOL3 )
        DATA PTUNITS / 'tons/yr', 'tons/day', '%', '%',
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

C...........   Unsorted pollutant/activity records
        INTEGER       INDEXA ( NIPPA_L )
        CHARACTER(LEN=IOVLEN3) ABRNAMA( NIPPA_L )    !  pollutant/activity names

C...........   Other local variables
        INTEGER         COD       !  tmp for pollutant/activity code
        INTEGER         I, J, K   !  counters and indices
        INTEGER         IDIF      !  abridged name length
        INTEGER         L, L2, LU !  length indices
        INTEGER         LCNT      !  same-abrigded-name counter

        CHARACTER(LEN=IOVLEN3)  LNAM  !  previous pollutant/activity name
        CHARACTER(LEN=IOVLEN3)  NAM   !  current pollutant/activity name
        CHARACTER*10    BUFFER        !  tmp string buffer
        CHARACTER*300   MESG          !  message buffer

        CHARACTER*16 :: PROGNAME = 'BLDENAMS' ! program name

C***********************************************************************
C   begin body of subroutine BLDENAMS

C.........  Truncate variable names based on original length and length of
C.........  fields to be inserted

        IDIF = IOVLEN3 - CPRTLEN3
        DO I = 1, NIPPA_L
           INDEXA ( I ) = I
           L = MIN( IDIF, LEN_TRIM( EANAM_L( I ) ) )        
           ABRNAMA( I ) = EANAM_L( I )( 1:L )
        END DO

C.........  Sort to set up for duplicates search
        CALL SORTIC( NIPPA_L, INDEXA, ABRNAMA )

C.........  Search for duplicates and rename them
        LNAM = EMCMISS3
        LCNT = 0
        DO I = 1, NIPPA_L

            J = INDEXA( I )
            NAM = ABRNAMA( J )

            IF( NAM .EQ. LNAM ) THEN
                LCNT = LCNT + 1
                WRITE( BUFFER, '(I10)' ) LCNT  ! turn in to string
                BUFFER = ADJUSTL( BUFFER )     ! left justify

                L  = LEN_TRIM( BUFFER )         ! get length of trailer
                L2 = MIN( IDIF -L -1, LEN_TRIM( NAM ) )  ! remaining name
                ABRNAMA( J ) = NAM( 1:L2 ) // '_' // BUFFER( 1:L ) ! store name
            ELSE
               LCNT = 0
               LNAM = NAM
            ENDIF

        ENDDO

C.........  Build output names
        DO I = 1, NIPPA_L

            OUTNAMES( I,1 ) = EANAM_L( I )  ! Annual emis name is same as pol/act

            L  = LEN_TRIM( ABRNAMA( I ) )
            L2 = LEN_TRIM( CATEGORY_L )

            SELECT CASE( CATEGORY_L )
            CASE( 'AREA' )
                OUTUNITS( I,1 ) = ARUNITS( 1 )          
                OUTTYPES( I,1 ) = ARTYPES( 1 )           
                OUTDESCS( I,1 ) = ARDESCS( 1 )
         
                DO J = 2, NPPOA
                    OUTNAMES( I,J ) = ARPREFIX( J )//ABRNAMA( I )( 1:L )
                    OUTUNITS( I,J ) = ARUNITS( J )          
                    OUTTYPES( I,J ) = ARTYPES( J )           
                    OUTDESCS( I,J ) = ARDESCS( J )          
                ENDDO

            CASE( 'MOBILE' )

                OUTUNITS( I,1 ) = MBUNITS ( 1 )          
                OUTTYPES( I,1 ) = MBTYPES ( 1 )           
                OUTDESCS( I,1 ) = MBDESCS ( 1 )
                DO J = 2, NPPOA
                    OUTNAMES( I,2 ) = MBPREFIX( 2 )//ABRNAMA( I )( 1:L )
                    OUTUNITS( I,2 ) = MBUNITS ( 2 )          
                    OUTTYPES( I,2 ) = MBTYPES ( 2 )           
                    OUTDESCS( I,2 ) = MBDESCS ( 2 )
                ENDDO

            CASE( 'POINT' )

                OUTUNITS( I,1 ) = PTUNITS( 1 )          
                OUTTYPES( I,1 ) = PTTYPES( 1 )           
                OUTDESCS( I,1 ) = PTDESCS( 1 )
         
                DO J = 2, NPPOA
                    OUTNAMES( I,J ) = PTPREFIX( J )//ABRNAMA( I )( 1:L )
                    OUTUNITS( I,J ) = PTUNITS( J )          
                    OUTTYPES( I,J ) = PTTYPES( J )           
                    OUTDESCS( I,J ) = PTDESCS( J )          
                ENDDO

            CASE DEFAULT
                MESG = 'INTERNAL ERROR: Do not know how to build ' //
     &                 'names for CATEGORY ' // CATEGORY_L( 1:L2 )
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END SELECT

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )   

        END
