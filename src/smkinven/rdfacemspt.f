
        SUBROUTINE RDFACEMSPT( FDEV, UTMZONE, NCNTY, FIPTOCSRC )

C***********************************************************************
C  subroutine body starts at line 232
C
C  DESCRIPTION:
C      This subroutine reads the EMS-95 point facility file.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Copied from rdemspt.f by C. Seppanen (2/03)
C
C****************************************************************************
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

C...........   MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC
        
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CSOURC, XLOCA, YLOCA, CPDESC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         FINDCFIRST
        INTEGER         FIND1
        INTEGER         STR2INT
        REAL            STR2REAL
        
        EXTERNAL        CRLF, FINDCFIRST, FIND1, STR2INT, STR2REAL

C...........   SUBROUTINE ARGUMENTS
        INTEGER,          INTENT (IN) :: FDEV   ! unit no. of inv file
        INTEGER,          INTENT(OUT) :: UTMZONE( NSRC )  ! UTM zone by source
        INTEGER,          INTENT (IN) :: NCNTY  ! no. counties in inventory
        INTEGER,          INTENT (IN) :: FIPTOCSRC( NCNTY+1,2 )  ! index into CSOURC

C...........   Other local variables
        INTEGER              I, K1                  ! counters and indices
        
        INTEGER              ACTCNT                 ! actual count of processed records
        INTEGER              IOS                    ! I/O status
        INTEGER              IREC                   ! no. records in file
        INTEGER              IFIP                   ! integer FIPS code
        INTEGER              STIDX, ENDIDX          ! start and end idx into CSOURC
        INTEGER              KEYLEN                 ! length of source key
        INTEGER              ZONE                   ! UTM zone

        REAL                 XVAL                   ! tmp X coordinate
        REAL                 YVAL                   ! tmp Y coordinate
        REAL                 XX                     ! Longitude
        REAL                 YY                     ! Latitude
        
        CHARACTER(LEN=FIPLEN3) CFIP      ! fip code
        CHARACTER(LEN=PLTLEN3) FCID      ! facility ID
        
        CHARACTER(LEN=ALLLEN3) SRCKEY                 ! source key
        CHARACTER(LEN=300)     LINE                   ! line from file
        CHARACTER(LEN=300)     MESG                   ! message buffer

        CHARACTER(LEN=16) :: PROGNAME = 'RDFACEMSPT' ! program name

C***********************************************************************
C   begin body of subroutine RDFACEMSPT

        ACTCNT = 0
        IREC = 0

C.........  Loop through lines in file
        DO
    
            READ( FDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1
            
C.............  Check for I/O errors
            IF( IOS > 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &              'reading facility file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
    
C.............  Check for end of file
            IF( IOS < 0 ) EXIT
    
C.............  Read FIPS code from file
            CFIP( 1:1 ) = '0'
            CFIP( 2:3 ) = ADJUSTR( LINE( 1:2 ) )
            CFIP( 4:6 ) = ADJUSTR( LINE( 3:5 ) )
                        
C.............  Replace blanks with zeros        
            DO I = 1,FIPLEN3
                IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
            END DO
            
C.............  Find county in FIPTOCSRC array
            IFIP = STR2INT( CFIP( 1:6 ) )
            K1 = FIND1( IFIP, NCNTY+1, FIPTOCSRC( :,1 ) )

C.............  Make sure county is in the inventory
            IF( K1 <= 0 ) CYCLE
            
C.............  Set starting and ending indices            
            STIDX = FIPTOCSRC( K1,2 )
            ENDIDX = FIPTOCSRC( K1+1,2 ) - 1
            
C.............  Build source key
            FCID = ADJUSTL( LINE(  6:20 ) )
            
            CALL BLDCSRC( CFIP, FCID, CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                    CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                    SRCKEY )
     
            KEYLEN = LEN_TRIM( SRCKEY )
            
C.............  Find source key in CSOURC array
            K1 = FINDCFIRST( SRCKEY, ENDIDX-STIDX+1, 
     &                       CSOURC( STIDX:ENDIDX )( 1:KEYLEN ) )
    
C.............  If key not found, go to next line
            IF( K1 <= 0 ) CYCLE

C.............  Shift index to account for only searching part of CSOURC array
            K1 = K1 + STIDX - 1
    
C.............  Read and check UTM zone
            ZONE = STR2INT( LINE( 43:44 ) )
            
            IF( ZONE == IMISS3 ) THEN
                WRITE( MESG,94010 )
     &                 'UTM zone is blank or badly formatted at ' //
     &                 'line', IREC, 'in facility file' // CRLF() // 
     &                 BLANK10 // 'Assuming lat/lon coordinates'
                CALL M3MESG( MESG )
                CYCLE
            END IF
    
C.............  Read and check coordinates
            XVAL = STR2REAL( LINE( 25:33 ) )
            YVAL = STR2REAL( LINE( 34:42 ) )
            
            IF( XVAL <= 0.0 .OR. YVAL <= 0.0 ) CYCLE
            
C.............  Increment count of actual records            
            ACTCNT = ACTCNT + 1
            
            IF( ZONE > 0 ) CALL UTM2LL( XVAL, YVAL, ZONE, XX, YY )

C.............  Loop through matching sources and set values
            DO
            
C.................  Check for end of array or end of matching sources
                IF( K1 > NSRC ) EXIT
                IF( CSOURC( K1 )( 1:KEYLEN ) /= SRCKEY ) EXIT

                UTMZONE( K1 ) = ZONE
                XLOCA  ( K1 ) = XX             ! in lat-lon
                YLOCA  ( K1 ) = YY             ! in lat-lon
                CPDESC ( K1 ) = LINE( 45:84 )
                
                K1 = K1 + 1
            END DO
    
        END DO

C.........  Write count of total records processed    
        WRITE( MESG,94010 )
     &         'FACILITY FILE processed' // CRLF() // BLANK10 //
     &         'Actual FACILITY record count ', ACTCNT
        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94125   FORMAT( I5 )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE RDFACEMSPT
