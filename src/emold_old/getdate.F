
        INTEGER FUNCTION GETDATE ( DEFAULT , PROMPT )

C********************************************************************
C       function body starts at line  72
C
C       CALLS:      GETNUM
C
C  FUNCTION:
C
C       Display the  PROMPT  for a date in either of two formats:
C       YYYYDDD or YYYYMMDD
C
C  REVISION HISTORY:
C
C       Create by M Houyoux: 2/97
C
C  ARGUMENT LIST DESCRIPTION:
C
C    Input arguments:
C
C        DEFAULT  Default return value
C        PROMPT   Prompt for user
C
C    Output arguments:  none
C
C  RETURNS   user response after checking its range; or default.
C
C********************************************************************
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

C.......   ARGUMENTS:

        INTEGER         DEFAULT
        CHARACTER*(*)   PROMPT

C.......   EXTERNAL FUNCTIONS:

        INTEGER         GETNUM, JULIAN
        EXTERNAL        GETNUM, JULIAN


C.......   LOCAL VARIABLES:

        INTEGER         JDATE, YEAR, MONTH, DAY

C......................................................................
C       begin GETDATE

        JDATE = GETNUM( 0, 99999999, DEFAULT, PROMPT )

        IF( JDATE .GT. 9999999 ) THEN
            YEAR = JDATE/10000
            MONTH = ( JDATE-YEAR*10000 ) / 100
            DAY   = MOD( JDATE-YEAR*10000, 100 )
            JDATE = YEAR*1000 + JULIAN( YEAR, MONTH, DAY )
        ENDIF

        GETDATE = JDATE

        RETURN

C................   end body of GETDATE  .......................................

        END

C................   end   GETDATE  ....................................

