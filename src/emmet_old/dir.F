C ****************************************************************************
C *
C * Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
C * File Version: @(#)$Id$
C *
C * COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C * All Rights Reserved
C *
C * See file COPYRIGHT for conditions of use.
C *
C * Environmental Programs Group
C * MCNC--North Carolina Supercomputing Center
C * P.O. Box 12889
C * Research Triangle Park, NC  27709-2889
C *
C * env_progs@mcnc.org
C *
C * REVISION  HISTORY:
C *     Prototype  7/96 by Kiran Alapaty
C * Pathname: $Source$
C * Last updated: $Date$ 
C *
C
****************************************************************************/

C......... Routine to calculate wind direction in degrees
C            given the horizontal wind components in m/s
       Subroutine direct(uu1,vv1,angle)

       pi=4.0*atan(1.0)
       todeg=180.0/pi
       IF (UU1.EQ.0.0) THEN
         IF (VV1.EQ.0.0) THEN
            angle = 0.0
         ELSE IF (VV1.GT.0.0) THEN
            angle = 180.0
         ELSE
            angle = 360.0
         END IF
       ELSE IF (VV1.EQ.0.0) THEN
         IF (UU1.GT.0.0) THEN
            angle = 270.0
         ELSE 
            angle = 90.0
         END IF
       ELSE IF(UU1.LT.0.0.AND.VV1.LT.0.0)THEN 
         ADI1=0.0
         uu1vv1=uu1/vv1
         angle=ABS(atan(UU1VV1))*todeg + ADI1
       ELSE IF(UU1.LT.0.0.AND.VV1.GT.0.0)  THEN
         ADI1=90.0
         uu1vv1=vv1/uu1
         angle=ABS(atan(UU1VV1))*todeg + ADI1
       ELSE IF(UU1.GT.0.0.AND.VV1.GT.0.0)  THEN 
         ADI1=270.0
         uu1vv1=vv1/uu1
         angle=ADI1 - ABS(atan(UU1VV1))*todeg 
       ELSE IF(UU1.GT.0.0.AND.VV1.LT.0.0)  THEN
         ADI1=360.0
         uu1vv1=uu1/vv1
         angle=ADI1 - ABS(atan(UU1VV1)) *todeg
       END IF
       IF (angle.EQ.0.0) THEN
          IF (UU1.NE.0.0 .OR. VV1.NE.0.0) angle=360.0
       END IF
       return 
       end
*.................... end of the routine .............................

