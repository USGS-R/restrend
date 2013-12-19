      SUBROUTINE porder(X,N,IX)

C===============================================================================
C
C        AUTHOR....TIM COHN
C        DATE......APRIL 1,  1986
C        REVISED...AUGUST 9, 1986      (TAC)
C
C===============================================================================

      REAL*8 X(*)
      INTEGER N,INDX,I,ICT,IX(*)
      INTEGER, allocatable :: L(:),R(:),P(:)
      allocate (L(0:N),R(0:N),P(0:N))

C===============================================================================
C
C    FIRST CHECK TO SEE IF WE HAVE AN ORDERED DATA VECTOR TO BEGIN WITH
C
      DO 50 I2=2,N
         IX(I2)    =  I2
         IF(X(I2) .LT. X(I2-1)) GOTO 1
 50   CONTINUE
      IX(1)     =  1
      RETURN

 1    CONTINUE

      L(1) =  0
      R(1) =  0
      P(1) =  0
      
      DO 10 I=2,N
         INDX =  1
         L(I)   =  0
         R(I)   =  0
         
 20      CONTINUE
         IF(X(I) .GE. X(INDX)) THEN
            IF(R(INDX) .EQ. 0) THEN
               R(INDX)   =  I
               P(I)      =  INDX
               GOTO 10
            ELSE
               INDX =  R(INDX)
               GOTO 20
            ENDIF
         ELSE
            IF(L(INDX) .EQ. 0) THEN
               L(INDX)   =  I
               P(I)      =  INDX
               GOTO 10
            ELSE
               INDX =  L(INDX)
               GOTO 20
            ENDIF
         ENDIF
 10   CONTINUE

      INDX =  1
      DO 40 ICT=1,N
         
 30      CONTINUE
         IF(L(INDX) .EQ. 0) THEN
            IX(ICT)     =  INDX
            P(R(INDX))   =  P(INDX)
            L(P(INDX))   =  R(INDX)
            INDX         =  P(INDX)
         ELSE
            INDX =  L(INDX)
            GOTO 30
         ENDIF
 40   CONTINUE
      RETURN
      END
