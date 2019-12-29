c ici, on copie toutes les routines gaN qui sont utiles au fonctionnement
c du code. Ne pas oublier les routines qui s'appelent...

c NIVEAU 0

c-----------------------------------------------------------------
c                    POUR LE CODE MATRIX3D.F
c-----------------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
      DOUBLE PRECISION X02CON
      DATA X02CON /1.11022302462516D-16 /
C     .. Executable Statements ..
      X02AJF = X02CON
      RETURN
      END
      SUBROUTINE F01AYF(N,TOL,A,IA,D,E,E2)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 15A REVISED. IER-907 (APR 1991).
C
C     TRED3
C     THIS SUBROUTINE REDUCES THE GIVEN LOWER TRIANGLE OF A
C     SYMMETRIC
C     MATRIX, A, STORED ROW BY ROW IN THE ARRAY A(N(N+1)/2), TO
C     TRIDIAGONAL FORM USING HOUSEHOLDERS REDUCTION. THE DIAGONAL
C     OF THE RESULT IS STORED IN THE ARRAY D(N) AND THE
C     SUB-DIAGONAL IN THE LAST N-1 STORES OF THE ARRAY E(N)
C     (WITH THE ADDITIONAL ELEMENT E(1) = 0). E2(I) IS SET TO EQUAL
C     E(I)**2. THE ARRAY A IS USED TO STORE SUFFICIENT INFORMATION
C     FOR THE DETAILS OF THE TRANSFORMATION TO BE RECOVERABLE IN
C     THE SUBROUTINE F01AZF
C     1ST. MARCH  1972
C
C     REVISED BY VINCE FERNANDO AT MARK 14 TO INTRODUCE SCALING INTO
C     THE GENERATION OF HOUSEHOLDER MATRICES AS PROPOSED BY
C     G.W. STEWART, INTRODUCTION TO MATRIX COMPUTATIONS, CHAPTER 7.
C     TOL IS NOW A DUMMY PARAMETER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA), D(N), E(N), E2(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F, G, H, HH, SCALE, SMALL
      INTEGER           I, II, IPOS, IZ, J, K, L
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AMF
      INTEGER           IDAMAX
      EXTERNAL          DDOT, IDAMAX, X02AMF
C     .. External Subroutines ..
      EXTERNAL          DSCAL, DSPMV, DSPR2
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      SMALL = X02AMF()
      DO 140 II = 1, N
         I = N - II + 1
         L = I - 1
         IZ = (I*L)/2
         H = 0.0D0
         SCALE = 0.0D0
         IF (L.EQ.0) GO TO 40
         DO 20 K = 1, L
            IPOS = IZ + K
            F = A(IPOS)
            D(K) = F
   20    CONTINUE
C        FIND THE ELEMENT OF LARGEST ABSOLUTE VALUE IN D
         K = IDAMAX(L,D,1)
         SCALE = ABS(D(K))
C        IF D IS A NULL VECTOR THEN SKIP THE TRANSFORMATION
         IF (SCALE.GE.SMALL) GO TO 60
   40    E(I) = 0.0D0
         E2(I) = 0.0D0
         H = 0.0D0
         GO TO 120
   60    CALL DSCAL(L,1.0D0/SCALE,D,1)
         H = DDOT(L,D,1,D,1)
         E2(I) = H*SCALE*SCALE
         F = D(L)
         G = SQRT(H)
         IF (F.GE.0.0D0) G = -G
         E(I) = G*SCALE
         H = H - F*G
         D(L) = F - G
         IPOS = IZ + L
         A(IPOS) = D(L)*SCALE
C        FORM A*U
         CALL DSPMV('U',L,1.0D0/H,A,D,1,0.0D0,E,1)
C        FORM P
         F = 0.0D0
         DO 80 J = 1, L
            F = F + E(J)*D(J)
   80    CONTINUE
C        FORM K
         HH = F/(H+H)
C        FORM Q
         DO 100 J = 1, L
            E(J) = E(J) - HH*D(J)
  100    CONTINUE
C        FORM REDUCED A
         CALL DSPR2('U',L,-1.0D0,D,1,E,1,A)
  120    IPOS = IZ + I
         D(I) = A(IPOS)
         A(IPOS) = H*SCALE*SCALE
  140 CONTINUE
      RETURN
      END
      SUBROUTINE F02BEF(N,C,ALB,UB,ACHEPS,EPS,B,BETA,M,MM,ROOT,VEC,IVEC,
     *                  ICOUNT,X,LOG,IFAIL)
C     NAG COPYRIGHT 1976. MARK  5 RELEASE.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 15A REVISED. IER-911 (APR 1991).
C     WRITTEN BY W.PHILLIPS     1ST OCTOBER 1975
C     OXFORD UNIVERSITY COMPUTING LABOTATORY.
C     THIS ROUTINE REPLACES F02ASF.
C
C     TRISTURM
C     C IS THE DIAGONAL, B THE SUB-DIAGONAL AND BETA THE SQUARED
C     SUBDIAGONAL OF A SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N. THE
C     EIGENVALUES WHICH ARE LESS THAN UB AND NOT LESS THAN ALB ARE
C     CALCULATED BY THE METHOD OF BISECTION AND STORED IN THE
C     VECTOR ROOT(M). THE SUBROUTINE FAILS IF M ON ENTRY IS LESS
C     THAN THE NUMBER OF EIGENVALUES REQUIRED AND ON EXIT MM GIVES
C     THE ACTUAL NUMBER OF EIGENVALUES FOUND. THE CORRESPONDING
C     EIGENVECTORS ARE CALCULATED BY INVERSE ITERATION AND ARE
C     STORED IN THE ARRAY VEC(N,M), NORMALISED SO THAT THE SUM
C     OF SQUARES IS 1, WITH THE NUMBER OF ITERATIONS STORED IN THE
C     VECTOR ICOUNT(M). THE SUBROUTINE FAILS IF ANY VECTOR HAS NOT
C     BEEN ACCEPTED AFTER 5 ITERATIONS. ELEMENTS OF B REGARDED AS
C     NEGLIGIBLE AND THE CORRESPONDING BETA ARE REPLACED BY ZERO.
C     ACHEPS IS THE RELATIVE MACHINE PRECISION AND EPS SHOULD BE
C     EQUAL TO THE ERROR TO BE TOLERATED IN THE SMALLEST
C     EIGENVALUE..
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02BEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACHEPS, ALB, EPS, UB
      INTEGER           IFAIL, IVEC, M, MM, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(N), BETA(N), C(N), ROOT(M), VEC(IVEC,M),
     *                  X(N,7)
      INTEGER           ICOUNT(M)
      LOGICAL           LOG(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANORM, BI, EPS1, EPS2, EPS3, EPS4, HALF, ONE,
     *                  THOU, TWO, U, V, X1, XM, XO, XU, ZERO
      INTEGER           I, IGROUP, II, IK, IP, IP1, IQ, IQ1, IR, IR1,
     *                  IRGROU, IS, ISAVE, ITS, J, K, M1, M2, N1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           F02BEZ, P01ABF
      EXTERNAL          F02BEZ, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT, DBLE
C     .. Data statements ..
      DATA              ONE/1.0D0/, ZERO/0.0D0/, HALF/0.5D0/,
     *                  TWO/2.0D0/, THOU/1.0D-3/
C     .. Executable Statements ..
      EPS1 = EPS
      ISAVE = IFAIL
      IF (N.LT.2) GO TO 40
C     LOOK FOR SMALL SUB-DIAGONAL ENTRIES
      DO 20 I = 2, N
         IF (ABS(B(I)).GT.ACHEPS*(ABS(C(I))+ABS(C(I-1)))) GO TO 20
         BETA(I) = ZERO
   20 CONTINUE
   40 MM = F02BEZ(1,N,C,BETA,N,UB,ACHEPS) -
     *     F02BEZ(1,N,C,BETA,N,ALB,ACHEPS)
      IF (MM.GT.M) GO TO 1120
      IQ = 0
      IR = 1
   60 IP = IQ + 1
      N1 = N - 1
      IF (N1.LT.IP) GO TO 100
      DO 80 IQ = IP, N1
         IF (BETA(IQ+1).EQ.ZERO) GO TO 120
   80 CONTINUE
  100 IQ = N
  120 IF (IP.NE.IQ) GO TO 160
      IF (ALB.GT.C(IP) .OR. C(IP).GE.UB) GO TO 1100
      DO 140 I = 1, N
         VEC(I,IR) = ZERO
  140 CONTINUE
      ROOT(IR) = C(IP)
      VEC(IP,IR) = ONE
      ICOUNT(IR) = 0
      IR = IR + 1
      GO TO 1100
  160 CONTINUE
      IF (EPS1.LE.ZERO) THEN
         XU = ABS(C(IP)) + ABS(B(IP+1))
         DO 170 I = IP + 1, IQ - 1
            XU = MAX(XU,(ABS(C(I))+ABS(B(I))+ABS(B(I+1))))
  170    CONTINUE
         XU = MAX(XU,(ABS(C(IQ))+ABS(B(IQ))))
         EPS1 = XU*EPS
      END IF
      M1 = F02BEZ(IP,IQ,C,BETA,N,ALB,ACHEPS) + 1
      M2 = F02BEZ(IP,IQ,C,BETA,N,UB,ACHEPS)
      IF (M1.GT.M2) GO TO 1100
C     FIND ROOTS BY BISECTION.
      XO = UB
      DO 180 I = M1, M2
         X(I,1) = UB
         X(I,2) = ALB
  180 CONTINUE
C     LOOP FOR K-TH EIGENVALUE.
      DO 320 IK = M1, M2
         K = M2 + M1 - IK
         XU = ALB
         IF (K.LT.M1) GO TO 220
         DO 200 II = M1, K
            I = M1 + K - II
            IF (XU.GE.X(I,2)) GO TO 200
            XU = X(I,2)
            GO TO 220
  200    CONTINUE
  220    IF (XO.GT.X(K,1)) XO = X(K,1)
  240    X1 = (XU+XO)*HALF
         XM = TWO*ACHEPS*(ABS(XU)+ABS(XO)) + EPS1
         IF (XM.GE.(XO-XU)) GO TO 300
         IS = F02BEZ(IP,IQ,C,BETA,N,X1,ACHEPS)
         IF (IS.GE.K) GO TO 280
         IF (IS.GE.M1) GO TO 260
         XU = X1
         X(M1,2) = X1
         GO TO 240
  260    XU = X1
         X(IS+1,2) = X1
         IF (X(IS,1).GT.X1) X(IS,1) = X1
         GO TO 240
  280    XO = X1
         GO TO 240
  300    X(K,1) = (XO+XU)*HALF
  320 CONTINUE
C     FIND VECTORS BY INVERSE ITERATION.
      ANORM = ABS(C(IP))
      IP1 = IP + 1
      IF (IQ.LT.IP1) GO TO 360
      DO 340 I = IP1, IQ
         ANORM = ANORM + ABS(C(I)) + ABS(B(I))
  340 CONTINUE
C     EPS2 IS THE CRITERION FOR GROUPING,
C     EPS3 REPLACES ZERO PIVOTS AND EQUAL ROOTS ARE
C     MODIFIED BY EPS3,
C     EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW,
  360 EPS2 = ANORM*THOU
      EPS3 = ACHEPS*ANORM
      EPS4 = EPS3*(DBLE(IQ)-DBLE(IP)+ONE)
      IGROUP = 0
      IS = IP
      IF (M2.LT.M1) GO TO 1100
      DO 1080 K = M1, M2
         ITS = 1
         ROOT(IR) = X(K,1)
         X1 = X(K,1)
C        LOOK FOR CLOSE OR COINCIDENT ROOTS.
         IF (K.EQ.M1) GO TO 420
         IF (X1-XO.GE.EPS2) GO TO 380
         IGROUP = IGROUP + 1
         GO TO 400
  380    IGROUP = 0
  400    IF (X1.LE.XO) X1 = XO + EPS3
  420    U = EPS4/SQRT(DBLE(IQ)-DBLE(IP)+ONE)
         IF (IQ.LT.IP) GO TO 460
         DO 440 I = IP, IQ
            X(I,7) = U
  440    CONTINUE
C        ELIMINATION WITH INTERCHANGES.
  460    U = C(IP) - X1
         V = B(IP+1)
         IP1 = IP + 1
         IF (IQ.LT.IP1) GO TO 560
         DO 540 I = IP1, IQ
            J = I - 1
            BI = B(I)
            LOG(I) = ABS(BI) .GE. ABS(U)
            IF (LOG(I)) GO TO 480
            X(I,6) = BI/U
            XU = BI/U
            X(J,3) = U
            X(J,4) = V
            X(J,5) = ZERO
            U = C(I) - X1 - XU*V
            IF (I.NE.IQ) V = B(I+1)
            GO TO 540
  480       X(I,6) = U/BI
            XU = U/BI
            X(J,3) = BI
            X(J,4) = C(I) - X1
            IF (I.EQ.IQ) GO TO 500
            X(J,5) = B(I+1)
            GO TO 520
  500       X(J,5) = ZERO
  520       U = V - XU*X(J,4)
            V = -XU*X(J,5)
  540    CONTINUE
  560    X(IQ,3) = U
         IF (U.EQ.ZERO) X(IQ,3) = EPS3
         X(IQ,4) = ZERO
         X(IQ,5) = ZERO
C        BACKSUBSTITUTION.
  580    IF (IQ.LT.IP) GO TO 620
         DO 600 II = IP, IQ
            I = IP + IQ - II
            X(I,7) = (X(I,7)-U*X(I,4)-V*X(I,5))/X(I,3)
            V = U
            U = X(I,7)
  600    CONTINUE
  620    IRGROU = IR - IGROUP
         IR1 = IR - 1
         IF (IRGROU.GT.IR1) GO TO 700
         DO 680 J = IRGROU, IR1
C           ORTHOGONALISE WITH RESPECT TO PREVIOUS MEMBERS OF GROUP
            XU = ZERO
            IF (IQ.LT.IP) GO TO 680
            DO 640 I = IP, IQ
               XU = XU + X(I,7)*VEC(I,J)
  640       CONTINUE
            DO 660 I = IP, IQ
               X(I,7) = X(I,7) - XU*VEC(I,J)
  660       CONTINUE
  680    CONTINUE
  700    ANORM = ZERO
         IF (IQ.LT.IP) GO TO 740
         DO 720 I = IP, IQ
            ANORM = ANORM + ABS(X(I,7))
  720    CONTINUE
C        FORWARD SUBSTITUTION.
  740    IF (ANORM.GE.ONE) GO TO 920
         IF (ITS.NE.5) GO TO 760
         ICOUNT(IR) = 6
         GO TO 1140
  760    IF (ANORM.NE.ZERO) GO TO 800
         X(IS,7) = EPS4
         IF (IS.EQ.IQ) GO TO 780
         IS = IS + 1
         GO TO 840
  780    IS = IP
         GO TO 840
  800    XU = EPS4/ANORM
         IF (IQ.LT.IP) GO TO 840
         DO 820 I = IP, IQ
            X(I,7) = X(I,7)*XU
  820    CONTINUE
  840    IP1 = IP + 1
         IF (IQ.LT.IP1) GO TO 900
         DO 880 I = IP1, IQ
            J = I - 1
            IF (LOG(I)) GO TO 860
            X(I,7) = X(I,7) - X(I,6)*X(J,7)
            GO TO 880
  860       U = X(J,7)
            X(J,7) = X(I,7)
            X(I,7) = U - X(I,6)*X(J,7)
  880    CONTINUE
  900    ITS = ITS + 1
         GO TO 580
C        NORMALISE SO THAT SUM OF SQUARES IS 1 AND EXPAND
C        TO FULL ORDER.
  920    U = ZERO
         IF (IQ.LT.IP) GO TO 980
         DO 940 I = IP, IQ
            U = U + X(I,7)**2
  940    CONTINUE
         XU = ONE/SQRT(U)
         DO 960 I = IP, IQ
            VEC(I,IR) = X(I,7)*XU
  960    CONTINUE
  980    IP1 = IP - 1
         IF (1.GT.IP1) GO TO 1020
         DO 1000 II = 1, IP1
            I = 1 + IP1 - II
            VEC(I,IR) = ZERO
 1000    CONTINUE
 1020    IQ1 = IQ + 1
         IF (IQ1.GT.N) GO TO 1060
         DO 1040 I = IQ1, N
            VEC(I,IR) = ZERO
 1040    CONTINUE
 1060    ICOUNT(IR) = ITS
         IR = IR + 1
         XO = X1
 1080 CONTINUE
 1100 IF (IQ.LT.N) GO TO 60
      IFAIL = 0
      RETURN
 1120 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
 1140 IFAIL = P01ABF(ISAVE,2,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F01AZF(N,M1,M2,A,IA,Z,IZI)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     TRBAK3
C     THIS SUBROUTINE PERFORMS, ON THE MATRIX OF EIGENVECTORS, Z,
C     STORED
C     IN COLUMNS M1 TO M2 OF THE ARRAY Z(N,N), A BACKTRANSFORMATION
C     TO
C     FORM THE EIGENVECTORS OF THE ORIGINAL SYMMETRIC MATRIX FROM
C     THE EIGENVECTORS OF THE TRIDIAGONAL MATRIX. THE NEW VECTORS
C     ARE OVERWRITTEN ON THE OLD ONES. THE DETAILS OF THE
C     HOUSEHOLDER REDUCTION MUST BE STORED IN THE ARRAY
C     A(N(N+1)/2), AS LEFT BY THE SUBROUTINE TRED3, F01AYF. IF Z
C     DENOTES ANY COLUMN OF THE RESULTANT MATRIX Z, THEN Z
C     SATISFIES ZT*Z = Z(INPUT)T * Z(INPUT).
C     IST. MARCH  1972
C
C     .. Scalar Arguments ..
      INTEGER           IA, IZI, M1, M2, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA), Z(IZI,M2)
C     .. Local Scalars ..
      DOUBLE PRECISION  H, S
      INTEGER           I, IPOS, IZ, J, K, L
C     .. Executable Statements ..
      IF (N.EQ.1) RETURN
      DO 80 I = 2, N
         L = I - 1
         IZ = I*L/2
         IPOS = IZ + I
         H = A(IPOS)
         IF (H.EQ.0.0D0) GO TO 80
         DO 60 J = M1, M2
            S = 0.0D0
            DO 20 K = 1, L
               IPOS = IZ + K
               S = S + A(IPOS)*Z(K,J)
   20       CONTINUE
            S = S/H
            DO 40 K = 1, L
               IPOS = IZ + K
               Z(K,J) = Z(K,J) - S*A(IPOS)
   40       CONTINUE
   60    CONTINUE
   80 CONTINUE
      RETURN
      END

c NIVEAU 1

      
      DOUBLE PRECISION FUNCTION F06EAF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DDOT
      ENTRY                     DDOT  ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER                           INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  X( * ), Y( * )
C     ..
C
C  F06EAF returns the value
C
C     F06EAF = x'y
C
C
C  Nag Fortran 77 version of the Blas routine DDOT.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ZERO
      PARAMETER           ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      SUM
      INTEGER               I, IX, IY
C     ..
C     .. Executable Statements ..
      SUM = ZERO
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               SUM = SUM + X( IX )*Y( IX )
   10       CONTINUE
         ELSE
            IF( INCY.GE.0 )THEN
               IY = 1
            ELSE
               IY = 1 - ( N - 1 )*INCY
            END IF
            IF( INCX.GT.0 )THEN
               DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  SUM = SUM + X( IX )*Y( IY )
                  IY  = IY  + INCY
   20          CONTINUE
            ELSE
               IX = 1 - ( N - 1 )*INCX
               DO 30, I = 1, N
                  SUM = SUM + X( IX )*Y( IY )
                  IX  = IX  + INCX
                  IY  = IY  + INCY
   30          CONTINUE
            END IF
         END IF
      END IF
C
      F06EAF = SUM
      RETURN
C
C     End of F06EAF. ( DDOT )
C
      END
      INTEGER FUNCTION F06JLF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      INTEGER          IDAMAX
      ENTRY            IDAMAX( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                  INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION         X( * )
C     ..
C
C  F06JLF returns the smallest value of i such that
C
C     abs( x( i ) ) = max( abs( x( j ) ) )
C                      j
C
C
C  Nag Fortran 77 version of the Blas routine IDAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION         XMAX
      INTEGER                  I, IMAX, IX
C     .. Intrinsic Functions ..
      INTRINSIC                ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IMAX = 1
         IF( N.GT.1 )THEN
            XMAX = ABS( X( 1 ) )
            IX   = 1
            DO 10, I = 2, N
               IX = IX + INCX
               IF( XMAX.LT.ABS( X( IX ) ) )THEN
                  XMAX = ABS( X( IX ) )
                  IMAX = I
               END IF
   10       CONTINUE
         END IF
      ELSE
         IMAX = 0
      END IF
C
      F06JLF = IMAX
      RETURN
C
C     End of F06JLF. ( IDAMAX )
C
      END
      DOUBLE PRECISION FUNCTION X02AMF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE 'SAFE RANGE' PARAMETER
C     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
C     FOR ANY X WHICH SATISFIES X.GE.Z AND X.LE.1/Z
C     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
C     ERROR
C
C        -X
C        1.0/X
C        SQRT(X)
C        LOG(X)
C        EXP(LOG(X))
C        Y**(LOG(X)/LOG(Y)) FOR ANY Y
C
      DOUBLE PRECISION X02CON
      DATA X02CON /2.22507385850721D-308 /
C     .. Executable Statements ..
      X02AMF = X02CON
      RETURN
      END
      SUBROUTINE F06EDF( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSCAL ( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06EDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine DSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   10       CONTINUE
         ELSE IF( ALPHA.EQ.( -ONE ) )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = -X( IX )
   20       CONTINUE
         ELSE IF( ALPHA.NE.ONE )THEN
            DO 30, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   30       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06EDF. ( DSCAL )
C
      END
      SUBROUTINE F06PEF( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DSPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DSPMV  performs the matrix-vector operation
C
C     y := alpha*A*x + beta*y,
C
C  where alpha and beta are scalars, x and y are n element vectors and
C  A is an n by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the matrix A is supplied in the packed
C           array AP as follows:
C
C              UPLO = 'U' or 'u'   The upper triangular part of A is
C                                  supplied in AP.
C
C              UPLO = 'L' or 'l'   The lower triangular part of A is
C                                  supplied in AP.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  AP     - DOUBLE PRECISION array of DIMENSION at least
C           ( ( n*( n + 1 ) )/2 ).
C           Before entry with UPLO = 'U' or 'u', the array AP must
C           contain the upper triangular part of the symmetric matrix
C           packed sequentially, column by column, so that AP( 1 )
C           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
C           and a( 2, 2 ) respectively, and so on.
C           Before entry with UPLO = 'L' or 'l', the array AP must
C           contain the lower triangular part of the symmetric matrix
C           packed sequentially, column by column, so that AP( 1 )
C           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
C           and a( 3, 1 ) respectively, and so on.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y. On exit, Y is overwritten by the updated
C           vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 6
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PEF/DSPMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set up the start points in  X  and  Y.
C
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of the array AP
C     are accessed sequentially with one pass through AP.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      KK = 1
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  y  when AP contains the upper triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               K     = KK
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*AP( K )
                  TEMP2  = TEMP2  + AP( K )*X( I )
                  K      = K      + 1
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*AP( KK + J - 1 ) + ALPHA*TEMP2
               KK     = KK     + J
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, K = KK, KK + J - 2
                  Y( IY ) = Y( IY ) + TEMP1*AP( K )
                  TEMP2   = TEMP2   + AP( K )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*AP( KK + J - 1 ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
               KK      = KK      + J
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y  when AP contains the lower triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*AP( KK )
               K      = KK           + 1
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*AP( K )
                  TEMP2  = TEMP2  + AP( K )*X( I )
                  K      = K      + 1
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
               KK     = KK     + ( N - J + 1 )
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*AP( KK )
               IX      = JX
               IY      = JY
               DO 110, K = KK + 1, KK + N - J
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*AP( K )
                  TEMP2   = TEMP2   + AP( K )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
               KK      = KK      + ( N - J + 1 )
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PEF (DSPMV ).
C
      END
      SUBROUTINE F06PSF( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DSPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DSPR2  performs the symmetric rank 2 operation
C
C     A := alpha*x*y' + alpha*y*x' + A,
C
C  where alpha is a scalar, x and y are n element vectors and A is an
C  n by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the matrix A is supplied in the packed
C           array AP as follows:
C
C              UPLO = 'U' or 'u'   The upper triangular part of A is
C                                  supplied in AP.
C
C              UPLO = 'L' or 'l'   The lower triangular part of A is
C                                  supplied in AP.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  AP     - DOUBLE PRECISION array of DIMENSION at least
C           ( ( n*( n + 1 ) )/2 ).
C           Before entry with  UPLO = 'U' or 'u', the array AP must
C           contain the upper triangular part of the symmetric matrix
C           packed sequentially, column by column, so that AP( 1 )
C           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
C           and a( 2, 2 ) respectively, and so on. On exit, the array
C           AP is overwritten by the upper triangular part of the
C           updated matrix.
C           Before entry with UPLO = 'L' or 'l', the array AP must
C           contain the lower triangular part of the symmetric matrix
C           packed sequentially, column by column, so that AP( 1 )
C           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
C           and a( 3, 1 ) respectively, and so on. On exit, the array
C           AP is overwritten by the lower triangular part of the
C           updated matrix.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PSF/DSPR2 ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Set up the start points in X and Y if the increments are not both
C     unity.
C
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
C
C     Start the operations. In this version the elements of the array AP
C     are accessed sequentially with one pass through AP.
C
      KK = 1
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  A  when upper triangle is stored in AP.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  K     = KK
                  DO 10, I = 1, J
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   10             CONTINUE
               END IF
               KK = KK + J
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = KX
                  IY    = KY
                  DO 30, K = KK, KK + J - 1
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
                     IX      = IX      + INCX
                     IY      = IY      + INCY
   30             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + J
   40       CONTINUE
         END IF
      ELSE
C
C        Form  A  when lower triangle is stored in AP.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  K     = KK
                  DO 50, I = J, N
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   50             CONTINUE
               END IF
               KK = KK + N - J + 1
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = JX
                  IY    = JY
                  DO 70, K = KK, KK + N - J
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
                     IX      = IX      + INCX
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + N - J + 1
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PSF (DSPR2 ).
C
      END
      INTEGER FUNCTION F02BEZ(IP,IQ,D,F,N,AMBDA,ACHEPS)
C     NAG COPYRIGHT 1976. MARK  5 RELEASE.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY W.PHILLIPS     1ST OCTOBER 1975
C     OXFORD UNIVERSITY COMPUTING LABORATORY.
C     THIS ROUTINE REPLACES F02ASZ.
C
C     STURMCNT
C     AUXILIARY ROUTINE CALLED BY F02BEF.
C     CALCULATES THE NUMBER OF EIGENVALUES LESS THAN AMBDA
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        ACHEPS, AMBDA
      INTEGER                 IP, IQ, N
C     .. Array Arguments ..
      DOUBLE PRECISION        D(N), F(N)
C     .. Local Scalars ..
      DOUBLE PRECISION        ONE, X, ZERO
      INTEGER                 I, ICOUNT
C     .. Intrinsic Functions ..
      INTRINSIC               SQRT
C     .. Data statements ..
      DATA                    ONE/1.0D0/, ZERO/0.0D0/
C     .. Executable Statements ..
      ICOUNT = 0
      X = ONE
      DO 60 I = IP, IQ
         IF (X.EQ.ZERO) GO TO 20
         X = D(I) - AMBDA - (F(I)/X)
         GO TO 40
   20    X = D(I) - AMBDA - (SQRT(F(I)))/ACHEPS
   40    IF (X.GE.ZERO) GO TO 60
         ICOUNT = ICOUNT + 1
   60 CONTINUE
      F02BEZ = ICOUNT
      RETURN
      END
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END


      SUBROUTINE C06FJF(NDIM,ND,N,X,Y,WORK,LWORK,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     MULTI-DIMENSIONAL DISCRETE FOURIER TRANSFORM OF A
C     MULTI-DIMENSIONAL SEQUENCE OF COMPLEX DATA VALUES
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FJF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LWORK, N, NDIM
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(LWORK), X(N), Y(N)
      INTEGER           ND(NDIM)
C     .. Local Scalars ..
      INTEGER           IFAIL1, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06FFF
C     .. Executable Statements ..
      IF (NDIM.LT.1) GO TO 40
      DO 20 L = 1, NDIM
         IFAIL1 = 1
         CALL C06FFF(NDIM,L,ND,N,X,Y,WORK,LWORK,IFAIL1)
         IF (IFAIL1.NE.0) GO TO 60
   20 CONTINUE
      IFAIL = 0
      RETURN
   40 IFAIL1 = 1
   60 IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE C06GCF(Y,PTS,IFAIL)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COMPLEX CONJUGATE
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06GCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(PTS)
C     .. Local Scalars ..
      INTEGER           IERROR, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      IF (PTS.LE.0) GO TO 40
      IERROR = 0
      DO 20 J = 1, PTS
         Y(J) = -Y(J)
   20 CONTINUE
      GO TO 60
C
   40 IERROR = 1
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE E01SAF(M,X,Y,F,TRIANG,GRADS,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Creates a Thiessen triangulation of (M,X,Y) and computes
C     derivatives required for interpolation.
C
C     Input arguments:
C     M is the number of data points (nodes).
C     X, Y, F are the scattered data to be interpolated, F = F(X,Y).
C
C     Output arguments:
C     TRIANG contains the data structure defining the triangulation;
C      used as parameters IADJ and IEND to subroutine E01SAZ.
C     GRADS contains the estimated partial derivatives at the nodes;
C      first row contains derivatives with respect to X, second row
C      with respect to Y. Used as parameter ZXZY to subroutine E01SAZ.
C
C     Parameters:
C     TOL is a convergence criterion for estimating gradients at nodes;
C      TOL .ge. 0.0; TOL = 0.01 is sufficient.
C     MAXIT is the maximum number of iterations allowed for computation
C      of estimated gradients; MAXIT .ge. 0.
C
C     M, X, Y, F, TRIANG and GRADS should be used as input to NAG
C     library routine E01SBF to compute interpolated values of F(X,Y).
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01SAF')
      DOUBLE PRECISION  TOL, ZERO
      PARAMETER         (TOL=0.01D0,ZERO=0.0D0)
      INTEGER           MAXIT
      PARAMETER         (MAXIT=6)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M
C     .. Array Arguments ..
      DOUBLE PRECISION  F(M), GRADS(2,M), X(M), Y(M)
      INTEGER           TRIANG(7*M)
C     .. Local Scalars ..
      INTEGER           COINC1, COINC2, I, IER, NIT, NREC, WK1, WK2
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E01SAY, E01SAZ
C     .. Executable Statements ..
      IER = 1
      NREC = 1
      IF (M.LT.3) THEN
         WRITE (REC,FMT=99999) M
      ELSE
C        Split the triangulation information array.
         WK1 = 1
         WK2 = 6*M + 1
C        Set up the Thiessen triangulation.
         IER = 0
         NREC = 0
         CALL E01SAY(M,X,Y,TRIANG(WK1),TRIANG(WK2),IER)
         IF (IER.EQ.2) THEN
            NREC = 1
            WRITE (REC,FMT=99998)
         ELSE
C           Initialise gradients to zero for E01SAZ.
            DO 20 I = 1, M
               GRADS(1,I) = ZERO
               GRADS(2,I) = ZERO
   20       CONTINUE
C           Estimate the gradients at the nodes.
            NIT = MAXIT
            CALL E01SAZ(M,X,Y,F,TRIANG(WK1),TRIANG(WK2),TOL,NIT,GRADS,
     *                  COINC1,COINC2,IER)
            IF (IER.EQ.3) THEN
               NREC = 2
               WRITE (REC,FMT=99997) COINC1, COINC2, X(COINC1),
     *           Y(COINC1)
            END IF
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, M .lt. 3: M =',I16,'.')
99998 FORMAT (1X,'** On entry, all (X,Y) pairs are collinear.')
99997 FORMAT (1X,'** On entry, (X(I),Y(I)) = (X(J),Y(J)) for I,J =',2I8,
     *       /4X,'X(I), Y(I) = ',1P,2D13.5,' .')
      END
      SUBROUTINE E01SBF(M,X,Y,F,TRIANG,GRADS,PX,PY,PF,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Takes a Thiessen triangulation of a set of points in the plane,
C     and estimated partial derivatives with respect to X and Y, as
C     returned by NAG Fortran Library routine E01SAF. Returns the
C     value of a C1 function F(X,Y) interpolating the points and their
C     partial derivatives, evaluated at the point (PX,PY).
C     If (PX,PY) lies outside the triangulation boundary, extrapolation
C     is performed. Interpolation is exact for quadratic data.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01SBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PF, PX, PY
      INTEGER           IFAIL, M
C     .. Array Arguments ..
      DOUBLE PRECISION  F(M), GRADS(2,M), X(M), Y(M)
      INTEGER           TRIANG(7*M)
C     .. Local Scalars ..
      DOUBLE PRECISION  DUM1, DUM2
      INTEGER           IER, IST, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E01SBZ
C     .. Save statement ..
      SAVE              IST
C     .. Data statements ..
      DATA              IST/1/
C     .. Executable Statements ..
      IER = 0
      NREC = 0
      CALL E01SBZ(M,PX,PY,X,Y,F,TRIANG(1),TRIANG(6*M+1),GRADS,IST,0,PF,
     *            DUM1,DUM2,IER)
      IF (IER.EQ.1) THEN
         NREC = 1
         WRITE (REC,FMT=99999) M
      ELSE IF (IER.EQ.2) THEN
         NREC = 2
         WRITE (REC,FMT=99998)
      ELSE IF (IER.EQ.3) THEN
         NREC = 2
         WRITE (REC,FMT=99997) PX, PY
      END IF
C
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, M .lt. 3: M =',I16,'.')
99998 FORMAT (1X,'** On entry, TRIANG does not contain a valid data po',
     *       'int triangulation;',/4X,'TRIANG may have been corrupted ',
     *       'since the call to E01SAF.')
99997 FORMAT (1X,'** Warning - the evaluation point (',1P,D12.4,',',
     *       D12.4,') lies outside the',/4X,'triangulation boundary. T',
     *       'he returned value was computed by extrapolation.')
      END
      SUBROUTINE C06FFF(NDIM,L,ND,N,X,Y,WORK,LWORK,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     DISCRETE FOURIER TRANSFORM OF ONE VARIABLE IN A
C     MULTI-VARIABLE SEQUENCE OF COMPLEX DATA VALUES
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, L, LWORK, N, NDIM
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(LWORK), X(N), Y(N)
      INTEGER           ND(NDIM)
C     .. Local Scalars ..
      INTEGER           I, IFAIL1, NI, NK, NL
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06FFZ
C     .. Executable Statements ..
      IF (NDIM.LT.1) GO TO 60
      IF (L.LT.1 .OR. L.GT.NDIM) GO TO 100
      NK = 1
      NI = 1
      DO 20 I = 1, NDIM
         IF (ND(I).LT.1) GO TO 120
         IF (I.GT.L) NK = NK*ND(I)
         IF (I.LT.L) NI = NI*ND(I)
   20 CONTINUE
      NL = ND(L)
      IF (NI*NL*NK.NE.N) GO TO 80
      IF (NL.EQ.1) GO TO 40
      IF (LWORK.LT.3*NL) GO TO 160
      CALL C06FFZ(X,Y,NI,NL,NK,WORK,WORK(NL+1),WORK(2*NL+1),IFAIL1)
      IF (IFAIL1.NE.0) GO TO 140
   40 IFAIL = 0
      RETURN
   60 IFAIL1 = 1
      GO TO 180
   80 IFAIL1 = 2
      GO TO 180
  100 IFAIL1 = 3
      GO TO 180
  120 IFAIL1 = 3 + 10*I
      GO TO 180
  140 IFAIL1 = IFAIL1 + 10*L
      GO TO 180
  160 IFAIL1 = 4 + 10*L
  180 IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE E01SAY(N,X,Y,IADJ,IEND,IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C     Original name: TRMESH.
C     This routine creates a Thiessen triangulation of N
C     arbitrarily spaced points in the plane referred to as
C     nodes.  The triangulation is optimal in the sense that it
C     is as nearly equiangular as possible.  E01SAY is part of
C     an interpolation package which also provides subroutines
C     to reorder the nodes, add a new node, delete an arc, plot
C     the mesh, and print the data structure.
C     Unless the nodes are already ordered in some reasonable
C     fashion, they should be sorted into ascending order for
C     increased efficiency before calling E01SAY.
C
C     Input Parameters -     N - number of nodes in the mesh.
C                            N .ge. 3.
C
C                      X,Y - N-vectors of coordinates.
C                            (X(I),Y(I)) defines node I.
C
C                     IADJ - vector of length .ge. 6*N-9.
C
C                     IEND - vector of length .ge. N.
C
C     N, X, and Y are not altered by this routine.
C
C     Output Parameters - IADJ - adjacency lists of neighbors in
C                            counterclockwise order.  The
C                            list for node I+1 follows that
C                            for node I where X and Y define
C                            the order.  The value 0 denotes
C                            the boundary (or a pseudo-node
C                            at infinity) and is always the
C                            last neighbor of a boundary
C                            node.  IADJ is unchanged if IER
C                            .ne. 0.
C
C                     IEND - pointers to the ends of
C                            adjacency lists (sets of
C                            neighbors) in IADJ.  The
C                            neighbors of node 1 begin in
C                            iadj(1).  For K .gt. 1, the
C                            neighbors of node K begin in
C                            IADJ(IEND(K-1)+1) and K has
C                            IEND(K) - IEND(K-1) neighbors
C                            including (possibly) the
C                            boundary.  IADJ(IEND(K)) .eq. 0
C                            iff node K is on the boundary.
C                            IEND is unchanged if IER = 1.
C                            If IER = 2 IEND contains the
C                            indices of a sequence of N
C                            nodes ordered from left to
C                            right where left and right are
C                            defined by assuming node 1 is
C                            to the left of node 2.
C
C                      IER - error indicator
C                            IER = 0 if no errors were
C                                    encountered.
C                            IER = 1 if N .lt. 3.
C                            IER = 2 if N .ge. 3 and all
C                                    nodes are collinear.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     NN =          local copy of N
C     K =           node (index) to be inserted into IEND
C     KM1 =         K-1 - (variable) length of IEND
C     NL,NR =       IEND(1), IEND(KM1) -- leftmost and rightmost
C                 nodes in IEND as viewed from the right of
C                 1-2 when IEND contains the initial ordered
C                 set of nodal indices
C     XL,YL,XR,YR = X and Y coordinates of NL and NR
C     DXR,DYR =     XR-XL, YR-YL
C     XK,YK =       X and Y coordinates of node K
C     DXK,DYK =     XK-XL, YK-YL
C     CPROD =       vector cross product of NL-NR and NL-K --
C                 used to determine the position of node K
C                 with respect to the line defined by the
C                 nodes in IEND
C     SPROD =       scalar product used to determine the
C                 interval containing node K when K is on
C                 the line defined by the nodes in IEND
C     IND,INDX =    indices for IEND and IADJ, respectively
C     N0,ITEMP =    temporary nodes (indices)
C     IERR =        dummy parameter for call to E01SAX
C     KM1D2,KMI,I = KM1/2, K-I, do-loop index -- used in IEND
C                 reordering loop
C     KMIN =        first node index sent to E01SAX
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           IER, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N)
      INTEGER           IADJ(6*N-9), IEND(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  CPROD, DXK, DXR, DYK, DYR, SPROD, XK, XL, XR,
     *                  YK, YL, YR
      INTEGER           I, IERR, IND, INDX, ITEMP, K, KM1, KM1D2, KMI,
     *                  KMIN, N0, NL, NN, NR
C     .. External Subroutines ..
      EXTERNAL          E01SAQ, E01SAX
C     .. Executable Statements ..
      NN = N
      IER = 1
      IF (NN.GE.3) THEN
         IER = 0
C
C        Initialize IEND, NL, NR, and K
C
         IEND(1) = 1
         IEND(2) = 2
         XL = X(1)
         YL = Y(1)
         XR = X(2)
         YR = Y(2)
         K = 2
   20    CONTINUE
C
C        Begin loop on nodes 3,4,...
C
         DXR = XR - XL
         DYR = YR - YL
   40    CONTINUE
C
C        Next loop begins here if NL and NR are unchanged
C
         IF (K.EQ.NN) THEN
            GO TO 260
         ELSE
            KM1 = K
            K = KM1 + 1
            XK = X(K)
            YK = Y(K)
            DXK = XK - XL
            DYK = YK - YL
            CPROD = DXR*DYK - DXK*DYR
            IF (CPROD.GT.ZERO) THEN
               GO TO 120
            ELSE IF (CPROD.LT.ZERO) THEN
               GO TO 160
            ELSE
C
C              Node K lies on the line containing nodes 1,2,...,K-1.
C              Set SPROD to (NL-NR,NL-K).
C
               SPROD = DXR*DXK + DYR*DYK
               IF (SPROD.GT.ZERO) THEN
C
C                 Node K is to the right of NL.  Find the leftmost node
C                 N0 which lies to the right of K.
C                 Set SPROD to (N0-NL,N0-K).
C
                  DO 60 IND = 2, KM1
                     N0 = IEND(IND)
                     SPROD = (XL-X(N0))*(XK-X(N0)) + (YL-Y(N0))
     *                       *(YK-Y(N0))
                     IF (SPROD.GE.ZERO) GO TO 80
   60             CONTINUE
                  GO TO 100
C
C                 Node K lies between IEND(IND-1) and IEND(IND).
C                 Insert K in IEND.
C
   80             CALL E01SAQ(IND,KM1,1,IEND)
                  IEND(IND) = K
                  GO TO 40
               END IF
            END IF
         END IF
C
C        Node K is to the left of NL.  Insert K as the first
C        (leftmost) node in IEND and set NL to K.
C
         CALL E01SAQ(1,KM1,1,IEND)
         IEND(1) = K
         XL = XK
         YL = YK
         GO TO 20
C
C        Node K is to the right of NR.  Insert K as the last
C        (rightmost) node in IEND and set NR to K.
C
  100    IEND(K) = K
         XR = XK
         YR = YK
         GO TO 20
C
C        Node K is to the left of NL-NR.  Reorder IEND so that NL
C        is the leftmost node as viewed from K.
C
  120    KM1D2 = KM1/2
         DO 140 I = 1, KM1D2
            KMI = K - I
            ITEMP = IEND(I)
            IEND(I) = IEND(KMI)
            IEND(KMI) = ITEMP
  140    CONTINUE
C
C        Node K is to the right of NL-NR.  Create a triangulation
C        consisting of nodes 1,2,...,K.
C
  160    NL = IEND(1)
         NR = IEND(KM1)
C
C        Create the adjacency lists for the first K-1 nodes.
C        Insert neighbors in reverse order.  Each node has four
C        neighbors except NL and NR which have three.
C
         DO 180 IND = 1, KM1
            N0 = IEND(IND)
            INDX = 4*N0
            IF (N0.GE.NL) INDX = INDX - 1
            IF (N0.GE.NR) INDX = INDX - 1
            IADJ(INDX) = 0
            INDX = INDX - 1
            IF (IND.LT.KM1) IADJ(INDX) = IEND(IND+1)
            IF (IND.LT.KM1) INDX = INDX - 1
            IADJ(INDX) = K
            IF (IND.NE.1) IADJ(INDX-1) = IEND(IND-1)
  180    CONTINUE
C
C        Create the adjacency list for node K
C
         INDX = 5*KM1 - 1
         IADJ(INDX) = 0
         DO 200 IND = 1, KM1
            INDX = INDX - 1
            IADJ(INDX) = IEND(IND)
  200    CONTINUE
C
C        Replace IEND elements with pointers to IADJ
C
         INDX = 0
         DO 220 IND = 1, KM1
            INDX = INDX + 4
            IF (IND.EQ.NL .OR. IND.EQ.NR) INDX = INDX - 1
            IEND(IND) = INDX
  220    CONTINUE
         INDX = INDX + K
         IEND(K) = INDX
C
C        Add the remaining nodes to the triangulation
C
         IF (K.NE.NN) THEN
            KMIN = K + 1
            DO 240 K = KMIN, NN
               CALL E01SAX(K,X,Y,IADJ,IEND,IERR)
  240       CONTINUE
         END IF
         RETURN
C
C        All nodes are collinear
C
  260    IER = 2
      END IF
      RETURN
      END
      SUBROUTINE E01SAZ(N,X,Y,Z,IADJ,IEND,EPS,NIT,ZXZY,COINC1,COINC2,
     *                  IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: GRADG.
C     Given a triangulation of N nodes in the plane with
C     associated data values, this routine uses a global method
C     to compute estimated gradients at the nodes.  The method
C     consists of minimizing a quadratic functional Q(G) over
C     the N-vector G of gradients where Q approximates the lin-
C     earized curvature of an interpolant F over the triangula-
C     tion.  The restriction of F to an arc of the triangulation
C     is taken to be the Hermite cubic interpolant of the data
C     values and tangential gradient components at the end-
C     points of the arc, and Q is the sum of the linearized
C     curvatures of F along the arcs -- the integrals over the
C     arcs of d2F(T)**2 where d2F(t) is the second derivative
C     of F with respect to distance T along the arc.  This min-
C     imization problem corresponds to an order 2N symmetric
C     positive-definite sparse linear system which is solved for
C     the X and Y partial derivatives by the block Gauss-Seidel
C     method with 2 by 2 blocks.
C
C     Input Parameters - N - number of nodes.  N .ge. 3.
C
C                  X,Y - cartesian coordinates of the nodes.
C
C                    Z - data values at the nodes.  Z(I) is
C                        associated with (X(I),Y(I)).
C
C            IADJ,IEND - data structure defining the trian-
C                        gulation.  See subroutine E01SAY.
C
C                  EPS - nonnegative convergence criterion.
C                        The method is terminated when the
C                        maximum change in a gradient compo-
C                        nent between iterations is at most
C                        EPS.  EPS = 1.E-2 is sufficient for
C                        effective convergence.
C
C     The above parameters are not altered by this routine.
C
C                  NIT - maximum number of Gauss-Seidel
C                        iterations to be applied.  This
C                        maximum will likely be achieved if
C                        EPS is smaller than the machine
C                        precision.  Optimal efficiency was
C                        achieved in testing with EPS = 0
C                        and NIT = 3 or 4.
C
C                 ZXZY - 2 by N array whose columns contain
C                        initial estimates of the partial
C                        derivatives (zero vectors are
C                        sufficient).
C
C     Output Parameters - NIT - number of Gauss-Seidel itera-
C                           tions employed.
C
C                    ZXZY - estimated X and Y partial deriv-
C                           atives at the nodes with X par-
C                           tials in the first row.  ZXZY is
C                           not changed if IER = 2.
C
C          COINC1, COINC2 - see description for IER = 3.
C
C                     IER - error indicator
C                           IER = 0 if the convergence cri-
C                                   terion was achieved.
C                           IER = 1 if convergence was not
C                                   achieved within nit
C                                   iterations.
C                           IER = 2 if N or EPS is out of
C                                   range or NIT .lt. 0 on
C                                   input.
C                           IER = 3 if two nodes are coincidental.
C                                   COINC1 and COINC2 are the
C                                   indices of the nodes.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     NN =          local copy of N
C     MAXIT =       input value of NIT
C     ITER =        number of iterations used
C     K =           do-loop and node index
C     INDF,INDL =   IADJ indices of the first and last neighbors
C                 of K
C     INDX =        IADJ index in the range INDF,...,INDL
C     NB =          neighbor of K
C     TOL =         local copy of EPS
C     DGMAX =       maximum change in a gradient component be-
C                 tween iterations
C     NORMG =       absolute value of largest ZXZY element.
C     XK,YK,ZK =    X(K), Y(K), Z(K)
C     ZXK,ZYK =     initial values of ZXZY(1,K) and ZXZY(2,K)
C     A11,A12,A22 = matrix components of the 2 by 2 block A*DG
C                 = R where A is symmetric, DG = (DZX,DZY)
C                 is the change in the gradient at K, and R
C                 is the residual
C     R1,R2 =       components of the residual -- derivatives of
C                 Q with respect to the components of the
C                 gradient at node K
C     DELX,DELY =   components of the arc NB-K
C     DELXS,DELYS = DELX**2, DELY**2
C     D =           the distance between K and NB
C     DCUB =        D**3
C     T =           factor of R1 and R2
C     DZX,DZY =     solution of the 2 by 2 system -- change in
C                 derivatives at K from the previous iterate
C
C     .. Parameters ..
      DOUBLE PRECISION  ONEHLF, TWO, ZERO
      PARAMETER         (ONEHLF=1.5D0,TWO=2.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           COINC1, COINC2, IER, N, NIT
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N), Z(N), ZXZY(2,N)
      INTEGER           IADJ(6*N), IEND(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A11, A12, A22, D, DCUB, DELX, DELXS, DELY,
     *                  DELYS, DGMAX, DZX, DZY, NORMG, R1, R2, T, TOL,
     *                  XK, YK, ZK, ZXK, ZYK
      INTEGER           INDF, INDL, INDX, ITER, K, MAXIT, NB, NN
C     .. External Functions ..
      DOUBLE PRECISION  F06BNF
      EXTERNAL          F06BNF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
      NN = N
      TOL = EPS
      MAXIT = NIT
      COINC1 = 0
      COINC2 = 0
      IF (NN.LT.3 .OR. TOL.LT.ZERO .OR. MAXIT.LT.0) THEN
C        Parameter out of range
         NIT = 0
         IER = 2
      ELSE
         NORMG = ZERO
         DO 20 K = 1, NN
            NORMG = MAX(NORMG,ABS(ZXZY(1,K)),ABS(ZXZY(2,K)))
   20    CONTINUE
         ITER = 0
C        Top of iteration loop
   40    IF (ITER.GT.MAXIT) THEN
C           Method has failed to converge within MAXIT iterations.
            IER = 0
         ELSE
            DGMAX = ZERO
            INDL = 0
            DO 80 K = 1, NN
               XK = X(K)
               YK = Y(K)
               ZK = Z(K)
               ZXK = ZXZY(1,K)
               ZYK = ZXZY(2,K)
C              Initialize components of the 2 by 2 system
               A11 = ZERO
               A12 = ZERO
               A22 = ZERO
               R1 = ZERO
               R2 = ZERO
C              Loop on neighbors NB of K
               INDF = INDL + 1
               INDL = IEND(K)
               DO 60 INDX = INDF, INDL
                  NB = IADJ(INDX)
                  IF (NB.NE.0) THEN
C                    Compute the components of arc NB-K
                     DELX = X(NB) - XK
                     DELY = Y(NB) - YK
                     DELXS = DELX*DELX
                     DELYS = DELY*DELY
                     D = F06BNF(DELX,DELY)
                     IF (D.EQ.ZERO) THEN
C                       Nodes K and NB are coincidental.
                        IER = 3
                        COINC1 = K
                        COINC2 = NB
                        RETURN
                     ELSE
                        DCUB = D*D*D
C                       Update the system components for node NB.
                        A11 = A11 + DELXS/DCUB
                        A12 = A12 + DELX*DELY/DCUB
                        A22 = A22 + DELYS/DCUB
                        T = (ONEHLF*(Z(NB)-ZK)-((ZXZY(1,NB)/TWO+ZXK)
     *                      *DELX+(ZXZY(2,NB)/TWO+ZYK)*DELY))/DCUB
                        R1 = R1 + T*DELX
                        R2 = R2 + T*DELY
                     END IF
                  END IF
   60          CONTINUE
C              Solve the 2 by 2 system and update DGMAX
               DZY = (A11*R2-A12*R1)/(A11*A22-A12*A12)
               DZX = (R1-A12*DZY)/A11
               DGMAX = MAX(DGMAX,ABS(DZX),ABS(DZY))
C              Update the partials at node K
               ZXZY(1,K) = ZXK + DZX
               ZXZY(2,K) = ZYK + DZY
   80       CONTINUE
C           Increment ITER and test for convergence
            ITER = ITER + 1
            IF (DGMAX.GT.TOL+TOL*NORMG) GO TO 40
C           Method has converged
            IER = 0
         END IF
      END IF
      NIT = ITER
      RETURN
      END
      SUBROUTINE E01SBZ(N,PX,PY,X,Y,Z,IADJ,IEND,ZXZY,IST,IFLAG,PZ,DZX,
     *                  DZY,IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14C REVISED. IER-875 (NOV 1990).
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: INTRC1.
C     Adapted for NAG by H. Scullion Leics. Univ.
C
C     Given a triangulation of a set of points in the plane,
C     this routine determines a piecewise cubic function F(X,Y)
C     which interpolates a set of data values and partial
C     derivatives at the vertices.  F has continuous first
C     derivatives over the mesh and extends beyond the mesh
C     boundary allowing extrapolation.  Interpolation is exact
C     for quadratic data.  The value of F at (PX,PY) is
C     returned.  E01SBZ is part of an interpolation package
C     which provides routines to generate, update and plot the
C     mesh.
C
C     Input Parameters -     N - number of nodes in the mesh.
C                            N .ge. 3.
C
C                    PX,PY - coordinates of a point at which
C                            F is to be evaluated.
C
C                      X,Y - vectors of coordinates of the
C                            nodes in the mesh.
C
C                        Z - vector of data values at the
C                            nodes.
C
C                     IADJ - set of adjacency lists of nodes
C                            in the mesh.
C
C                     IEND - pointers to the ends of
C                            adjacency lists in IADJ for
C                            each node in the mesh.
C
C                     ZXZY - 2 by N array whose columns
C                            contain estimated partial der-
C                            ivatives at the nodes (X par-
C                            tials in the first row)
C
C                      IST - index of the starting node in
C                            the search for a triangle con-
C                            taining (PX,PY).  1 .le. IST
C                            .le. N.  The output value of
C                            IST from a previous call may
C                            be a good choice.
C
C     IADJ and IEND may be created by E01SAY and derivative
C     estimates are computed by E01SAZ.
C
C     Input parameters other than IST are not altered by this
C     routine.
C
C     Output Parameters - IST - index of one of the vertices of
C                           the triangle containing (PX,PY)
C                           unless IER .gt. 0.
C
C                      PZ - value of F at (PX,PY), or 0 if
C                           IER .gt. 0.
C
C                     IER - error indicator
C                           IER = 0 if no errors were
C                                   encountered.
C                           IER = 1 if N or IST is
C                                    out of range.
C                           IER = 2 if the nodes are col-
C                                    linear.
C                           IER = 3 if extrapolation was used.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     NN =                      local copy of N
C     I1,I2,I3 =                vertices determined by E01SAW
C     IERR =                    error flag for calls to E01SBY
C     N1,N2 =                   endpoints of the closest bound-
C                             ary edge to P when P is out-
C                             side of the mesh boundary
C     INDX =                    IADJ index of N1 as a neighbor
C                             of N2
C     XP,YP =                   local copies of the coordinates
C                             of P=(PX,PY)
C     ZX1,ZY1,ZX2,ZY2,ZX3,ZY3 = X and Y derivatives at the
C                             vertices of a triangle T which
C                             contains P or at N1 and N2
C     X1,Y1,X2,Y2,X3,Y3 =       X,Y coordinates of the vertices
C                             of T or of N1 and N2
C     Z1,Z2,Z3 =                data values at the vertices of T
C     DP =                      inner product of N1-N2 and P-N2
C     U,V =                     X,Y coordinates of the vector
C                             N2-N1
C     XQ,YQ =                   X,Y coordinates of the closest
C                             boundary point to P when P is
C                             outside of the mesh boundary
C     R1,R2 =                   barycentric coordinates of Q
C                             with respect to the line seg-
C                             ment N2-N1 containing Q
C     A1,A2,B1,B2,C1,C2 =       cardinal functions for evaluat-
C                             ing the interpolatory surface
C                             at Q
C     F1,F2 =                   cubic factors used to compute
C                             the cardinal functions
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DZX, DZY, PX, PY, PZ
      INTEGER           IER, IFLAG, IST, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N), Z(N), ZXZY(2,N)
      INTEGER           IADJ(6*N), IEND(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, B1, B2, C1, C2, DP, F1, F2, R1, R2, U,
     *                  V, X1, X2, X3, XP, XQ, Y1, Y2, Y3, YP, YQ, Z1,
     *                  Z2, Z3, ZX1, ZX2, ZX3, ZY1, ZY2, ZY3
      INTEGER           I1, I2, I3, IERR, INDX, N1, N2, NN
C     .. External Subroutines ..
      EXTERNAL          E01SAW, E01SBY
C     .. Executable Statements ..
      IER = 0
      NN = N
      PZ = ZERO
      DZX = ZERO
      DZY = ZERO
      IF (NN.LT.3) THEN
C
C        N or IST out of range.
C
         IER = 1
      ELSE
         IF (IST.LT.1 .OR. IST.GT.NN) IST = 1
         XP = PX
         YP = PY
C
C        Find a triangle containing P if P is within the mesh
C        boundary
C
         CALL E01SAW(IST,XP,YP,X,Y,IADJ,IEND,I1,I2,I3)
         IF (I1.NE.0) THEN
            IST = I1
            IF (I3.LE.0) THEN
               IF (I3.EQ.0) IER = 3
C
C              P is outside of the mesh boundary.  Extrapolate to P by
C              passing a linear function of one variable through the
C              value and directional derivative (in the direction
C              P-Q) of the interpolatory surface (E01SBY) at Q where
C              Q is the closest boundary point to P.
C
C              Determine Q by traversing the boundary starting from
C              the rightmost visible node I1.
C
               N2 = I1
   20          CONTINUE
C
C              Set N1 to the last nonzero neighbor of N2 and compute DP
C
               INDX = IEND(N2) - 1
               N1 = IADJ(INDX)
               X1 = X(N1)
               Y1 = Y(N1)
               X2 = X(N2)
               Y2 = Y(N2)
               DP = (X1-X2)*(XP-X2) + (Y1-Y2)*(YP-Y2)
               IF (DP.LE.ZERO) THEN
                  GO TO 40
               ELSE IF ((XP-X1)*(X2-X1)+(YP-Y1)*(Y2-Y1).LE.ZERO) THEN
                  N2 = N1
                  GO TO 20
               END IF
C
C              The closest boundary point Q lies on N2-N1.  Compute
C              partials at N1 and N2.
C
               ZX1 = ZXZY(1,N1)
               ZY1 = ZXZY(2,N1)
               ZX2 = ZXZY(1,N2)
               ZY2 = ZXZY(2,N2)
C
C              Compute Q, its barycentric coordinates, and the cardinal
C              functions for extrapolation
C
               U = X2 - X1
               V = Y2 - Y1
               R1 = DP/(U**2+V**2)
               R2 = ONE - R1
               XQ = R1*X1 + R2*X2
               YQ = R1*Y1 + R2*Y2
               F1 = R1*R1*R2
               F2 = R1*R2*R2
               A1 = R1 + (F1-F2)
               A2 = R2 - (F1-F2)
               B1 = U*F1
               B2 = -U*F2
               C1 = V*F1
               C2 = -V*F2
C
C              Compute the value of the interpolatory surface (E01SBY)
C              at Q
C
               PZ = A1*Z(N1) + A2*Z(N2) + B1*ZX1 + B2*ZX2 + C1*ZY1 +
     *              C2*ZY2
C
C              Compute the extrapolated value at P
C
               PZ = PZ + (R1*ZX1+R2*ZX2)*(XP-XQ) + (R1*ZY1+R2*ZY2)
     *              *(YP-YQ)
               RETURN
C
C              N2 is the closest boundary point to P.  Compute partial
C              derivatives at N2.
C
   40          CONTINUE
               ZX2 = ZXZY(1,N2)
               ZY2 = ZXZY(2,N2)
C
C              Compute extrapolated value at P
C
               PZ = Z(N2) + ZX2*(XP-X2) + ZY2*(YP-Y2)
               RETURN
            ELSE
C
C              Derivatives are user provided
C
               ZX1 = ZXZY(1,I1)
               ZX2 = ZXZY(1,I2)
               ZX3 = ZXZY(1,I3)
               ZY1 = ZXZY(2,I1)
               ZY2 = ZXZY(2,I2)
               ZY3 = ZXZY(2,I3)
C
C              Set local parameters for call to E01SBY
C
               X1 = X(I1)
               Y1 = Y(I1)
               X2 = X(I2)
               Y2 = Y(I2)
               X3 = X(I3)
               Y3 = Y(I3)
               Z1 = Z(I1)
               Z2 = Z(I2)
               Z3 = Z(I3)
               CALL E01SBY(XP,YP,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,ZX1,ZX2,ZX3,
     *                     ZY1,ZY2,ZY3,IFLAG,PZ,DZX,DZY,IERR)
               IF (IERR.EQ.0) RETURN
            END IF
         END IF
C
C        Nodes are collinear
C
         IER = 2
      END IF
      RETURN
      END



c NIVEAU 2

      SUBROUTINE F06AAZ ( SRNAME, INFO )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 15 REVISED. IER-915 (APR 1991).
C     .. Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*13       SRNAME
C     ..
C
C  Purpose
C  =======
C
C  F06AAZ  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*13.
C           On entry, SRNAME specifies the name of the routine which
C           called F06AAZ.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Local Scalars ..
      INTEGER            IERR, IFAIL
      CHARACTER*4        VARBNM
C     .. Local Arrays ..
      CHARACTER*80       REC (1)
C     .. External Functions ..
      INTEGER            P01ACF
      EXTERNAL           P01ACF
C     ..
C     .. Executable Statements ..
      WRITE (REC (1),99999) SRNAME, INFO
      IF (SRNAME(1:3).EQ.'F06') THEN
         IERR = -1
         VARBNM = '    '
      ELSE
         IERR = -INFO
         VARBNM = 'INFO'
      END IF
      IFAIL = 0
      IFAIL = P01ACF (IFAIL, IERR, SRNAME(1:6), VARBNM, 1, REC)
C
      RETURN
C
99999 FORMAT ( ' ** On entry to ', A13, ' parameter number ', I2,
     $         ' had an illegal value' )
C
C     End of F06AAZ.
C
      END
      SUBROUTINE C06FFZ(X,Y,NI,NJ,NK,W1,W2,W3,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     DISCRETE FOURIER TRANSFORM OF THE 2ND VARIABLE IN A
C     3-DIMENSIONAL SEQUENCE OF COMPLEX DATA VALUES
C
C     .. Scalar Arguments ..
      INTEGER           IFAIL, NI, NJ, NK
C     .. Array Arguments ..
      DOUBLE PRECISION  W1(NJ), W2(NJ), W3(NJ), X(NI,NJ,NK), Y(NI,NJ,NK)
C     .. Local Scalars ..
      INTEGER           I, J, K
C     .. External Subroutines ..
      EXTERNAL          C06FCF
C     .. Executable Statements ..
      DO 80 K = 1, NK
         DO 60 I = 1, NI
            DO 20 J = 1, NJ
               W1(J) = X(I,J,K)
               W2(J) = Y(I,J,K)
   20       CONTINUE
            IFAIL = 1
            CALL C06FCF(W1,W2,NJ,W3,IFAIL)
            DO 40 J = 1, NJ
               X(I,J,K) = W1(J)
               Y(I,J,K) = W2(J)
   40       CONTINUE
   60    CONTINUE
   80 CONTINUE
      RETURN
      END
      SUBROUTINE E01SAQ(NFRST,NLAST,KK,IARR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: SHIFTD.
C     This routine shifts a set of contiguous elements of an
C     integer array KK positions downward (upward if KK .lt. 0).
C     The loops are unrolled in order to increase efficiency.
C
C     Input parameters - NFRST,NLAST - bounds on the portion of
C                                  IARR to be shifted.  All
C                                  elements between and
C                                  including the bounds are
C                                  shifted unless NFRST .gt.
C                                  NLAST, in which case no
C                                  shift occurs.
C
C                             KK - number of positions each
C                                  element is to be shifted.
C                                  if KK .lt. 0 shift up.
C                                  if KK .gt. 0 shift down.
C
C                           IARR - integer array of length
C                                  .ge. NLAST + max(KK,0).
C
C     NFRST, NLAST, and KK are not altered by this routine.
C
C     Output Parameter -        IARR - shifted array.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     INC =  do-loop increment (unrolling factor) -- if INC is
C          changed, statements must be added to or deleted
C          from the do-loops
C     K =    local copy of KK
C     NF =   local copy of NFRST
C     NL =   local copy of NLAST
C     NLP1 = NL + 1
C     NS =   number of shifts
C     NSL =  number of shifts done in unrolled do-loop (multiple
C          of INC)
C     I =    do-loop index and index for IARR
C     IBAK = index for downward shift of IARR
C     INDX = index for IARR
C     IMAX = bound on do-loop index
C
C     .. Scalar Arguments ..
      INTEGER           KK, NFRST, NLAST
C     .. Array Arguments ..
      INTEGER           IARR(*)
C     .. Local Scalars ..
      INTEGER           I, IBAK, IMAX, INC, INDX, K, NF, NL, NLP1, NS,
     *                  NSL
C     .. Data statements ..
      DATA              INC/5/
C     .. Executable Statements ..
      K = KK
      NF = NFRST
      NL = NLAST
      IF (NF.LE.NL .AND. K.NE.0) THEN
         NLP1 = NL + 1
         NS = NLP1 - NF
         NSL = INC*(NS/INC)
         IF (K.LT.0) THEN
C
C           Shift upward starting from the top
C
            IF (NSL.GT.0) THEN
               IMAX = NLP1 - INC
               DO 20 I = NF, IMAX, INC
                  INDX = I + K
                  IARR(INDX) = IARR(I)
                  IARR(INDX+1) = IARR(I+1)
                  IARR(INDX+2) = IARR(I+2)
                  IARR(INDX+3) = IARR(I+3)
                  IARR(INDX+4) = IARR(I+4)
   20          CONTINUE
            END IF
C
C           Perform the remaining NS-NSL shifts one at a time
C
            I = NSL + NF
   40       CONTINUE
            IF (I.LE.NL) THEN
               INDX = I + K
               IARR(INDX) = IARR(I)
               I = I + 1
               GO TO 40
            END IF
         ELSE
C
C           Shift downward starting from the bottom
C
            IF (NSL.GT.0) THEN
               DO 60 I = 1, NSL, INC
                  IBAK = NLP1 - I
                  INDX = IBAK + K
                  IARR(INDX) = IARR(IBAK)
                  IARR(INDX-1) = IARR(IBAK-1)
                  IARR(INDX-2) = IARR(IBAK-2)
                  IARR(INDX-3) = IARR(IBAK-3)
                  IARR(INDX-4) = IARR(IBAK-4)
   60          CONTINUE
            END IF
C
C           Perform the remaining NS-NSL shifts one at a time
C
            IBAK = NLP1 - NSL
   80       CONTINUE
            IF (IBAK.GT.NF) THEN
               IBAK = IBAK - 1
               INDX = IBAK + K
               IARR(INDX) = IARR(IBAK)
               GO TO 80
            END IF
         END IF
      END IF
      RETURN
      END
      SUBROUTINE E01SAX(KK,X,Y,IADJ,IEND,IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: ADNODE.
C     This routine adds node KK to a triangulation of a set
C     of points in the plane producing a new triangulation.  A
C     sequence of edge swaps is then applied to the mesh,
C     resulting in an optimal triangulation.  E01SAX is part
C     of an interpolation package which also provides routines
C     to initialize the data structure, plot the mesh, and
C     delete arcs.
C
C     Input Parameters -   KK - index of the node to be added
C                           to the mesh.  KK .ge. 4.
C
C                     X,Y - vectors of coordinates of the
C                           nodes in the mesh.  (X(I),Y(I))
C                           defines node I for I = 1,..,KK.
C
C                    IADJ - set of adjacency lists of nodes
C                           1,..,KK-1.
C
C                    IEND - pointers to the ends of
C                           adjacency lists in IADJ for
C                           each node in the mesh.
C
C     IADJ and IEND may be created by E01SAY.
C
C     KK, X, and Y are not altered by this routine.
C
C     Output Parameters - IADJ,IEND - updated with the addition
C                                 of node KK as the last
C                                 entry.
C
C                           IER - error indicator
C                                 IER = 0 if no errors
C                                         were encountered.
C                                 IER = 1 if all nodes
C                                         (including KK) are
C                                         collinear.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     K =        local copy of KK
C     KM1 =      K - 1
C     I1,I2,I3 = vertices of a triangle containing K
C     INDKF =    IADJ index of the first neighbor of K
C     INDKL =    IADJ index of the last neighbor of K
C     NABOR1 =   first neighbor of K before any swaps occur
C     IO1,IO2 =  adjacent neighbors of K defining an arc to
C              be tested for a swap
C     IN1 =      vertex opposite K -- first neighbor of IO2
C              which precedes IO1.  IN1,IO1,IO2 are in
C              counterclockwise order.
C     INDK1 =    index of IO1 in the adjacency list for K
C     IND2F =    index of the first neighbor of IO2
C     IND21 =    index of IO1 in the adjacency list for IO2
C     XK,YK =    X(K), Y(K)
C
C     .. Scalar Arguments ..
      INTEGER           IER, KK
C     .. Array Arguments ..
      DOUBLE PRECISION  X(KK), Y(KK)
      INTEGER           IADJ(*), IEND(KK)
C     .. Local Scalars ..
      DOUBLE PRECISION  XK, YK
      INTEGER           I1, I2, I3, IN1, IND21, IND2F, INDK1, INDKF,
     *                  INDKL, IO1, IO2, K, KM1, NABOR1
C     .. External Functions ..
      INTEGER           E01SAR
      LOGICAL           E01SAT
      EXTERNAL          E01SAR, E01SAT
C     .. External Subroutines ..
      EXTERNAL          E01SAS, E01SAU, E01SAV, E01SAW
C     .. Executable Statements ..
      IER = 0
      K = KK
C
C     Initialization
C
      KM1 = K - 1
      XK = X(K)
      YK = Y(K)
C
C     Add node K to the mesh
C
      CALL E01SAW(KM1,XK,YK,X,Y,IADJ,IEND,I1,I2,I3)
      IF (I1.EQ.0) THEN
C
C        All nodes are collinear
C
         IER = 1
      ELSE
         IF (I3.LE.0) THEN
            CALL E01SAU(K,I1,I2,IADJ,IEND)
         ELSE
            CALL E01SAV(K,I1,I2,I3,IADJ,IEND)
         END IF
C
C        Initialize variables for optimization of the mesh
C
         INDKF = IEND(KM1) + 1
         INDKL = IEND(K)
         NABOR1 = IADJ(INDKF)
         IO2 = NABOR1
         INDK1 = INDKF + 1
         IO1 = IADJ(INDK1)
   20    CONTINUE
C
C        Begin loop -- find the vertex opposite K
C
         IND2F = 1
         IF (IO2.NE.1) IND2F = IEND(IO2-1) + 1
         IND21 = E01SAR(IO2,IO1,IADJ,IEND)
         IF (IND2F.EQ.IND21) THEN
C
C           IN1 is the last neighbor of IO2
C
            IND21 = IEND(IO2)
            IN1 = IADJ(IND21)
            IF (IN1.EQ.0) GO TO 40
         ELSE
            IN1 = IADJ(IND21-1)
         END IF
C
C        Swap test -- if a swap occurs, two new arcs are opposite K
C              and must be tested.  INDK1 and INDKF must be
C              decremented.
C
         IF (E01SAT(IN1,K,IO1,IO2,X,Y)) THEN
            CALL E01SAS(IN1,K,IO1,IO2,IADJ,IEND)
            IO1 = IN1
            INDK1 = INDK1 - 1
            INDKF = INDKF - 1
            GO TO 20
         END IF
C
C        No swap occurred.  Reset IO2 and IO1, and test for
C        termination.
C
   40    IF (IO1.NE.NABOR1) THEN
            IO2 = IO1
            INDK1 = INDK1 + 1
            IF (INDK1.GT.INDKL) INDK1 = INDKF
            IO1 = IADJ(INDK1)
            IF (IO1.NE.0) GO TO 20
         END IF
      END IF
      RETURN
      END
      DOUBLE PRECISION FUNCTION F06BNF( A, B )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  A, B
C     ..
C
C  F06BNF returns the value
C
C     p = sqrt( a*a + b*b )
C
C  via the function name.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 17-January-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ZERO
      PARAMETER           ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      P
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
C     ..
C     .. Executable Statements ..
      IF( A.EQ.ZERO )THEN
         P = ABS( B )
      ELSE IF( B.EQ.ZERO )THEN
         P = ABS( A )
      ELSE IF( ABS( A ).GE.ABS( B ) )THEN
         P = ABS( A )*SQRT( 1 + ( B/A )**2 )
      ELSE
         P = ABS( B )*SQRT( 1 + ( A/B )**2 )
      END IF
C
      F06BNF = P
      RETURN
C
C     End of F06BNF. ( SPYTH )
C
      END
      SUBROUTINE E01SAW(NST,PX,PY,X,Y,IADJ,IEND,I1,I2,I3)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 13B REVISED. IER-656 (AUG 1988).
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: TRFIND.
C     This routine locates a point P in a Thiessen triangu-
C     lation, returning the vertex indices of a triangle which
C     contains P.  E01SAW is part of an interpolation package
C     which provides subroutines for creating the mesh.
C
C     Input Parameters -    NST - index of node at which E01SAW
C                             begins search.  Search time
C                             depends on the proximity of
C                             NST to P.
C
C                     PX,PY - X and Y-coordinates of the
C                             point to be located.
C
C                       X,Y - vectors of coordinates of
C                             nodes in the mesh.  (X(I),Y(I))
C                             defines node I for I = 1,...,N
C                             where N .ge. 3.
C
C                      IADJ - set of adjacency lists of
C                             nodes in the mesh.
C
C                      IEND - pointers to the ends of
C                             adjacency lists in IADJ for
C                             each node in the mesh.
C
C     IADJ and IEND may be created by E01SAY.
C
C     Input parameters are not altered by this routine.
C
C     Output Parameters - I1,I2,I3 - vertex indices in counter-
C                                clockwise order - vertices
C                                of a triangle containing P
C                                if P is an interior node.
C                                If P is outside of the
C                                boundary of the mesh, I1
C                                and I2 are the first (right
C                                -most) and last (leftmost)
C                                nodes which are visible
C                                from P, and I3 = 0.  If P
C                                is on the mesh boundary, I1
C                                and I2 are the first (right
C                                -most) and last (leftmost)
C                                nodes which are visible
C                                from P, and I3 = -1. If P
C                                and all of the nodes lie on
C                                a single line then I1 = I2
C                                = I3 = 0.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     XP,YP =     local variables containing PX and PY
C     N0,N1,N2 =  nodes in counterclockwise order defining a
C               cone (with vertex N0) containing P
C     N3,N4 =     nodes opposite N1-N2 and N2-N1, respectively
C     INDX,IND =  indices for IADJ
C     NF,NL =     first and last neighbors of N0 in IADJ, or
C               first (rightmost) and last (leftmost) nodes
C               visible from P when P is outside the
C               boundary
C     NEXT =      candidate for I1 or I2 when P is outside of
C               the boundary
C     LEFT =      statement function which computes the sign of
C               a cross product (Z-component).  LEFT(X1,...,
C               Y0) = .TRUE. iff (X0,Y0) is on or to the
C               left of the vector from (X1,Y1) to (X2,Y2).
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PX, PY
      INTEGER           I1, I2, I3, NST
C     .. Array Arguments ..
      DOUBLE PRECISION  X(*), Y(*)
      INTEGER           IADJ(*), IEND(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  X0, X1, X2, XP, Y0, Y1, Y2, YP
      INTEGER           IND, INDX, N0, N1, N2, N3, N4, NEXT, NF, NL
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Statement Functions ..
      LOGICAL           LEFT
C     .. Statement Function definitions ..
      LEFT(X1,Y1,X2,Y2,X0,Y0) = (X2-X1)*(Y0-Y1) .GE. (X0-X1)*(Y2-Y1)
C     .. Executable Statements ..
      XP = PX
      YP = PY
C
C     Initialize variables and find a cone containing P
C
      N0 = MAX(NST,1)
   20 CONTINUE
      INDX = IEND(N0)
      NL = IADJ(INDX)
      INDX = 1
      IF (N0.NE.1) INDX = IEND(N0-1) + 1
      NF = IADJ(INDX)
      N1 = NF
      IF (NL.NE.0) THEN
   40    CONTINUE
C
C        N0 is an interior node.  Find N1.
C
         IF ( .NOT. LEFT(X(N0),Y(N0),X(N1),Y(N1),XP,YP)) THEN
            INDX = INDX + 1
            N1 = IADJ(INDX)
            IF (N1.EQ.NL) THEN
               GO TO 160
            ELSE
               GO TO 40
            END IF
         END IF
      ELSE
C
C        N0 is a boundary node.  Set NL to the last nonzero
C        neighbor of N0.
C
         IND = IEND(N0) - 1
         NL = IADJ(IND)
         IF ( .NOT. LEFT(X(N0),Y(N0),X(NF),Y(NF),XP,YP)) THEN
            GO TO 120
         ELSE IF ( .NOT. LEFT(X(NL),Y(NL),X(N0),Y(N0),XP,YP)) THEN
            GO TO 100
         END IF
      END IF
   60 CONTINUE
C
C     P is to the left of arc N0-N1.  Initialize N2 to the next
C     neighbor of N0.
C
      INDX = INDX + 1
      N2 = IADJ(INDX)
      IF (LEFT(X(N0),Y(N0),X(N2),Y(N2),XP,YP)) THEN
         N1 = N2
         IF (N1.NE.NL) GO TO 60
      ELSE
         GO TO 180
      END IF
      IF (LEFT(X(N0),Y(N0),X(NF),Y(NF),XP,YP)) THEN
         IF (XP.NE.X(N0) .OR. YP.NE.Y(N0)) THEN
   80       CONTINUE
C
C           P is left of or on arcs N0-NB for all neighbors NB
C           of N0.
C           All points are collinear iff P is left of NB-N0 for
C           all neighbors NB of N0.  Search the neighbors of N0
C           in reverse order.  Note -- N1 = NL and INDX points to
C           NL.
C
            IF (LEFT(X(N1),Y(N1),X(N0),Y(N0),XP,YP)) THEN
               IF (N1.EQ.NF) THEN
                  GO TO 140
               ELSE
                  INDX = INDX - 1
                  N1 = IADJ(INDX)
                  GO TO 80
               END IF
            END IF
         END IF
C
C        P is to the right of N1-N0, or P=N0.  Set N0 to N1 and
C        start over.
C
         N0 = N1
         GO TO 20
      ELSE
         GO TO 160
      END IF
C
C     P is outside the boundary and N0 is the rightmost
C     visible boundary node
C
  100 I1 = N0
      GO TO 300
C
C     P is outside the boundary
C
  120 NL = N0
      GO TO 280
C
C     All points are collinear
C
  140 I1 = 0
      I2 = 0
      I3 = 0
      RETURN
C
C     P is between arcs N0-N1 and N0-NF
C
  160 N2 = NF
C
C     P is contained in a cone defined by line segments N0-N1
C     and N0-N2 where N1 is adjacent to N2
C
  180 N3 = N0
  200 CONTINUE
      IF (LEFT(X(N1),Y(N1),X(N2),Y(N2),XP,YP)) THEN
C        The algorithm believes that P is contained in
C        triangle (N1,N2,N3). However, if N1, N2, N3 and P
C        are almost collinear then P may in fact be far
C        outside the triangle, the algorithm being fooled
C        by rounding errors. The following code tests
C        whether P lies inside a rectangle containing the
C        triangle. If not, the triangle is rejected, and
C        the algorithm continues.
         IF (XP.GE.MIN(X(N1),X(N2),X(N3)) .AND.
     *       XP.LE.MAX(X(N1),X(N2),X(N3)) .AND.
     *       YP.GE.MIN(Y(N1),Y(N2),Y(N3)) .AND.
     *       YP.LE.MAX(Y(N1),Y(N2),Y(N3))) GO TO 260
      END IF
C
C     Set N4 to the first neighbor of N2 following N1
C
      INDX = IEND(N2)
      IF (IADJ(INDX).NE.N1) THEN
  220    CONTINUE
C
C        N1 is not the last neighbor of N2
C
         INDX = INDX - 1
         IF (IADJ(INDX).NE.N1) GO TO 220
         N4 = IADJ(INDX+1)
         IF (N4.EQ.0) GO TO 240
      ELSE
C
C        N1 is the last neighbor of N2.
C        Set N4 to the first neighbor.
C
         INDX = 1
         IF (N2.NE.1) INDX = IEND(N2-1) + 1
         N4 = IADJ(INDX)
      END IF
C
C     Define a new arc N1-N2 which intersects the line
C     segment N0-P
C
      IF (LEFT(X(N0),Y(N0),X(N4),Y(N4),XP,YP)) THEN
         N3 = N1
         N1 = N4
      ELSE
         N3 = N2
         N2 = N4
      END IF
      GO TO 200
C
C     P is outside the boundary
C
  240 NF = N2
      NL = N1
      GO TO 280
C
C     P is in the triangle (N1,N2,N3) and not on N2-N3.  If
C     N3-N1 or N1-N2 is a boundary arc containing P, treat P
C     as exterior.
C
  260 INDX = IEND(N1)
      IF (IADJ(INDX).EQ.0) THEN
C
C        N1 is a boundary node.  N3-N1 is a boundary arc iff N3
C        is the last nonzero neighbor of N1.
C
         IF (N3.EQ.IADJ(INDX-1)) THEN
C
C           N3-N1 is a boundary arc
C
            IF (LEFT(X(N1),Y(N1),X(N3),Y(N3),XP,YP)) THEN
C
C              P lies on N1-N3
C
               I1 = N1
               I2 = N3
               I3 = -1
               RETURN
            END IF
         END IF
C
C        N3-N1 is not a boundary arc containing P.  N1-N2 is a
C        boundary arc iff N2 is the first neighbor of N1.
C
         INDX = 1
         IF (N1.NE.1) INDX = IEND(N1-1) + 1
         IF (N2.EQ.IADJ(INDX)) THEN
C
C           N1-N2 is a boundary arc
C
            IF (LEFT(X(N2),Y(N2),X(N1),Y(N1),XP,YP)) THEN
C
C              P lies on N1-N2
C
               I1 = N2
               I2 = N1
               I3 = -1
               RETURN
            END IF
         END IF
      END IF
C
C     P does not lie on a boundary arc.
C
      I1 = N1
      I2 = N2
      I3 = N3
      RETURN
  280 CONTINUE
C
C     NF and NL are adjacent boundary nodes which are visible
C     from P.  Find the first visible boundary node.
C     Set next to the first neighbor of NF.
C
      INDX = 1
      IF (NF.NE.1) INDX = IEND(NF-1) + 1
      NEXT = IADJ(INDX)
      IF ( .NOT. LEFT(X(NF),Y(NF),X(NEXT),Y(NEXT),XP,YP)) THEN
         NF = NEXT
         GO TO 280
      END IF
C
C     NF is the first (rightmost) visible boundary node
C
      I1 = NF
  300 CONTINUE
C
C     Find the last visible boundary node.  NL is the first
C     candidate for I2.
C     Set next to the last neighbor of NL.
C
      INDX = IEND(NL) - 1
      NEXT = IADJ(INDX)
      IF ( .NOT. LEFT(X(NEXT),Y(NEXT),X(NL),Y(NL),XP,YP)) THEN
         NL = NEXT
         GO TO 300
      END IF
C
C     NL is the last (leftmost) visible boundary node
C
      I2 = NL
      I3 = 0
      RETURN
      END
      SUBROUTINE E01SBY(X,Y,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,ZX1,ZX2,ZX3,ZY1,
     *                  ZY2,ZY3,IFLAG,W,WX,WY,IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Given function values and first partial derivatives at
C     the three vertices of a triangle, this routine determines
C     a function W which agrees with the given data, returning
C     the value and (optionally) first partial derivatives of W
C     at a point (X,Y) in the triangle.  The interpolation
C     method is exact for quadratic polynomial data.  The
C     triangle is partitioned into three subtriangles with
C     equal areas.  W is cubic in each subtriangle and along
C     the edges, but has only one continuous derivative across
C     edges.  The normal derivative of W varies linearly along
C     each outer edge.  The values and partial derivatives of W
C     along a triangle edge depend only on the data values at
C     the endpoints of the edge.  Thus the method yields C-1
C     continuity when used to interpolate over a triangular
C     grid.  This algorithm is due to C. L. Lawson.
C
C     Input Parameters -   X,Y - coordinates of a point at which
C                            W is to be evaluated.
C
C        X1,X2,X3,Y1,Y2,Y3 - coordinates of the vertices of
C                            a triangle containing (X,Y).
C
C                 Z1,Z2,Z3 - function values at the vertices
C                            to be interpolated.
C
C              ZX1,ZX2,ZX3 - X-derivative values at the
C                            vertices.
C
C              ZY1,ZY2,ZY3 - Y-derivative values at the
C                            vertices.
C
C                    IFLAG - option indicator
C                            IFLAG = 0 if only W is to be
C                                      computed.
C                            IFLAG = 1 if W, WX, and WY are
C                                      to be returned.
C
C     Input parameters are not altered by this routine.
C
C     Output Parameters -   W - estimated value of the interp-
C                           olatory function at (X,Y) if
C                           IER = 0.  Otherwise W = 0.
C
C                   WX,WY - partial derivatives of W at
C                           (X,Y) if IER = 0 and IFLAG = 1,
C                           unchanged if IFLAG .ne. 1, zero
C                           if IER .ne. 0 and IFLAG = 1.
C
C                     IER - error indicator
C                           IER = 0 if no errors were
C                                   encountered.
C                           IER = 1 if the vertices of the
C                                   triangle are collinear.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     I =               do-loop index
C     IP1,IP2,IP3 =     permuted indices for computing RO, ROX,
C                       and ROY
C     U(K) =            X-component of the vector representing
C                       the side opposite vertex K
C     V(K) =            Y-component of the vector representing
C                       the side opposite vertex K
C     SL(K) =           square of the length of the side
C                       opposite vertex K
C     AREA =            twice the area of the triangle
C     XP,YP =           X-X1, Y-Y1
C     R(K) =            K-th barycentric coordinate
C     RX(K),RY(K) =     X,Y partial derivatives of R(K)
C     PHI(K)            R(K-1)*R(K+1) -- quadratic
C     PHIX(K),PHIY(K) = X,Y partials of PHI(K)
C     RMIN =            min(R1,R2,R3)
C     C1,C2 =           factors for computing RO
C     RO(K) =           factors for computing G -- cubic
C                       correction terms
C     ROX(K),ROY(K) =   X,Y partials of RO(K)
C     F(K) =            factors for computing G, GX, and GY --
C                       constant
C     G(K) =            factors for computing the cardinal
C                       functions -- cubic
C     GX(K),GY(K) =     X,Y partials of G(K)
C     P(K) =            G(K) + PHI(K)
C     PX(K),PY(K) =     X,Y partials of P(K)
C     Q(K) =            G(K) - PHI(K)
C     QX(K),QY(K) =     X,Y partials of Q(K)
C     A(K) =            cardinal function whose coefficient is
C                       Z(K)
C     AX(K),AY(K) =     X,Y partials of A(K) -- cardinal
C                       functions for WX and WY
C     B(K) =            twice the cardinal function whose
C                       coefficient is ZX(K)
C     BX(K),BY(K) =     X,Y partials of B(K)
C     C(K) =            twice the cardinal function whose
C                       coefficient is ZY(K)
C     CX(K),CY(K) =     X,Y partials of C(K)
C
C     .. Parameters ..
      DOUBLE PRECISION  FIVE, THREE, TWO, ZERO
      PARAMETER         (FIVE=5.0D0,THREE=3.0D0,TWO=2.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  W, WX, WY, X, X1, X2, X3, Y, Y1, Y2, Y3, Z1, Z2,
     *                  Z3, ZX1, ZX2, ZX3, ZY1, ZY2, ZY3
      INTEGER           IER, IFLAG
C     .. Local Scalars ..
      DOUBLE PRECISION  AREA, C1, C2, RMIN, XP, YP
      INTEGER           I, IP1, IP2, IP3
C     .. Local Arrays ..
      DOUBLE PRECISION  A(3), AX(3), AY(3), B(3), BX(3), BY(3), C(3),
     *                  CX(3), CY(3), F(3), G(3), GX(3), GY(3), P(3),
     *                  PHI(3), PHIX(3), PHIY(3), PX(3), PY(3), Q(3),
     *                  QX(3), QY(3), R(3), RO(3), ROX(3), ROY(3),
     *                  RX(3), RY(3), SL(3), U(3), V(3)
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      U(1) = X3 - X2
      U(2) = X1 - X3
      U(3) = X2 - X1
C
      V(1) = Y3 - Y2
      V(2) = Y1 - Y3
      V(3) = Y2 - Y1
C
      DO 20 I = 1, 3
         SL(I) = U(I)*U(I) + V(I)*V(I)
   20 CONTINUE
C
C     Area = 3-1 x 3-2
C
      AREA = U(1)*V(2) - U(2)*V(1)
      IF (AREA.EQ.ZERO) THEN
C
C        Vertices are collinear
C
         IER = 1
         W = ZERO
         IF (IFLAG.EQ.1) THEN
            WX = ZERO
            WY = ZERO
         END IF
      ELSE
C
C        R(1) = (2-3 X 2-(X,Y))/AREA, R(2) = (1-(X,Y) X 1-3)/AREA,
C        R(3) = (1-2 X 1-(X,Y))/AREA
C
         R(1) = (U(1)*(Y-Y2)-V(1)*(X-X2))/AREA
         XP = X - X1
         YP = Y - Y1
         R(2) = (U(2)*YP-V(2)*XP)/AREA
         R(3) = (U(3)*YP-V(3)*XP)/AREA
         IER = 0
C
         PHI(1) = R(2)*R(3)
         PHI(2) = R(3)*R(1)
         PHI(3) = R(1)*R(2)
C
         RMIN = MIN(R(1),R(2),R(3))
         IF (RMIN.EQ.R(1)) THEN
            IP1 = 1
            IP2 = 2
            IP3 = 3
         ELSE IF (RMIN.NE.R(2)) THEN
            IP1 = 3
            IP2 = 1
            IP3 = 2
         ELSE
            IP1 = 2
            IP2 = 3
            IP3 = 1
         END IF
C
         C1 = RMIN*RMIN/TWO
         C2 = RMIN/THREE
         RO(IP1) = (PHI(IP1)+FIVE*C1/THREE)*R(IP1) - C1
         RO(IP2) = C1*(R(IP3)-C2)
         RO(IP3) = C1*(R(IP2)-C2)
C
         F(1) = THREE*(SL(2)-SL(3))/SL(1)
         F(2) = THREE*(SL(3)-SL(1))/SL(2)
         F(3) = THREE*(SL(1)-SL(2))/SL(3)
C
         G(1) = (R(2)-R(3))*PHI(1) + F(1)*RO(1) - RO(2) + RO(3)
         G(2) = (R(3)-R(1))*PHI(2) + F(2)*RO(2) - RO(3) + RO(1)
         G(3) = (R(1)-R(2))*PHI(3) + F(3)*RO(3) - RO(1) + RO(2)
C
         DO 40 I = 1, 3
            P(I) = G(I) + PHI(I)
            Q(I) = G(I) - PHI(I)
   40    CONTINUE
C
         A(1) = R(1) + G(3) - G(2)
         A(2) = R(2) + G(1) - G(3)
         A(3) = R(3) + G(2) - G(1)
C
         B(1) = U(3)*P(3) + U(2)*Q(2)
         B(2) = U(1)*P(1) + U(3)*Q(3)
         B(3) = U(2)*P(2) + U(1)*Q(1)
C
         C(1) = V(3)*P(3) + V(2)*Q(2)
         C(2) = V(1)*P(1) + V(3)*Q(3)
         C(3) = V(2)*P(2) + V(1)*Q(1)
C
C        W is a linear combination of the cardinal functions
C
         W = A(1)*Z1 + A(2)*Z2 + A(3)*Z3 + (B(1)*ZX1+B(2)*ZX2+B(3)
     *       *ZX3+C(1)*ZY1+C(2)*ZY2+C(3)*ZY3)/TWO
         IF (IFLAG.EQ.1) THEN
C
C           Compute WX and WY
C
            DO 60 I = 1, 3
               RX(I) = -V(I)/AREA
               RY(I) = U(I)/AREA
   60       CONTINUE
            PHIX(1) = R(2)*RX(3) + RX(2)*R(3)
            PHIY(1) = R(2)*RY(3) + RY(2)*R(3)
            PHIX(2) = R(3)*RX(1) + RX(3)*R(1)
            PHIY(2) = R(3)*RY(1) + RY(3)*R(1)
            PHIX(3) = R(1)*RX(2) + RX(1)*R(2)
            PHIY(3) = R(1)*RY(2) + RY(1)*R(2)
C
            ROX(IP1) = RX(IP1)*(PHI(IP1)+FIVE*C1) + R(IP1)*(PHIX(IP1)
     *                 -RX(IP1))
            ROY(IP1) = RY(IP1)*(PHI(IP1)+FIVE*C1) + R(IP1)*(PHIY(IP1)
     *                 -RY(IP1))
            ROX(IP2) = RX(IP1)*(PHI(IP2)-C1) + C1*RX(IP3)
            ROY(IP2) = RY(IP1)*(PHI(IP2)-C1) + C1*RY(IP3)
            ROX(IP3) = RX(IP1)*(PHI(IP3)-C1) + C1*RX(IP2)
            ROY(IP3) = RY(IP1)*(PHI(IP3)-C1) + C1*RY(IP2)
C
            GX(1) = (RX(2)-RX(3))*PHI(1) + (R(2)-R(3))*PHIX(1) + F(1)
     *              *ROX(1) - ROX(2) + ROX(3)
            GY(1) = (RY(2)-RY(3))*PHI(1) + (R(2)-R(3))*PHIY(1) + F(1)
     *              *ROY(1) - ROY(2) + ROY(3)
            GX(2) = (RX(3)-RX(1))*PHI(2) + (R(3)-R(1))*PHIX(2) + F(2)
     *              *ROX(2) - ROX(3) + ROX(1)
            GY(2) = (RY(3)-RY(1))*PHI(2) + (R(3)-R(1))*PHIY(2) + F(2)
     *              *ROY(2) - ROY(3) + ROY(1)
            GX(3) = (RX(1)-RX(2))*PHI(3) + (R(1)-R(2))*PHIX(3) + F(3)
     *              *ROX(3) - ROX(1) + ROX(2)
            GY(3) = (RY(1)-RY(2))*PHI(3) + (R(1)-R(2))*PHIY(3) + F(3)
     *              *ROY(3) - ROY(1) + ROY(2)
C
            DO 80 I = 1, 3
               PX(I) = GX(I) + PHIX(I)
               PY(I) = GY(I) + PHIY(I)
               QX(I) = GX(I) - PHIX(I)
               QY(I) = GY(I) - PHIY(I)
   80       CONTINUE
C
            AX(1) = RX(1) + GX(3) - GX(2)
            AY(1) = RY(1) + GY(3) - GY(2)
            AX(2) = RX(2) + GX(1) - GX(3)
            AY(2) = RY(2) + GY(1) - GY(3)
            AX(3) = RX(3) + GX(2) - GX(1)
            AY(3) = RY(3) + GY(2) - GY(1)
C
            BX(1) = U(3)*PX(3) + U(2)*QX(2)
            BY(1) = U(3)*PY(3) + U(2)*QY(2)
            BX(2) = U(1)*PX(1) + U(3)*QX(3)
            BY(2) = U(1)*PY(1) + U(3)*QY(3)
            BX(3) = U(2)*PX(2) + U(1)*QX(1)
            BY(3) = U(2)*PY(2) + U(1)*QY(1)
C
            CX(1) = V(3)*PX(3) + V(2)*QX(2)
            CY(1) = V(3)*PY(3) + V(2)*QY(2)
            CX(2) = V(1)*PX(1) + V(3)*QX(3)
            CY(2) = V(1)*PY(1) + V(3)*QY(3)
            CX(3) = V(2)*PX(2) + V(1)*QX(1)
            CY(3) = V(2)*PY(2) + V(1)*QY(1)
C
C           WX and WY are linear combinations of the cardinal
C           functions
C
            WX = AX(1)*Z1 + AX(2)*Z2 + AX(3)*Z3 + (BX(1)*ZX1+BX(2)
     *           *ZX2+BX(3)*ZX3+CX(1)*ZY1+CX(2)*ZY2+CX(3)*ZY3)/TWO
            WY = AY(1)*Z1 + AY(2)*Z2 + AY(3)*Z3 + (BY(1)*ZX1+BY(2)
     *           *ZX2+BY(3)*ZX3+CY(1)*ZY1+CY(2)*ZY2+CY(3)*ZY3)/TWO
         END IF
      END IF
      RETURN
      END


c NIVEAU 3

      
      INTEGER FUNCTION P01ACF(IFAIL,IERROR,SRNAME,VARBNM,NREC,REC)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     P01ACF is the error-handling routine for the F06 AND F07
C     Chapters of the NAG Fortran Library. It is a slightly modified
C     version of P01ABF.
C
C     P01ACF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ACF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ACF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME, VARBNM
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR, VARLEN
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, LEN, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
         VARLEN = 0
         DO 20 I = LEN(VARBNM), 1, -1
            IF (VARBNM(I:I).NE.' ') THEN
               VARLEN = I
               GO TO 40
            END IF
   20    CONTINUE
   40    CONTINUE
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 60 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   60       CONTINUE
            IF (IFAIL.NE.-13) THEN
               IF (VARLEN.NE.0) THEN
                  WRITE (MESS,FMT=99999) SRNAME, VARBNM(1:VARLEN),
     *              IERROR
               ELSE
                  WRITE (MESS,FMT=99998) SRNAME
               END IF
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ACF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': ',A,
     *       ' =',I6)
99998 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A)
      END
      SUBROUTINE C06FCF(X,Y,PTS,WORK,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COMPLEX FOURIER TRANSFORM
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(PTS), X(PTS), Y(PTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  SQPTS
      INTEGER           IERROR, IPTS, PMAX, TWOGRP
C     .. Local Arrays ..
      INTEGER           RFACT(21), TFACT(21)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06ECW, C06FAY, C06FAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Data statements ..
      DATA              PMAX/19/
      DATA              TWOGRP/8/
C     .. Executable Statements ..
      IF (PTS.LE.1) GO TO 40
      CALL C06FAZ(PTS,PMAX,TWOGRP,TFACT,RFACT,IERROR)
      IF (IERROR.NE.0) GO TO 60
      CALL C06ECW(X,Y,PTS,TFACT)
      CALL C06FAY(X,PTS,RFACT,WORK)
      CALL C06FAY(Y,PTS,RFACT,WORK)
      SQPTS = SQRT(DBLE(PTS))
      DO 20 IPTS = 1, PTS
         X(IPTS) = X(IPTS)/SQPTS
         Y(IPTS) = Y(IPTS)/SQPTS
   20 CONTINUE
      IFAIL = 0
      GO TO 80
C
   40 IERROR = 3
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
   80 RETURN
      END
      INTEGER FUNCTION E01SAR(NVERTX,NABOR,IADJ,IEND)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: INDEX.
C     This function returns the index of NABOR in the
C     adjacency list for NVERTX.
C
C     Input Parameters - NVERTX - node whose adjacency list is
C                             to be searched.
C
C                     NABOR - node whose index is to be
C                             returned.  NABOR must be
C                             connected to NVERTX.
C
C                      IADJ - set of adjacency lists.
C
C                      IEND - pointers to the ends of
C                             adjacency lists in IADJ.
C
C     Input parameters are not altered by this function.
C
C     Output Parameter -  E01SAR - IADJ(INDEX) = NABOR.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     NB =   local copy of NABOR
C     INDX = index for IADJ
C
C     .. Scalar Arguments ..
      INTEGER                 NABOR, NVERTX
C     .. Array Arguments ..
      INTEGER                 IADJ(*), IEND(*)
C     .. Local Scalars ..
      INTEGER                 INDX, NB
C     .. Executable Statements ..
      NB = NABOR
C
C     Initialization
C
      INDX = IEND(NVERTX) + 1
   20 CONTINUE
C
C     Search the list of NVERTX neighbors for NB
C
      INDX = INDX - 1
      IF (IADJ(INDX).NE.NB) GO TO 20
C
      E01SAR = INDX
      RETURN
      END
      LOGICAL FUNCTION E01SAT(IN1,IN2,IO1,IO2,X,Y)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: SWPTST.
C     This function decides whether or not to replace a
C     diagonal arc in a quadrilateral with the other diagonal.
C     The determination is based on the sizes of the angles
C     contained in the 2 triangles defined by the diagonal.
C     The diagonal is chosen to maximize the smallest of the
C     six angles over the two pairs of triangles.
C
C     Input Parameters -  IN1,IN2,IO1,IO2 - node indices of the
C                              four points defining the
C                              quadrilateral.  IO1 and IO2
C                              are currently connected by a
C                              diagonal arc.  This arc
C                              should be replaced by an arc
C                              connecting IN1, IN2 if the
C                              decision is made to swap.
C                              IN1,IO1,IO2 must be in
C                              counterclockwise order.
C
C                        X,Y - vectors of nodal coordinates.
C                              (X(I),Y(I)) are the coord-
C                              inates of node I for I = IN1,
C                              IN2, IO1, or IO2.
C
C     None of the input parameters are altered by this routine.
C
C     Output Parameter -  E01SAT - .TRUE. iff the arc connecting
C                              IO1 and IO2 is to be replaced
C
C     ***********************************************************
C
C     Local Parameters -
C
C     DX11,DY11 = X,Y coordinates of the vector IN1-IO1
C     DX12,DY12 = X,Y coordinates of the vector IN1-IO2
C     DX22,DY22 = X,Y coordinates of the vector IN2-IO2
C     DX21,DY21 = X,Y coordinates of the vector IN2-IO1
C     SIN1 =      cross product of the vectors IN1-IO1 and
C               IN1-IO2 -- proportional to sin(T1) where T1
C               is the angle at IN1 formed by the vectors
C     COS1 =      inner product of the vectors IN1-IO1 and
C               IN1-IO2 -- proportional to cos(T1)
C     SIN2 =      cross product of the vectors IN2-IO2 and
C               IN2-IO1 -- proportional to sin(T2) where T2
C               is the angle at IN2 formed by the vectors
C     COS2 =      inner product of the vectors IN2-IO2 and
C               IN2-IO1 -- proportional to cos(T2)
C     SIN12 =     SIN1*COS2 + COS1*SIN2 -- proportional to
C               sin(T1+T2)
C
C     .. Parameters ..
      DOUBLE PRECISION        ZERO
      PARAMETER               (ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER                 IN1, IN2, IO1, IO2
C     .. Array Arguments ..
      DOUBLE PRECISION        X(*), Y(*)
C     .. Local Scalars ..
      DOUBLE PRECISION        COS1, COS2, DX11, DX12, DX21, DX22, DY11,
     *                        DY12, DY21, DY22, SIN1, SIN12, SIN2
C     .. Executable Statements ..
      E01SAT = .FALSE.
C
C     Compute the vectors containing the angles T1, T2
C
      DX11 = X(IO1) - X(IN1)
      DX12 = X(IO2) - X(IN1)
      DX22 = X(IO2) - X(IN2)
      DX21 = X(IO1) - X(IN2)
C
      DY11 = Y(IO1) - Y(IN1)
      DY12 = Y(IO2) - Y(IN1)
      DY22 = Y(IO2) - Y(IN2)
      DY21 = Y(IO1) - Y(IN2)
C
C     Compute inner products
C
      COS1 = DX11*DX12 + DY11*DY12
      COS2 = DX22*DX21 + DY22*DY21
C
C     The diagonals should be swapped iff (T1+T2) .gt. 180
C     degrees.  The following two tests insure numerical
C     stability.
C
      IF (COS1.LT.ZERO .OR. COS2.LT.ZERO) THEN
         IF (COS1.GE.ZERO .OR. COS2.GE.ZERO) THEN
C
C           Compute vector cross products
C
            SIN1 = DX11*DY12 - DX12*DY11
            SIN2 = DX22*DY21 - DX21*DY22
            SIN12 = SIN1*COS2 + COS1*SIN2
            IF (SIN12.GE.ZERO) RETURN
         END IF
         E01SAT = .TRUE.
      END IF
      RETURN
      END
      SUBROUTINE E01SAS(NIN1,NIN2,NOUT1,NOUT2,IADJ,IEND)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: SWAP.
C     This subroutine swaps the diagonals in a convex quadri-
C     lateral.
C
C     Input Parameters -  NIN1,NIN2,NOUT1,NOUT2 - nodal indices
C                            of a pair of adjacent triangles
C                            which form a convex quadrilat-
C                            eral.  NOUT1 and NOUT2 are con-
C                            nected by an arc which is to be
C                            replaced by the arc NIN1-NIN2.
C                            (NIN1,NOUT1,NOUT2) must be tri-
C                            angle vertices in counterclock-
C                            wise order.
C
C     The above parameters are not altered by this routine.
C
C                iadj,iend - triangulation data structure
C                            (see subroutine E01SAY).
C
C     output parameters - IADJ,IEND - updated with the arc
C                                 replacement.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     IN =        NIN1 and NIN2 ordered by increasing magnitude
C               (the neighbors of IN(1) precede those of
C               IN(2) in IADJ)
C     IO =        NOUT1 and NOUT2 in increasing order
C     IP1,IP2 =   permutation of (1,2) such that IO(IP1)
C               precedes IO(IP2) as a neighbor of IN(1)
C     J,K =       permutation of (1,2) used as indices of in
C               and IO
C     NF,NL =     IADJ indices boundary a portion of the array
C               to be shifted
C     I =         IEND index
C     IMIN,IMAX = bounds on the portion of IEND to be incre-
C               mented or decremented
C
C     .. Scalar Arguments ..
      INTEGER           NIN1, NIN2, NOUT1, NOUT2
C     .. Array Arguments ..
      INTEGER           IADJ(*), IEND(*)
C     .. Local Scalars ..
      INTEGER           I, IMAX, IMIN, IP1, IP2, J, K, NF, NL
C     .. Local Arrays ..
      INTEGER           IN(2), IO(2)
C     .. External Functions ..
      INTEGER           E01SAR
      EXTERNAL          E01SAR
C     .. External Subroutines ..
      EXTERNAL          E01SAQ
C     .. Executable Statements ..
      IN(1) = NIN1
      IN(2) = NIN2
      IO(1) = NOUT1
      IO(2) = NOUT2
      IP1 = 1
C
C     Order the indices so that IN(1) .lt. IN(2) and IO(1) .lt.
C     IO(2), and choose IP1 and IP2 such that (IN(1),IO(IP1),
C     IO(IP2)) forms a triangle.
C
      IF (IN(1).GE.IN(2)) THEN
         IN(1) = IN(2)
         IN(2) = NIN1
         IP1 = 2
      END IF
      IF (IO(1).GE.IO(2)) THEN
         IO(1) = IO(2)
         IO(2) = NOUT1
         IP1 = 3 - IP1
      END IF
      IP2 = 3 - IP1
      IF (IO(2).LT.IN(1)) THEN
C
C        The vertices are ordered (IO(1),IO(2),IN(1),IN(2)).
C        Delete IO(2) by shifting up by 1
C
         NF = 1 + E01SAR(IO(1),IO(2),IADJ,IEND)
         NL = -1 + E01SAR(IO(2),IO(1),IADJ,IEND)
         IF (NF.LE.NL) CALL E01SAQ(NF,NL,-1,IADJ)
         IMIN = IO(1)
         IMAX = IO(2) - 1
         DO 20 I = IMIN, IMAX
            IEND(I) = IEND(I) - 1
   20    CONTINUE
C
C        Delete IO(1) by shifting up by 2 and insert IN(2)
C
         NF = NL + 2
         NL = -1 + E01SAR(IN(1),IO(IP2),IADJ,IEND)
         IF (NF.LE.NL) CALL E01SAQ(NF,NL,-2,IADJ)
         IADJ(NL-1) = IN(2)
         IMIN = IO(2)
         IMAX = IN(1) - 1
         DO 40 I = IMIN, IMAX
            IEND(I) = IEND(I) - 2
   40    CONTINUE
C
C        Shift up by 1 and insert IN(1)
C
         NF = NL + 1
         NL = -1 + E01SAR(IN(2),IO(IP1),IADJ,IEND)
         CALL E01SAQ(NF,NL,-1,IADJ)
         IADJ(NL) = IN(1)
         IMIN = IN(1)
         IMAX = IN(2) - 1
         DO 60 I = IMIN, IMAX
            IEND(I) = IEND(I) - 1
   60    CONTINUE
      ELSE IF (IN(2).LT.IO(1)) THEN
C
C        The vertices are ordered (IN(1),IN(2),IO(1),IO(2)).
C        Delete IO(1) by shifting down by 1
C
         NF = 1 + E01SAR(IO(1),IO(2),IADJ,IEND)
         NL = -1 + E01SAR(IO(2),IO(1),IADJ,IEND)
         IF (NF.LE.NL) CALL E01SAQ(NF,NL,1,IADJ)
         IMIN = IO(1)
         IMAX = IO(2) - 1
         DO 80 I = IMIN, IMAX
            IEND(I) = IEND(I) + 1
   80    CONTINUE
C
C        Delete IO(2) by shifting down by 2 and insert IN(1)
C
         NL = NF - 2
         NF = 1 + E01SAR(IN(2),IO(IP2),IADJ,IEND)
         IF (NF.LE.NL) CALL E01SAQ(NF,NL,2,IADJ)
         IADJ(NF+1) = IN(1)
         IMIN = IN(2)
         IMAX = IO(1) - 1
         DO 100 I = IMIN, IMAX
            IEND(I) = IEND(I) + 2
  100    CONTINUE
C
C        Shift down by 1 and insert IN(2)
C
         NL = NF - 1
         NF = 1 + E01SAR(IN(1),IO(IP1),IADJ,IEND)
         CALL E01SAQ(NF,NL,1,IADJ)
         IADJ(NF) = IN(2)
         IMIN = IN(1)
         IMAX = IN(2) - 1
         DO 120 I = IMIN, IMAX
            IEND(I) = IEND(I) + 1
  120    CONTINUE
      ELSE
C
C        IN(1) and IO(1) precede IN(2) and IO(2).  For (J,K) =
C        (1,2) and (2,1), delete IO(K) as a neighbor of IO(J)
C        by shifting a portion of IADJ either up or down and
C        and insert IN(K) as a neighbor of IN(J).
C
         DO 180 J = 1, 2
            K = 3 - J
            IF (IN(J).GT.IO(J)) THEN
C
C              The neighbors of IO(J) precede those of IN(J) -- shift
C              up by 1
C
               NF = 1 + E01SAR(IO(J),IO(K),IADJ,IEND)
               NL = -1 + E01SAR(IN(J),IO(IP2),IADJ,IEND)
               IF (NF.LE.NL) CALL E01SAQ(NF,NL,-1,IADJ)
               IADJ(NL) = IN(K)
               IMIN = IO(J)
               IMAX = IN(J) - 1
               DO 140 I = IMIN, IMAX
                  IEND(I) = IEND(I) - 1
  140          CONTINUE
            ELSE
C
C              The neighbors of IN(J) precede those of IO(J) -- shift
C              down by 1
C
               NF = 1 + E01SAR(IN(J),IO(IP1),IADJ,IEND)
               NL = -1 + E01SAR(IO(J),IO(K),IADJ,IEND)
               IF (NF.LE.NL) CALL E01SAQ(NF,NL,1,IADJ)
               IADJ(NF) = IN(K)
               IMIN = IN(J)
               IMAX = IO(J) - 1
               DO 160 I = IMIN, IMAX
                  IEND(I) = IEND(I) + 1
  160          CONTINUE
            END IF
C
C           Reverse (IP1,IP2) for (J,K) = (2,1)
C
            IP1 = IP2
            IP2 = 3 - IP1
  180    CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE E01SAU(KK,I1,I2,IADJ,IEND)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: BDYADD.
C     This routine adds a boundary node to a triangulation
C     of a set of KK-1 points in the plane.  IADJ and IEND are
C     updated with the insertion of node KK.
C
C     Input Parameters -   KK - index of an exterior node to be
C                           added.  KK .ge. 4.
C
C                      I1 - first (rightmost as viewed from
C                           KK) boundary node in the mesh
C                           which is visible from KK - the
C                           line segment KK-I1 intersects
C                           no arcs.
C
C                      I2 - last (leftmost) boundary node
C                           which is visible from KK.
C
C                    IADJ - set of adjacency lists of nodes
C                           in the mesh.
C
C                    IEND - pointers to the ends of
C                           adjacency lists in IADJ for
C                           each node in the mesh.
C
C     IADJ and IEND may be created by E01SAY and must contain
C     the vertices I1 and I2.  I1 and I2 may be determined by
C     E01SAW.
C
C     KK, I1, and I2 are not altered by this routine.
C
C     Output Parameters - IADJ,IEND - updated with the addition
C                                 of node KK as the last
C                                 entry.  Node KK will be
C                                 connected to I1, I2, and
C                                 all boundary nodes between
C                                 them.  No optimization of
C                                 the mesh is performed.
C
C     Module referenced by E01SAU - E01SAQ
C
C     ***********************************************************
C
C     Local Parameters -
C
C     K =            local copy of KK
C     KM1 =          K - 1
C     NRIGHT,NLEFT = local copies of I1, I2
C     NF,NL =        indices of IADJ bounding the portion of the
C                  array to be shifted
C     N1 =           IADJ index of the first neighbor of NLEFT
C     N2 =           IADJ index of the last neighbor of NRIGHT
C     I =            do-loop index
C     IMIN,IMAX =    bounds on do-loop index -- first and last
C                  elements of IEND to be incremented
C     KEND =         pointer to the last neighbor of K in IADJ
C     NEXT =         next boundary node to be connected to KK
C     INDX =         index for IADJ
C
C     .. Scalar Arguments ..
      INTEGER           I1, I2, KK
C     .. Array Arguments ..
      INTEGER           IADJ(*), IEND(KK)
C     .. Local Scalars ..
      INTEGER           I, IMAX, IMIN, INDX, K, KEND, KM1, N1, N2, NEXT,
     *                  NF, NL, NLEFT, NRIGHT
C     .. External Subroutines ..
      EXTERNAL          E01SAQ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
      K = KK
      KM1 = K - 1
      NRIGHT = I1
      NLEFT = I2
C
C     Initialize variables
C
      NL = IEND(KM1)
      N1 = 1
      IF (NLEFT.NE.1) N1 = IEND(NLEFT-1) + 1
      N2 = IEND(NRIGHT)
      NF = MAX(N1,N2)
C
C     Insert K as a neighbor of max(NRIGHT,NLEFT)
C
      CALL E01SAQ(NF,NL,2,IADJ)
      IADJ(NF+1) = K
      IMIN = MAX(NRIGHT,NLEFT)
      DO 20 I = IMIN, KM1
         IEND(I) = IEND(I) + 2
   20 CONTINUE
C
C     Initialize KEND and insert K as a neighbor of
C     min(NRIGHT,NLEFT)
C
      KEND = NL + 3
      NL = NF - 1
      NF = MIN(N1,N2)
      CALL E01SAQ(NF,NL,1,IADJ)
      IADJ(NF) = K
      IMAX = IMIN - 1
      IMIN = MIN(NRIGHT,NLEFT)
      DO 40 I = IMIN, IMAX
         IEND(I) = IEND(I) + 1
   40 CONTINUE
C
C     INSERT NRIGHT AS THE FIRST NEIGHBOR OF K
C
      IADJ(KEND) = NRIGHT
C
C     Initialize INDX for loop on boundary nodes between nright
C     and NLEFT
C
      INDX = IEND(NRIGHT) - 2
   60 CONTINUE
      NEXT = IADJ(INDX)
      IF (NEXT.NE.NLEFT) THEN
C
C        Connect NEXT and K
C
         KEND = KEND + 1
         IADJ(KEND) = NEXT
         INDX = IEND(NEXT)
         IADJ(INDX) = K
         INDX = INDX - 1
         GO TO 60
      END IF
C
C     Insert NLEFT and 0 as the last neighbors of K
C
      IADJ(KEND+1) = NLEFT
      KEND = KEND + 2
      IADJ(KEND) = 0
      IEND(K) = KEND
      RETURN
      END
      SUBROUTINE E01SAV(KK,I1,I2,I3,IADJ,IEND)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: INTADD.
C     This routine adds an interior node to a triangulation
C     of a set of KK-1 points in the plane.  IADJ and IEND are
C     updated with the insertion of node KK in the triangle
C     whose vertices are I1, I2, and I3.
C
C     Input Parameters -        KK - index of node to be
C                                inserted.  KK .ge. 4.
C
C                     I1,I2,I3 - indices of the vertices of
C                                a triangle containing node
C                                KK -- in counterclockwise
C                                order.
C
C                         IADJ - set of adjacency lists
C                                of nodes in the mesh.
C
C                         IEND - pointers to the ends of
C                                adjacency lists in IADJ for
C                                each node in the mesh.
C
C     IADJ and IEND may be created by E01SAY and must contain
C     the vertices I1, I2, and I3.  I1,I2,I3 may be determined
C     by E01SAW.
C
C     KK, I1, I2, and I3 are not altered by this routine.
C
C     Output Parameters - IADJ,IEND - updated with the addition
C                                 of node KK as the last
C                                 entry.  Node KK will be
C                                 connected to nodes I1, I2,
C                                 and I3.  No optimization
C                                 of the mesh is performed.
C
C     Module referenced by E01SAV - E01SAQ
C
C     ***********************************************************
C
C     Local Parameters -
C
C     K =           local copy of KK
C     KM1 =         K - 1
C     N =           vector containing I1, I2, I3
C     NFT =         pointers to the tops of the 3 sets of IADJ
C                 elements to be shifted downward
C     IP1,IP2,IP3 = permutation indices for N and NFT
C     INDX =        index for IADJ and N
C     NF,NL =       indices of first and last entries in IADJ
C                 to be shifted down
C     N1,N2 =       first 2 vertices of a new triangle --
C                 (N1,N2,KK)
C     IMIN,IMAX =   bounds on do-loop index -- first and last
C                 elements of IEND to be incremented
C     I =           do-loop index
C     ITEMP =       temporary storage location
C
C     .. Scalar Arguments ..
      INTEGER           I1, I2, I3, KK
C     .. Array Arguments ..
      INTEGER           IADJ(*), IEND(KK)
C     .. Local Scalars ..
      INTEGER           I, IMAX, IMIN, INDX, IP1, IP2, IP3, ITEMP, K,
     *                  KM1, N1, N2, NF, NL
C     .. Local Arrays ..
      INTEGER           N(3), NFT(3)
C     .. External Subroutines ..
      EXTERNAL          E01SAQ
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      K = KK
C
C     Initialization
C
      N(1) = I1
      N(2) = I2
      N(3) = I3
C
C     set up NFT
C
      DO 40 I = 1, 3
         N1 = N(I)
         INDX = MOD(I,3) + 1
         N2 = N(INDX)
         INDX = IEND(N1) + 1
   20    CONTINUE
C
C        Find the index of N2 as a neighbor of N1
C
         INDX = INDX - 1
         IF (IADJ(INDX).NE.N2) GO TO 20
         NFT(I) = INDX + 1
   40 CONTINUE
C
C     Order the vertices by decreasing magnitude.
C     N(IP(I+1)) precedes N(IP(I)) in IEND for
C     I = 1,2.
C
      IP1 = 1
      IP2 = 2
      IP3 = 3
      IF (N(2).GT.N(1)) THEN
         IP1 = 2
         IP2 = 1
      END IF
      IF (N(3).GT.N(IP1)) THEN
         IP3 = IP1
         IP1 = 3
      END IF
      IF (N(IP3).GT.N(IP2)) THEN
         ITEMP = IP2
         IP2 = IP3
         IP3 = ITEMP
      END IF
C
C     Add node K to the adjacency lists of each vertex and
C     update IEND.  For each vertex, a set of IADJ elements
C     is shifted downward and K is inserted.  Shifting starts
C     at the end of the array.
C
      KM1 = K - 1
      NL = IEND(KM1)
      NF = NFT(IP1)
      IF (NF.LE.NL) CALL E01SAQ(NF,NL,3,IADJ)
      IADJ(NF+2) = K
      IMIN = N(IP1)
      IMAX = KM1
      DO 60 I = IMIN, IMAX
         IEND(I) = IEND(I) + 3
   60 CONTINUE
C
      NL = NF - 1
      NF = NFT(IP2)
      CALL E01SAQ(NF,NL,2,IADJ)
      IADJ(NF+1) = K
      IMAX = IMIN - 1
      IMIN = N(IP2)
      DO 80 I = IMIN, IMAX
         IEND(I) = IEND(I) + 2
   80 CONTINUE
C
      NL = NF - 1
      NF = NFT(IP3)
      CALL E01SAQ(NF,NL,1,IADJ)
      IADJ(NF) = K
      IMAX = IMIN - 1
      IMIN = N(IP3)
      DO 100 I = IMIN, IMAX
         IEND(I) = IEND(I) + 1
  100 CONTINUE
C
C     Add node K to IEND and its neighbors to IADJ
C
      INDX = IEND(KM1)
      IEND(K) = INDX + 3
      DO 120 I = 1, 3
         INDX = INDX + 1
         IADJ(INDX) = N(I)
  120 CONTINUE
      RETURN
      END


c niveau 4
      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      STOP
      END
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/6/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END
      SUBROUTINE C06ECW(X,Y,PTS,FACTOR)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COMPLEX FOURIER TRANSFORM KERNEL DRIVER
C     .. Scalar Arguments ..
      INTEGER           PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS), Y(PTS)
      INTEGER           FACTOR(21)
C     .. Local Scalars ..
      INTEGER           F, M, M1, M2, M3, M4, M5, M6, M7, P
C     .. External Subroutines ..
      EXTERNAL          C06ECQ, C06ECR, C06ECS, C06ECT, C06ECU, C06ECV
C     .. Executable Statements ..
      F = 0
      M = PTS
   20 CONTINUE
      F = F + 1
      P = FACTOR(F)
      IF (P.EQ.0) RETURN
      IF (P.EQ.1) GO TO 20
      M = M/P
      M1 = PTS - M
      M2 = M1 - M
      M3 = M2 - M
      M4 = M3 - M
      M5 = M4 - M
      M6 = M5 - M
      M7 = M6 - M
      IF (P.EQ.2) GO TO 40
      IF (P.EQ.3) GO TO 60
      IF (P.EQ.4) GO TO 80
      IF (P.EQ.5) GO TO 100
      IF (P.EQ.8) GO TO 120
      GO TO 140
C
   40 CONTINUE
      CALL C06ECV(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,M)
      GO TO 20
C
   60 CONTINUE
      CALL C06ECU(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1),M2,M)
      GO TO 20
C
   80 CONTINUE
      CALL C06ECT(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1)
     *            ,M2,X(3*M+1),Y(3*M+1),M3,M)
      GO TO 20
C
  100 CONTINUE
      CALL C06ECS(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1)
     *            ,M2,X(3*M+1),Y(3*M+1),M3,X(4*M+1),Y(4*M+1),M4,M)
      GO TO 20
C
  120 CONTINUE
      CALL C06ECR(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1)
     *            ,M2,X(3*M+1),Y(3*M+1),M3,X(4*M+1),Y(4*M+1),M4,X(5*M+1)
     *            ,Y(5*M+1),M5,X(6*M+1),Y(6*M+1),M6,X(7*M+1),Y(7*M+1)
     *            ,M7,M)
      GO TO 20
C
  140 CONTINUE
      CALL C06ECQ(X,Y,PTS,M,P)
      GO TO 20
C
      END
      SUBROUTINE C06FAY(X,PTS,FACTOR,WORK)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11 REVISED. IER-444 (FEB 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     SINGLE COPY REORDERING PROGRAMME
C     EQUIVALENCE (I1,I(1)), (K1,K(1)), (L1,L(1))
C     EQUIVALENCE (I2,I(2)), (K2,K(2)), (L2,L(2))
C     EQUIVALENCE (I3,I(3)), (K3,K(3)), (L3,L(3))
C     EQUIVALENCE (I4,I(4)), (K4,K(4)), (L4,L(4))
C     EQUIVALENCE (I5,I(5)), (K5,K(5)), (L5,L(5))
C     EQUIVALENCE (I6,I(6)), (K6,K(6)), (L6,L(6))
C     EQUIVALENCE (I7,I(7)), (K7,K(7)), (L7,L(7))
C     EQUIVALENCE (I8,I(8)), (K8,K(8)), (L8,L(8))
C     EQUIVALENCE (I9,I(9)), (K9,K(9)), (L9,L(9))
C     EQUIVALENCE (I10,I(10)), (K10,K(10)), (L10,L(10))
C     EQUIVALENCE (K11,K(11))
C     .. Scalar Arguments ..
      INTEGER           PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(PTS), X(PTS)
      INTEGER           FACTOR(21)
C     .. Local Scalars ..
      INTEGER           I1, I10, I2, I3, I4, I5, I6, I7, I8, I9, J, JJ,
     *                  K1, K10, K11, K2, K3, K4, K5, K6, K7, K8, K9,
     *                  KK, L1, L10, L2, L3, L4, L5, L6, L7, L8, L9,
     *                  LEVEL, LOOP, NEST
C     .. Local Arrays ..
      INTEGER           I(20), K(20), L(20)
C     .. Data statements ..
      DATA              NEST/20/
      DATA              LOOP/10/
C     .. Executable Statements ..
      IF (FACTOR(1).EQ.0) GO TO 380
      DO 20 J = 1, NEST
         L(J) = 1
         I(J) = 1
   20 CONTINUE
      KK = PTS
      DO 40 J = 1, NEST
         IF (FACTOR(J).EQ.0) GO TO 60
         L(J) = KK
         I(J) = KK/FACTOR(J)
         KK = KK/FACTOR(J)
   40 CONTINUE
   60 CONTINUE
C
      L1 = L(1)
      L2 = L(2)
      L3 = L(3)
      L4 = L(4)
      L5 = L(5)
      L6 = L(6)
      L7 = L(7)
      L8 = L(8)
      L9 = L(9)
      L10 = L(10)
      I1 = I(1)
      I2 = I(2)
      I3 = I(3)
      I4 = I(4)
      I5 = I(5)
      I6 = I(6)
      I7 = I(7)
      I8 = I(8)
      I9 = I(9)
      I10 = I(10)
C
      KK = 0
      LEVEL = NEST
      K(LEVEL) = 1
      GO TO 100
   80 CONTINUE
      IF (LEVEL.GE.NEST) GO TO 340
      LEVEL = LEVEL + 1
      K(LEVEL) = K(LEVEL) + I(LEVEL)
      IF (K(LEVEL).GT.L(LEVEL)) GO TO 80
  100 CONTINUE
      LEVEL = LEVEL - 1
      DO 120 J = LOOP, LEVEL
         JJ = LEVEL + LOOP - J
         K(JJ) = K(JJ+1)
  120 CONTINUE
      K11 = K(11)
      DO 320 K10 = K11, L10, I10
         DO 300 K9 = K10, L9, I9
            DO 280 K8 = K9, L8, I8
               DO 260 K7 = K8, L7, I7
                  DO 240 K6 = K7, L6, I6
                     DO 220 K5 = K6, L5, I5
                        DO 200 K4 = K5, L4, I4
                           DO 180 K3 = K4, L3, I3
                              DO 160 K2 = K3, L2, I2
                                 DO 140 K1 = K2, L1, I1
                                    KK = KK + 1
                                    WORK(KK) = X(K1)
  140                            CONTINUE
  160                         CONTINUE
  180                      CONTINUE
  200                   CONTINUE
  220                CONTINUE
  240             CONTINUE
  260          CONTINUE
  280       CONTINUE
  300    CONTINUE
  320 CONTINUE
      LEVEL = LOOP
      GO TO 80
  340 CONTINUE
      DO 360 J = 1, PTS
         X(J) = WORK(J)
  360 CONTINUE
C
  380 CONTINUE
      RETURN
      END
      SUBROUTINE C06FAZ(PTS,PMAX,TWOGRP,TFACT,RFACT,IERROR)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COPY REORDERING FACTORING PROGRAMME
C     .. Scalar Arguments ..
      INTEGER           IERROR, PMAX, PTS, TWOGRP
C     .. Array Arguments ..
      INTEGER           RFACT(21), TFACT(21)
C     .. Local Scalars ..
      INTEGER           F, J, JJ, N, NEST, P, PTWO, Q
C     .. Local Arrays ..
      INTEGER           PP(10), QQ(20)
C     .. Data statements ..
      DATA              NEST/20/
C     .. Executable Statements ..
      N = PTS
      F = 2
      P = 0
      Q = 0
   20 CONTINUE
      IF (N.LE.1) GO TO 100
      DO 40 J = F, PMAX
         IF (N.EQ.(N/J)*J) GO TO 60
   40 CONTINUE
      GO TO 280
   60 CONTINUE
      IF (2*P+Q.GE.NEST) GO TO 300
      F = J
      N = N/F
      IF (N.EQ.(N/F)*F) GO TO 80
      Q = Q + 1
      QQ(Q) = F
      GO TO 20
   80 CONTINUE
      N = N/F
      P = P + 1
      PP(P) = F
      GO TO 20
C
  100 CONTINUE
      IF (P.LT.1) GO TO 140
      DO 120 J = 1, P
         JJ = P + 1 - J
         TFACT(J) = PP(JJ)
         JJ = P + Q + J
         TFACT(JJ) = PP(J)
  120 CONTINUE
  140 CONTINUE
      IF (Q.LT.1) GO TO 180
      DO 160 J = 1, Q
         JJ = P + J
         TFACT(JJ) = QQ(J)
  160 CONTINUE
  180 CONTINUE
      JJ = 2*P + Q
      TFACT(JJ+1) = 0
      RFACT(JJ+1) = 0
      DO 200 J = 1, JJ
         RFACT(J) = TFACT(J)
  200 CONTINUE
      IF (JJ.EQ.1) RFACT(1) = 0
      PTWO = 1
      J = 0
  220 CONTINUE
      J = J + 1
      IF (TFACT(J).EQ.0) GO TO 260
      IF (TFACT(J).NE.2) GO TO 220
      PTWO = PTWO*2
      TFACT(J) = 1
      IF (PTWO.GE.TWOGRP) GO TO 240
      IF (TFACT(J+1).EQ.2) GO TO 220
  240 CONTINUE
      TFACT(J) = PTWO
      PTWO = 1
      GO TO 220
  260 CONTINUE
      IERROR = 0
      GO TO 320
  280 IERROR = 1
      GO TO 320
C
  300 IERROR = 2
  320 RETURN
      END



c niveau 5
      

      SUBROUTINE C06ECQ(X,Y,PTS,M,P)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX PRIME COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, P, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS), Y(PTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, IS, IU, RS, RU, T, TWOPI, XT, YT
      INTEGER           J, JJ, K, K0, KS1, KS2, MOVER2, MP, PM, PP, U, V
      LOGICAL           FOLD, ZERO
C     .. Local Arrays ..
      DOUBLE PRECISION  A(18), AA(9,9), B(18), BB(9,9), C(18), IA(9),
     *                  IB(9), RA(9), RB(9), S(18)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN
C     .. Executable Statements ..
      TWOPI = 2.0D0*X01AAF(0.0D0)
      MOVER2 = M/2 + 1
      MP = M*P
      PP = P/2
      PM = P - 1
      DO 20 U = 1, PP
         JJ = P - U
         ANGLE = TWOPI*DBLE(U)/DBLE(P)
         A(U) = COS(ANGLE)
         B(U) = SIN(ANGLE)
         A(JJ) = A(U)
         B(JJ) = -B(U)
   20 CONTINUE
      DO 60 U = 1, PP
         DO 40 V = 1, PP
            JJ = U*V - ((U*V)/P)*P
            AA(V,U) = A(JJ)
            BB(V,U) = B(JJ)
   40    CONTINUE
   60 CONTINUE
C
      DO 300 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(MP)
         ZERO = ANGLE .EQ. 0.0D0
         C(1) = COS(ANGLE)
         S(1) = SIN(ANGLE)
         DO 80 U = 2, PM
            C(U) = C(U-1)*C(1) - S(U-1)*S(1)
            S(U) = S(U-1)*C(1) + C(U-1)*S(1)
   80    CONTINUE
         GO TO 140
  100    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         DO 120 U = 1, PM
            T = C(U)*A(U) + S(U)*B(U)
            S(U) = -S(U)*A(U) + C(U)*B(U)
            C(U) = T
  120    CONTINUE
  140    CONTINUE
C
         DO 280 K = K0, PTS, MP
            XT = X(K)
            YT = Y(K)
            KS1 = M + K
            KS2 = (P-1)*M + K
            RS = X(KS1) + X(KS2)
            IS = Y(KS1) + Y(KS2)
            RU = X(KS1) - X(KS2)
            IU = Y(KS1) - Y(KS2)
            DO 160 U = 1, PP
               RA(U) = XT + RS*AA(U,1)
               IA(U) = YT + IS*AA(U,1)
               RB(U) = RU*BB(U,1)
               IB(U) = IU*BB(U,1)
  160       CONTINUE
            XT = XT + RS
            YT = YT + IS
            DO 200 U = 2, PP
               JJ = P - U
               KS1 = U*M + K
               KS2 = JJ*M + K
               RS = X(KS1) + X(KS2)
               IS = Y(KS1) + Y(KS2)
               RU = X(KS1) - X(KS2)
               IU = Y(KS1) - Y(KS2)
               XT = XT + RS
               YT = YT + IS
               DO 180 V = 1, PP
                  RA(V) = RA(V) + RS*AA(V,U)
                  IA(V) = IA(V) + IS*AA(V,U)
                  RB(V) = RB(V) + RU*BB(V,U)
                  IB(V) = IB(V) + IU*BB(V,U)
  180          CONTINUE
  200       CONTINUE
            X(K) = XT
            Y(K) = YT
            DO 260 U = 1, PP
               JJ = P - U
               IF (ZERO) GO TO 220
               XT = RA(U) + IB(U)
               YT = IA(U) - RB(U)
               KS1 = U*M + K
               X(KS1) = XT*C(U) + YT*S(U)
               Y(KS1) = YT*C(U) - XT*S(U)
               XT = RA(U) - IB(U)
               YT = IA(U) + RB(U)
               KS1 = JJ*M + K
               X(KS1) = XT*C(JJ) + YT*S(JJ)
               Y(KS1) = YT*C(JJ) - XT*S(JJ)
               GO TO 240
  220          CONTINUE
               KS1 = U*M + K
               X(KS1) = RA(U) + IB(U)
               Y(KS1) = IA(U) - RB(U)
               KS1 = JJ*M + K
               X(KS1) = RA(U) - IB(U)
               Y(KS1) = IA(U) + RB(U)
  240          CONTINUE
  260       CONTINUE
  280    CONTINUE
         IF (FOLD) GO TO 100
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE C06ECR(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,X3,Y3,M3,X4,Y4,M4,
     *                  X5,Y5,M5,X6,Y6,M6,X7,Y7,M7,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX EIGHT COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, M3, M4, M5, M6, M7, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), X3(M3), X4(M4), X5(M5),
     *                  X6(M6), X7(M7), Y0(PTS), Y1(M1), Y2(M2), Y3(M3),
     *                  Y4(M4), Y5(M5), Y6(M6), Y7(M7)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, C1, C2, C3, C4, C5, C6, C7, E, I1, I2,
     *                  I3, I4, I5, I6, I7, IS0, IS1, IS2, IS3, ISS0,
     *                  ISS1, ISU0, ISU1, IU0, IU1, IU2, IU3, IUS0,
     *                  IUS1, IUU0, IUU1, R1, R2, R3, R4, R5, R6, R7,
     *                  RS0, RS1, RS2, RS3, RSS0, RSS1, RSU0, RSU1, RU0,
     *                  RU1, RU2, RU3, RUS0, RUS1, RUU0, RUU1, S1, S2,
     *                  S3, S4, S5, S6, S7, T, TWOPI
      INTEGER           J, K, K0, M8, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN
C     .. Executable Statements ..
      M8 = M*8
      MOVER2 = M/2 + 1
      TWOPI = 2.0D0*X01AAF(0.0D0)
      E = COS(TWOPI/8.0D0)
C
      DO 120 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M8)
         ZERO = ANGLE .EQ. 0.0D0
         C1 = COS(ANGLE)
         S1 = SIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         C3 = C2*C1 - S2*S1
         S3 = S2*C1 + C2*S1
         C4 = C2*C2 - S2*S2
         S4 = S2*C2 + C2*S2
         C5 = C4*C1 - S4*S1
         S5 = S4*C1 + C4*S1
         C6 = C4*C2 - S4*S2
         S6 = S4*C2 + C4*S2
         C7 = C4*C3 - S4*S3
         S7 = S4*C3 + C4*S3
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         T = (C1+S1)*E
         S1 = (C1-S1)*E
         C1 = T
         T = S2
         S2 = C2
         C2 = T
         T = (-C3+S3)*E
         S3 = (C3+S3)*E
         C3 = T
         C4 = -C4
         T = -(C5+S5)*E
         S5 = (-C5+S5)*E
         C5 = T
         T = -S6
         S6 = -C6
         C6 = T
         T = (C7-S7)*E
         S7 = -(C7+S7)*E
         C7 = T
   40    CONTINUE
C
         DO 100 K = K0, PTS, M8
            RS0 = X0(K) + X4(K)
            IS0 = Y0(K) + Y4(K)
            RU0 = X0(K) - X4(K)
            IU0 = Y0(K) - Y4(K)
            RS1 = X1(K) + X5(K)
            IS1 = Y1(K) + Y5(K)
            RU1 = X1(K) - X5(K)
            IU1 = Y1(K) - Y5(K)
            RS2 = X2(K) + X6(K)
            IS2 = Y2(K) + Y6(K)
            RU2 = X2(K) - X6(K)
            IU2 = Y2(K) - Y6(K)
            RS3 = X3(K) + X7(K)
            IS3 = Y3(K) + Y7(K)
            RU3 = X3(K) - X7(K)
            IU3 = Y3(K) - Y7(K)
            RSS0 = RS0 + RS2
            ISS0 = IS0 + IS2
            RSU0 = RS0 - RS2
            ISU0 = IS0 - IS2
            RSS1 = RS1 + RS3
            ISS1 = IS1 + IS3
            RSU1 = RS1 - RS3
            ISU1 = IS1 - IS3
            RUS0 = RU0 - IU2
            IUS0 = IU0 + RU2
            RUU0 = RU0 + IU2
            IUU0 = IU0 - RU2
            RUS1 = RU1 - IU3
            IUS1 = IU1 + RU3
            RUU1 = RU1 + IU3
            IUU1 = IU1 - RU3
            T = (RUS1+IUS1)*E
            IUS1 = (IUS1-RUS1)*E
            RUS1 = T
            T = (RUU1+IUU1)*E
            IUU1 = (IUU1-RUU1)*E
            RUU1 = T
            X0(K) = RSS0 + RSS1
            Y0(K) = ISS0 + ISS1
            IF (ZERO) GO TO 60
            R1 = RUU0 + RUU1
            I1 = IUU0 + IUU1
            R2 = RSU0 + ISU1
            I2 = ISU0 - RSU1
            R3 = RUS0 + IUS1
            I3 = IUS0 - RUS1
            R4 = RSS0 - RSS1
            I4 = ISS0 - ISS1
            R5 = RUU0 - RUU1
            I5 = IUU0 - IUU1
            R6 = RSU0 - ISU1
            I6 = ISU0 + RSU1
            R7 = RUS0 - IUS1
            I7 = IUS0 + RUS1
            X4(K) = R1*C1 + I1*S1
            Y4(K) = I1*C1 - R1*S1
            X2(K) = R2*C2 + I2*S2
            Y2(K) = I2*C2 - R2*S2
            X6(K) = R3*C3 + I3*S3
            Y6(K) = I3*C3 - R3*S3
            X1(K) = R4*C4 + I4*S4
            Y1(K) = I4*C4 - R4*S4
            X5(K) = R5*C5 + I5*S5
            Y5(K) = I5*C5 - R5*S5
            X3(K) = R6*C6 + I6*S6
            Y3(K) = I6*C6 - R6*S6
            X7(K) = R7*C7 + I7*S7
            Y7(K) = I7*C7 - R7*S7
            GO TO 80
   60       CONTINUE
            X4(K) = RUU0 + RUU1
            Y4(K) = IUU0 + IUU1
            X2(K) = RSU0 + ISU1
            Y2(K) = ISU0 - RSU1
            X6(K) = RUS0 + IUS1
            Y6(K) = IUS0 - RUS1
            X1(K) = RSS0 - RSS1
            Y1(K) = ISS0 - ISS1
            X5(K) = RUU0 - RUU1
            Y5(K) = IUU0 - IUU1
            X3(K) = RSU0 - ISU1
            Y3(K) = ISU0 + RSU1
            X7(K) = RUS0 - IUS1
            Y7(K) = IUS0 + RUS1
   80       CONTINUE
  100    CONTINUE
         IF (FOLD) GO TO 20
  120 CONTINUE
C
      RETURN
      END
      SUBROUTINE C06ECS(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,X3,Y3,M3,X4,Y4,M4,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX FIVE COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, M3, M4, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), X3(M3), X4(M4),
     *                  Y0(PTS), Y1(M1), Y2(M2), Y3(M3), Y4(M4)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, ANGLE, AS, AU, B1, B2, C1, C2, C3, C4,
     *                  I0, I1, I2, I3, I4, IA1, IA2, IAS, IAU, IB1,
     *                  IB2, IS1, IS2, ISS, IU1, IU2, R0, R1, R2, R3,
     *                  R4, RA1, RA2, RAS, RAU, RB1, RB2, RS1, RS2, RSS,
     *                  RU1, RU2, S1, S2, S3, S4, T, TWOPI
      INTEGER           J, K, K0, M5, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN, SQRT
C     .. Executable Statements ..
      M5 = M*5
      MOVER2 = M/2 + 1
      TWOPI = 2.0D0*X01AAF(0.0D0)
      A1 = COS(TWOPI/5.0D0)
      B1 = SIN(TWOPI/5.0D0)
      A2 = COS(2.0D0*TWOPI/5.0D0)
      B2 = SIN(2.0D0*TWOPI/5.0D0)
      AS = -1.0D0/4.0D0
      AU = SQRT(5.0D0)/4.0D0
C
      DO 120 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M5)
         ZERO = ANGLE .EQ. 0.0D0
         C1 = COS(ANGLE)
         S1 = SIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         C3 = C2*C1 - S2*S1
         S3 = S2*C1 + C2*S1
         C4 = C2*C2 - S2*S2
         S4 = S2*C2 + C2*S2
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         T = C1*A1 + S1*B1
         S1 = C1*B1 - S1*A1
         C1 = T
         T = C2*A2 + S2*B2
         S2 = C2*B2 - S2*A2
         C2 = T
         T = C3*A2 - S3*B2
         S3 = -C3*B2 - S3*A2
         C3 = T
         T = C4*A1 - S4*B1
         S4 = -C4*B1 - S4*A1
         C4 = T
   40    CONTINUE
C
         DO 100 K = K0, PTS, M5
            R0 = X0(K)
            I0 = Y0(K)
            RS1 = X1(K) + X4(K)
            IS1 = Y1(K) + Y4(K)
            RU1 = X1(K) - X4(K)
            IU1 = Y1(K) - Y4(K)
            RS2 = X2(K) + X3(K)
            IS2 = Y2(K) + Y3(K)
            RU2 = X2(K) - X3(K)
            IU2 = Y2(K) - Y3(K)
            RSS = RS1 + RS2
            ISS = IS1 + IS2
            RAS = R0 + RSS*AS
            IAS = I0 + ISS*AS
            RAU = (RS1-RS2)*AU
            IAU = (IS1-IS2)*AU
            RA1 = RAS + RAU
            IA1 = IAS + IAU
            RA2 = RAS - RAU
            IA2 = IAS - IAU
            RB1 = RU1*B1 + RU2*B2
            IB1 = IU1*B1 + IU2*B2
            RB2 = RU1*B2 - RU2*B1
            IB2 = IU1*B2 - IU2*B1
            X0(K) = R0 + RSS
            Y0(K) = I0 + ISS
            IF (ZERO) GO TO 60
            R1 = RA1 + IB1
            I1 = IA1 - RB1
            R2 = RA2 + IB2
            I2 = IA2 - RB2
            R3 = RA2 - IB2
            I3 = IA2 + RB2
            R4 = RA1 - IB1
            I4 = IA1 + RB1
            X1(K) = R1*C1 + I1*S1
            Y1(K) = I1*C1 - R1*S1
            X2(K) = R2*C2 + I2*S2
            Y2(K) = I2*C2 - R2*S2
            X3(K) = R3*C3 + I3*S3
            Y3(K) = I3*C3 - R3*S3
            X4(K) = R4*C4 + I4*S4
            Y4(K) = I4*C4 - R4*S4
            GO TO 80
   60       CONTINUE
            X1(K) = RA1 + IB1
            Y1(K) = IA1 - RB1
            X2(K) = RA2 + IB2
            Y2(K) = IA2 - RB2
            X3(K) = RA2 - IB2
            Y3(K) = IA2 + RB2
            X4(K) = RA1 - IB1
            Y4(K) = IA1 + RB1
   80       CONTINUE
  100    CONTINUE
         IF (FOLD) GO TO 20
  120 CONTINUE
C
      RETURN
      END
      SUBROUTINE C06ECT(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,X3,Y3,M3,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX FOUR COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, M3, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), X3(M3), Y0(PTS),
     *                  Y1(M1), Y2(M2), Y3(M3)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, C1, C2, C3, I1, I2, I3, IS0, IS1, IU0,
     *                  IU1, R1, R2, R3, RS0, RS1, RU0, RU1, S1, S2, S3,
     *                  T, TWOPI
      INTEGER           J, K, K0, M4, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN
C     .. Executable Statements ..
      M4 = M*4
      MOVER2 = M/2 + 1
      TWOPI = 2.0D0*X01AAF(0.0D0)
C
      DO 120 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M4)
         ZERO = ANGLE .EQ. 0.0D0
         C1 = COS(ANGLE)
         S1 = SIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         C3 = C2*C1 - S2*S1
         S3 = S2*C1 + C2*S1
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         T = C1
         C1 = S1
         S1 = T
         C2 = -C2
         T = C3
         C3 = -S3
         S3 = -T
   40    CONTINUE
C
         DO 100 K = K0, PTS, M4
            RS0 = X0(K) + X2(K)
            IS0 = Y0(K) + Y2(K)
            RU0 = X0(K) - X2(K)
            IU0 = Y0(K) - Y2(K)
            RS1 = X1(K) + X3(K)
            IS1 = Y1(K) + Y3(K)
            RU1 = X1(K) - X3(K)
            IU1 = Y1(K) - Y3(K)
            X0(K) = RS0 + RS1
            Y0(K) = IS0 + IS1
            IF (ZERO) GO TO 60
            R1 = RU0 + IU1
            I1 = IU0 - RU1
            R2 = RS0 - RS1
            I2 = IS0 - IS1
            R3 = RU0 - IU1
            I3 = IU0 + RU1
            X2(K) = R1*C1 + I1*S1
            Y2(K) = I1*C1 - R1*S1
            X1(K) = R2*C2 + I2*S2
            Y1(K) = I2*C2 - R2*S2
            X3(K) = R3*C3 + I3*S3
            Y3(K) = I3*C3 - R3*S3
            GO TO 80
   60       CONTINUE
            X2(K) = RU0 + IU1
            Y2(K) = IU0 - RU1
            X1(K) = RS0 - RS1
            Y1(K) = IS0 - IS1
            X3(K) = RU0 - IU1
            Y3(K) = IU0 + RU1
   80       CONTINUE
  100    CONTINUE
         IF (FOLD) GO TO 20
  120 CONTINUE
C
      RETURN
      END
      SUBROUTINE C06ECU(X0,Y0,PTS,X1,Y1,M1,X2,Y2,M2,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX THREE COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, M2, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), X2(M2), Y0(PTS), Y1(M1), Y2(M2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ANGLE, B, C1, C2, I0, I1, I2, IA, IB, IS, R0,
     *                  R1, R2, RA, RB, RS, S1, S2, T, TWOPI
      INTEGER           J, K, K0, M3, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN, SQRT
C     .. Executable Statements ..
      M3 = M*3
      MOVER2 = M/2 + 1
      TWOPI = 2.0D0*X01AAF(0.0D0)
C     A = COS(TWOPI/3.0)
      A = -0.5D0
C     B = SIN(TWOPI/3.0)
      B = SQRT(0.75D0)
C
      DO 120 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M3)
         ZERO = ANGLE .EQ. 0.0D0
         C1 = COS(ANGLE)
         S1 = SIN(ANGLE)
         C2 = C1*C1 - S1*S1
         S2 = S1*C1 + C1*S1
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         T = C1*A + S1*B
         S1 = C1*B - S1*A
         C1 = T
         T = C2*A - S2*B
         S2 = -C2*B - S2*A
         C2 = T
   40    CONTINUE
C
         DO 100 K = K0, PTS, M3
            R0 = X0(K)
            I0 = Y0(K)
            RS = X1(K) + X2(K)
            IS = Y1(K) + Y2(K)
            X0(K) = R0 + RS
            Y0(K) = I0 + IS
            RA = R0 + RS*A
            IA = I0 + IS*A
            RB = (X1(K)-X2(K))*B
            IB = (Y1(K)-Y2(K))*B
            IF (ZERO) GO TO 60
            R1 = RA + IB
            I1 = IA - RB
            R2 = RA - IB
            I2 = IA + RB
            X1(K) = R1*C1 + I1*S1
            Y1(K) = I1*C1 - R1*S1
            X2(K) = R2*C2 + I2*S2
            Y2(K) = I2*C2 - R2*S2
            GO TO 80
   60       CONTINUE
            X1(K) = RA + IB
            Y1(K) = IA - RB
            X2(K) = RA - IB
            Y2(K) = IA + RB
   80       CONTINUE
  100    CONTINUE
         IF (FOLD) GO TO 20
  120 CONTINUE
C
      RETURN
      END
      SUBROUTINE C06ECV(X0,Y0,PTS,X1,Y1,M1,M)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX TWO COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, M1, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X0(PTS), X1(M1), Y0(PTS), Y1(M1)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, C, IS, IU, RS, RU, S, TWOPI
      INTEGER           J, K, K0, M2, MOVER2
      LOGICAL           FOLD, ZERO
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN
C     .. Executable Statements ..
      M2 = M*2
      MOVER2 = M/2 + 1
      TWOPI = 2.0D0*X01AAF(0.0D0)
C
      DO 120 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(M2)
         ZERO = ANGLE .EQ. 0.0D0
         C = COS(ANGLE)
         S = SIN(ANGLE)
         GO TO 40
   20    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         C = -C
   40    CONTINUE
C
         DO 100 K = K0, PTS, M2
            RS = X0(K) + X1(K)
            IS = Y0(K) + Y1(K)
            RU = X0(K) - X1(K)
            IU = Y0(K) - Y1(K)
            X0(K) = RS
            Y0(K) = IS
            IF (ZERO) GO TO 60
            X1(K) = RU*C + IU*S
            Y1(K) = IU*C - RU*S
            GO TO 80
   60       CONTINUE
            X1(K) = RU
            Y1(K) = IU
   80       CONTINUE
  100    CONTINUE
         IF (FOLD) GO TO 20
  120 CONTINUE
C
      RETURN
      END


c niveau 6

      
      DOUBLE PRECISION FUNCTION X01AAF(X)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS THE VALUE OF THE MATHEMATICAL CONSTANT PI.
C
C     X IS A DUMMY ARGUMENT
C
C     IT MAY BE NECESSARY TO ROUND THE REAL CONSTANT IN THE
C     ASSIGNMENT STATEMENT TO A SMALLER NUMBER OF SIGNIFICANT
C     DIGITS IN ORDER TO AVOID COMPILATION PROBLEMS.  IF SO, THEN
C     THE NUMBER OF DIGITS RETAINED SHOULD NOT BE LESS THAN
C     .     2 + INT(FLOAT(IT)*ALOG10(IB))
C     WHERE  IB  IS THE BASE FOR THE REPRESENTATION OF FLOATING-
C     .             POINT NUMBERS
C     . AND  IT  IS THE NUMBER OF IB-ARY DIGITS IN THE MANTISSA OF
C     .             A FLOATING-POINT NUMBER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Executable Statements ..
      X01AAF = 3.14159265358979323846264338328D0
      RETURN
      END
      SUBROUTINE F02ABF(A,IA,N,R,V,IV,E,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-735 (DEC 1989).
C
C     EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIX MATRIX
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02ABF')
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, IV, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), E(N), R(N), V(IV,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  XXXX
      INTEGER           ISAVE
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01AJF, F02AMF
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      CALL F01AJF(N,XXXX,A,IA,R,E,V,IV)
      CALL F02AMF(N,X02AJF(),R,E,V,IV,IFAIL)
      IF (IFAIL.NE.0) IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F01AJF(N,ATOL,A,IA,D,E,Z,IZ)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 15A REVISED. IER-906 (APR 1991).
C
C     TRED2
C     THIS SUBROUTINE REDUCES THE GIVEN LOWER TRIANGLE OF A
C     SYMMETRIC MATRIX, A, STORED IN THE ARRAY A(N,N), TO
C     TRIDIAGONAL FORM USING HOUSEHOLDERS REDUCTION. THE DIAGONAL
C     OF THE RESULT IS STORED IN THE ARRAY D(N) AND THE
C     SUB-DIAGONAL IN THE LAST N - 1 STORES OF THE ARRAY E(N)
C     (WITH THE ADDITIONAL ELEMENT E(1) = 0). THE TRANSFORMATION
C     MATRICES ARE ACCUMULATED IN THE ARRAY Z(N,N). THE ARRAY
C     A IS LEFT UNALTERED UNLESS THE ACTUAL PARAMETERS
C     CORRESPONDING TO A AND Z ARE IDENTICAL.
C     1ST AUGUST 1971
C
C     REVISED BY VINCE FERNANDO AT MARK 14 TO INTRODUCE SCALING INTO
C     THE GENERATION OF HOUSEHOLDER MATRICES AS PROPOSED BY
C     G.W. STEWART, INTRODUCTION TO MATRIX COMPUTATIONS, CHAPTER 7.
C     ATOL IS NOW A DUMMY PARAMETER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ATOL
      INTEGER           IA, IZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), D(N), E(N), Z(IZ,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F, G, H, HH, SCALE, SMALL
      INTEGER           I, II, J, K, L
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AMF
      INTEGER           IDAMAX
      EXTERNAL          DDOT, IDAMAX, X02AMF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER, DSCAL, DSYMV, DSYR2
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      SMALL = X02AMF()
      DO 40 I = 1, N
         DO 20 J = I, N
            Z(J,I) = A(J,I)
   20    CONTINUE
         D(I) = A(N,I)
   40 CONTINUE
      IF (N.EQ.1) GO TO 340
      DO 240 II = 2, N
         I = N - II + 2
         L = I - 1
         IF (L.EQ.1) GO TO 60
C        FIND THE ELEMENT OF LARGEST ABSOLUTE VALUE IN D
         K = IDAMAX(L,D,1)
         SCALE = ABS(D(K))
C        IF D IS A NULL VECTOR THEN SKIP THE TRANSFORMATION
         IF (SCALE.GE.SMALL) GO TO 120
   60    E(I) = D(L)
         H = 0.0D0
         DO 80 J = 1, L
            Z(J,I) = 0.0D0
   80    CONTINUE
         DO 100 J = 1, L
            Z(I,J) = 0.0D0
            D(J) = Z(I-1,J)
  100    CONTINUE
         GO TO 220
  120    CALL DSCAL(L,1.0D0/SCALE,D,1)
         H = DDOT(L,D,1,D,1)
         F = D(I-1)
         G = SQRT(H)
         IF (F.GE.0.0D0) G = -G
         E(I) = G*SCALE
         H = H - F*G
         D(I-1) = F - G
C        COPY U
         DO 140 J = 1, L
            Z(J,I) = D(J)
  140    CONTINUE
C        FORM A*U
         CALL DSYMV('L',L,1.0D0/H,Z,IZ,D,1,0.0D0,E,1)
C        FORM P
         F = 0.0D0
         DO 160 J = 1, L
            F = F + E(J)*D(J)
  160    CONTINUE
C        FORM K
         HH = F/(H+H)
C        FORM Q
         DO 180 J = 1, L
            E(J) = E(J) - HH*D(J)
  180    CONTINUE
C        FORM REDUCED A
         CALL DSYR2('L',L,-1.0D0,D,1,E,1,Z,IZ)
         DO 200 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
  200    CONTINUE
  220    D(I) = H
  240 CONTINUE
C     ACCUMULATION OF TRANSFORMATION MATRICES
      DO 300 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0D0
         H = D(I)
         IF (H.EQ.0.0D0) GO TO 260
         CALL DGEMV('T',L,L,1.0D0/H,Z,IZ,Z(1,I),1,0.0D0,D,1)
         CALL DGER(L,L,-1.0D0,Z(1,I),1,D,1,Z,IZ)
  260    DO 280 J = 1, L
            Z(J,I) = 0.0D0
  280    CONTINUE
  300 CONTINUE
      DO 320 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0D0
  320 CONTINUE
  340 Z(N,N) = 1.0D0
      E(1) = 0.0D0
      RETURN
      END
      SUBROUTINE F02AMF(N,EPS,D,E,Z,IZ,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     TQL2
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS OF A
C     TRIDIAGONAL MATRIX, T, GIVEN WITH ITS DIAGONAL ELEMENTS IN
C     THE ARRAY D(N) AND ITS SUB-DIAGONAL ELEMENTS IN THE LAST N
C     - 1 STORES OF THE ARRAY E(N), USING QL TRANSFORMATIONS. THE
C     EIGENVALUES ARE OVERWRITTEN ON THE DIAGONAL ELEMENTS IN THE
C     ARRAY D IN ASCENDING ORDER. THE EIGENVECTORS ARE FORMED IN
C     THE ARRAY Z(N,N), OVERWRITING THE ACCUMULATED
C     TRANSFORMATIONS AS SUPPLIED BY THE SUBROUTINE F01AJF. THE
C     SUBROUTINE WILL FAIL IF ALL EIGENVALUES TAKE MORE THAN 30*N
C     ITERATIONS.
C     1ST APRIL 1972
C
C     .. Parameters ..
      INTEGER           VLEN
      PARAMETER         (VLEN=128)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AMF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           IFAIL, IZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), E(N), Z(IZ,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, C, F, G, H, P, R, S
      INTEGER           I, I1, IPOS, ISAVE, ISEG, J, K, L, M
C     .. Local Arrays ..
      DOUBLE PRECISION  CC(VLEN), SS(VLEN)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06QXF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
      IF (N.EQ.1) GO TO 40
      DO 20 I = 2, N
         E(I-1) = E(I)
   20 CONTINUE
   40 E(N) = 0.0D0
      B = 0.0D0
      F = 0.0D0
      J = 30*N
      DO 300 L = 1, N
         H = EPS*(ABS(D(L))+ABS(E(L)))
         IF (B.LT.H) B = H
C        LOOK FOR SMALL SUB-DIAG ELEMENT
         DO 60 M = L, N
            IF (ABS(E(M)).LE.B) GO TO 80
   60    CONTINUE
   80    IF (M.EQ.L) GO TO 280
  100    IF (J.LE.0) GO TO 400
         J = J - 1
C        FORM SHIFT
         G = D(L)
         H = D(L+1) - G
         IF (ABS(H).GE.ABS(E(L))) GO TO 120
         P = H*0.5D0/E(L)
         R = SQRT(P*P+1.0D0)
         H = P + R
         IF (P.LT.0.0D0) H = P - R
         D(L) = E(L)/H
         GO TO 140
  120    P = 2.0D0*E(L)/H
         R = SQRT(P*P+1.0D0)
         D(L) = E(L)*P/(1.0D0+R)
  140    H = G - D(L)
         I1 = L + 1
         IF (I1.GT.N) GO TO 180
         DO 160 I = I1, N
            D(I) = D(I) - H
  160    CONTINUE
  180    F = F + H
C        QL TRANSFORMATION
         P = D(M)
         C = 1.0D0
         S = 0.0D0
         DO 260 K = M - 1, L, -VLEN
            ISEG = MAX(K-VLEN+1,L)
            DO 240 I = K, ISEG, -1
               G = C*E(I)
               H = C*P
               IF (ABS(P).LT.ABS(E(I))) GO TO 200
               C = E(I)/P
               R = SQRT(C*C+1.0D0)
               E(I+1) = S*P*R
               S = C/R
               C = 1.0D0/R
               GO TO 220
  200          C = P/E(I)
               R = SQRT(C*C+1.0D0)
               E(I+1) = S*E(I)*R
               S = 1.0D0/R
               C = C/R
  220          P = C*D(I) - S*G
               D(I+1) = H + S*(C*G+S*D(I))
C           STORE ROTATIONS
               CC(VLEN-K+I) = C
               SS(VLEN-K+I) = -S
  240       CONTINUE
C        UPDATE VECTORS
            IPOS = VLEN - K + ISEG
            CALL F06QXF('Right','Variable','Backward',N,K-ISEG+2,1,
     *                  K-ISEG+2,CC(IPOS),SS(IPOS),Z(1,ISEG),IZ)
  260    CONTINUE
         E(L) = S*P
         D(L) = C*P
         IF (ABS(E(L)).GT.B) GO TO 100
  280    D(L) = D(L) + F
  300 CONTINUE
C     ORDER EIGENVALUES AND EIGENVECTORS
      DO 380 I = 1, N
         K = I
         P = D(I)
         I1 = I + 1
         IF (I1.GT.N) GO TO 340
         DO 320 J = I1, N
            IF (D(J).GE.P) GO TO 320
            K = J
            P = D(J)
  320    CONTINUE
  340    IF (K.EQ.I) GO TO 380
         D(K) = D(I)
         D(I) = P
         DO 360 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  360    CONTINUE
  380 CONTINUE
      IFAIL = 0
      RETURN
  400 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F06QXF( SIDE, PIVOT, DIRECT, M, N, K1, K2, C, S, A,
     $                   LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, M, N
      CHARACTER*1        DIRECT, PIVOT, SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QXF  performs the transformation
C
C     A := P*A,   when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C     A := A*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where A is an m by n matrix and P is an orthogonal matrix, consisting
C  of a  sequence  of  plane  rotations,  applied  in  planes  k1 to k2,
C  determined by the parameters PIVOT and DIRECT as follows:
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C  c( k ) and s( k )  must contain the  cosine and sine  that define the
C  matrix  P( k ).  The  two by two  plane rotation  part of the  matrix
C  P( k ), R( k ), is assumed to be of the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  If m, n or k1 are less than unity,  or k2 is not greater than k1,  or
C  SIDE = 'L' or 'l'  and  k2  is greater than  m, or  SIDE = 'R' or 'r'
C  and  k2  is greater than  n,  then an  immediate return  is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 20-November-1986.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, TEMP
      INTEGER            I, J
      LOGICAL            LEFT, RIGHT
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      IF( ( MIN( M, N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $    ( ( LEFT ).AND.( K2.GT.M ) ).OR.
     $    ( ( RIGHT ).AND.( K2.GT.N ) ) )RETURN
      IF( LEFT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 20 J = 1, N
                  AIJ = A( K1, J )
                  DO 10 I = K1, K2 - 1
                     TEMP = A( I + 1, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     AIJ = C( I )*TEMP - S( I )*AIJ
   10             CONTINUE
                  A( K2, J ) = AIJ
   20          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 40 J = 1, N
                  AIJ = A( K2, J )
                  DO 30 I = K2 - 1, K1, -1
                     TEMP = A( I, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     AIJ = S( I )*AIJ + C( I )*TEMP
   30             CONTINUE
                  A( K1, J ) = AIJ
   40          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 60 J = 1, N
                  TEMP = A( K1, J )
                  DO 50 I = K1, K2 - 1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   50             CONTINUE
                  A( K1, J ) = TEMP
   60          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 80 J = 1, N
                  TEMP = A( K1, J )
                  DO 70 I = K2 - 1, K1, -1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   70             CONTINUE
                  A( K1, J ) = TEMP
   80          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 100 J = 1, N
                  TEMP = A( K2, J )
                  DO 90 I = K1, K2 - 1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
   90             CONTINUE
                  A( K2, J ) = TEMP
  100          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 120 J = 1, N
                  TEMP = A( K2, J )
                  DO 110 I = K2 - 1, K1, -1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
  110             CONTINUE
                  A( K2, J ) = TEMP
  120          CONTINUE
            END IF
         END IF
      ELSE IF( RIGHT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 140 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 130 I = 1, M
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 160 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 150 I = M, 1, -1
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 180 J = K1 + 1, K2
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 200 J = K2, K1 + 1, -1
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 190 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 220 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 240 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 230 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06QXF. ( SGESRC )
C
      END
      SUBROUTINE F06PCF( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DSYMV  performs the matrix-vector  operation
C
C     y := alpha*A*x + beta*y,
C
C  where alpha and beta are scalars, x and y are n element vectors and
C  A is an n by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y. On exit, Y is overwritten by the updated
C           vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PCF/DSYMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set up the start points in  X  and  Y.
C
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  y  when A is stored in upper triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( J, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, I = 1, J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( J, J ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y  when A is stored in lower triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*A( J, J )
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*A( J, J )
               IX      = JX
               IY      = JY
               DO 110, I = J + 1, N
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PCF (DSYMV ).
C
      END

      SUBROUTINE F06PRF( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DSYR2  performs the symmetric rank 2 operation
C
C     A := alpha*x*y' + alpha*y*x' + A,
C
C  where alpha is a scalar, x and y are n element vectors and A is an n
C  by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced. On exit, the
C           upper triangular part of the array A is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced. On exit, the
C           lower triangular part of the array A is overwritten by the
C           lower triangular part of the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PRF/DSYR2 ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Set up the start points in X and Y if the increments are not both
C     unity.
C
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  A  when A is stored in the upper triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = KX
                  IY    = KY
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   30             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
   40       CONTINUE
         END IF
      ELSE
C
C        Form  A  when A is stored in the lower triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = JX
                  IY    = JY
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PRF (DSYR2 ).
C
      END
      SUBROUTINE F06PAF( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PAF/DGEMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  y := alpha*A*x + y.
C
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y := alpha*A'*x + y.
C
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PAF (DGEMV ).
C
      END
      SUBROUTINE F06PMF( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGER   performs the rank 1 operation
C
C     A := alpha*x*y' + A,
C
C  where alpha is a scalar, x is an m element vector, y is an n element
C  vector and A is an m by n matrix.
C
C  Parameters
C  ==========
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( m - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the m
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients. On exit, A is
C           overwritten by the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PMF/DGER  ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06PMF (DGER  ).
C
      END
