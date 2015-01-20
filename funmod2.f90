module funmod2

contains

subroutine matout(A,m,n)
integer :: i,m,n
double precision :: A(m,n)
do i = 1,m
   print "(10f16.8)", A(i,:); call flush(6)
end do
end subroutine matout


subroutine bytes_in_rec( bytes )
  implicit none
  integer     bytes

  character*8 string
  integer     i
  integer     ierr
  double precision d;
  real r;

  d = dble(1.0);
  r = 1.0;
  bytes = 0
  do i = 1,8
     !open(unit=15,file='/home/jason/tmpbytetst',access='direct',recl=i)
     open(unit=15,status='scratch',access='direct',recl=i)
     write(15,rec=1,iostat=ierr) r
     close(15, status='delete')
     !print *, 'i = ', i, ' real ierr = ', ierr
     if (ierr == 0) then
           bytes = i
           exit
     end if 
     !open(unit=16,file='/home/jason/tmpbytetst2',access='direct',recl=i)
     !write(16,rec=1,iostat=ierr) d
     !close(16, status='delete')
     !print *, 'i = ', i, ' double ierr = ', ierr 
  end do
  print *, 'bytes in real = ', bytes
  !open( 10, status = 'scratch', access = 'direct', recl = 1 )
  !write( 10,rec=1,iostat=ierr) d

  !do i = 1,8
  !   write( 10, rec = 1, iostat = ierr ) string(1:i)
  !   print *, 'i = ', i, ' ierr = ', ierr
  !   if ( ierr /= 0 ) exit
  !   bytes = i
  !end do
  !close( 10, status = 'delete' )
  
end subroutine bytes_in_rec


!--------------------------------------------------------------------------------------
! The following gamma/psi code was obtained from Alan Miller's Fortran Software website.

FUNCTION gamln (a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!          NAVAL SURFACE WARFARE CENTER
!          DAHLGREN, VIRGINIA
!--------------------------
!     D = 0.5*(LN(2*PI) - 1)
!--------------------------
IMPLICIT NONE
double precision, INTENT(IN) :: a
double precision             :: fn_val

double precision :: c0 = .833333333333333D-01, c1 = -.277777777760991D-02,  &
             c2 = .793650666825390D-03, c3 = -.595202931351870D-03,  &
             c4 = .837308034031215D-03, c5 = -.165322962780713D-02,  &
             d = .418938533204673D0, t, w
INTEGER   :: i, n
!--------------------------
IF (a > 0.8D0) GO TO 10
fn_val = gamln1(a) - LOG(a)
RETURN
10 IF (a > 2.25D0) GO TO 20
t = (a - 0.5D0) - 0.5D0
fn_val = gamln1(t)
RETURN

20 IF (a >= 10.0D0) GO TO 30
n = a - 1.25D0
t = a
w = 1.0D0
DO i = 1, n
  t = t - 1.0D0
  w = t*w
END DO
fn_val = gamln1(t - 1.0D0) + LOG(w)
RETURN

30 t = (1.0D0/a)**2
w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a
fn_val = (d + w) + (a - 0.5D0)*(LOG(a) - 1.0D0)

RETURN
END FUNCTION gamln



FUNCTION gamln1 (a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
!-----------------------------------------------------------------------
IMPLICIT NONE
double precision, INTENT(IN) :: a
double precision             :: fn_val

double precision :: w, x, &
             p0 =  .577215664901533D+00, p1 =  .844203922187225D+00,  &
             p2 = -.168860593646662D+00, p3 = -.780427615533591D+00,  &
             p4 = -.402055799310489D+00, p5 = -.673562214325671D-01,  &
             p6 = -.271935708322958D-02,   &
             q1 =  .288743195473681D+01, q2 =  .312755088914843D+01,  &
             q3 =  .156875193295039D+01, q4 =  .361951990101499D+00,  &
             q5 =  .325038868253937D-01, q6 =  .667465618796164D-03,  &
             r0 = .422784335098467D+00,  r1 = .848044614534529D+00,  &
             r2 = .565221050691933D+00,  r3 = .156513060486551D+00,  &
             r4 = .170502484022650D-01,  r5 = .497958207639485D-03,  &
             s1 = .124313399877507D+01,  s2 = .548042109832463D+00,  &
             s3 = .101552187439830D+00,  s4 = .713309612391000D-02,  &
             s5 = .116165475989616D-03
!----------------------
IF (a >= 0.6D0) GO TO 10
w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0)/  &
    ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0D0)
fn_val = -a*w
RETURN

10 x = (a - 0.5D0) - 0.5D0
w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0)/  &
    (((((s5*x + s4)*x + s3)*x + s2)*x + s1)*x + 1.0D0)
fn_val = x*w
RETURN
END FUNCTION gamln1

!-----------------------------------------------------------------------

FUNCTION psifun(xx) RESULT(fn_val)
!---------------------------------------------------------------------

!                 EVALUATION OF THE DIGAMMA FUNCTION

!                           -----------

!     PSIFUN(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT
!     BE COMPUTED.

!     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV
!     APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY
!     CODY, STRECOK AND THACHER.

!---------------------------------------------------------------------
!     PSIFUN WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK
!     PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSIFUN WAS MODIFIED BY
!     A.H. MORRIS (NSWC).
!---------------------------------------------------------------------
IMPLICIT NONE
double precision, INTENT(IN) :: xx
double precision             :: fn_val

double precision :: dx0 = 1.461632144968362341262659542325721325D0
!---------------------------------------------------------------------

!     PIOV4 = PI/4
!     DX0 = ZERO OF PSIFUN TO EXTENDED PRECISION

!---------------------------------------------------------------------
double precision :: aug, den, piov4 = .785398163397448D0, sgn, upper,  &
             w, x, xmax1, xmx0, xsmall, z
INTEGER   :: i, m, n, nq
!---------------------------------------------------------------------

!     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
!     PSIFUN(X) / (X - X0),  0.5 <= X <= 3.0

!---------------------------------------------------------------------
double precision :: p1(7) = (/ .895385022981970D-02, .477762828042627D+01,  &
                        .142441585084029D+03, .118645200713425D+04,  &
                        .363351846806499D+04, .413810161269013D+04,  &
                        .130560269827897D+04 /),   &
             q1(6) = (/ .448452573429826D+02, .520752771467162D+03,  &
                        .221000799247830D+04, .364127349079381D+04,  &
                        .190831076596300D+04, .691091682714533D-05 /)
!---------------------------------------------------------------------

!     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
!     PSIFUN(X) - LN(X) + 1 / (2*X),  X > 3.0

!---------------------------------------------------------------------
double precision :: p2(4) = (/ -.212940445131011D+01, -.701677227766759D+01,  &
                        -.448616543918019D+01, -.648157123766197D+00 /), &
             q2(4) = (/  .322703493791143D+02,  .892920700481861D+02,  &
                         .546117738103215D+02,  .777788548522962D+01 /)
!---------------------------------------------------------------------

!     MACHINE DEPENDENT CONSTANTS ...

!        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
!                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
!                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
!                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
!                 PSIFUN MAY BE REPRESENTED AS ALOG(X).

!        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
!                 MAY BE REPRESENTED BY 1/X.

!---------------------------------------------------------------------
xmax1 = ipmpar(3)
xmax1 = MIN(xmax1, 1.0D0/dpmpar(1))
xsmall = 1.d-9
!---------------------------------------------------------------------
x = xx
aug = 0.0D0
IF (x >= 0.5D0) GO TO 200
!---------------------------------------------------------------------
!     X .LT. 0.5,  USE REFLECTION FORMULA
!     PSIFUN(1-X) = PSI(X) + PI * COTAN(PI*X)
!---------------------------------------------------------------------
IF (ABS(x) > xsmall) GO TO 100
IF (x == 0.0D0) GO TO 400
!---------------------------------------------------------------------
!     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
!     FOR  PI*COTAN(PI*X)
!---------------------------------------------------------------------
aug = -1.0D0 / x
GO TO 150
!---------------------------------------------------------------------
!     REDUCTION OF ARGUMENT FOR COTAN
!---------------------------------------------------------------------
100 w = - x
sgn = piov4
IF (w > 0.0D0) GO TO 120
w = - w
sgn = -sgn
!---------------------------------------------------------------------
!     MAKE AN ERROR EXIT IF X .LE. -XMAX1
!---------------------------------------------------------------------
120 IF (w >= xmax1) GO TO 400
nq = INT(w)
w = w - nq
nq = INT(w*4.0D0)
w = 4.0D0 * (w - nq * .25D0)
!---------------------------------------------------------------------
!     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
!     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
!     QUADRANT AND DETERMINE SIGN
!---------------------------------------------------------------------
n = nq / 2
IF ((n+n) /= nq) w = 1.0D0 - w
z = piov4 * w
m = n / 2
IF ((m+m) /= n) sgn = - sgn
!---------------------------------------------------------------------
!     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
!---------------------------------------------------------------------
n = (nq + 1) / 2
m = n / 2
m = m + m
IF (m /= n) GO TO 140
!---------------------------------------------------------------------
!     CHECK FOR SINGULARITY
!---------------------------------------------------------------------
IF (z == 0.0D0) GO TO 400
!---------------------------------------------------------------------
!     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
!     SIN/COS AS A SUBSTITUTE FOR TAN
!---------------------------------------------------------------------
aug = sgn * ((COS(z) / SIN(z)) * 4.0D0)
GO TO 150
140 aug = sgn * ((SIN(z) / COS(z)) * 4.0D0)
150 x = 1.0D0 - x
200 IF (x > 3.0D0) GO TO 300
!---------------------------------------------------------------------
!     0.5 .LE. X .LE. 3.0
!---------------------------------------------------------------------
den = x
upper = p1(1) * x

DO i = 1, 5
  den = (den + q1(i)) * x
  upper = (upper + p1(i+1)) * x
END DO

den = (upper + p1(7)) / (den + q1(6))
xmx0 = x - dx0
fn_val = den * xmx0 + aug
RETURN
!---------------------------------------------------------------------
!     IF X .GE. XMAX1, PSIFUN = LN(X)
!---------------------------------------------------------------------
300 IF (x >= xmax1) GO TO 350
!---------------------------------------------------------------------
!     3.0 .LT. X .LT. XMAX1
!---------------------------------------------------------------------
w = 1.0D0 / (x * x)
den = w
upper = p2(1) * w

DO i = 1, 3
  den = (den + q2(i)) * w
  upper = (upper + p2(i+1)) * w
END DO

aug = upper / (den + q2(4)) - 0.5D0 / x + aug
350 fn_val = aug + LOG(x)
RETURN
!---------------------------------------------------------------------
!     ERROR RETURN
!---------------------------------------------------------------------
400 fn_val = 0.0D0
RETURN
END FUNCTION psifun


FUNCTION ipmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
!     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
!     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...

!  INTEGERS.

!     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM

!               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )

!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.

!     IPMPAR(1) = A, THE BASE (radix).

!     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS (digits).

!     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE (huge).

!  FLOATING-POINT NUMBERS.

!     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
!     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
!     NONZERO NUMBERS ARE REPRESENTED IN THE FORM

!               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)

!               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
!               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.

!     IPMPAR(4) = B, THE BASE.

!  SINGLE-PRECISION

!     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.

!  DOUBLE-PRECISION

!     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.

!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
INTEGER             :: fn_val

SELECT CASE(i)
  CASE( 1)
    fn_val = RADIX(i)
  CASE( 2)
    fn_val = DIGITS(i)
  CASE( 3)
    fn_val = HUGE(i)
  CASE( 4)
    fn_val = RADIX(1.0)
  CASE( 5)
    fn_val = DIGITS(1.0)
  CASE( 6)
    fn_val = MINEXPONENT(1.0)
  CASE( 7)
    fn_val = MAXEXPONENT(1.0)
  CASE( 8)
    fn_val = DIGITS(1.0D0)
  CASE( 9)
    fn_val = MINEXPONENT(1.0D0)
  CASE(10)
    fn_val = MAXEXPONENT(1.0D0)
  CASE DEFAULT
    RETURN
END SELECT

RETURN
END FUNCTION ipmpar

FUNCTION spmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
double precision    :: fn_val

! Local variable
double precision                :: one = 1.0

SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
END SELECT

RETURN
END FUNCTION spmpar


FUNCTION dpmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     DPMPAR PROVIDES THE DOUBLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     DOUBLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
double precision           :: fn_val

! Local variable
double precision    :: one = 1.0

SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
END SELECT

RETURN
END FUNCTION dpmpar


FUNCTION epsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.0 + EPS .GT. 1.0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
IMPLICIT NONE
double precision                :: fn_val

! Local variable
double precision                :: one = 1.0

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION epsln


FUNCTION exparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     EXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
double precision                :: fn_val

! Local variable
double precision                :: one = 1.0

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION exparg


FUNCTION depsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.D0 + EPS .GT. 1.D0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
IMPLICIT NONE
double precision           :: fn_val

! Local variable
double precision    :: one = 1.0

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION depsln


FUNCTION dxparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  DXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     DEXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  DXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF DEXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR DXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
double precision           :: fn_val

! Local variable
double precision    :: one = 1.0

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION dxparg





!----------------------------------------------------------------------
end module funmod2
