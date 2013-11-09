module ddassl_mod

  type, abstract :: problem_ddassl
     contains
       procedure(res_i), deferred :: res
       procedure(jac_i), deferred :: jac
  end type problem_ddassl

  abstract interface

     subroutine res_i(eqn, t, y, ydot, res, ires, rwk, iwk)
       import problem_ddassl
       class(problem_ddassl) :: eqn
       integer :: ires
       integer, dimension(*) :: iwk
       real(8) :: t
       real(8), dimension(*) :: y, ydot, res, rwk
     end subroutine res_i

     subroutine jac_i(eqn, t, y, ydot, pd,  cj,   rwk, iwk)
       import problem_ddassl
       class(problem_ddassl) :: eqn
       integer, dimension(*) :: iwk
       real(8) :: t, cj
       real(8), dimension(*) :: y, ydot, rwk
       real(8), dimension(10,*) :: pd
     end subroutine jac_i

  end interface

contains

!DECK DDASSL
!-----------------------------------------------------------------------
! NOTE:  Users of this solver, DDASSL, are encouraged to use the
! solver DDASPK instead.  DDASPK has a much improved initial condition
! calculation algorithm.  In addition, DDASPK includes iterative
! (Krylov) methods for the linear systems that arise, in addition to
! the direct (dense/banded) methods in DDASSL.
!-----------------------------------------------------------------------
!***BEGIN PROLOGUE  DDASSL
!***PURPOSE  This code solves a system of differential/algebraic
!            equations of the form G(T,Y,YPRIME) = 0.
!***LIBRARY   SLATEC (DASSL)
!***CATEGORY  I1A2
!***TYPE      REAL(8) (SDASSL-S, DDASSL-D)
!***KEYWORDS  BACKWARD DIFFERENTIATION FORMULAS, DASSL,
!             DIFFERENTIAL/ALGEBRAIC, IMPLICIT DIFFERENTIAL SYSTEMS
!***AUTHOR  Petzold, Linda R., (LLNL)
!             Computing and Mathematics Research Division
!             Lawrence Livermore National Laboratory
!             L - 316, P.O. Box 808,
!             Livermore, CA.    94550
!***DESCRIPTION
!
! *Usage:
!
!      EXTERNAL RES, JAC
!      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR
!      REAL(8) T, Y(NEQ), YPRIME(NEQ), TOUT, RTOL, ATOL,
!     *   RWORK(LRW), RPAR
!
!      CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
!     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
!
!
! *Arguments:
!  (In the following, all real arrays should be type REAL(8).)
!
!  RES:EXT     This is a subroutine which you provide to define the
!              differential/algebraic system.
!
!  NEQ:IN      This is the number of equations to be solved.
!
!  T:INOUT     This is the current value of the independent variable.
!
!  Y(*):INOUT  This array contains the solution components at T.
!
!  YPRIME(*):INOUT  This array contains the derivatives of the solution
!              components at T.
!
!  TOUT:IN     This is a point at which a solution is desired.
!
!  INFO(N):IN  The basic task of the code is to solve the system from T
!              to TOUT and return an answer at TOUT.  INFO is an integer
!              array which is used to communicate exactly how you want
!              this task to be carried out.  (See below for details.)
!              N must be greater than or equal to 15.
!
!  RTOL,ATOL:INOUT  These quantities represent relative and absolute
!              error tolerances which you provide to indicate how
!              accurately you wish the solution to be computed.  You
!              may choose them to be both scalars or else both vectors.
!              Caution:  In Fortran 77, a scalar is not the same as an
!                        array of length 1.  Some compilers may object
!                        to using scalars for RTOL,ATOL.
!
!  IDID:OUT    This scalar quantity is an indicator reporting what the
!              code did.  You must monitor this integer variable to
!              decide  what action to take next.
!
!  RWORK:WORK  A real work array of length LRW which provides the
!              code with needed storage space.
!
!  LRW:IN      The length of RWORK.  (See below for required length.)
!
!  IWORK:WORK  An integer work array of length LIW which provides the
!              code with needed storage space.
!
!  LIW:IN      The length of IWORK.  (See below for required length.)
!
!  RPAR,IPAR:IN  These are real and integer parameter arrays which
!              you can use for communication between your calling
!              program and the RES subroutine (and the JAC subroutine)
!
!  JAC:EXT     This is the name of a subroutine which you may choose
!              to provide for defining a matrix of partial derivatives
!              described below.
!
!  Quantities which may be altered by DDASSL are:
!     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL,
!     IDID, RWORK(*) AND IWORK(*)
!
! *Description
!
!  Subroutine DDASSL uses the backward differentiation formulas of
!  orders one through five to solve a system of the above form for Y and
!  YPRIME.  Values for Y and YPRIME at the initial time must be given as
!  input.  These values must be consistent, (that is, if T,Y,YPRIME are
!  the given initial values, they must satisfy G(T,Y,YPRIME) = 0.).  The
!  subroutine solves the system from T to TOUT.  It is easy to continue
!  the solution to get results at additional TOUT.  This is the interval
!  mode of operation.  Intermediate results can also be obtained easily
!  by using the intermediate-output capability.
!
!  The following detailed description is divided into subsections:
!    1. Input required for the first call to DDASSL.
!    2. Output after any return from DDASSL.
!    3. What to do to continue the integration.
!    4. Error messages.
!
!
!  -------- INPUT -- WHAT TO DO ON THE FIRST CALL TO DDASSL ------------
!
!  The first call of the code is defined to be the start of each new
!  problem. Read through the descriptions of all the following items,
!  provide sufficient storage space for designated arrays, set
!  appropriate variables for the initialization of the problem, and
!  give information about how you want the problem to be solved.
!
!
!  RES -- Provide a subroutine of the form
!             SUBROUTINE RES(T,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
!         to define the system of differential/algebraic
!         equations which is to be solved. For the given values
!         of T,Y and YPRIME, the subroutine should
!         return the residual of the differential/algebraic
!         system
!             DELTA = G(T,Y,YPRIME)
!         (DELTA(*) is a vector of length NEQ which is
!         output for RES.)
!
!         Subroutine RES must not alter T,Y or YPRIME.
!         You must declare the name RES in an external
!         statement in your program that calls DDASSL.
!         You must dimension Y,YPRIME and DELTA in RES.
!
!         IRES is an integer flag which is always equal to
!         zero on input. Subroutine RES should alter IRES
!         only if it encounters an illegal value of Y or
!         a stop condition. Set IRES = -1 if an input value
!         is illegal, and DDASSL will try to solve the problem
!         without getting IRES = -1. If IRES = -2, DDASSL
!         will return control to the calling program
!         with IDID = -11.
!
!         RPAR and IPAR are real and integer parameter arrays which
!         you can use for communication between your calling program
!         and subroutine RES. They are not altered by DDASSL. If you
!         do not need RPAR or IPAR, ignore these parameters by treat-
!         ing them as dummy arguments. If you do choose to use them,
!         dimension them in your calling program and in RES as arrays
!         of appropriate length.
!
!  NEQ -- Set it to the number of differential equations.
!         (NEQ .GE. 1)
!
!  T -- Set it to the initial point of the integration.
!         T must be defined as a variable.
!
!  Y(*) -- Set this vector to the initial values of the NEQ solution
!         components at the initial point. You must dimension Y of
!         length at least NEQ in your calling program.
!
!  YPRIME(*) -- Set this vector to the initial values of the NEQ
!         first derivatives of the solution components at the initial
!         point.  You must dimension YPRIME at least NEQ in your
!         calling program. If you do not know initial values of some
!         of the solution components, see the explanation of INFO(11).
!
!  TOUT -- Set it to the first point at which a solution
!         is desired. You can not take TOUT = T.
!         integration either forward in T (TOUT .GT. T) or
!         backward in T (TOUT .LT. T) is permitted.
!
!         The code advances the solution from T to TOUT using
!         step sizes which are automatically selected so as to
!         achieve the desired accuracy. If you wish, the code will
!         return with the solution and its derivative at
!         intermediate steps (intermediate-output mode) so that
!         you can monitor them, but you still must provide TOUT in
!         accord with the basic aim of the code.
!
!         The first step taken by the code is a critical one
!         because it must reflect how fast the solution changes near
!         the initial point. The code automatically selects an
!         initial step size which is practically always suitable for
!         the problem. By using the fact that the code will not step
!         past TOUT in the first step, you could, if necessary,
!         restrict the length of the initial step size.
!
!         For some problems it may not be permissible to integrate
!         past a point TSTOP because a discontinuity occurs there
!         or the solution or its derivative is not defined beyond
!         TSTOP. When you have declared a TSTOP point (SEE INFO(4)
!         and RWORK(1)), you have told the code not to integrate
!         past TSTOP. In this case any TOUT beyond TSTOP is invalid
!         input.
!
!  INFO(*) -- Use the INFO array to give the code more details about
!         how you want your problem solved.  This array should be
!         dimensioned of length 15, though DDASSL uses only the first
!         eleven entries.  You must respond to all of the following
!         items, which are arranged as questions.  The simplest use
!         of the code corresponds to answering all questions as yes,
!         i.e. setting all entries of INFO to 0.
!
!       INFO(1) - This parameter enables the code to initialize
!              itself. You must set it to indicate the start of every
!              new problem.
!
!          **** Is this the first call for this problem ...
!                Yes - Set INFO(1) = 0
!                 No - Not applicable here.
!                      See below for continuation calls.  ****
!
!       INFO(2) - How much accuracy you want of your solution
!              is specified by the error tolerances RTOL and ATOL.
!              The simplest use is to take them both to be scalars.
!              To obtain more flexibility, they can both be vectors.
!              The code must be told your choice.
!
!          **** Are both error tolerances RTOL, ATOL scalars ...
!                Yes - Set INFO(2) = 0
!                      and input scalars for both RTOL and ATOL
!                 No - Set INFO(2) = 1
!                      and input arrays for both RTOL and ATOL ****
!
!       INFO(3) - The code integrates from T in the direction
!              of TOUT by steps. If you wish, it will return the
!              computed solution and derivative at the next
!              intermediate step (the intermediate-output mode) or
!              TOUT, whichever comes first. This is a good way to
!              proceed if you want to see the behavior of the solution.
!              If you must have solutions at a great many specific
!              TOUT points, this code will compute them efficiently.
!
!          **** Do you want the solution only at
!                TOUT (and not at the next intermediate step) ...
!                 Yes - Set INFO(3) = 0
!                  No - Set INFO(3) = 1 ****
!
!       INFO(4) - To handle solutions at a great many specific
!              values TOUT efficiently, this code may integrate past
!              TOUT and interpolate to obtain the result at TOUT.
!              Sometimes it is not possible to integrate beyond some
!              point TSTOP because the equation changes there or it is
!              not defined past TSTOP. Then you must tell the code
!              not to go past.
!
!           **** Can the integration be carried out without any
!                restrictions on the independent variable T ...
!                 Yes - Set INFO(4)=0
!                  No - Set INFO(4)=1
!                       and define the stopping point TSTOP by
!                       setting RWORK(1)=TSTOP ****
!
!       INFO(5) - To solve differential/algebraic problems it is
!              necessary to use a matrix of partial derivatives of the
!              system of differential equations. If you do not
!              provide a subroutine to evaluate it analytically (see
!              description of the item JAC in the call list), it will
!              be approximated by numerical differencing in this code.
!              although it is less trouble for you to have the code
!              compute partial derivatives by numerical differencing,
!              the solution will be more reliable if you provide the
!              derivatives via JAC. Sometimes numerical differencing
!              is cheaper than evaluating derivatives in JAC and
!              sometimes it is not - this depends on your problem.
!
!           **** Do you want the code to evaluate the partial
!                derivatives automatically by numerical differences ...
!                   Yes - Set INFO(5)=0
!                    No - Set INFO(5)=1
!                  and provide subroutine JAC for evaluating the
!                  matrix of partial derivatives ****
!
!       INFO(6) - DDASSL will perform much better if the matrix of
!              partial derivatives, DG/DY + CJ*DG/DYPRIME,
!              (here CJ is a scalar determined by DDASSL)
!              is banded and the code is told this. In this
!              case, the storage needed will be greatly reduced,
!              numerical differencing will be performed much cheaper,
!              and a number of important algorithms will execute much
!              faster. The differential equation is said to have
!              half-bandwidths ML (lower) and MU (upper) if equation i
!              involves only unknowns Y(J) with
!                             I-ML .LE. J .LE. I+MU
!              for all I=1,2,...,NEQ. Thus, ML and MU are the widths
!              of the lower and upper parts of the band, respectively,
!              with the main diagonal being excluded. If you do not
!              indicate that the equation has a banded matrix of partial
!              derivatives, the code works with a full matrix of NEQ**2
!              elements (stored in the conventional way). Computations
!              with banded matrices cost less time and storage than with
!              full matrices if 2*ML+MU .LT. NEQ. If you tell the
!              code that the matrix of partial derivatives has a banded
!              structure and you want to provide subroutine JAC to
!              compute the partial derivatives, then you must be careful
!              to store the elements of the matrix in the special form
!              indicated in the description of JAC.
!
!          **** Do you want to solve the problem using a full
!               (dense) matrix (and not a special banded
!               structure) ...
!                Yes - Set INFO(6)=0
!                 No - Set INFO(6)=1
!                       and provide the lower (ML) and upper (MU)
!                       bandwidths by setting
!                       IWORK(1)=ML
!                       IWORK(2)=MU ****
!
!
!        INFO(7) -- You can specify a maximum (absolute value of)
!              stepsize, so that the code
!              will avoid passing over very
!              large regions.
!
!          ****  Do you want the code to decide
!                on its own maximum stepsize?
!                Yes - Set INFO(7)=0
!                 No - Set INFO(7)=1
!                      and define HMAX by setting
!                      RWORK(2)=HMAX ****
!
!        INFO(8) -- Differential/algebraic problems
!              may occasionally suffer from
!              severe scaling difficulties on the
!              first step. If you know a great deal
!              about the scaling of your problem, you can
!              help to alleviate this problem by
!              specifying an initial stepsize HO.
!
!          ****  Do you want the code to define
!                its own initial stepsize?
!                Yes - Set INFO(8)=0
!                 No - Set INFO(8)=1
!                      and define HO by setting
!                      RWORK(3)=HO ****
!
!        INFO(9) -- If storage is a severe problem,
!              you can save some locations by
!              restricting the maximum order MAXORD.
!              the default value is 5. for each
!              order decrease below 5, the code
!              requires NEQ fewer locations, however
!              it is likely to be slower. In any
!              case, you must have 1 .LE. MAXORD .LE. 5
!          ****  Do you want the maximum order to
!                default to 5?
!                Yes - Set INFO(9)=0
!                 No - Set INFO(9)=1
!                      and define MAXORD by setting
!                      IWORK(3)=MAXORD ****
!
!        INFO(10) --If you know that the solutions to your equations
!               will always be nonnegative, it may help to set this
!               parameter. However, it is probably best to
!               try the code without using this option first,
!               and only to use this option if that doesn't
!               work very well.
!           ****  Do you want the code to solve the problem without
!                 invoking any special nonnegativity constraints?
!                  Yes - Set INFO(10)=0
!                   No - Set INFO(10)=1
!
!        INFO(11) --DDASSL normally requires the initial T,
!               Y, and YPRIME to be consistent. That is,
!               you must have G(T,Y,YPRIME) = 0 at the initial
!               time. If you do not know the initial
!               derivative precisely, you can let DDASSL try
!               to compute it.
!          ****   Are the initial T, Y, YPRIME consistent?
!                 Yes - Set INFO(11) = 0
!                  No - Set INFO(11) = 1,
!                       and set YPRIME to an initial approximation
!                       to YPRIME.  (If you have no idea what
!                       YPRIME should be, set it to zero. Note
!                       that the initial Y should be such
!                       that there must exist a YPRIME so that
!                       G(T,Y,YPRIME) = 0.)
!
!  RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL
!         error tolerances to tell the code how accurately you
!         want the solution to be computed.  They must be defined
!         as variables because the code may change them.  You
!         have two choices --
!               Both RTOL and ATOL are scalars. (INFO(2)=0)
!               Both RTOL and ATOL are vectors. (INFO(2)=1)
!         in either case all components must be non-negative.
!
!         The tolerances are used by the code in a local error
!         test at each step which requires roughly that
!               ABS(LOCAL ERROR) .LE. RTOL*ABS(Y)+ATOL
!         for each vector component.
!         (More specifically, a root-mean-square norm is used to
!         measure the size of vectors, and the error test uses the
!         magnitude of the solution at the beginning of the step.)
!
!         The true (global) error is the difference between the
!         true solution of the initial value problem and the
!         computed approximation.  Practically all present day
!         codes, including this one, control the local error at
!         each step and do not even attempt to control the global
!         error directly.
!         Usually, but not always, the true accuracy of the
!         computed Y is comparable to the error tolerances. This
!         code will usually, but not always, deliver a more
!         accurate solution if you reduce the tolerances and
!         integrate again.  By comparing two such solutions you
!         can get a fairly reliable idea of the true error in the
!         solution at the bigger tolerances.
!
!         Setting ATOL=0. results in a pure relative error test on
!         that component.  Setting RTOL=0. results in a pure
!         absolute error test on that component.  A mixed test
!         with non-zero RTOL and ATOL corresponds roughly to a
!         relative error test when the solution component is much
!         bigger than ATOL and to an absolute error test when the
!         solution component is smaller than the threshhold ATOL.
!
!         The code will not attempt to compute a solution at an
!         accuracy unreasonable for the machine being used.  It will
!         advise you if you ask for too much accuracy and inform
!         you as to the maximum accuracy it believes possible.
!
!  RWORK(*) --  Dimension this real work array of length LRW in your
!         calling program.
!
!  LRW -- Set it to the declared length of the RWORK array.
!               You must have
!                    LRW .GE. 40+(MAXORD+4)*NEQ+NEQ**2
!               for the full (dense) JACOBIAN case (when INFO(6)=0), or
!                    LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
!               for the banded user-defined JACOBIAN case
!               (when INFO(5)=1 and INFO(6)=1), or
!                     LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
!                           +2*(NEQ/(ML+MU+1)+1)
!               for the banded finite-difference-generated JACOBIAN case
!               (when INFO(5)=0 and INFO(6)=1)
!
!  IWORK(*) --  Dimension this integer work array of length LIW in
!         your calling program.
!
!  LIW -- Set it to the declared length of the IWORK array.
!               You must have LIW .GE. 20+NEQ
!
!  RPAR, IPAR -- These are parameter arrays, of real and integer
!         type, respectively.  You can use them for communication
!         between your program that calls DDASSL and the
!         RES subroutine (and the JAC subroutine).  They are not
!         altered by DDASSL.  If you do not need RPAR or IPAR,
!         ignore these parameters by treating them as dummy
!         arguments.  If you do choose to use them, dimension
!         them in your calling program and in RES (and in JAC)
!         as arrays of appropriate length.
!
!  JAC -- If you have set INFO(5)=0, you can ignore this parameter
!         by treating it as a dummy argument.  Otherwise, you must
!         provide a subroutine of the form
!             SUBROUTINE JAC(T,Y,YPRIME,PD,CJ,RPAR,IPAR)
!         to define the matrix of partial derivatives
!             PD=DG/DY+CJ*DG/DYPRIME
!         CJ is a scalar which is input to JAC.
!         For the given values of T,Y,YPRIME, the
!         subroutine must evaluate the non-zero partial
!         derivatives for each equation and each solution
!         component, and store these values in the
!         matrix PD.  The elements of PD are set to zero
!         before each call to JAC so only non-zero elements
!         need to be defined.
!
!         Subroutine JAC must not alter T,Y,(*),YPRIME(*), or CJ.
!         You must declare the name JAC in an EXTERNAL statement in
!         your program that calls DDASSL.  You must dimension Y,
!         YPRIME and PD in JAC.
!
!         The way you must store the elements into the PD matrix
!         depends on the structure of the matrix which you
!         indicated by INFO(6).
!               *** INFO(6)=0 -- Full (dense) matrix ***
!                   Give PD a first dimension of NEQ.
!                   When you evaluate the (non-zero) partial derivative
!                   of equation I with respect to variable J, you must
!                   store it in PD according to
!                   PD(I,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
!               *** INFO(6)=1 -- Banded JACOBIAN with ML lower and MU
!                   upper diagonal bands (refer to INFO(6) description
!                   of ML and MU) ***
!                   Give PD a first dimension of 2*ML+MU+1.
!                   when you evaluate the (non-zero) partial derivative
!                   of equation I with respect to variable J, you must
!                   store it in PD according to
!                   IROW = I - J + ML + MU + 1
!                   PD(IROW,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
!
!         RPAR and IPAR are real and integer parameter arrays
!         which you can use for communication between your calling
!         program and your JACOBIAN subroutine JAC. They are not
!         altered by DDASSL. If you do not need RPAR or IPAR,
!         ignore these parameters by treating them as dummy
!         arguments. If you do choose to use them, dimension
!         them in your calling program and in JAC as arrays of
!         appropriate length.
!
!
!  OPTIONALLY REPLACEABLE NORM ROUTINE:
!
!     DDASSL uses a weighted norm DDANRM to measure the size
!     of vectors such as the estimated error in each step.
!     A FUNCTION subprogram
!       REAL(8) FUNCTION DDANRM(NEQ,V,WT,RPAR,IPAR)
!       DIMENSION V(NEQ),WT(NEQ)
!     is used to define this norm. Here, V is the vector
!     whose norm is to be computed, and WT is a vector of
!     weights.  A DDANRM routine has been included with DDASSL
!     which computes the weighted root-mean-square norm
!     given by
!       DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
!     this norm is suitable for most problems. In some
!     special cases, it may be more convenient and/or
!     efficient to define your own norm by writing a function
!     subprogram to be called instead of DDANRM. This should,
!     however, be attempted only after careful thought and
!     consideration.
!
!
!  -------- OUTPUT -- AFTER ANY RETURN FROM DDASSL ---------------------
!
!  The principal aim of the code is to return a computed solution at
!  TOUT, although it is also possible to obtain intermediate results
!  along the way. To find out whether the code achieved its goal
!  or if the integration process was interrupted before the task was
!  completed, you must check the IDID parameter.
!
!
!  T -- The solution was successfully advanced to the
!               output value of T.
!
!  Y(*) -- Contains the computed solution approximation at T.
!
!  YPRIME(*) -- Contains the computed derivative
!               approximation at T.
!
!  IDID -- Reports what the code did.
!
!                     *** Task completed ***
!                Reported by positive values of IDID
!
!           IDID = 1 -- A step was successfully taken in the
!                   intermediate-output mode. The code has not
!                   yet reached TOUT.
!
!           IDID = 2 -- The integration to TSTOP was successfully
!                   completed (T=TSTOP) by stepping exactly to TSTOP.
!
!           IDID = 3 -- The integration to TOUT was successfully
!                   completed (T=TOUT) by stepping past TOUT.
!                   Y(*) is obtained by interpolation.
!                   YPRIME(*) is obtained by interpolation.
!
!                    *** Task interrupted ***
!                Reported by negative values of IDID
!
!           IDID = -1 -- A large amount of work has been expended.
!                   (About 500 steps)
!
!           IDID = -2 -- The error tolerances are too stringent.
!
!           IDID = -3 -- The local error test cannot be satisfied
!                   because you specified a zero component in ATOL
!                   and the corresponding computed solution
!                   component is zero. Thus, a pure relative error
!                   test is impossible for this component.
!
!           IDID = -6 -- DDASSL had repeated error test
!                   failures on the last attempted step.
!
!           IDID = -7 -- The corrector could not converge.
!
!           IDID = -8 -- The matrix of partial derivatives
!                   is singular.
!
!           IDID = -9 -- The corrector could not converge.
!                   there were repeated error test failures
!                   in this step.
!
!           IDID =-10 -- The corrector could not converge
!                   because IRES was equal to minus one.
!
!           IDID =-11 -- IRES equal to -2 was encountered
!                   and control is being returned to the
!                   calling program.
!
!           IDID =-12 -- DDASSL failed to compute the initial
!                   YPRIME.
!
!
!
!           IDID = -13,..,-32 -- Not applicable for this code
!
!                    *** Task terminated ***
!                Reported by the value of IDID=-33
!
!           IDID = -33 -- The code has encountered trouble from which
!                   it cannot recover. A message is printed
!                   explaining the trouble and control is returned
!                   to the calling program. For example, this occurs
!                   when invalid input is detected.
!
!  RTOL, ATOL -- These quantities remain unchanged except when
!               IDID = -2. In this case, the error tolerances have been
!               increased by the code to values which are estimated to
!               be appropriate for continuing the integration. However,
!               the reported solution at T was obtained using the input
!               values of RTOL and ATOL.
!
!  RWORK, IWORK -- Contain information which is usually of no
!               interest to the user but necessary for subsequent calls.
!               However, you may find use for
!
!               RWORK(3)--Which contains the step size H to be
!                       attempted on the next step.
!
!               RWORK(4)--Which contains the current value of the
!                       independent variable, i.e., the farthest point
!                       integration has reached. This will be different
!                       from T only when interpolation has been
!                       performed (IDID=3).
!
!               RWORK(7)--Which contains the stepsize used
!                       on the last successful step.
!
!               IWORK(7)--Which contains the order of the method to
!                       be attempted on the next step.
!
!               IWORK(8)--Which contains the order of the method used
!                       on the last step.
!
!               IWORK(11)--Which contains the number of steps taken so
!                        far.
!
!               IWORK(12)--Which contains the number of calls to RES
!                        so far.
!
!               IWORK(13)--Which contains the number of evaluations of
!                        the matrix of partial derivatives needed so
!                        far.
!
!               IWORK(14)--Which contains the total number
!                        of error test failures so far.
!
!               IWORK(15)--Which contains the total number
!                        of convergence test failures so far.
!                        (includes singular iteration matrix
!                        failures.)
!
!
!  -------- INPUT -- WHAT TO DO TO CONTINUE THE INTEGRATION ------------
!                    (CALLS AFTER THE FIRST)
!
!  This code is organized so that subsequent calls to continue the
!  integration involve little (if any) additional effort on your
!  part. You must monitor the IDID parameter in order to determine
!  what to do next.
!
!  Recalling that the principal task of the code is to integrate
!  from T to TOUT (the interval mode), usually all you will need
!  to do is specify a new TOUT upon reaching the current TOUT.
!
!  Do not alter any quantity not specifically permitted below,
!  in particular do not alter NEQ,T,Y(*),YPRIME(*),RWORK(*),IWORK(*)
!  or the differential equation in subroutine RES. Any such
!  alteration constitutes a new problem and must be treated as such,
!  i.e., you must start afresh.
!
!  You cannot change from vector to scalar error control or vice
!  versa (INFO(2)), but you can change the size of the entries of
!  RTOL, ATOL. Increasing a tolerance makes the equation easier
!  to integrate. Decreasing a tolerance will make the equation
!  harder to integrate and should generally be avoided.
!
!  You can switch from the intermediate-output mode to the
!  interval mode (INFO(3)) or vice versa at any time.
!
!  If it has been necessary to prevent the integration from going
!  past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
!  code will not integrate to any TOUT beyond the currently
!  specified TSTOP. Once TSTOP has been reached you must change
!  the value of TSTOP or set INFO(4)=0. You may change INFO(4)
!  or TSTOP at any time but you must supply the value of TSTOP in
!  RWORK(1) whenever you set INFO(4)=1.
!
!  Do not change INFO(5), INFO(6), IWORK(1), or IWORK(2)
!  unless you are going to restart the code.
!
!                 *** Following a completed task ***
!  If
!     IDID = 1, call the code again to continue the integration
!                  another step in the direction of TOUT.
!
!     IDID = 2 or 3, define a new TOUT and call the code again.
!                  TOUT must be different from T. You cannot change
!                  the direction of integration without restarting.
!
!                 *** Following an interrupted task ***
!               To show the code that you realize the task was
!               interrupted and that you want to continue, you
!               must take appropriate action and set INFO(1) = 1
!  If
!    IDID = -1, The code has taken about 500 steps.
!                  If you want to continue, set INFO(1) = 1 and
!                  call the code again. An additional 500 steps
!                  will be allowed.
!
!    IDID = -2, The error tolerances RTOL, ATOL have been
!                  increased to values the code estimates appropriate
!                  for continuing. You may want to change them
!                  yourself. If you are sure you want to continue
!                  with relaxed error tolerances, set INFO(1)=1 and
!                  call the code again.
!
!    IDID = -3, A solution component is zero and you set the
!                  corresponding component of ATOL to zero. If you
!                  are sure you want to continue, you must first
!                  alter the error criterion to use positive values
!                  for those components of ATOL corresponding to zero
!                  solution components, then set INFO(1)=1 and call
!                  the code again.
!
!    IDID = -4,-5  --- Cannot occur with this code.
!
!    IDID = -6, Repeated error test failures occurred on the
!                  last attempted step in DDASSL. A singularity in the
!                  solution may be present. If you are absolutely
!                  certain you want to continue, you should restart
!                  the integration. (Provide initial values of Y and
!                  YPRIME which are consistent)
!
!    IDID = -7, Repeated convergence test failures occurred
!                  on the last attempted step in DDASSL. An inaccurate
!                  or ill-conditioned JACOBIAN may be the problem. If
!                  you are absolutely certain you want to continue, you
!                  should restart the integration.
!
!    IDID = -8, The matrix of partial derivatives is singular.
!                  Some of your equations may be redundant.
!                  DDASSL cannot solve the problem as stated.
!                  It is possible that the redundant equations
!                  could be removed, and then DDASSL could
!                  solve the problem. It is also possible
!                  that a solution to your problem either
!                  does not exist or is not unique.
!
!    IDID = -9, DDASSL had multiple convergence test
!                  failures, preceded by multiple error
!                  test failures, on the last attempted step.
!                  It is possible that your problem
!                  is ill-posed, and cannot be solved
!                  using this code. Or, there may be a
!                  discontinuity or a singularity in the
!                  solution. If you are absolutely certain
!                  you want to continue, you should restart
!                  the integration.
!
!    IDID =-10, DDASSL had multiple convergence test failures
!                  because IRES was equal to minus one.
!                  If you are absolutely certain you want
!                  to continue, you should restart the
!                  integration.
!
!    IDID =-11, IRES=-2 was encountered, and control is being
!                  returned to the calling program.
!
!    IDID =-12, DDASSL failed to compute the initial YPRIME.
!                  This could happen because the initial
!                  approximation to YPRIME was not very good, or
!                  if a YPRIME consistent with the initial Y
!                  does not exist. The problem could also be caused
!                  by an inaccurate or singular iteration matrix.
!
!    IDID = -13,..,-32  --- Cannot occur with this code.
!
!
!                 *** Following a terminated task ***
!
!  If IDID= -33, you cannot continue the solution of this problem.
!                  An attempt to do so will result in your
!                  run being terminated.
!
!
!  -------- ERROR MESSAGES ---------------------------------------------
!
!      The SLATEC error print routine XERMSG is called in the event of
!   unsuccessful completion of a task.  Most of these are treated as
!   "recoverable errors", which means that (unless the user has directed
!   otherwise) control will be returned to the calling program for
!   possible action after the message has been printed.
!
!   In the event of a negative value of IDID other than -33, an appro-
!   priate message is printed and the "error number" printed by XERMSG
!   is the value of IDID.  There are quite a number of illegal input
!   errors that can lead to a returned value IDID=-33.  The conditions
!   and their printed "error numbers" are as follows:
!
!   Error number       Condition
!
      SUBROUTINE DDASSL (problem, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,&
     &   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR)
!        1       Some element of INFO vector is not zero or one.
!        2       NEQ .le. 0
!        3       MAXORD not in range.
!        4       LRW is less than the required length for RWORK.
!        5       LIW is less than the required length for IWORK.
!        6       Some element of RTOL is .lt. 0
!        7       Some element of ATOL is .lt. 0
!        8       All elements of RTOL and ATOL are zero.
!        9       INFO(4)=1 and TSTOP is behind TOUT.
!       10       HMAX .lt. 0.0
!       11       TOUT is behind T.
!       12       INFO(8)=1 and H0=0.0
!       13       Some element of WT is .le. 0.0
!       14       TOUT is too close to T to start integration.
!       15       INFO(4)=1 and TSTOP is behind T.
!       16       --( Not used in this version )--
!       17       ML illegal.  Either .lt. 0 or .gt. NEQ
!       18       MU illegal.  Either .lt. 0 or .gt. NEQ
!       19       TOUT = T.
!
!   If DDASSL is called again without any action taken to remove the
!   cause of an unsuccessful return, XERMSG will be called with a fatal
!   error flag, which will cause unconditional termination of the
!   program.  There are two such fatal errors:
!
!   Error number -998:  The last step was terminated with a negative
!       value of IDID other than -33, and no appropriate action was
!       taken.
!
!   Error number -999:  The previous call was terminated because of
!       illegal input (IDID=-33) and there is illegal input in the
!       present call, as well.  (Suspect infinite loop.)
!
!  ---------------------------------------------------------------------
!
!***REFERENCES  A DESCRIPTION OF DASSL: A DIFFERENTIAL/ALGEBRAIC
!                 SYSTEM SOLVER, L. R. PETZOLD, SAND82-8637,
!                 SANDIA NATIONAL LABORATORIES, SEPTEMBER 1982.
!***ROUTINES CALLED  D1MACH, DDAINI, DDANRM, DDASTP, DDATRP, DDAWTS,
!                    XERMSG
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   880387  Code changes made.  All common statements have been
!           replaced by a DATA statement, which defines pointers into
!           RWORK, and PARAMETER statements which define pointers
!           into IWORK.  As well the documentation has gone through
!           grammatical changes.
!   881005  The prologue has been changed to mixed case.
!           The subordinate routines had revision dates changed to
!           this date, although the documentation for these routines
!           is all upper case.  No code changes.
!   890511  Code changes made.  The DATA statement in the declaration
!           section of DDASSL was replaced with a PARAMETER
!           statement.  Also the statement S = 100.D0 was removed
!           from the top of the Newton iteration in DDASTP.
!           The subordinate routines had revision dates changed to
!           this date.
!   890517  The revision date syntax was replaced with the revision
!           history syntax.  Also the "DECK" comment was added to
!           the top of all subroutines.  These changes are consistent
!           with new SLATEC guidelines.
!           The subordinate routines had revision dates changed to
!           this date.  No code changes.
!   891013  Code changes made.
!           Removed all occurrences of FLOAT or DBLE.  All operations
!           are now performed with "mixed-mode" arithmetic.
!           Also, specific function names were replaced with generic
!           function names to be consistent with new SLATEC guidelines.
!           In particular:
!              Replaced DSQRT with SQRT everywhere.
!              Replaced DABS with ABS everywhere.
!              Replaced DMIN1 with MIN everywhere.
!              Replaced MIN0 with MIN everywhere.
!              Replaced DMAX1 with MAX everywhere.
!              Replaced MAX0 with MAX everywhere.
!              Replaced DSIGN with SIGN everywhere.
!           Also replaced REVISION DATE with REVISION HISTORY in all
!           subordinate routines.
!   901004  Miscellaneous changes to prologue to complete conversion
!           to SLATEC 4.0 format.  No code changes.  (F.N.Fritsch)
!   901009  Corrected GAMS classification code and converted subsidiary
!           routines to 4.0 format.  No code changes.  (F.N.Fritsch)
!   901010  Converted XERRWV calls to XERMSG calls.  (R.Clemens, AFWL)
!   901019  Code changes made.
!           Merged SLATEC 4.0 changes with previous changes made
!           by C. Ulrich.  Below is a history of the changes made by
!           C. Ulrich. (Changes in subsidiary routines are implied
!           by this history)
!           891228  Bug was found and repaired inside the DDASSL
!                   and DDAINI routines.  DDAINI was incorrectly
!                   returning the initial T with Y and YPRIME
!                   computed at T+H.  The routine now returns T+H
!                   rather than the initial T.
!                   Cosmetic changes made to DDASTP.
!           900904  Three modifications were made to fix a bug (inside
!                   DDASSL) re interpolation for continuation calls and
!                   cases where TN is very close to TSTOP:
!
!                   1) In testing for whether H is too large, just
!                      compare H to (TSTOP - TN), rather than
!                      (TSTOP - TN) * (1-4*UROUND), and set H to
!                      TSTOP - TN.  This will force DDASTP to step
!                      exactly to TSTOP under certain situations
!                      (i.e. when H returned from DDASTP would otherwise
!                      take TN beyond TSTOP).
!
!                   2) Inside the DDASTP loop, interpolate exactly to
!                      TSTOP if TN is very close to TSTOP (rather than
!                      interpolating to within roundoff of TSTOP).
!
!                   3) Modified IDID description for IDID = 2 to say
!                      that the solution is returned by stepping exactly
!                      to TSTOP, rather than TOUT.  (In some cases the
!                      solution is actually obtained by extrapolating
!                      over a distance near unit roundoff to TSTOP,
!                      but this small distance is deemed acceptable in
!                      these circumstances.)
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue, removed unreferenced labels,
!           and improved XERMSG calls.  (FNF)
!   901030  Added ERROR MESSAGES section and reworked other sections to
!           be of more uniform format.  (FNF)
!   910624  Fixed minor bug related to HMAX (six lines after label
!           525).  (LRP)
!   000711  Fixed tests on (TN - TOUT) at 420 and 440 (ACH)
!***END PROLOGUE  DDASSL
!
!**End
!
!     Declare arguments.
!
      INTEGER  NEQ, INFO(15), IDID, LRW, IWORK(*), LIW, IPAR(*)
      REAL(8)                                                  &
     &   T, Y(*), YPRIME(*), TOUT, RTOL(*), ATOL(*), RWORK(*),          &
     &   RPAR(*)
      class(problem_ddassl), target :: problem
!
!     Declare externals.
!
      EXTERNAL  D1MACH
      REAL(8)  D1MACH
!
!     Declare local variables.
!
      INTEGER  I, ITEMP, LALPHA, LBETA, LCJ, LCJOLD, LCTF, LDELTA,      &
     &   LENIW, LENPD, LENRW, LE, LETF, LGAMMA, LH, LHMAX, LHOLD, LIPVT,&
     &   LJCALC, LK, LKOLD, LIWM, LML, LMTYPE, LMU, LMXORD, LNJE, LNPD, &
     &   LNRE, LNS, LNST, LNSTL, LPD, LPHASE, LPHI, LPSI, LROUND, LS,   &
     &   LSIGMA, LTN, LTSTOP, LWM, LWT, MBAND, MSAVE, MXORD, NPD, NTEMP,&
     &   NZFLG
      REAL(8)                                                  &
     &   ATOLI, H, HMAX, HMIN, HO, R, RH, RTOLI, TDIST, TN, TNEXT,      &
     &   TSTOP, UROUND, YPNORM
      LOGICAL  DONE
!       Auxiliary variables for conversion of values to be included in
!       error messages.
      CHARACTER(len=8)  XERN1, XERN2
      CHARACTER(len=16) XERN3, XERN4
!
!     SET POINTERS INTO IWORK
      PARAMETER (LML=1, LMU=2, LMXORD=3, LMTYPE=4, LNST=11,             &
     &  LNRE=12, LNJE=13, LETF=14, LCTF=15, LNPD=16,                    &
     &  LIPVT=21, LJCALC=5, LPHASE=6, LK=7, LKOLD=8,                    &
     &  LNS=9, LNSTL=10, LIWM=1)
!
!     SET RELATIVE OFFSET INTO RWORK
      PARAMETER (NPD=1)
!
!     SET POINTERS INTO RWORK
      PARAMETER (LTSTOP=1, LHMAX=2, LH=3, LTN=4,                        &
     &  LCJ=5, LCJOLD=6, LHOLD=7, LS=8, LROUND=9,                       &
     &  LALPHA=11, LBETA=17, LGAMMA=23,                                 &
     &  LPSI=29, LSIGMA=35, LDELTA=41)

!
!***FIRST EXECUTABLE STATEMENT  DDASSL
      IF(INFO(1).NE.0)GO TO 100
!
!-----------------------------------------------------------------------
!     THIS BLOCK IS EXECUTED FOR THE INITIAL CALL ONLY.
!     IT CONTAINS CHECKING OF INPUTS AND INITIALIZATIONS.
!-----------------------------------------------------------------------
!
!     FIRST CHECK INFO ARRAY TO MAKE SURE ALL ELEMENTS OF INFO
!     ARE EITHER ZERO OR ONE.
      DO 10 I=2,11
         IF(INFO(I).NE.0.AND.INFO(I).NE.1)GO TO 701
   10    CONTINUE
!
      IF(NEQ.LE.0)GO TO 702
!
!     CHECK AND COMPUTE MAXIMUM ORDER
      MXORD=5
      IF(INFO(9).EQ.0)GO TO 20
         MXORD=IWORK(LMXORD)
         IF(MXORD.LT.1.OR.MXORD.GT.5)GO TO 703
   20    IWORK(LMXORD)=MXORD
!
!     COMPUTE MTYPE,LENPD,LENRW.CHECK ML AND MU.
      IF(INFO(6).NE.0)GO TO 40
         LENPD=NEQ**2
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
         IF(INFO(5).NE.0)GO TO 30
            IWORK(LMTYPE)=2
            GO TO 60
   30       IWORK(LMTYPE)=1
            GO TO 60
   40 IF(IWORK(LML).LT.0.OR.IWORK(LML).GE.NEQ)GO TO 717
      IF(IWORK(LMU).LT.0.OR.IWORK(LMU).GE.NEQ)GO TO 718
      LENPD=(2*IWORK(LML)+IWORK(LMU)+1)*NEQ
      IF(INFO(5).NE.0)GO TO 50
         IWORK(LMTYPE)=5
         MBAND=IWORK(LML)+IWORK(LMU)+1
         MSAVE=(NEQ/MBAND)+1
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD+2*MSAVE
         GO TO 60
   50    IWORK(LMTYPE)=4
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
!
!     CHECK LENGTHS OF RWORK AND IWORK
   60 LENIW=20+NEQ
      IWORK(LNPD)=LENPD
      IF(LRW.LT.LENRW)GO TO 704
      IF(LIW.LT.LENIW)GO TO 705
!
!     CHECK TO SEE THAT TOUT IS DIFFERENT FROM T
      IF(TOUT .EQ. T)GO TO 719
!
!     CHECK HMAX
      IF(INFO(7).EQ.0)GO TO 70
         HMAX=RWORK(LHMAX)
         IF(HMAX.LE.0.0D0)GO TO 710
   70 CONTINUE
!
!     INITIALIZE COUNTERS
      IWORK(LNST)=0
      IWORK(LNRE)=0
      IWORK(LNJE)=0
!
      IWORK(LNSTL)=0
      IDID=1
      GO TO 200
!
!-----------------------------------------------------------------------
!     THIS BLOCK IS FOR CONTINUATION CALLS
!     ONLY. HERE WE CHECK INFO(1), AND IF THE
!     LAST STEP WAS INTERRUPTED WE CHECK WHETHER
!     APPROPRIATE ACTION WAS TAKEN.
!-----------------------------------------------------------------------
!
  100 CONTINUE
      IF(INFO(1).EQ.1)GO TO 110
      IF(INFO(1).NE.-1)GO TO 701
!
!     IF WE ARE HERE, THE LAST STEP WAS INTERRUPTED
!     BY AN ERROR CONDITION FROM DDASTP, AND
!     APPROPRIATE ACTION WAS NOT TAKEN. THIS
!     IS A FATAL ERROR.
      WRITE (XERN1, '(I8)') IDID
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'THE LAST STEP TERMINATED WITH A NEGATIVE VALUE OF IDID = ' // &
     &   XERN1 // ' AND NO APPROPRIATE ACTION WAS TAKEN.  ' //          &
     &   'RUN TERMINATED', -998, 2)
      RETURN
  110 CONTINUE
      IWORK(LNSTL)=IWORK(LNST)
!
!-----------------------------------------------------------------------
!     THIS BLOCK IS EXECUTED ON ALL CALLS.
!     THE ERROR TOLERANCE PARAMETERS ARE
!     CHECKED, AND THE WORK ARRAY POINTERS
!     ARE SET.
!-----------------------------------------------------------------------
!
  200 CONTINUE
!     CHECK RTOL,ATOL
      NZFLG=0
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 210 I=1,NEQ
         IF(INFO(2).EQ.1)RTOLI=RTOL(I)
         IF(INFO(2).EQ.1)ATOLI=ATOL(I)
         IF(RTOLI.GT.0.0D0.OR.ATOLI.GT.0.0D0)NZFLG=1
         IF(RTOLI.LT.0.0D0)GO TO 706
         IF(ATOLI.LT.0.0D0)GO TO 707
  210    CONTINUE
      IF(NZFLG.EQ.0)GO TO 708
!
!     SET UP RWORK STORAGE.IWORK STORAGE IS FIXED
!     IN DATA STATEMENT.
      LE=LDELTA+NEQ
      LWT=LE+NEQ
      LPHI=LWT+NEQ
      LPD=LPHI+(IWORK(LMXORD)+1)*NEQ
      LWM=LPD
      NTEMP=NPD+IWORK(LNPD)
      IF(INFO(1).EQ.1)GO TO 400
!
!-----------------------------------------------------------------------
!     THIS BLOCK IS EXECUTED ON THE INITIAL CALL
!     ONLY. SET THE INITIAL STEP SIZE, AND
!     THE ERROR WEIGHT VECTOR, AND PHI.
!     COMPUTE INITIAL YPRIME, IF NECESSARY.
!-----------------------------------------------------------------------
!
      TN=T
      IDID=1
!
!     SET ERROR WEIGHT VECTOR WT
      CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      DO 305 I = 1,NEQ
         IF(RWORK(LWT+I-1).LE.0.0D0) GO TO 713
  305    CONTINUE
!
!     COMPUTE UNIT ROUNDOFF AND HMIN
      UROUND = D1MACH(4)
      RWORK(LROUND) = UROUND
      HMIN = 4.0D0*UROUND*MAX(ABS(T),ABS(TOUT))
!
!     CHECK INITIAL INTERVAL TO SEE THAT IT IS LONG ENOUGH
      TDIST = ABS(TOUT - T)
      IF(TDIST .LT. HMIN) GO TO 714
!
!     CHECK HO, IF THIS WAS INPUT
      IF (INFO(8) .EQ. 0) GO TO 310
         HO = RWORK(LH)
         IF ((TOUT - T)*HO .LT. 0.0D0) GO TO 711
         IF (HO .EQ. 0.0D0) GO TO 712
         GO TO 320
  310  CONTINUE
!
!     COMPUTE INITIAL STEPSIZE, TO BE USED BY EITHER
!     DDASTP OR DDAINI, DEPENDING ON INFO(11)
      HO = 0.001D0*TDIST
      YPNORM = DDANRM(NEQ,YPRIME,RWORK(LWT),RPAR,IPAR)
      IF (YPNORM .GT. 0.5D0/HO) HO = 0.5D0/YPNORM
      HO = SIGN(HO,TOUT-T)
!     ADJUST HO IF NECESSARY TO MEET HMAX BOUND
  320 IF (INFO(7) .EQ. 0) GO TO 330
         RH = ABS(HO)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) HO = HO/RH
!     COMPUTE TSTOP, IF APPLICABLE
  330 IF (INFO(4) .EQ. 0) GO TO 340
         TSTOP = RWORK(LTSTOP)
         IF ((TSTOP - T)*HO .LT. 0.0D0) GO TO 715
         IF ((T + HO - TSTOP)*HO .GT. 0.0D0) HO = TSTOP - T
         IF ((TSTOP - TOUT)*HO .LT. 0.0D0) GO TO 709
!
!     COMPUTE INITIAL DERIVATIVE, UPDATING TN AND Y, IF APPLICABLE
  340 IF (INFO(11) .EQ. 0) GO TO 350
      CALL DDAINI(problem,TN,Y,YPRIME,NEQ,                                  &
     &  HO,RWORK(LWT),IDID,RPAR,IPAR,                                   &
     &  RWORK(LPHI),RWORK(LDELTA),RWORK(LE),                            &
     &  RWORK(LWM),IWORK(LIWM),HMIN,RWORK(LROUND),                      &
     &  INFO(10),NTEMP)
      IF (IDID .LT. 0) GO TO 390
!
!     LOAD H WITH HO.  STORE H IN RWORK(LH)
  350 H = HO
      RWORK(LH) = H
!
!     LOAD Y AND H*YPRIME INTO PHI(*,1) AND PHI(*,2)
      ITEMP = LPHI + NEQ
      DO I = 1,NEQ
         RWORK(LPHI + I - 1) = Y(I)
         RWORK(ITEMP + I - 1) = H*YPRIME(I)
      end DO
!
  390 GO TO 500
!
!-------------------------------------------------------
!     THIS BLOCK IS FOR CONTINUATION CALLS ONLY. ITS
!     PURPOSE IS TO CHECK STOP CONDITIONS BEFORE
!     TAKING A STEP.
!     ADJUST H IF NECESSARY TO MEET HMAX BOUND
!-------------------------------------------------------
!
  400 CONTINUE
      UROUND=RWORK(LROUND)
      DONE = .FALSE.
      TN=RWORK(LTN)
      H=RWORK(LH)
      IF(INFO(7) .EQ. 0) GO TO 410
         RH = ABS(H)/RWORK(LHMAX)
         IF(RH .GT. 1.0D0) H = H/RH
  410 CONTINUE
      IF(T .EQ. TOUT) GO TO 719
      IF((T - TOUT)*H .GT. 0.0D0) GO TO 711
      IF(INFO(4) .EQ. 1) GO TO 430
      IF(INFO(3) .EQ. 1) GO TO 420
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 490
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),                    &
     &  RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
  420 IF((TN-T)*H .LE. 0.0D0) GO TO 490
      IF((TN - TOUT)*H .GE. 0.0D0) GO TO 425
      CALL DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),                      &
     &  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
  425 CONTINUE
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),                    &
     &  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
  430 IF(INFO(3) .EQ. 1) GO TO 440
      TSTOP=RWORK(LTSTOP)
      IF((TN-TSTOP)*H.GT.0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H.LT.0.0D0)GO TO 709
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 450
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),                    &
     &   RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
  440 TSTOP = RWORK(LTSTOP)
      IF((TN-TSTOP)*H .GT. 0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H .LT. 0.0D0) GO TO 709
      IF((TN-T)*H .LE. 0.0D0) GO TO 450
      IF((TN - TOUT)*H .GE. 0.0D0) GO TO 445
      CALL DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),                      &
     &  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
  445 CONTINUE
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),                    &
     &  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
  450 CONTINUE
!     CHECK WHETHER WE ARE WITHIN ROUNDOFF OF TSTOP
      IF(ABS(TN-TSTOP).GT.100.0D0*UROUND*                               &
     &   (ABS(TN)+ABS(H)))GO TO 460
      CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,IWORK(LKOLD),                   &
     &  RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      DONE = .TRUE.
      GO TO 490
  460 TNEXT=TN+H
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 490
      H=TSTOP-TN
      RWORK(LH)=H
!
  490 IF (DONE) GO TO 580
!
!-------------------------------------------------------
!     THE NEXT BLOCK CONTAINS THE CALL TO THE
!     ONE-STEP INTEGRATOR DDASTP.
!     THIS IS A LOOPING POINT FOR THE INTEGRATION STEPS.
!     CHECK FOR TOO MANY STEPS.
!     UPDATE WT.
!     CHECK FOR TOO MUCH ACCURACY REQUESTED.
!     COMPUTE MINIMUM STEPSIZE.
!-------------------------------------------------------
!
  500 CONTINUE
!     CHECK FOR FAILURE TO COMPUTE INITIAL YPRIME
      IF (IDID .EQ. -12) GO TO 527
!
!     CHECK FOR TOO MANY STEPS
      IF((IWORK(LNST)-IWORK(LNSTL)).LT.500)                             &
     &   GO TO 510
           IDID=-1
           GO TO 527
!
!     UPDATE WT
  510 CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,RWORK(LPHI),                    &
     &  RWORK(LWT),RPAR,IPAR)
      DO 520 I=1,NEQ
         IF(RWORK(I+LWT-1).GT.0.0D0)GO TO 520
           IDID=-3
           GO TO 527
  520 END DO
!
!     TEST FOR TOO MUCH ACCURACY REQUESTED.
      R=DDANRM(NEQ,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)*                   &
     &   100.0D0*UROUND
      IF(R.LE.1.0D0)GO TO 525
!     MULTIPLY RTOL AND ATOL BY R AND RETURN
      IF(INFO(2).EQ.1)GO TO 523
           RTOL(1)=R*RTOL(1)
           ATOL(1)=R*ATOL(1)
           IDID=-2
           GO TO 527
  523 DO  I=1,NEQ
           RTOL(I)=R*RTOL(I)
           ATOL(I)=R*ATOL(I)
      end DO

      IDID=-2
      GO TO 527
  525 CONTINUE
!
!     COMPUTE MINIMUM STEPSIZE
      HMIN=4.0D0*UROUND*MAX(ABS(TN),ABS(TOUT))
!
!     TEST H VS. HMAX
      IF (INFO(7) .NE. 0) THEN
         RH = ABS(H)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) H = H/RH
      ENDIF
!
      CALL DDASTP(problem,TN,Y,YPRIME,NEQ,                                  &
     &   H,RWORK(LWT),INFO(1),IDID,RPAR,IPAR,                           &
     &   RWORK(LPHI),RWORK(LDELTA),RWORK(LE),                           &
     &   RWORK(LWM),IWORK(LIWM),                                        &
     &   RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),                      &
     &   RWORK(LPSI),RWORK(LSIGMA),                                     &
     &   RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),                         &
     &   RWORK(LS),HMIN,RWORK(LROUND),                                  &
     &   IWORK(LPHASE),IWORK(LJCALC),IWORK(LK),                         &
     &   IWORK(LKOLD),IWORK(LNS),INFO(10),NTEMP)
  527 IF(IDID.LT.0)GO TO 600
!
!--------------------------------------------------------
!     THIS BLOCK HANDLES THE CASE OF A SUCCESSFUL RETURN
!     FROM DDASTP (IDID=1).  TEST FOR STOP CONDITIONS.
!--------------------------------------------------------
!
      IF(INFO(4).NE.0)GO TO 540
           IF(INFO(3).NE.0)GO TO 530
             IF((TN-TOUT)*H.LT.0.0D0)GO TO 500
             CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,                          &
     &         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
  530        IF((TN-TOUT)*H.GE.0.0D0)GO TO 535
             T=TN
             IDID=1
             GO TO 580
  535        CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,                          &
     &         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
  540 IF(INFO(3).NE.0)GO TO 550
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 542
         CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,                              &
     &     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
         T=TOUT
         IDID=3
         GO TO 580
  542 IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*                               &
     &   (ABS(TN)+ABS(H)))GO TO 545
      TNEXT=TN+H
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 500
      H=TSTOP-TN
      GO TO 500
  545 CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,                                &
     &  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      GO TO 580
  550 IF((TN-TOUT)*H.GE.0.0D0)GO TO 555
      IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*(ABS(TN)+ABS(H)))GO TO 552
      T=TN
      IDID=1
      GO TO 580
  552 CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,                                &
     &  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      GO TO 580
  555 CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,                                 &
     &   IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID=3
      GO TO 580
!
!--------------------------------------------------------
!     ALL SUCCESSFUL RETURNS FROM DDASSL ARE MADE FROM
!     THIS BLOCK.
!--------------------------------------------------------
!
  580 CONTINUE
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
!
!-----------------------------------------------------------------------
!     THIS BLOCK HANDLES ALL UNSUCCESSFUL
!     RETURNS OTHER THAN FOR ILLEGAL INPUT.
!-----------------------------------------------------------------------
!
  600 CONTINUE
      ITEMP=-IDID
      GO TO (610,620,630,690,690,640,650,660,670,675,                   &
     &  680,685), ITEMP
!
!     THE MAXIMUM NUMBER OF STEPS WAS TAKEN BEFORE
!     REACHING TOUT
  610 WRITE (XERN3, '(1P,D15.6)') TN
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'AT CURRENT T = ' // XERN3 // ' 500 STEPS TAKEN ON THIS ' //   &
     &   'CALL BEFORE REACHING TOUT', IDID, 1)
      GO TO 690
!
!     TOO MUCH ACCURACY FOR MACHINE PRECISION
  620 WRITE (XERN3, '(1P,D15.6)') TN
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'AT T = ' // XERN3 // ' TOO MUCH ACCURACY REQUESTED FOR ' //   &
     &   'PRECISION OF MACHINE. RTOL AND ATOL WERE INCREASED TO ' //    &
     &   'APPROPRIATE VALUES', IDID, 1)
      GO TO 690
!
!     WT(I) .LE. 0.0 FOR SOME I (NOT AT START OF PROBLEM)
  630 WRITE (XERN3, '(1P,D15.6)') TN
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'AT T = ' // XERN3 // ' SOME ELEMENT OF WT HAS BECOME .LE. ' //&
     &   '0.0', IDID, 1)
      GO TO 690
!
!     ERROR TEST FAILED REPEATEDLY OR WITH H=HMIN
  640 WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //         &
     &   ' THE ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN',       &
     &   IDID, 1)
      GO TO 690
!
!     CORRECTOR CONVERGENCE FAILED REPEATEDLY OR WITH H=HMIN
  650 WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //         &
     &   ' THE CORRECTOR FAILED TO CONVERGE REPEATEDLY OR WITH ' //     &
     &   'ABS(H)=HMIN', IDID, 1)
      GO TO 690
!
!     THE ITERATION MATRIX IS SINGULAR
  660 WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //         &
     &   ' THE ITERATION MATRIX IS SINGULAR', IDID, 1)
      GO TO 690
!
!     CORRECTOR FAILURE PRECEDED BY ERROR TEST FAILURES.
  670 WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //         &
     &   ' THE CORRECTOR COULD NOT CONVERGE.  ALSO, THE ERROR TEST ' // &
     &   'FAILED REPEATEDLY.', IDID, 1)
      GO TO 690
!
!     CORRECTOR FAILURE BECAUSE IRES = -1
  675 WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //         &
     &   ' THE CORRECTOR COULD NOT CONVERGE BECAUSE IRES WAS EQUAL ' // &
     &   'TO MINUS ONE', IDID, 1)
      GO TO 690
!
!     FAILURE BECAUSE IRES = -2
  680 WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //         &
     &   ' IRES WAS EQUAL TO MINUS TWO', IDID, 1)
      GO TO 690
!
!     FAILED TO COMPUTE INITIAL YPRIME
  685 WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') HO
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //         &
     &   ' THE INITIAL YPRIME COULD NOT BE COMPUTED', IDID, 1)
      GO TO 690
!
  690 CONTINUE
      INFO(1)=-1
      T=TN
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
!
!-----------------------------------------------------------------------
!     THIS BLOCK HANDLES ALL ERROR RETURNS DUE
!     TO ILLEGAL INPUT, AS DETECTED BEFORE CALLING
!     DDASTP. FIRST THE ERROR MESSAGE ROUTINE IS
!     CALLED. IF THIS HAPPENS TWICE IN
!     SUCCESSION, EXECUTION IS TERMINATED
!
!-----------------------------------------------------------------------
  701 CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR ONE', 1, 1)
      GO TO 750
!
  702 WRITE (XERN1, '(I8)') NEQ
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'NEQ = ' // XERN1 // ' .LE. 0', 2, 1)
      GO TO 750
!
  703 WRITE (XERN1, '(I8)') MXORD
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'MAXORD = ' // XERN1 // ' NOT IN RANGE', 3, 1)
      GO TO 750
!
  704 WRITE (XERN1, '(I8)') LENRW
      WRITE (XERN2, '(I8)') LRW
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'RWORK LENGTH NEEDED, LENRW = ' // XERN1 //                    &
     &   ', EXCEEDS LRW = ' // XERN2, 4, 1)
      GO TO 750
!
  705 WRITE (XERN1, '(I8)') LENIW
      WRITE (XERN2, '(I8)') LIW
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'IWORK LENGTH NEEDED, LENIW = ' // XERN1 //                    &
     &   ', EXCEEDS LIW = ' // XERN2, 5, 1)
      GO TO 750
!
  706 CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'SOME ELEMENT OF RTOL IS .LT. 0', 6, 1)
      GO TO 750
!
  707 CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'SOME ELEMENT OF ATOL IS .LT. 0', 7, 1)
      GO TO 750
!
  708 CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'ALL ELEMENTS OF RTOL AND ATOL ARE ZERO', 8, 1)
      GO TO 750
!
  709 WRITE (XERN3, '(1P,D15.6)') TSTOP
      WRITE (XERN4, '(1P,D15.6)') TOUT
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'INFO(4) = 1 AND TSTOP = ' // XERN3 // ' BEHIND TOUT = ' //    &
     &   XERN4, 9, 1)
      GO TO 750
!
  710 WRITE (XERN3, '(1P,D15.6)') HMAX
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'HMAX = ' // XERN3 // ' .LT. 0.0', 10, 1)
      GO TO 750
!
  711 WRITE (XERN3, '(1P,D15.6)') TOUT
      WRITE (XERN4, '(1P,D15.6)') T
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'TOUT = ' // XERN3 // ' BEHIND T = ' // XERN4, 11, 1)
      GO TO 750
!
  712 CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'INFO(8)=1 AND H0=0.0', 12, 1)
      GO TO 750
!
  713 CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'SOME ELEMENT OF WT IS .LE. 0.0', 13, 1)
      GO TO 750
!
  714 WRITE (XERN3, '(1P,D15.6)') TOUT
      WRITE (XERN4, '(1P,D15.6)') T
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'TOUT = ' // XERN3 // ' TOO CLOSE TO T = ' // XERN4 //         &
     &   ' TO START INTEGRATION', 14, 1)
      GO TO 750
!
  715 WRITE (XERN3, '(1P,D15.6)') TSTOP
      WRITE (XERN4, '(1P,D15.6)') T
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'INFO(4)=1 AND TSTOP = ' // XERN3 // ' BEHIND T = ' // XERN4,  &
     &   15, 1)
      GO TO 750
!
  717 WRITE (XERN1, '(I8)') IWORK(LML)
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'ML = ' // XERN1 // ' ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ',    &
     &   17, 1)
      GO TO 750
!
  718 WRITE (XERN1, '(I8)') IWORK(LMU)
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &   'MU = ' // XERN1 // ' ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ',    &
     &   18, 1)
      GO TO 750
!
  719 WRITE (XERN3, '(1P,D15.6)') TOUT
      CALL XERMSG ('SLATEC', 'DDASSL',                                  &
     &  'TOUT = T = ' // XERN3, 19, 1)
      GO TO 750
!
  750 IDID=-33
      IF(INFO(1).EQ.-1) THEN
         CALL XERMSG ('SLATEC', 'DDASSL',                               &
     &      'REPEATED OCCURRENCES OF ILLEGAL INPUT$$' //                &
     &      'RUN TERMINATED. APPARENT INFINITE LOOP', -999, 2)
      ENDIF
!
      INFO(1)=-1
      RETURN
!-----------END OF SUBROUTINE DDASSL------------------------------------
      END subroutine ddassl
!DECK DDAINI
      SUBROUTINE DDAINI (problem, X, Y, YPRIME, NEQ, H, WT, IDID, RPAR,&
     &   IPAR, PHI, DELTA, E, WM, IWM, HMIN, UROUND, NONNEG, NTEMP)
!***BEGIN PROLOGUE  DDAINI
!***SUBSIDIARY
!***PURPOSE  Initialization routine for DDASSL.
!***LIBRARY   SLATEC (DASSL)
!***TYPE      REAL(8) (SDAINI-S, DDAINI-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------
!     DDAINI TAKES ONE STEP OF SIZE H OR SMALLER
!     WITH THE BACKWARD EULER METHOD, TO
!     FIND YPRIME.  X AND Y ARE UPDATED TO BE CONSISTENT WITH THE
!     NEW STEP.  A MODIFIED DAMPED NEWTON ITERATION IS USED TO
!     SOLVE THE CORRECTOR ITERATION.
!
!     THE INITIAL GUESS FOR YPRIME IS USED IN THE
!     PREDICTION, AND IN FORMING THE ITERATION
!     MATRIX, BUT IS NOT INVOLVED IN THE
!     ERROR TEST. THIS MAY HAVE TROUBLE
!     CONVERGING IF THE INITIAL GUESS IS NO
!     GOOD, OR IF G(X,Y,YPRIME) DEPENDS
!     NONLINEARLY ON YPRIME.
!
!     THE PARAMETERS REPRESENT:
!     X --         INDEPENDENT VARIABLE
!     Y --         SOLUTION VECTOR AT X
!     YPRIME --    DERIVATIVE OF SOLUTION VECTOR
!     NEQ --       NUMBER OF EQUATIONS
!     H --         STEPSIZE. IMDER MAY USE A STEPSIZE
!                  SMALLER THAN H.
!     WT --        VECTOR OF WEIGHTS FOR ERROR
!                  CRITERION
!     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS
!                  IDID= 1 -- YPRIME WAS FOUND SUCCESSFULLY
!                  IDID=-12 -- DDAINI FAILED TO FIND YPRIME
!     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS
!                  THAT ARE NOT ALTERED BY DDAINI
!     PHI --       WORK SPACE FOR DDAINI
!     DELTA,E --   WORK SPACE FOR DDAINI
!     WM,IWM --    REAL AND INTEGER ARRAYS STORING
!                  MATRIX INFORMATION
!
!-----------------------------------------------------------------
!***ROUTINES CALLED  DDAJAC, DDANRM, DDASLV
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!   901030  Minor corrections to declarations.  (FNF)
!***END PROLOGUE  DDAINI
!
      INTEGER  NEQ, IDID, IPAR(*), IWM(*), NONNEG, NTEMP
      REAL(8)                                                  &
     &   X, Y(*), YPRIME(*), H, WT(*), RPAR(*), PHI(NEQ,*), DELTA(*),   &
     &   E(*), WM(*), HMIN, UROUND
      class(problem_ddassl), target :: problem
!
      INTEGER  I, IER, IRES, JCALC, LNJE, LNRE, M, MAXIT, MJAC, NCF,    &
     &   NEF, NSF
      REAL(8)                                                  &
     &   CJ, DAMP, DELNRM, ERR, OLDNRM, R, RATE, S, XOLD, YNORM
      LOGICAL  CONVGD
!
      PARAMETER (LNRE=12)
      PARAMETER (LNJE=13)
!
      DATA MAXIT/10/,MJAC/5/
      DATA DAMP/0.75D0/
!
!
!---------------------------------------------------
!     BLOCK 1.
!     INITIALIZATIONS.
!---------------------------------------------------
!
!***FIRST EXECUTABLE STATEMENT  DDAINI
      IDID=1
      NEF=0
      NCF=0
      NSF=0
      XOLD=X
      YNORM=DDANRM(NEQ,Y,WT,RPAR,IPAR)
!
!     SAVE Y AND YPRIME IN PHI
      DO I=1,NEQ
         PHI(I,1)=Y(I)
         PHI(I,2)=YPRIME(I)
      end DO

!
!
!----------------------------------------------------
!     BLOCK 2.
!     DO ONE BACKWARD EULER STEP.
!----------------------------------------------------
!
!     SET UP FOR START OF CORRECTOR ITERATION
  200 CJ=1.0D0/H
      X=X+H
!
!     PREDICT SOLUTION AND DERIVATIVE
      DO I=1,NEQ
         Y(I)=Y(I)+H*YPRIME(I)
      end DO
!
      JCALC=-1
      M=0
      CONVGD=.TRUE.
!
!
!     CORRECTOR LOOP.
  300 IWM(LNRE)=IWM(LNRE)+1
      IRES=0
!
      CALL problem%res(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES.LT.0) GO TO 430
!
!
!     EVALUATE THE ITERATION MATRIX
      IF (JCALC.NE.-1) GO TO 310
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      CALL DDAJAC(problem,NEQ,X,Y,YPRIME,DELTA,CJ,H,                        &
     &   IER,WT,E,WM,IWM,IRES,                                          &
     &   UROUND,RPAR,IPAR,NTEMP)
!
      S=1000000.D0
      IF (IRES.LT.0) GO TO 430
      IF (IER.NE.0) GO TO 430
      NSF=0
!
!
!
!     MULTIPLY RESIDUAL BY DAMPING FACTOR
  310 CONTINUE
      DO I=1,NEQ
         DELTA(I)=DELTA(I)*DAMP
      end DO
!
!     COMPUTE A NEW ITERATE (BACK SUBSTITUTION)
!     STORE THE CORRECTION IN DELTA
!
      CALL DDASLV(NEQ,DELTA,WM,IWM)
!
!     UPDATE Y AND YPRIME
      DO I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
         YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
      end DO
!
!     TEST FOR CONVERGENCE OF THE ITERATION.
!
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.LE.100.D0*UROUND*YNORM)                                &
     &   GO TO 400
!
      IF (M.GT.0) GO TO 340
         OLDNRM=DELNRM
         GO TO 350
!
  340 RATE=(DELNRM/OLDNRM)**(1.0D0/M)
      IF (RATE.GT.0.90D0) GO TO 430
      S=RATE/(1.0D0-RATE)
!
  350 IF (S*DELNRM .LE. 0.33D0) GO TO 400
!
!
!     THE CORRECTOR HAS NOT YET CONVERGED. UPDATE
!     M AND AND TEST WHETHER THE MAXIMUM
!     NUMBER OF ITERATIONS HAVE BEEN TRIED.
!     EVERY MJAC ITERATIONS, GET A NEW
!     ITERATION MATRIX.
!
      M=M+1
      IF (M.GE.MAXIT) GO TO 430
!
      IF ((M/MJAC)*MJAC.EQ.M) JCALC=-1
      GO TO 300
!
!
!     THE ITERATION HAS CONVERGED.
!     CHECK NONNEGATIVITY CONSTRAINTS
  400 IF (NONNEG.EQ.0) GO TO 450
      DO I=1,NEQ
         DELTA(I)=MIN(Y(I),0.0D0)
      end DO
!
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.GT.0.33D0) GO TO 430
!
      DO I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
         YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
      end DO
      GO TO 450
!
!
!     EXITS FROM CORRECTOR LOOP.
  430 CONVGD=.FALSE.
  450 IF (.NOT.CONVGD) GO TO 600
!
!
!
!-----------------------------------------------------
!     BLOCK 3.
!     THE CORRECTOR ITERATION CONVERGED.
!     DO ERROR TEST.
!-----------------------------------------------------
!
      DO I=1,NEQ
         E(I)=Y(I)-PHI(I,1)
      end DO
      ERR=DDANRM(NEQ,E,WT,RPAR,IPAR)
!
      IF (ERR.LE.1.0D0) RETURN
!
!
!
!--------------------------------------------------------
!     BLOCK 4.
!     THE BACKWARD EULER STEP FAILED. RESTORE X, Y
!     AND YPRIME TO THEIR ORIGINAL VALUES.
!     REDUCE STEPSIZE AND TRY AGAIN, IF
!     POSSIBLE.
!---------------------------------------------------------
!
  600 CONTINUE
      X = XOLD
      DO 610 I=1,NEQ
         Y(I)=PHI(I,1)
  610    YPRIME(I)=PHI(I,2)
!
      IF (CONVGD) GO TO 640
      IF (IER.EQ.0) GO TO 620
         NSF=NSF+1
         H=H*0.25D0
         IF (NSF.LT.3.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
  620 IF (IRES.GT.-2) GO TO 630
         IDID=-12
         RETURN
  630 NCF=NCF+1
      H=H*0.25D0
      IF (NCF.LT.10.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
!
  640 NEF=NEF+1
      R=0.90D0/(2.0D0*ERR+0.0001D0)
      R=MAX(0.1D0,MIN(0.5D0,R))
      H=H*R
      IF (ABS(H).GE.HMIN.AND.NEF.LT.10) GO TO 690
         IDID=-12
         RETURN
  690    GO TO 200
!
!-------------END OF SUBROUTINE DDAINI----------------------
      END subroutine ddaini
!DECK DDAJAC
      SUBROUTINE DDAJAC (problem, NEQ, X, Y, YPRIME, DELTA, CJ, H, &
           &IER, WT, E, WM, IWM, IRES, UROUND, RPAR, IPAR, NTEMP)
!***BEGIN PROLOGUE  DDAJAC
!***SUBSIDIARY
!***PURPOSE  Compute the iteration matrix for DDASSL and form the
!            LU-decomposition.
!***LIBRARY   SLATEC (DASSL)
!***TYPE      REAL(8) (SDAJAC-S, DDAJAC-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------------
!     THIS ROUTINE COMPUTES THE ITERATION MATRIX
!     PD=DG/DY+CJ*DG/DYPRIME (WHERE G(X,Y,YPRIME)=0).
!     HERE PD IS COMPUTED BY THE USER-SUPPLIED
!     ROUTINE JAC IF IWM(MTYPE) IS 1 OR 4, AND
!     IT IS COMPUTED BY NUMERICAL FINITE DIFFERENCING
!     IF IWM(MTYPE)IS 2 OR 5
!     THE PARAMETERS HAVE THE FOLLOWING MEANINGS.
!     Y        = ARRAY CONTAINING PREDICTED VALUES
!     YPRIME   = ARRAY CONTAINING PREDICTED DERIVATIVES
!     DELTA    = RESIDUAL EVALUATED AT (X,Y,YPRIME)
!                (USED ONLY IF IWM(MTYPE)=2 OR 5)
!     CJ       = SCALAR PARAMETER DEFINING ITERATION MATRIX
!     H        = CURRENT STEPSIZE IN INTEGRATION
!     IER      = VARIABLE WHICH IS .NE. 0
!                IF ITERATION MATRIX IS SINGULAR,
!                AND 0 OTHERWISE.
!     WT       = VECTOR OF WEIGHTS FOR COMPUTING NORMS
!     E        = WORK SPACE (TEMPORARY) OF LENGTH NEQ
!     WM       = REAL WORK SPACE FOR MATRICES. ON
!                OUTPUT IT CONTAINS THE LU DECOMPOSITION
!                OF THE ITERATION MATRIX.
!     IWM      = INTEGER WORK SPACE CONTAINING
!                MATRIX INFORMATION
!     RES      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
!                TO EVALUATE THE RESIDUAL FUNCTION G(X,Y,YPRIME)
!     IRES     = FLAG WHICH IS EQUAL TO ZERO IF NO ILLEGAL VALUES
!                IN RES, AND LESS THAN ZERO OTHERWISE.  (IF IRES
!                IS LESS THAN ZERO, THE MATRIX WAS NOT COMPLETED)
!                IN THIS CASE (IF IRES .LT. 0), THEN IER = 0.
!     UROUND   = THE UNIT ROUNDOFF ERROR OF THE MACHINE BEING USED.
!     JAC      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
!                TO EVALUATE THE ITERATION MATRIX (THIS ROUTINE
!                IS ONLY USED IF IWM(MTYPE) IS 1 OR 4)
!-----------------------------------------------------------------------
!***ROUTINES CALLED  DGBFA, DGEFA
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901010  Modified three MAX calls to be all on one line.  (FNF)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!   901101  Corrected PURPOSE.  (FNF)
!***END PROLOGUE  DDAJAC
!
      INTEGER  NEQ, IER, IWM(*), IRES, IPAR(*), NTEMP
      REAL(8)                                                  &
     &   X, Y(*), YPRIME(*), DELTA(*), CJ, H, WT(*), E(*), WM(*),       &
     &   UROUND, RPAR(*)
      class(problem_ddassl), target :: problem
!
      EXTERNAL  DGBFA, DGEFA
!
      INTEGER  I, I1, I2, II, IPSAVE, ISAVE, J, K, L, LENPD, LIPVT,     &
     &   LML, LMTYPE, LMU, MBA, MBAND, MEB1, MEBAND, MSAVE, MTYPE, N,   &
     &   NPD, NPDM1, NROW
      REAL(8)  DEL, DELINV, SQUR, YPSAVE, YSAVE
!
      PARAMETER (NPD=1)
      PARAMETER (LML=1)
      PARAMETER (LMU=2)
      PARAMETER (LMTYPE=4)
      PARAMETER (LIPVT=21)
!
!***FIRST EXECUTABLE STATEMENT  DDAJAC
      IER = 0
      NPDM1=NPD-1
      MTYPE=IWM(LMTYPE)
      GO TO (100,200,300,400,500),MTYPE
!
!
!     DENSE USER-SUPPLIED MATRIX
  100 LENPD=NEQ*NEQ
      DO 110 I=1,LENPD
  110    WM(NPDM1+I)=0.0D0
      CALL PROBLEM%JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      GO TO 230
!
!
!     DENSE FINITE-DIFFERENCE-GENERATED MATRIX
  200 IRES=0
      NROW=NPDM1
      SQUR = SQRT(UROUND)
      DO 210 I=1,NEQ
         DEL=SQUR*MAX(ABS(Y(I)),ABS(H*YPRIME(I)),ABS(WT(I)))
         DEL=SIGN(DEL,H*YPRIME(I))
         DEL=(Y(I)+DEL)-Y(I)
         YSAVE=Y(I)
         YPSAVE=YPRIME(I)
         Y(I)=Y(I)+DEL
         YPRIME(I)=YPRIME(I)+CJ*DEL
         CALL problem%res(X,Y,YPRIME,E,IRES,RPAR,IPAR)
         IF (IRES .LT. 0) RETURN
         DELINV=1.0D0/DEL
         DO 220 L=1,NEQ
  220    WM(NROW+L)=(E(L)-DELTA(L))*DELINV
      NROW=NROW+NEQ
      Y(I)=YSAVE
      YPRIME(I)=YPSAVE
  210 END DO
!
!
!     DO DENSE-MATRIX LU DECOMPOSITION ON PD
  230    CALL DGEFA(WM(NPD),NEQ,NEQ,IWM(LIPVT),IER)
      RETURN
!
!
!     DUMMY SECTION FOR IWM(MTYPE)=3
  300 RETURN
!
!
!     BANDED USER-SUPPLIED MATRIX
  400 LENPD=(2*IWM(LML)+IWM(LMU)+1)*NEQ
      DO 410 I=1,LENPD
  410    WM(NPDM1+I)=0.0D0
      CALL PROBLEM%JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 550
!
!
!     BANDED FINITE-DIFFERENCE-GENERATED MATRIX
  500 MBAND=IWM(LML)+IWM(LMU)+1
      MBA=MIN(MBAND,NEQ)
      MEBAND=MBAND+IWM(LML)
      MEB1=MEBAND-1
      MSAVE=(NEQ/MBAND)+1
      ISAVE=NTEMP-1
      IPSAVE=ISAVE+MSAVE
      IRES=0
      SQUR=SQRT(UROUND)
      DO 540 J=1,MBA
         DO 510 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          WM(ISAVE+K)=Y(N)
          WM(IPSAVE+K)=YPRIME(N)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          Y(N)=Y(N)+DEL
  510     YPRIME(N)=YPRIME(N)+CJ*DEL
      CALL problem%res(X,Y,YPRIME,E,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) RETURN
      DO 530 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          Y(N)=WM(ISAVE+K)
          YPRIME(N)=WM(IPSAVE+K)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          DELINV=1.0D0/DEL
          I1=MAX(1,(N-IWM(LMU)))
          I2=MIN(NEQ,(N+IWM(LML)))
          II=N*MEB1-IWM(LML)+NPDM1
          DO 520 I=I1,I2
  520       WM(II+I)=(E(I)-DELTA(I))*DELINV
  530    CONTINUE
  540 END DO
!
!
!     DO LU DECOMPOSITION OF BANDED PD
  550 CALL DGBFA(WM(NPD),MEBAND,NEQ,                                    &
     &    IWM(LML),IWM(LMU),IWM(LIPVT),IER)
      RETURN
!------END OF SUBROUTINE DDAJAC------
      END subroutine ddajac
!DECK DDANRM
      REAL(8) FUNCTION DDANRM (NEQ, V, WT, RPAR, IPAR)
!***BEGIN PROLOGUE  DDANRM
!***SUBSIDIARY
!***PURPOSE  Compute vector norm for DDASSL.
!***LIBRARY   SLATEC (DASSL)
!***TYPE      REAL(8) (SDANRM-S, DDANRM-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------------
!     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED
!     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH
!     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS
!     CONTAINED IN THE ARRAY WT OF LENGTH NEQ.
!        DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
!-----------------------------------------------------------------------
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!***END PROLOGUE  DDANRM
!
      INTEGER  NEQ, IPAR(*)
      REAL(8)  V(NEQ), WT(NEQ), RPAR(*)
!
      INTEGER  I
      REAL(8)  SUM, VMAX
!
!***FIRST EXECUTABLE STATEMENT  DDANRM
      DDANRM = 0.0D0
      VMAX = 0.0D0
      DO 10 I = 1,NEQ
        IF(ABS(V(I)/WT(I)) .GT. VMAX) VMAX = ABS(V(I)/WT(I))
   10   CONTINUE
      IF(VMAX .LE. 0.0D0) GO TO 30
      SUM = 0.0D0
      DO 20 I = 1,NEQ
   20   SUM = SUM + ((V(I)/WT(I))/VMAX)**2
      DDANRM = VMAX*SQRT(SUM/NEQ)
   30 CONTINUE
      RETURN
!------END OF FUNCTION DDANRM------
      END function ddanrm
!DECK DDASLV
      SUBROUTINE DDASLV (NEQ, DELTA, WM, IWM)
!***BEGIN PROLOGUE  DDASLV
!***SUBSIDIARY
!***PURPOSE  Linear system solver for DDASSL.
!***LIBRARY   SLATEC (DASSL)
!***TYPE      REAL(8) (SDASLV-S, DDASLV-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------------
!     THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR
!     SYSTEM ARISING IN THE NEWTON ITERATION.
!     MATRICES AND REAL TEMPORARY STORAGE AND
!     REAL INFORMATION ARE STORED IN THE ARRAY WM.
!     INTEGER MATRIX INFORMATION IS STORED IN
!     THE ARRAY IWM.
!     FOR A DENSE MATRIX, THE LINPACK ROUTINE
!     DGESL IS CALLED.
!     FOR A BANDED MATRIX,THE LINPACK ROUTINE
!     DGBSL IS CALLED.
!-----------------------------------------------------------------------
!***ROUTINES CALLED  DGBSL, DGESL
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!***END PROLOGUE  DDASLV
!
      INTEGER  NEQ, IWM(*)
      REAL(8)  DELTA(*), WM(*)
!
      EXTERNAL  DGBSL, DGESL
!
      INTEGER  LIPVT, LML, LMU, LMTYPE, MEBAND, MTYPE, NPD
      PARAMETER (NPD=1)
      PARAMETER (LML=1)
      PARAMETER (LMU=2)
      PARAMETER (LMTYPE=4)
      PARAMETER (LIPVT=21)
!
!***FIRST EXECUTABLE STATEMENT  DDASLV
      MTYPE=IWM(LMTYPE)
      GO TO(100,100,300,400,400),MTYPE
!
!     DENSE MATRIX
  100 CALL DGESL(WM(NPD),NEQ,NEQ,IWM(LIPVT),DELTA,0)
      RETURN
!
!     DUMMY SECTION FOR MTYPE=3
  300 CONTINUE
      RETURN
!
!     BANDED MATRIX
  400 MEBAND=2*IWM(LML)+IWM(LMU)+1
      CALL DGBSL(WM(NPD),MEBAND,NEQ,IWM(LML),                           &
     &  IWM(LMU),IWM(LIPVT),DELTA,0)
      RETURN
!------END OF SUBROUTINE DDASLV------
      END subroutine ddaslv
!DECK DDASTP
      SUBROUTINE DDASTP (problem, X, Y, YPRIME, NEQ, H, WT, JSTART,         &
     &   IDID, RPAR, IPAR, PHI, DELTA, E, WM, IWM, ALPHA, BETA, GAMMA,  &
     &   PSI, SIGMA, CJ, CJOLD, HOLD, S, HMIN, UROUND, IPHASE, JCALC, K,&
     &   KOLD, NS, NONNEG, NTEMP)
!***BEGIN PROLOGUE  DDASTP
!***SUBSIDIARY
!***PURPOSE  Perform one step of the DDASSL integration.
!***LIBRARY   SLATEC (DASSL)
!***TYPE      REAL(8) (SDASTP-S, DDASTP-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------------
!     DDASTP SOLVES A SYSTEM OF DIFFERENTIAL/
!     ALGEBRAIC EQUATIONS OF THE FORM
!     G(X,Y,YPRIME) = 0,  FOR ONE STEP (NORMALLY
!     FROM X TO X+H).
!
!     THE METHODS USED ARE MODIFIED DIVIDED
!     DIFFERENCE,FIXED LEADING COEFFICIENT
!     FORMS OF BACKWARD DIFFERENTIATION
!     FORMULAS. THE CODE ADJUSTS THE STEPSIZE
!     AND ORDER TO CONTROL THE LOCAL ERROR PER
!     STEP.
!
!
!     THE PARAMETERS REPRESENT
!     X  --        INDEPENDENT VARIABLE
!     Y  --        SOLUTION VECTOR AT X
!     YPRIME --    DERIVATIVE OF SOLUTION VECTOR
!                  AFTER SUCCESSFUL STEP
!     NEQ --       NUMBER OF EQUATIONS TO BE INTEGRATED
!     RES --       EXTERNAL USER-SUPPLIED SUBROUTINE
!                  TO EVALUATE THE RESIDUAL.  THE CALL IS
!                  CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
!                  X,Y,YPRIME ARE INPUT.  DELTA IS OUTPUT.
!                  ON INPUT, IRES=0.  RES SHOULD ALTER IRES ONLY
!                  IF IT ENCOUNTERS AN ILLEGAL VALUE OF Y OR A
!                  STOP CONDITION.  SET IRES=-1 IF AN INPUT VALUE
!                  OF Y IS ILLEGAL, AND DDASTP WILL TRY TO SOLVE
!                  THE PROBLEM WITHOUT GETTING IRES = -1.  IF
!                  IRES=-2, DDASTP RETURNS CONTROL TO THE CALLING
!                  PROGRAM WITH IDID = -11.
!     JAC --       EXTERNAL USER-SUPPLIED ROUTINE TO EVALUATE
!                  THE ITERATION MATRIX (THIS IS OPTIONAL)
!                  THE CALL IS OF THE FORM
!                  CALL JAC(X,Y,YPRIME,PD,CJ,RPAR,IPAR)
!                  PD IS THE MATRIX OF PARTIAL DERIVATIVES,
!                  PD=DG/DY+CJ*DG/DYPRIME
!     H --         APPROPRIATE STEP SIZE FOR NEXT STEP.
!                  NORMALLY DETERMINED BY THE CODE
!     WT --        VECTOR OF WEIGHTS FOR ERROR CRITERION.
!     JSTART --    INTEGER VARIABLE SET 0 FOR
!                  FIRST STEP, 1 OTHERWISE.
!     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS:
!                  IDID= 1 -- THE STEP WAS COMPLETED SUCCESSFULLY
!                  IDID=-6 -- THE ERROR TEST FAILED REPEATEDLY
!                  IDID=-7 -- THE CORRECTOR COULD NOT CONVERGE
!                  IDID=-8 -- THE ITERATION MATRIX IS SINGULAR
!                  IDID=-9 -- THE CORRECTOR COULD NOT CONVERGE.
!                             THERE WERE REPEATED ERROR TEST
!                             FAILURES ON THIS STEP.
!                  IDID=-10-- THE CORRECTOR COULD NOT CONVERGE
!                             BECAUSE IRES WAS EQUAL TO MINUS ONE
!                  IDID=-11-- IRES EQUAL TO -2 WAS ENCOUNTERED,
!                             AND CONTROL IS BEING RETURNED TO
!                             THE CALLING PROGRAM
!     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS THAT
!                  ARE USED FOR COMMUNICATION BETWEEN THE
!                  CALLING PROGRAM AND EXTERNAL USER ROUTINES
!                  THEY ARE NOT ALTERED BY DDASTP
!     PHI --       ARRAY OF DIVIDED DIFFERENCES USED BY
!                  DDASTP. THE LENGTH IS NEQ*(K+1),WHERE
!                  K IS THE MAXIMUM ORDER
!     DELTA,E --   WORK VECTORS FOR DDASTP OF LENGTH NEQ
!     WM,IWM --    REAL AND INTEGER ARRAYS STORING
!                  MATRIX INFORMATION SUCH AS THE MATRIX
!                  OF PARTIAL DERIVATIVES,PERMUTATION
!                  VECTOR, AND VARIOUS OTHER INFORMATION.
!
!     THE OTHER PARAMETERS ARE INFORMATION
!     WHICH IS NEEDED INTERNALLY BY DDASTP TO
!     CONTINUE FROM STEP TO STEP.
!
!-----------------------------------------------------------------------
!***ROUTINES CALLED  DDAJAC, DDANRM, DDASLV, DDATRP
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!   951030  Reset PSI(1), PHI(*,2) at 690. (ACH)
!   000711  Fixed Newton convergence test below 360 (ACH)
!***END PROLOGUE  DDASTP
!
      INTEGER  NEQ, JSTART, IDID, IPAR(*), IWM(*), IPHASE, JCALC, K,    &
     &   KOLD, NS, NONNEG, NTEMP
      REAL(8)                                                  &
     &   X, Y(*), YPRIME(*), H, WT(*), RPAR(*), PHI(NEQ,*), DELTA(*),   &
     &   E(*), WM(*), ALPHA(*), BETA(*), GAMMA(*), PSI(*), SIGMA(*), CJ,&
     &   CJOLD, HOLD, S, HMIN, UROUND
      class(problem_ddassl), target :: problem
!
      INTEGER  I, IER, IRES, J, J1, KDIFF, KM1, KNEW, KP1, KP2, LCTF,   &
     &   LETF, LMXORD, LNJE, LNRE, LNST, M, MAXIT, NCF, NEF, NSF, NSP1
      REAL(8)                                                  &
     &   ALPHA0, ALPHAS, CJLAST, CK, DELNRM, ENORM, ERK, ERKM1,         &
     &   ERKM2, ERKP1, ERR, EST, HNEW, OLDNRM, PNORM, R, RATE, TEMP1,   &
     &   TEMP2, TERK, TERKM1, TERKM2, TERKP1, XOLD, XRATE
      LOGICAL  CONVGD
!
      PARAMETER (LMXORD=3)
      PARAMETER (LNST=11)
      PARAMETER (LNRE=12)
      PARAMETER (LNJE=13)
      PARAMETER (LETF=14)
      PARAMETER (LCTF=15)
!
      DATA MAXIT/4/
      DATA XRATE/0.25D0/
!
!
!
!
!
!-----------------------------------------------------------------------
!     BLOCK 1.
!     INITIALIZE. ON THE FIRST CALL,SET
!     THE ORDER TO 1 AND INITIALIZE
!     OTHER VARIABLES.
!-----------------------------------------------------------------------
!
!     INITIALIZATIONS FOR ALL CALLS
!***FIRST EXECUTABLE STATEMENT  DDASTP
      IDID=1
      XOLD=X
      NCF=0
      NSF=0
      NEF=0
      IF(JSTART .NE. 0) GO TO 120
!
!     IF THIS IS THE FIRST STEP,PERFORM
!     OTHER INITIALIZATIONS
      IWM(LETF) = 0
      IWM(LCTF) = 0
      K=1
      KOLD=0
      HOLD=0.0D0
      JSTART=1
      PSI(1)=H
      CJOLD = 1.0D0/H
      CJ = CJOLD
      S = 100.D0
      JCALC = -1
      DELNRM=1.0D0
      IPHASE = 0
      NS=0
  120 CONTINUE
!
!
!
!
!
!-----------------------------------------------------------------------
!     BLOCK 2
!     COMPUTE COEFFICIENTS OF FORMULAS FOR
!     THIS STEP.
!-----------------------------------------------------------------------
  200 CONTINUE
      KP1=K+1
      KP2=K+2
      KM1=K-1
      XOLD=X
      IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0
      NS=MIN(NS+1,KOLD+2)
      NSP1=NS+1
      IF(KP1 .LT. NS)GO TO 230
!
      BETA(1)=1.0D0
      ALPHA(1)=1.0D0
      TEMP1=H
      GAMMA(1)=0.0D0
      SIGMA(1)=1.0D0
      DO 210 I=2,KP1
         TEMP2=PSI(I-1)
         PSI(I-1)=TEMP1
         BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
         TEMP1=TEMP2+H
         ALPHA(I)=H/TEMP1
         SIGMA(I)=(I-1)*SIGMA(I-1)*ALPHA(I)
         GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
  210    CONTINUE
      PSI(KP1)=TEMP1
  230 CONTINUE
!
!     COMPUTE ALPHAS, ALPHA0
      ALPHAS = 0.0D0
      ALPHA0 = 0.0D0
      DO 240 I = 1,K
        ALPHAS = ALPHAS - 1.0D0/I
        ALPHA0 = ALPHA0 - ALPHA(I)
  240   CONTINUE
!
!     COMPUTE LEADING COEFFICIENT CJ
      CJLAST = CJ
      CJ = -ALPHAS/H
!
!     COMPUTE VARIABLE STEPSIZE ERROR COEFFICIENT CK
      CK = ABS(ALPHA(KP1) + ALPHAS - ALPHA0)
      CK = MAX(CK,ALPHA(KP1))
!
!     DECIDE WHETHER NEW JACOBIAN IS NEEDED
      TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)
      TEMP2 = 1.0D0/TEMP1
      IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
      IF (CJ .NE. CJLAST) S = 100.D0
!
!     CHANGE PHI TO PHI STAR
      IF(KP1 .LT. NSP1) GO TO 280
      DO 270 J=NSP1,KP1
         DO 260 I=1,NEQ
  260       PHI(I,J)=BETA(J)*PHI(I,J)
  270    CONTINUE
  280 CONTINUE
!
!     UPDATE TIME
      X=X+H
!
!
!
!
!
!-----------------------------------------------------------------------
!     BLOCK 3
!     PREDICT THE SOLUTION AND DERIVATIVE,
!     AND SOLVE THE CORRECTOR EQUATION
!-----------------------------------------------------------------------
!
!     FIRST,PREDICT THE SOLUTION AND DERIVATIVE
  300 CONTINUE
      DO 310 I=1,NEQ
         Y(I)=PHI(I,1)
  310    YPRIME(I)=0.0D0
      DO 330 J=2,KP1
         DO 320 I=1,NEQ
            Y(I)=Y(I)+PHI(I,J)
  320       YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
  330 END DO
      PNORM = DDANRM (NEQ,Y,WT,RPAR,IPAR)
!
!
!
!     SOLVE THE CORRECTOR EQUATION USING A
!     MODIFIED NEWTON SCHEME.
      CONVGD= .TRUE.
      M=0
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL problem%res(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
!
!
!     IF INDICATED,REEVALUATE THE
!     ITERATION MATRIX PD = DG/DY + CJ*DG/DYPRIME
!     (WHERE G(X,Y,YPRIME)=0). SET
!     JCALC TO 0 AS AN INDICATOR THAT
!     THIS HAS BEEN DONE.
      IF(JCALC .NE. -1)GO TO 340
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      CALL DDAJAC(PROBLEM,NEQ,X,Y,YPRIME,DELTA,CJ,H,                    &
     & IER,WT,E,WM,IWM,IRES,UROUND,RPAR,                                &
     & IPAR,NTEMP)
      CJOLD=CJ
      S = 100.D0
      IF (IRES .LT. 0) GO TO 380
      IF(IER .NE. 0)GO TO 380
      NSF=0
!
!
!     INITIALIZE THE ERROR ACCUMULATION VECTOR E.
  340 CONTINUE
      DO 345 I=1,NEQ
  345    E(I)=0.0D0
!
!
!     CORRECTOR LOOP.
  350 CONTINUE
!
!     MULTIPLY RESIDUAL BY TEMP1 TO ACCELERATE CONVERGENCE
      TEMP1 = 2.0D0/(1.0D0 + CJ/CJOLD)
      DO 355 I = 1,NEQ
  355   DELTA(I) = DELTA(I) * TEMP1
!
!     COMPUTE A NEW ITERATE (BACK-SUBSTITUTION).
!     STORE THE CORRECTION IN DELTA.
      CALL DDASLV(NEQ,DELTA,WM,IWM)
!
!     UPDATE Y, E, AND YPRIME
      DO 360 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
         E(I)=E(I)-DELTA(I)
  360    YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
!
!     TEST FOR CONVERGENCE OF THE ITERATION
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (M .GT. 0) GO TO 365
         IF (DELNRM .LE. 100.D0*UROUND*PNORM) GO TO 375
         OLDNRM = DELNRM
         GO TO 367
  365 RATE = (DELNRM/OLDNRM)**(1.0D0/M)
      IF (RATE .GT. 0.90D0) GO TO 370
      S = RATE/(1.0D0 - RATE)
  367 IF (S*DELNRM .LE. 0.33D0) GO TO 375
!
!     THE CORRECTOR HAS NOT YET CONVERGED.
!     UPDATE M AND TEST WHETHER THE
!     MAXIMUM NUMBER OF ITERATIONS HAVE
!     BEEN TRIED.
      M=M+1
      IF(M.GE.MAXIT)GO TO 370
!
!     EVALUATE THE RESIDUAL
!     AND GO BACK TO DO ANOTHER ITERATION
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL problem%res(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
      GO TO 350
!
!
!     THE CORRECTOR FAILED TO CONVERGE IN MAXIT
!     ITERATIONS. IF THE ITERATION MATRIX
!     IS NOT CURRENT,RE-DO THE STEP WITH
!     A NEW ITERATION MATRIX.
  370 CONTINUE
      IF(JCALC.EQ.0)GO TO 380
      JCALC=-1
      GO TO 300
!
!
!     THE ITERATION HAS CONVERGED.  IF NONNEGATIVITY OF SOLUTION IS
!     REQUIRED, SET THE SOLUTION NONNEGATIVE, IF THE PERTURBATION
!     TO DO IT IS SMALL ENOUGH.  IF THE CHANGE IS TOO LARGE, THEN
!     CONSIDER THE CORRECTOR ITERATION TO HAVE FAILED.
  375 IF(NONNEG .EQ. 0) GO TO 390
      DO 377 I = 1,NEQ
  377    DELTA(I) = MIN(Y(I),0.0D0)
      DELNRM = DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF(DELNRM .GT. 0.33D0) GO TO 380
      DO 378 I = 1,NEQ
  378    E(I) = E(I) - DELTA(I)
      GO TO 390
!
!
!     EXITS FROM BLOCK 3
!     NO CONVERGENCE WITH CURRENT ITERATION
!     MATRIX,OR SINGULAR ITERATION MATRIX
  380 CONVGD= .FALSE.
  390 JCALC = 1
      IF(.NOT.CONVGD)GO TO 600
!
!
!
!
!
!-----------------------------------------------------------------------
!     BLOCK 4
!     ESTIMATE THE ERRORS AT ORDERS K,K-1,K-2
!     AS IF CONSTANT STEPSIZE WAS USED. ESTIMATE
!     THE LOCAL ERROR AT ORDER K AND TEST
!     WHETHER THE CURRENT STEP IS SUCCESSFUL.
!-----------------------------------------------------------------------
!
!     ESTIMATE ERRORS AT ORDERS K,K-1,K-2
      ENORM = DDANRM(NEQ,E,WT,RPAR,IPAR)
      ERK = SIGMA(K+1)*ENORM
      TERK = (K+1)*ERK
      EST = ERK
      KNEW=K
      IF(K .EQ. 1)GO TO 430
      DO 405 I = 1,NEQ
  405   DELTA(I) = PHI(I,KP1) + E(I)
      ERKM1=SIGMA(K)*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM1 = K*ERKM1
      IF(K .GT. 2)GO TO 410
      IF(TERKM1 .LE. 0.5D0*TERK)GO TO 420
      GO TO 430
  410 CONTINUE
      DO 415 I = 1,NEQ
  415   DELTA(I) = PHI(I,K) + DELTA(I)
      ERKM2=SIGMA(K-1)*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM2 = (K-1)*ERKM2
      IF(MAX(TERKM1,TERKM2).GT.TERK)GO TO 430
!     LOWER THE ORDER
  420 CONTINUE
      KNEW=K-1
      EST = ERKM1
!
!
!     CALCULATE THE LOCAL ERROR FOR THE CURRENT STEP
!     TO SEE IF THE STEP WAS SUCCESSFUL
  430 CONTINUE
      ERR = CK * ENORM
      IF(ERR .GT. 1.0D0)GO TO 600
!
!
!
!
!
!-----------------------------------------------------------------------
!     BLOCK 5
!     THE STEP IS SUCCESSFUL. DETERMINE
!     THE BEST ORDER AND STEPSIZE FOR
!     THE NEXT STEP. UPDATE THE DIFFERENCES
!     FOR THE NEXT STEP.
!-----------------------------------------------------------------------
      IDID=1
      IWM(LNST)=IWM(LNST)+1
      KDIFF=K-KOLD
      KOLD=K
      HOLD=H
!
!
!     ESTIMATE THE ERROR AT ORDER K+1 UNLESS:
!        ALREADY DECIDED TO LOWER ORDER, OR
!        ALREADY USING MAXIMUM ORDER, OR
!        STEPSIZE NOT CONSTANT, OR
!        ORDER RAISED IN PREVIOUS STEP
      IF(KNEW.EQ.KM1.OR.K.EQ.IWM(LMXORD))IPHASE=1
      IF(IPHASE .EQ. 0)GO TO 545
      IF(KNEW.EQ.KM1)GO TO 540
      IF(K.EQ.IWM(LMXORD)) GO TO 550
      IF(KP1.GE.NS.OR.KDIFF.EQ.1)GO TO 550
      DO 510 I=1,NEQ
  510    DELTA(I)=E(I)-PHI(I,KP2)
      ERKP1 = (1.0D0/(K+2))*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKP1 = (K+2)*ERKP1
      IF(K.GT.1)GO TO 520
      IF(TERKP1.GE.0.5D0*TERK)GO TO 550
      GO TO 530
  520 IF(TERKM1.LE.MIN(TERK,TERKP1))GO TO 540
      IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD))GO TO 550
!
!     RAISE ORDER
  530 K=KP1
      EST = ERKP1
      GO TO 550
!
!     LOWER ORDER
  540 K=KM1
      EST = ERKM1
      GO TO 550
!
!     IF IPHASE = 0, INCREASE ORDER BY ONE AND MULTIPLY STEPSIZE BY
!     FACTOR TWO
  545 K = KP1
      HNEW = H*2.0D0
      H = HNEW
      GO TO 575
!
!
!     DETERMINE THE APPROPRIATE STEPSIZE FOR
!     THE NEXT STEP.
  550 HNEW=H
      TEMP2=K+1
      R=(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      IF(R .LT. 2.0D0) GO TO 555
      HNEW = 2.0D0*H
      GO TO 560
  555 IF(R .GT. 1.0D0) GO TO 560
      R = MAX(0.5D0,MIN(0.9D0,R))
      HNEW = H*R
  560 H=HNEW
!
!
!     UPDATE DIFFERENCES FOR NEXT STEP
  575 CONTINUE
      IF(KOLD.EQ.IWM(LMXORD))GO TO 585
      DO 580 I=1,NEQ
  580    PHI(I,KP2)=E(I)
  585 CONTINUE
      DO 590 I=1,NEQ
  590    PHI(I,KP1)=PHI(I,KP1)+E(I)
      DO J1=2,KP1
         J=KP1-J1+1
         DO I=1,NEQ
            PHI(I,J)=PHI(I,J)+PHI(I,J+1)
         end DO
      end DO

      RETURN
!
!
!
!
!
!-----------------------------------------------------------------------
!     BLOCK 6
!     THE STEP IS UNSUCCESSFUL. RESTORE X,PSI,PHI
!     DETERMINE APPROPRIATE STEPSIZE FOR
!     CONTINUING THE INTEGRATION, OR EXIT WITH
!     AN ERROR FLAG IF THERE HAVE BEEN MANY
!     FAILURES.
!-----------------------------------------------------------------------
  600 IPHASE = 1
!
!     RESTORE X,PHI,PSI
      X=XOLD
      IF(KP1.LT.NSP1)GO TO 630
      DO 620 J=NSP1,KP1
         TEMP1=1.0D0/BETA(J)
         DO 610 I=1,NEQ
  610       PHI(I,J)=TEMP1*PHI(I,J)
  620    CONTINUE
  630 CONTINUE
      DO 640 I=2,KP1
  640    PSI(I-1)=PSI(I)-H
!
!
!     TEST WHETHER FAILURE IS DUE TO CORRECTOR ITERATION
!     OR ERROR TEST
      IF(CONVGD)GO TO 660
      IWM(LCTF)=IWM(LCTF)+1
!
!
!     THE NEWTON ITERATION FAILED TO CONVERGE WITH
!     A CURRENT ITERATION MATRIX.  DETERMINE THE CAUSE
!     OF THE FAILURE AND TAKE APPROPRIATE ACTION.
      IF(IER.EQ.0)GO TO 650
!
!     THE ITERATION MATRIX IS SINGULAR. REDUCE
!     THE STEPSIZE BY A FACTOR OF 4. IF
!     THIS HAPPENS THREE TIMES IN A ROW ON
!     THE SAME STEP, RETURN WITH AN ERROR FLAG
      NSF=NSF+1
      R = 0.25D0
      H=H*R
      IF (NSF .LT. 3 .AND. ABS(H) .GE. HMIN) GO TO 690
      IDID=-8
      GO TO 675
!
!
!     THE NEWTON ITERATION FAILED TO CONVERGE FOR A REASON
!     OTHER THAN A SINGULAR ITERATION MATRIX.  IF IRES = -2, THEN
!     RETURN.  OTHERWISE, REDUCE THE STEPSIZE AND TRY AGAIN, UNLESS
!     TOO MANY FAILURES HAVE OCCURRED.
  650 CONTINUE
      IF (IRES .GT. -2) GO TO 655
      IDID = -11
      GO TO 675
  655 NCF = NCF + 1
      R = 0.25D0
      H = H*R
      IF (NCF .LT. 10 .AND. ABS(H) .GE. HMIN) GO TO 690
      IDID = -7
      IF (IRES .LT. 0) IDID = -10
      IF (NEF .GE. 3) IDID = -9
      GO TO 675
!
!
!     THE NEWTON SCHEME CONVERGED, AND THE CAUSE
!     OF THE FAILURE WAS THE ERROR ESTIMATE
!     EXCEEDING THE TOLERANCE.
  660 NEF=NEF+1
      IWM(LETF)=IWM(LETF)+1
      IF (NEF .GT. 1) GO TO 665
!
!     ON FIRST ERROR TEST FAILURE, KEEP CURRENT ORDER OR LOWER
!     ORDER BY ONE.  COMPUTE NEW STEPSIZE BASED ON DIFFERENCES
!     OF THE SOLUTION.
      K = KNEW
      TEMP2 = K + 1
      R = 0.90D0*(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      R = MAX(0.25D0,MIN(0.9D0,R))
      H = H*R
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
!
!     ON SECOND ERROR TEST FAILURE, USE THE CURRENT ORDER OR
!     DECREASE ORDER BY ONE.  REDUCE THE STEPSIZE BY A FACTOR OF
!     FOUR.
  665 IF (NEF .GT. 2) GO TO 670
      K = KNEW
      R = 0.25D0
      H = R*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
!
!     ON THIRD AND SUBSEQUENT ERROR TEST FAILURES, SET THE ORDER TO
!     ONE AND REDUCE THE STEPSIZE BY A FACTOR OF FOUR.
  670 K = 1
      R = 0.25D0
      H = R*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
!
!
!     FOR ALL CRASHES, RESTORE Y TO ITS LAST VALUE,
!     INTERPOLATE TO FIND YPRIME AT LAST X, AND RETURN
  675 CONTINUE
      CALL DDATRP(X,X,Y,YPRIME,NEQ,K,PHI,PSI)
      RETURN
!
!
!     GO BACK AND TRY THIS STEP AGAIN.
!     IF THIS IS THE FIRST STEP, RESET PSI(1) AND RESCALE PHI(*,2).
  690 IF (KOLD .EQ. 0) THEN
        PSI(1) = H
        DO 695 I = 1,NEQ
  695     PHI(I,2) = R*PHI(I,2)
        ENDIF
      GO TO 200
!
!------END OF SUBROUTINE DDASTP------
      END subroutine ddastp
!DECK DDATRP
      SUBROUTINE DDATRP (X, XOUT, YOUT, YPOUT, NEQ, KOLD, PHI, PSI)
!***BEGIN PROLOGUE  DDATRP
!***SUBSIDIARY
!***PURPOSE  Interpolation routine for DDASSL.
!***LIBRARY   SLATEC (DASSL)
!***TYPE      REAL(8) (SDATRP-S, DDATRP-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------------
!     THE METHODS IN SUBROUTINE DDASTP USE POLYNOMIALS
!     TO APPROXIMATE THE SOLUTION. DDATRP APPROXIMATES THE
!     SOLUTION AND ITS DERIVATIVE AT TIME XOUT BY EVALUATING
!     ONE OF THESE POLYNOMIALS, AND ITS DERIVATIVE,THERE.
!     INFORMATION DEFINING THIS POLYNOMIAL IS PASSED FROM
!     DDASTP, SO DDATRP CANNOT BE USED ALONE.
!
!     THE PARAMETERS ARE:
!     X     THE CURRENT TIME IN THE INTEGRATION.
!     XOUT  THE TIME AT WHICH THE SOLUTION IS DESIRED
!     YOUT  THE INTERPOLATED APPROXIMATION TO Y AT XOUT
!           (THIS IS OUTPUT)
!     YPOUT THE INTERPOLATED APPROXIMATION TO YPRIME AT XOUT
!           (THIS IS OUTPUT)
!     NEQ   NUMBER OF EQUATIONS
!     KOLD  ORDER USED ON LAST SUCCESSFUL STEP
!     PHI   ARRAY OF SCALED DIVIDED DIFFERENCES OF Y
!     PSI   ARRAY OF PAST STEPSIZE HISTORY
!-----------------------------------------------------------------------
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!***END PROLOGUE  DDATRP
!
      INTEGER  NEQ, KOLD
      REAL(8)  X, XOUT, YOUT(*), YPOUT(*), PHI(NEQ,*), PSI(*)
!
      INTEGER  I, J, KOLDP1
      REAL(8)  C, D, GAMMA, TEMP1
!
!***FIRST EXECUTABLE STATEMENT  DDATRP
      KOLDP1=KOLD+1
      TEMP1=XOUT-X
      DO 10 I=1,NEQ
         YOUT(I)=PHI(I,1)
   10    YPOUT(I)=0.0D0
      C=1.0D0
      D=0.0D0
      GAMMA=TEMP1/PSI(1)
      DO 30 J=2,KOLDP1
         D=D*GAMMA+C/PSI(J-1)
         C=C*GAMMA
         GAMMA=(TEMP1+PSI(J-1))/PSI(J)
         DO 20 I=1,NEQ
            YOUT(I)=YOUT(I)+C*PHI(I,J)
   20       YPOUT(I)=YPOUT(I)+D*PHI(I,J)
   30    CONTINUE
      RETURN
!
!------END OF SUBROUTINE DDATRP------
      END subroutine ddatrp
!DECK DDAWTS
      SUBROUTINE DDAWTS (NEQ, IWT, RTOL, ATOL, Y, WT, RPAR, IPAR)
!***BEGIN PROLOGUE  DDAWTS
!***SUBSIDIARY
!***PURPOSE  Set error weight vector for DDASSL.
!***LIBRARY   SLATEC (DASSL)
!***TYPE      REAL(8) (SDAWTS-S, DDAWTS-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------------
!     THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR
!     WT ACCORDING TO WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
!     I=1,-,N.
!     RTOL AND ATOL ARE SCALARS IF IWT = 0,
!     AND VECTORS IF IWT = 1.
!-----------------------------------------------------------------------
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!***END PROLOGUE  DDAWTS
!
      INTEGER  NEQ, IWT, IPAR(*)
      REAL(8)  RTOL(*), ATOL(*), Y(*), WT(*), RPAR(*)
!
      INTEGER  I
      REAL(8)  ATOLI, RTOLI
!
!***FIRST EXECUTABLE STATEMENT  DDAWTS
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 20 I=1,NEQ
         IF (IWT .EQ.0) GO TO 10
           RTOLI=RTOL(I)
           ATOLI=ATOL(I)
   10      WT(I)=RTOLI*ABS(Y(I))+ATOLI
   20      CONTINUE
      RETURN
!-----------END OF SUBROUTINE DDAWTS------------------------------------
      END subroutine ddawts
!DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
!-----------------------------------------------------------------------
! Subroutines XERMSG, XSETF, XSETUN, and the function routine IXSAV, as
! given here, constitute a simplified version of the SLATEC error
! handling package.  Written by A. C. Hindmarsh, 18 November 1992.
!
! All arguments are input arguments.
! LIBRAR = Library name (character array).  Prefixed to message.
! SUBROU = Routine name (character array).  Prefixed to message.
! MESSG  = The message (character array).
! NERR   = Integer error number.  Prefixed to message.
! LEVEL  = The error level..
!          0 or 1 means recoverable (control returns to caller).
!          2 means fatal (run is aborted--see note below).
!
! Note..  This routine has been simplified in the following ways..
! 1. A single prefix line is printed with NERR, SUBROU, and LIBRAR.
! 2. The message in MESSG is printed, unaltered, on lines of up to 72
!    characters each using a format of (A).
! 3. If LEVEL = 2, control passes to the statement   STOP
!    to abort the run.  This statement may be machine-dependent.
!
! For a different default logical unit number, change the data
! statement in function routine IXSAV.
! For a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
! Subroutines called by XERMSG.. None
! Function routines called by XERMSG.. IXSAV
! Intrinsic function used by XERMSG.. LEN
!-----------------------------------------------------------------------
      CHARACTER(len=*) :: LIBRAR, SUBROU, MESSG
      INTEGER NERR, LEVEL
      INTEGER I1, I2, IL, LENMSG, LLEN, LUNIT, MESFLG, NLINES
      PARAMETER (LLEN = 72)
!
! Get message print flag and logical unit number. ----------------------
      LUNIT = IXSAV (1, 0, .FALSE.)
      MESFLG = IXSAV (2, 0, .FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
! Write NERR, SUBROU, and LIBRAR. --------------------------------------
      I1 = LEN(SUBROU)
      I2 = LEN(LIBRAR)
      WRITE (LUNIT, 10) NERR, SUBROU(1:I1), LIBRAR(1:I2)
   10 FORMAT(/,'***Error number ',I6,' from ',A,' in library ',A,'***')
! Write the message. ---------------------------------------------------
      LENMSG = LEN(MESSG)
      NLINES = ( (LENMSG - 1)/LLEN ) + 1
      DO 20 IL = 1,NLINES
        I1 = 1 + (IL - 1)*LLEN
        I2 = MIN(IL*LLEN,LENMSG)
        WRITE (LUNIT,'(A)') MESSG(I1:I2)
   20   CONTINUE
! Abort the run if LEVEL = 2. ------------------------------------------
  100 IF (LEVEL .NE. 2) RETURN
      STOP
!----------------------- End of Subroutine XERMSG ----------------------
      END subroutine xermsg
!DECK XSETUN
      SUBROUTINE XSETUN (LUN)
!-----------------------------------------------------------------------
! This routine resets the logical unit number for messages.
!
! Subroutines called by XSETUN.. None
! Function routine called by XSETUN.. IXSAV
!-----------------------------------------------------------------------
      INTEGER LUN, JUNK
!
      IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
      RETURN
!----------------------- End of Subroutine XSETUN ----------------------
      END subroutine xsetun
!DECK XSETF
      SUBROUTINE XSETF (MFLAG)
!-----------------------------------------------------------------------
! This routine resets the print control flag MFLAG.
!
! Subroutines called by XSETF.. None
! Function routine called by XSETF.. IXSAV
!-----------------------------------------------------------------------
      INTEGER MFLAG, JUNK
!
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = IXSAV (2,MFLAG,.TRUE.)
      RETURN
!----------------------- End of Subroutine XSETF -----------------------
      END subroutine xsetf
!DECK IXSAV
      INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
      LOGICAL ISET
      INTEGER IPAR, IVALUE
!-----------------------------------------------------------------------
! IXSAV saves and recalls one of two error message parameters:
!   LUNIT, the logical unit number to which messages are printed, and
!   MESFLG, the message print flag.
! This is a modification of the SLATEC library routine J4SAVE.
!
! Saved local variables..
!  LUNIT  = Logical unit number for messages.
!           The default is 6 (machine-dependent).
!  MESFLG = Print control flag..
!           1 means print all messages (the default).
!           0 means no printing.
!
! On input..
!   IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
!   IVALUE = The value to be set for the parameter, if ISET = .TRUE.
!   ISET   = Logical flag to indicate whether to read or write.
!            If ISET = .TRUE., the parameter will be given
!            the value IVALUE.  If ISET = .FALSE., the parameter
!            will be unchanged, and IVALUE is a dummy argument.
!
! On return..
!   IXSAV = The (old) value of the parameter.
!
! Subroutines/functions called by IXSAV.. None
!-----------------------------------------------------------------------
      INTEGER LUNIT, MESFLG
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this routine.
!-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/6/, MESFLG/1/
!
      IF (IPAR .EQ. 1) THEN
        IXSAV = LUNIT
        IF (ISET) LUNIT = IVALUE
        ENDIF
!
      IF (IPAR .EQ. 2) THEN
        IXSAV = MESFLG
        IF (ISET) MESFLG = IVALUE
        ENDIF
!
      RETURN
!----------------------- End of Function IXSAV -------------------------
      END function ixsav
end module ddassl_mod
