module ddassl_mod

  type, abstract :: problem_ddassl
   contains
     procedure(res_i), deferred :: res
     procedure(jac_i), deferred :: jac
     procedure :: ddajac
     procedure :: ddaini
     procedure :: ddastp
     procedure :: ddassl
  end type problem_ddassl


  abstract interface

     subroutine res_i (eqn, t, y, ydot, res, ires, rwk, iwk)
       import problem_ddassl
       class(problem_ddassl) :: eqn
       integer :: ires
       integer, dimension(*) :: iwk
       real(8) :: t
       real(8), dimension(*) :: y, ydot, res, rwk
     end subroutine res_i

     subroutine jac_i (eqn, t, y, ydot, pd, cj, rwk, iwk)
       import problem_ddassl
       class(problem_ddassl) :: eqn
       integer :: ires
       integer, dimension(*) :: iwk
       real(8) :: t, cj
       real(8), dimension(*) :: y, ydot, rwk
       real(8), dimension(10,*) :: pd
     end subroutine jac_i

  end interface

contains
!deck ddassl
!-----------------------------------------------------------------------
! note:  users of this solver, ddassl, are encouraged to use the
! solver ddaspk instead.  ddaspk has a much improved initial condition
! calculation algorithm.  in addition, ddaspk includes iterative
! (krylov) methods for the linear systems that arise, in addition to
! the direct (dense/banded) methods in ddassl.
!-----------------------------------------------------------------------
!***begin prologue  ddassl
!***purpose  this code solves a system of differential/algebraic
!            equations of the form g(t,y,yprime) = 0.
!***library   slatec (dassl)
!***category  i1a2
!***type      double precision (sdassl-s, ddassl-d)
!***keywords  backward differentiation formulas, dassl,
!             differential/algebraic, implicit differential systems
!***author  petzold, linda r., (llnl)
!             computing and mathematics research division
!             lawrence livermore national laboratory
!             l - 316, p.o. box 808,
!             livermore, ca.    94550
!***description
!
! *usage:
!
!      external res, jac
!      integer neq, info(n), idid, lrw, liw, iwork(liw), ipar
!      double precision t, y(neq), yprime(neq), tout, rtol, atol,
!     *   rwork(lrw), rpar
!
!      call ddassl (res, neq, t, y, yprime, tout, info, rtol, atol,
!     *   idid, rwork, lrw, iwork, liw, rpar, ipar, jac)
!
!
! *arguments:
!  (in the following, all real arrays should be type double precision.)
!
!  res:ext     this is a subroutine which you provide to define the
!              differential/algebraic system.
!
!  neq:in      this is the number of equations to be solved.
!
!  t:inout     this is the current value of the independent variable.
!
!  y(*):inout  this array contains the solution components at t.
!
!  yprime(*):inout  this array contains the derivatives of the solution
!              components at t.
!
!  tout:in     this is a point at which a solution is desired.
!
!  info(n):in  the basic task of the code is to solve the system from t
!              to tout and return an answer at tout.  info is an integer
!              array which is used to communicate exactly how you want
!              this task to be carried out.  (see below for details.)
!              n must be greater than or equal to 15.
!
!  rtol,atol:inout  these quantities represent relative and absolute
!              error tolerances which you provide to indicate how
!              accurately you wish the solution to be computed.  you
!              may choose them to be both scalars or else both vectors.
!              caution:  in fortran 77, a scalar is not the same as an
!                        array of length 1.  some compilers may object
!                        to using scalars for rtol,atol.
!
!  idid:out    this scalar quantity is an indicator reporting what the
!              code did.  you must monitor this integer variable to
!              decide  what action to take next.
!
!  rwork:work  a real work array of length lrw which provides the
!              code with needed storage space.
!
!  lrw:in      the length of rwork.  (see below for required length.)
!
!  iwork:work  an integer work array of length liw which provides the
!              code with needed storage space.
!
!  liw:in      the length of iwork.  (see below for required length.)
!
!  rpar,ipar:in  these are real and integer parameter arrays which
!              you can use for communication between your calling
!              program and the res subroutine (and the jac subroutine)
!
!  jac:ext     this is the name of a subroutine which you may choose
!              to provide for defining a matrix of partial derivatives
!              described below.
!
!  quantities which may be altered by ddassl are:
!     t, y(*), yprime(*), info(1), rtol, atol,
!     idid, rwork(*) and iwork(*)
!
! *description
!
!  subroutine ddassl uses the backward differentiation formulas of
!  orders one through five to solve a system of the above form for y and
!  yprime.  values for y and yprime at the initial time must be given as
!  input.  these values must be consistent, (that is, if t,y,yprime are
!  the given initial values, they must satisfy g(t,y,yprime) = 0.).  the
!  subroutine solves the system from t to tout.  it is easy to continue
!  the solution to get results at additional tout.  this is the interval
!  mode of operation.  intermediate results can also be obtained easily
!  by using the intermediate-output capability.
!
!  the following detailed description is divided into subsections:
!    1. input required for the first call to ddassl.
!    2. output after any return from ddassl.
!    3. what to do to continue the integration.
!    4. error messages.
!
!
!  -------- input -- what to do on the first call to ddassl ------------
!
!  the first call of the code is defined to be the start of each new
!  problem. read through the descriptions of all the following items,
!  provide sufficient storage space for designated arrays, set
!  appropriate variables for the initialization of the problem, and
!  give information about how you want the problem to be solved.
!
!
!  res -- provide a subroutine of the form
!             subroutine res(t,y,yprime,delta,ires,rpar,ipar)
!         to define the system of differential/algebraic
!         equations which is to be solved. for the given values
!         of t,y and yprime, the subroutine should
!         return the residual of the differential/algebraic
!         system
!             delta = g(t,y,yprime)
!         (delta(*) is a vector of length neq which is
!         output for res.)
!
!         subroutine res must not alter t,y or yprime.
!         you must declare the name res in an external
!         statement in your program that calls ddassl.
!         you must dimension y,yprime and delta in res.
!
!         ires is an integer flag which is always equal to
!         zero on input. subroutine res should alter ires
!         only if it encounters an illegal value of y or
!         a stop condition. set ires = -1 if an input value
!         is illegal, and ddassl will try to solve the problem
!         without getting ires = -1. if ires = -2, ddassl
!         will return control to the calling program
!         with idid = -11.
!
!         rpar and ipar are real and integer parameter arrays which
!         you can use for communication between your calling program
!         and subroutine res. they are not altered by ddassl. if you
!         do not need rpar or ipar, ignore these parameters by treat-
!         ing them as dummy arguments. if you do choose to use them,
!         dimension them in your calling program and in res as arrays
!         of appropriate length.
!
!  neq -- set it to the number of differential equations.
!         (neq .ge. 1)
!
!  t -- set it to the initial point of the integration.
!         t must be defined as a variable.
!
!  y(*) -- set this vector to the initial values of the neq solution
!         components at the initial point. you must dimension y of
!         length at least neq in your calling program.
!
!  yprime(*) -- set this vector to the initial values of the neq
!         first derivatives of the solution components at the initial
!         point.  you must dimension yprime at least neq in your
!         calling program. if you do not know initial values of some
!         of the solution components, see the explanation of info(11).
!
!  tout -- set it to the first point at which a solution
!         is desired. you can not take tout = t.
!         integration either forward in t (tout .gt. t) or
!         backward in t (tout .lt. t) is permitted.
!
!         the code advances the solution from t to tout using
!         step sizes which are automatically selected so as to
!         achieve the desired accuracy. if you wish, the code will
!         return with the solution and its derivative at
!         intermediate steps (intermediate-output mode) so that
!         you can monitor them, but you still must provide tout in
!         accord with the basic aim of the code.
!
!         the first step taken by the code is a critical one
!         because it must reflect how fast the solution changes near
!         the initial point. the code automatically selects an
!         initial step size which is practically always suitable for
!         the problem. by using the fact that the code will not step
!         past tout in the first step, you could, if necessary,
!         restrict the length of the initial step size.
!
!         for some problems it may not be permissible to integrate
!         past a point tstop because a discontinuity occurs there
!         or the solution or its derivative is not defined beyond
!         tstop. when you have declared a tstop point (see info(4)
!         and rwork(1)), you have told the code not to integrate
!         past tstop. in this case any tout beyond tstop is invalid
!         input.
!
!  info(*) -- use the info array to give the code more details about
!         how you want your problem solved.  this array should be
!         dimensioned of length 15, though ddassl uses only the first
!         eleven entries.  you must respond to all of the following
!         items, which are arranged as questions.  the simplest use
!         of the code corresponds to answering all questions as yes,
!         i.e. setting all entries of info to 0.
!
!       info(1) - this parameter enables the code to initialize
!              itself. you must set it to indicate the start of every
!              new problem.
!
!          **** is this the first call for this problem ...
!                yes - set info(1) = 0
!                 no - not applicable here.
!                      see below for continuation calls.  ****
!
!       info(2) - how much accuracy you want of your solution
!              is specified by the error tolerances rtol and atol.
!              the simplest use is to take them both to be scalars.
!              to obtain more flexibility, they can both be vectors.
!              the code must be told your choice.
!
!          **** are both error tolerances rtol, atol scalars ...
!                yes - set info(2) = 0
!                      and input scalars for both rtol and atol
!                 no - set info(2) = 1
!                      and input arrays for both rtol and atol ****
!
!       info(3) - the code integrates from t in the direction
!              of tout by steps. if you wish, it will return the
!              computed solution and derivative at the next
!              intermediate step (the intermediate-output mode) or
!              tout, whichever comes first. this is a good way to
!              proceed if you want to see the behavior of the solution.
!              if you must have solutions at a great many specific
!              tout points, this code will compute them efficiently.
!
!          **** do you want the solution only at
!                tout (and not at the next intermediate step) ...
!                 yes - set info(3) = 0
!                  no - set info(3) = 1 ****
!
!       info(4) - to handle solutions at a great many specific
!              values tout efficiently, this code may integrate past
!              tout and interpolate to obtain the result at tout.
!              sometimes it is not possible to integrate beyond some
!              point tstop because the equation changes there or it is
!              not defined past tstop. then you must tell the code
!              not to go past.
!
!           **** can the integration be carried out without any
!                restrictions on the independent variable t ...
!                 yes - set info(4)=0
!                  no - set info(4)=1
!                       and define the stopping point tstop by
!                       setting rwork(1)=tstop ****
!
!       info(5) - to solve differential/algebraic problems it is
!              necessary to use a matrix of partial derivatives of the
!              system of differential equations. if you do not
!              provide a subroutine to evaluate it analytically (see
!              description of the item jac in the call list), it will
!              be approximated by numerical differencing in this code.
!              although it is less trouble for you to have the code
!              compute partial derivatives by numerical differencing,
!              the solution will be more reliable if you provide the
!              derivatives via jac. sometimes numerical differencing
!              is cheaper than evaluating derivatives in jac and
!              sometimes it is not - this depends on your problem.
!
!           **** do you want the code to evaluate the partial
!                derivatives automatically by numerical differences ...
!                   yes - set info(5)=0
!                    no - set info(5)=1
!                  and provide subroutine jac for evaluating the
!                  matrix of partial derivatives ****
!
!       info(6) - ddassl will perform much better if the matrix of
!              partial derivatives, dg/dy + cj*dg/dyprime,
!              (here cj is a scalar determined by ddassl)
!              is banded and the code is told this. in this
!              case, the storage needed will be greatly reduced,
!              numerical differencing will be performed much cheaper,
!              and a number of important algorithms will execute much
!              faster. the differential equation is said to have
!              half-bandwidths ml (lower) and mu (upper) if equation i
!              involves only unknowns y(j) with
!                             i-ml .le. j .le. i+mu
!              for all i=1,2,...,neq. thus, ml and mu are the widths
!              of the lower and upper parts of the band, respectively,
!              with the main diagonal being excluded. if you do not
!              indicate that the equation has a banded matrix of partial
!              derivatives, the code works with a full matrix of neq**2
!              elements (stored in the conventional way). computations
!              with banded matrices cost less time and storage than with
!              full matrices if 2*ml+mu .lt. neq. if you tell the
!              code that the matrix of partial derivatives has a banded
!              structure and you want to provide subroutine jac to
!              compute the partial derivatives, then you must be careful
!              to store the elements of the matrix in the special form
!              indicated in the description of jac.
!
!          **** do you want to solve the problem using a full
!               (dense) matrix (and not a special banded
!               structure) ...
!                yes - set info(6)=0
!                 no - set info(6)=1
!                       and provide the lower (ml) and upper (mu)
!                       bandwidths by setting
!                       iwork(1)=ml
!                       iwork(2)=mu ****
!
!
!        info(7) -- you can specify a maximum (absolute value of)
!              stepsize, so that the code
!              will avoid passing over very
!              large regions.
!
!          ****  do you want the code to decide
!                on its own maximum stepsize?
!                yes - set info(7)=0
!                 no - set info(7)=1
!                      and define hmax by setting
!                      rwork(2)=hmax ****
!
!        info(8) -- differential/algebraic problems
!              may occasionally suffer from
!              severe scaling difficulties on the
!              first step. if you know a great deal
!              about the scaling of your problem, you can
!              help to alleviate this problem by
!              specifying an initial stepsize ho.
!
!          ****  do you want the code to define
!                its own initial stepsize?
!                yes - set info(8)=0
!                 no - set info(8)=1
!                      and define ho by setting
!                      rwork(3)=ho ****
!
!        info(9) -- if storage is a severe problem,
!              you can save some locations by
!              restricting the maximum order maxord.
!              the default value is 5. for each
!              order decrease below 5, the code
!              requires neq fewer locations, however
!              it is likely to be slower. in any
!              case, you must have 1 .le. maxord .le. 5
!          ****  do you want the maximum order to
!                default to 5?
!                yes - set info(9)=0
!                 no - set info(9)=1
!                      and define maxord by setting
!                      iwork(3)=maxord ****
!
!        info(10) --if you know that the solutions to your equations
!               will always be nonnegative, it may help to set this
!               parameter. however, it is probably best to
!               try the code without using this option first,
!               and only to use this option if that doesn't
!               work very well.
!           ****  do you want the code to solve the problem without
!                 invoking any special nonnegativity constraints?
!                  yes - set info(10)=0
!                   no - set info(10)=1
!
!        info(11) --ddassl normally requires the initial t,
!               y, and yprime to be consistent. that is,
!               you must have g(t,y,yprime) = 0 at the initial
!               time. if you do not know the initial
!               derivative precisely, you can let ddassl try
!               to compute it.
!          ****   are the initial t, y, yprime consistent?
!                 yes - set info(11) = 0
!                  no - set info(11) = 1,
!                       and set yprime to an initial approximation
!                       to yprime.  (if you have no idea what
!                       yprime should be, set it to zero. note
!                       that the initial y should be such
!                       that there must exist a yprime so that
!                       g(t,y,yprime) = 0.)
!
!  rtol, atol -- you must assign relative (rtol) and absolute (atol
!         error tolerances to tell the code how accurately you
!         want the solution to be computed.  they must be defined
!         as variables because the code may change them.  you
!         have two choices --
!               both rtol and atol are scalars. (info(2)=0)
!               both rtol and atol are vectors. (info(2)=1)
!         in either case all components must be non-negative.
!
!         the tolerances are used by the code in a local error
!         test at each step which requires roughly that
!               abs(local error) .le. rtol*abs(y)+atol
!         for each vector component.
!         (more specifically, a root-mean-square norm is used to
!         measure the size of vectors, and the error test uses the
!         magnitude of the solution at the beginning of the step.)
!
!         the true (global) error is the difference between the
!         true solution of the initial value problem and the
!         computed approximation.  practically all present day
!         codes, including this one, control the local error at
!         each step and do not even attempt to control the global
!         error directly.
!         usually, but not always, the true accuracy of the
!         computed y is comparable to the error tolerances. this
!         code will usually, but not always, deliver a more
!         accurate solution if you reduce the tolerances and
!         integrate again.  by comparing two such solutions you
!         can get a fairly reliable idea of the true error in the
!         solution at the bigger tolerances.
!
!         setting atol=0. results in a pure relative error test on
!         that component.  setting rtol=0. results in a pure
!         absolute error test on that component.  a mixed test
!         with non-zero rtol and atol corresponds roughly to a
!         relative error test when the solution component is much
!         bigger than atol and to an absolute error test when the
!         solution component is smaller than the threshhold atol.
!
!         the code will not attempt to compute a solution at an
!         accuracy unreasonable for the machine being used.  it will
!         advise you if you ask for too much accuracy and inform
!         you as to the maximum accuracy it believes possible.
!
!  rwork(*) --  dimension this real work array of length lrw in your
!         calling program.
!
!  lrw -- set it to the declared length of the rwork array.
!               you must have
!                    lrw .ge. 40+(maxord+4)*neq+neq**2
!               for the full (dense) jacobian case (when info(6)=0), or
!                    lrw .ge. 40+(maxord+4)*neq+(2*ml+mu+1)*neq
!               for the banded user-defined jacobian case
!               (when info(5)=1 and info(6)=1), or
!                     lrw .ge. 40+(maxord+4)*neq+(2*ml+mu+1)*neq
!                           +2*(neq/(ml+mu+1)+1)
!               for the banded finite-difference-generated jacobian case
!               (when info(5)=0 and info(6)=1)
!
!  iwork(*) --  dimension this integer work array of length liw in
!         your calling program.
!
!  liw -- set it to the declared length of the iwork array.
!               you must have liw .ge. 20+neq
!
!  rpar, ipar -- these are parameter arrays, of real and integer
!         type, respectively.  you can use them for communication
!         between your program that calls ddassl and the
!         res subroutine (and the jac subroutine).  they are not
!         altered by ddassl.  if you do not need rpar or ipar,
!         ignore these parameters by treating them as dummy
!         arguments.  if you do choose to use them, dimension
!         them in your calling program and in res (and in jac)
!         as arrays of appropriate length.
!
!  jac -- if you have set info(5)=0, you can ignore this parameter
!         by treating it as a dummy argument.  otherwise, you must
!         provide a subroutine of the form
!             subroutine jac(t,y,yprime,pd,cj,rpar,ipar)
!         to define the matrix of partial derivatives
!             pd=dg/dy+cj*dg/dyprime
!         cj is a scalar which is input to jac.
!         for the given values of t,y,yprime, the
!         subroutine must evaluate the non-zero partial
!         derivatives for each equation and each solution
!         component, and store these values in the
!         matrix pd.  the elements of pd are set to zero
!         before each call to jac so only non-zero elements
!         need to be defined.
!
!         subroutine jac must not alter t,y,(*),yprime(*), or cj.
!         you must declare the name jac in an external statement in
!         your program that calls ddassl.  you must dimension y,
!         yprime and pd in jac.
!
!         the way you must store the elements into the pd matrix
!         depends on the structure of the matrix which you
!         indicated by info(6).
!               *** info(6)=0 -- full (dense) matrix ***
!                   give pd a first dimension of neq.
!                   when you evaluate the (non-zero) partial derivative
!                   of equation i with respect to variable j, you must
!                   store it in pd according to
!                   pd(i,j) = "dg(i)/dy(j)+cj*dg(i)/dyprime(j)"
!               *** info(6)=1 -- banded jacobian with ml lower and mu
!                   upper diagonal bands (refer to info(6) description
!                   of ml and mu) ***
!                   give pd a first dimension of 2*ml+mu+1.
!                   when you evaluate the (non-zero) partial derivative
!                   of equation i with respect to variable j, you must
!                   store it in pd according to
!                   irow = i - j + ml + mu + 1
!                   pd(irow,j) = "dg(i)/dy(j)+cj*dg(i)/dyprime(j)"
!
!         rpar and ipar are real and integer parameter arrays
!         which you can use for communication between your calling
!         program and your jacobian subroutine jac. they are not
!         altered by ddassl. if you do not need rpar or ipar,
!         ignore these parameters by treating them as dummy
!         arguments. if you do choose to use them, dimension
!         them in your calling program and in jac as arrays of
!         appropriate length.
!
!
!  optionally replaceable norm routine:
!
!     ddassl uses a weighted norm ddanrm to measure the size
!     of vectors such as the estimated error in each step.
!     a function subprogram
!       double precision function ddanrm(neq,v,wt,rpar,ipar)
!       dimension v(neq),wt(neq)
!     is used to define this norm. here, v is the vector
!     whose norm is to be computed, and wt is a vector of
!     weights.  a ddanrm routine has been included with ddassl
!     which computes the weighted root-mean-square norm
!     given by
!       ddanrm=sqrt((1/neq)*sum(v(i)/wt(i))**2)
!     this norm is suitable for most problems. in some
!     special cases, it may be more convenient and/or
!     efficient to define your own norm by writing a function
!     subprogram to be called instead of ddanrm. this should,
!     however, be attempted only after careful thought and
!     consideration.
!
!
!  -------- output -- after any return from ddassl ---------------------
!
!  the principal aim of the code is to return a computed solution at
!  tout, although it is also possible to obtain intermediate results
!  along the way. to find out whether the code achieved its goal
!  or if the integration process was interrupted before the task was
!  completed, you must check the idid parameter.
!
!
!  t -- the solution was successfully advanced to the
!               output value of t.
!
!  y(*) -- contains the computed solution approximation at t.
!
!  yprime(*) -- contains the computed derivative
!               approximation at t.
!
!  idid -- reports what the code did.
!
!                     *** task completed ***
!                reported by positive values of idid
!
!           idid = 1 -- a step was successfully taken in the
!                   intermediate-output mode. the code has not
!                   yet reached tout.
!
!           idid = 2 -- the integration to tstop was successfully
!                   completed (t=tstop) by stepping exactly to tstop.
!
!           idid = 3 -- the integration to tout was successfully
!                   completed (t=tout) by stepping past tout.
!                   y(*) is obtained by interpolation.
!                   yprime(*) is obtained by interpolation.
!
!                    *** task interrupted ***
!                reported by negative values of idid
!
!           idid = -1 -- a large amount of work has been expended.
!                   (about 500 steps)
!
!           idid = -2 -- the error tolerances are too stringent.
!
!           idid = -3 -- the local error test cannot be satisfied
!                   because you specified a zero component in atol
!                   and the corresponding computed solution
!                   component is zero. thus, a pure relative error
!                   test is impossible for this component.
!
!           idid = -6 -- ddassl had repeated error test
!                   failures on the last attempted step.
!
!           idid = -7 -- the corrector could not converge.
!
!           idid = -8 -- the matrix of partial derivatives
!                   is singular.
!
!           idid = -9 -- the corrector could not converge.
!                   there were repeated error test failures
!                   in this step.
!
!           idid =-10 -- the corrector could not converge
!                   because ires was equal to minus one.
!
!           idid =-11 -- ires equal to -2 was encountered
!                   and control is being returned to the
!                   calling program.
!
!           idid =-12 -- ddassl failed to compute the initial
!                   yprime.
!
!
!
!           idid = -13,..,-32 -- not applicable for this code
!
!                    *** task terminated ***
!                reported by the value of idid=-33
!
!           idid = -33 -- the code has encountered trouble from which
!                   it cannot recover. a message is printed
!                   explaining the trouble and control is returned
!                   to the calling program. for example, this occurs
!                   when invalid input is detected.
!
!  rtol, atol -- these quantities remain unchanged except when
!               idid = -2. in this case, the error tolerances have been
!               increased by the code to values which are estimated to
!               be appropriate for continuing the integration. however,
!               the reported solution at t was obtained using the input
!               values of rtol and atol.
!
!  rwork, iwork -- contain information which is usually of no
!               interest to the user but necessary for subsequent calls.
!               however, you may find use for
!
!               rwork(3)--which contains the step size h to be
!                       attempted on the next step.
!
!               rwork(4)--which contains the current value of the
!                       independent variable, i.e., the farthest point
!                       integration has reached. this will be different
!                       from t only when interpolation has been
!                       performed (idid=3).
!
!               rwork(7)--which contains the stepsize used
!                       on the last successful step.
!
!               iwork(7)--which contains the order of the method to
!                       be attempted on the next step.
!
!               iwork(8)--which contains the order of the method used
!                       on the last step.
!
!               iwork(11)--which contains the number of steps taken so
!                        far.
!
!               iwork(12)--which contains the number of calls to res
!                        so far.
!
!               iwork(13)--which contains the number of evaluations of
!                        the matrix of partial derivatives needed so
!                        far.
!
!               iwork(14)--which contains the total number
!                        of error test failures so far.
!
!               iwork(15)--which contains the total number
!                        of convergence test failures so far.
!                        (includes singular iteration matrix
!                        failures.)
!
!
!  -------- input -- what to do to continue the integration ------------
!                    (calls after the first)
!
!  this code is organized so that subsequent calls to continue the
!  integration involve little (if any) additional effort on your
!  part. you must monitor the idid parameter in order to determine
!  what to do next.
!
!  recalling that the principal task of the code is to integrate
!  from t to tout (the interval mode), usually all you will need
!  to do is specify a new tout upon reaching the current tout.
!
!  do not alter any quantity not specifically permitted below,
!  in particular do not alter neq,t,y(*),yprime(*),rwork(*),iwork(*)
!  or the differential equation in subroutine res. any such
!  alteration constitutes a new problem and must be treated as such,
!  i.e., you must start afresh.
!
!  you cannot change from vector to scalar error control or vice
!  versa (info(2)), but you can change the size of the entries of
!  rtol, atol. increasing a tolerance makes the equation easier
!  to integrate. decreasing a tolerance will make the equation
!  harder to integrate and should generally be avoided.
!
!  you can switch from the intermediate-output mode to the
!  interval mode (info(3)) or vice versa at any time.
!
!  if it has been necessary to prevent the integration from going
!  past a point tstop (info(4), rwork(1)), keep in mind that the
!  code will not integrate to any tout beyond the currently
!  specified tstop. once tstop has been reached you must change
!  the value of tstop or set info(4)=0. you may change info(4)
!  or tstop at any time but you must supply the value of tstop in
!  rwork(1) whenever you set info(4)=1.
!
!  do not change info(5), info(6), iwork(1), or iwork(2)
!  unless you are going to restart the code.
!
!                 *** following a completed task ***
!  if
!     idid = 1, call the code again to continue the integration
!                  another step in the direction of tout.
!
!     idid = 2 or 3, define a new tout and call the code again.
!                  tout must be different from t. you cannot change
!                  the direction of integration without restarting.
!
!                 *** following an interrupted task ***
!               to show the code that you realize the task was
!               interrupted and that you want to continue, you
!               must take appropriate action and set info(1) = 1
!  if
!    idid = -1, the code has taken about 500 steps.
!                  if you want to continue, set info(1) = 1 and
!                  call the code again. an additional 500 steps
!                  will be allowed.
!
!    idid = -2, the error tolerances rtol, atol have been
!                  increased to values the code estimates appropriate
!                  for continuing. you may want to change them
!                  yourself. if you are sure you want to continue
!                  with relaxed error tolerances, set info(1)=1 and
!                  call the code again.
!
!    idid = -3, a solution component is zero and you set the
!                  corresponding component of atol to zero. if you
!                  are sure you want to continue, you must first
!                  alter the error criterion to use positive values
!                  for those components of atol corresponding to zero
!                  solution components, then set info(1)=1 and call
!                  the code again.
!
!    idid = -4,-5  --- cannot occur with this code.
!
!    idid = -6, repeated error test failures occurred on the
!                  last attempted step in ddassl. a singularity in the
!                  solution may be present. if you are absolutely
!                  certain you want to continue, you should restart
!                  the integration. (provide initial values of y and
!                  yprime which are consistent)
!
!    idid = -7, repeated convergence test failures occurred
!                  on the last attempted step in ddassl. an inaccurate
!                  or ill-conditioned jacobian may be the problem. if
!                  you are absolutely certain you want to continue, you
!                  should restart the integration.
!
!    idid = -8, the matrix of partial derivatives is singular.
!                  some of your equations may be redundant.
!                  ddassl cannot solve the problem as stated.
!                  it is possible that the redundant equations
!                  could be removed, and then ddassl could
!                  solve the problem. it is also possible
!                  that a solution to your problem either
!                  does not exist or is not unique.
!
!    idid = -9, ddassl had multiple convergence test
!                  failures, preceded by multiple error
!                  test failures, on the last attempted step.
!                  it is possible that your problem
!                  is ill-posed, and cannot be solved
!                  using this code. or, there may be a
!                  discontinuity or a singularity in the
!                  solution. if you are absolutely certain
!                  you want to continue, you should restart
!                  the integration.
!
!    idid =-10, ddassl had multiple convergence test failures
!                  because ires was equal to minus one.
!                  if you are absolutely certain you want
!                  to continue, you should restart the
!                  integration.
!
!    idid =-11, ires=-2 was encountered, and control is being
!                  returned to the calling program.
!
!    idid =-12, ddassl failed to compute the initial yprime.
!                  this could happen because the initial
!                  approximation to yprime was not very good, or
!                  if a yprime consistent with the initial y
!                  does not exist. the problem could also be caused
!                  by an inaccurate or singular iteration matrix.
!
!    idid = -13,..,-32  --- cannot occur with this code.
!
!
!                 *** following a terminated task ***
!
!  if idid= -33, you cannot continue the solution of this problem.
!                  an attempt to do so will result in your
!                  run being terminated.
!
!
!  -------- error messages ---------------------------------------------
!
!      the slatec error print routine xermsg is called in the event of
!   unsuccessful completion of a task.  most of these are treated as
!   "recoverable errors", which means that (unless the user has directed
!   otherwise) control will be returned to the calling program for
!   possible action after the message has been printed.
!
!   in the event of a negative value of idid other than -33, an appro-
!   priate message is printed and the "error number" printed by xermsg
!   is the value of idid.  there are quite a number of illegal input
!   errors that can lead to a returned value idid=-33.  the conditions
!   and their printed "error numbers" are as follows:
!
!   error number       condition
!
      subroutine ddassl (eqn, neq, t, y, yprime, tout, info, rtol, atol,&
     &   idid, rwork, lrw, iwork, liw, rpar, ipar)
!        1       some element of info vector is not zero or one.
!        2       neq .le. 0
!        3       maxord not in range.
!        4       lrw is less than the required length for rwork.
!        5       liw is less than the required length for iwork.
!        6       some element of rtol is .lt. 0
!        7       some element of atol is .lt. 0
!        8       all elements of rtol and atol are zero.
!        9       info(4)=1 and tstop is behind tout.
!       10       hmax .lt. 0.0
!       11       tout is behind t.
!       12       info(8)=1 and h0=0.0
!       13       some element of wt is .le. 0.0
!       14       tout is too close to t to start integration.
!       15       info(4)=1 and tstop is behind t.
!       16       --( not used in this version )--
!       17       ml illegal.  either .lt. 0 or .gt. neq
!       18       mu illegal.  either .lt. 0 or .gt. neq
!       19       tout = t.
!
!   if ddassl is called again without any action taken to remove the
!   cause of an unsuccessful return, xermsg will be called with a fatal
!   error flag, which will cause unconditional termination of the
!   program.  there are two such fatal errors:
!
!   error number -998:  the last step was terminated with a negative
!       value of idid other than -33, and no appropriate action was
!       taken.
!
!   error number -999:  the previous call was terminated because of
!       illegal input (idid=-33) and there is illegal input in the
!       present call, as well.  (suspect infinite loop.)
!
!  ---------------------------------------------------------------------
!
!***references  a description of dassl: a differential/algebraic
!                 system solver, l. r. petzold, sand82-8637,
!                 sandia national laboratories, september 1982.
!***routines called  d1mach, ddaini, ddanrm, ddastp, ddatrp, ddawts,
!                    xermsg
!***revision history  (yymmdd)
!   830315  date written
!   880387  code changes made.  all common statements have been
!           replaced by a data statement, which defines pointers into
!           rwork, and parameter statements which define pointers
!           into iwork.  as well the documentation has gone through
!           grammatical changes.
!   881005  the prologue has been changed to mixed case.
!           the subordinate routines had revision dates changed to
!           this date, although the documentation for these routines
!           is all upper case.  no code changes.
!   890511  code changes made.  the data statement in the declaration
!           section of ddassl was replaced with a parameter
!           statement.  also the statement s = 100.d0 was removed
!           from the top of the newton iteration in ddastp.
!           the subordinate routines had revision dates changed to
!           this date.
!   890517  the revision date syntax was replaced with the revision
!           history syntax.  also the "deck" comment was added to
!           the top of all subroutines.  these changes are consistent
!           with new slatec guidelines.
!           the subordinate routines had revision dates changed to
!           this date.  no code changes.
!   891013  code changes made.
!           removed all occurrences of float or dble.  all operations
!           are now performed with "mixed-mode" arithmetic.
!           also, specific function names were replaced with generic
!           function names to be consistent with new slatec guidelines.
!           in particular:
!              replaced dsqrt with sqrt everywhere.
!              replaced dabs with abs everywhere.
!              replaced dmin1 with min everywhere.
!              replaced min0 with min everywhere.
!              replaced dmax1 with max everywhere.
!              replaced max0 with max everywhere.
!              replaced dsign with sign everywhere.
!           also replaced revision date with revision history in all
!           subordinate routines.
!   901004  miscellaneous changes to prologue to complete conversion
!           to slatec 4.0 format.  no code changes.  (f.n.fritsch)
!   901009  corrected gams classification code and converted subsidiary
!           routines to 4.0 format.  no code changes.  (f.n.fritsch)
!   901010  converted xerrwv calls to xermsg calls.  (r.clemens, afwl)
!   901019  code changes made.
!           merged slatec 4.0 changes with previous changes made
!           by c. ulrich.  below is a history of the changes made by
!           c. ulrich. (changes in subsidiary routines are implied
!           by this history)
!           891228  bug was found and repaired inside the ddassl
!                   and ddaini routines.  ddaini was incorrectly
!                   returning the initial t with y and yprime
!                   computed at t+h.  the routine now returns t+h
!                   rather than the initial t.
!                   cosmetic changes made to ddastp.
!           900904  three modifications were made to fix a bug (inside
!                   ddassl) re interpolation for continuation calls and
!                   cases where tn is very close to tstop:
!
!                   1) in testing for whether h is too large, just
!                      compare h to (tstop - tn), rather than
!                      (tstop - tn) * (1-4*uround), and set h to
!                      tstop - tn.  this will force ddastp to step
!                      exactly to tstop under certain situations
!                      (i.e. when h returned from ddastp would otherwise
!                      take tn beyond tstop).
!
!                   2) inside the ddastp loop, interpolate exactly to
!                      tstop if tn is very close to tstop (rather than
!                      interpolating to within roundoff of tstop).
!
!                   3) modified idid description for idid = 2 to say
!                      that the solution is returned by stepping exactly
!                      to tstop, rather than tout.  (in some cases the
!                      solution is actually obtained by extrapolating
!                      over a distance near unit roundoff to tstop,
!                      but this small distance is deemed acceptable in
!                      these circumstances.)
!   901026  added explicit declarations for all variables and minor
!           cosmetic changes to prologue, removed unreferenced labels,
!           and improved xermsg calls.  (fnf)
!   901030  added error messages section and reworked other sections to
!           be of more uniform format.  (fnf)
!   910624  fixed minor bug related to hmax (six lines after label
!           525).  (lrp)
!   000711  fixed tests on (tn - tout) at 420 and 440 (ach)
!***end prologue  ddassl
!
!**end
!
!     declare arguments.
!
      integer  neq, info(15), idid, lrw, iwork(*), liw, ipar(*)
      double precision                                                  &
     &   t, y(*), yprime(*), tout, rtol(*), atol(*), rwork(*),          &
     &   rpar(*)
      class(problem_ddassl) :: eqn
!
!     declare local variables.
!
      integer  i, itemp, lalpha, lbeta, lcj, lcjold, lctf, ldelta,      &
     &   leniw, lenpd, lenrw, le, letf, lgamma, lh, lhmax, lhold, lipvt,&
     &   ljcalc, lk, lkold, liwm, lml, lmtype, lmu, lmxord, lnje, lnpd, &
     &   lnre, lns, lnst, lnstl, lpd, lphase, lphi, lpsi, lround, ls,   &
     &   lsigma, ltn, ltstop, lwm, lwt, mband, msave, mxord, npd, ntemp,&
     &   nzflg
      double precision                                                  &
     &   atoli, h, hmax, hmin, ho, r, rh, rtoli, tdist, tn, tnext,      &
     &   tstop, uround, ypnorm
      logical  done
!       auxiliary variables for conversion of values to be included in
!       error messages.
      character*8  xern1, xern2
      character*16 xern3, xern4
!
!     set pointers into iwork
      parameter (lml=1, lmu=2, lmxord=3, lmtype=4, lnst=11,             &
     &  lnre=12, lnje=13, letf=14, lctf=15, lnpd=16,                    &
     &  lipvt=21, ljcalc=5, lphase=6, lk=7, lkold=8,                    &
     &  lns=9, lnstl=10, liwm=1)
!
!     set relative offset into rwork
      parameter (npd=1)
!
!     set pointers into rwork
      parameter (ltstop=1, lhmax=2, lh=3, ltn=4,                        &
     &  lcj=5, lcjold=6, lhold=7, ls=8, lround=9,                       &
     &  lalpha=11, lbeta=17, lgamma=23,                                 &
     &  lpsi=29, lsigma=35, ldelta=41)
!
!***first executable statement  ddassl
      if(info(1).ne.0)go to 100
!
!-----------------------------------------------------------------------
!     this block is executed for the initial call only.
!     it contains checking of inputs and initializations.
!-----------------------------------------------------------------------
!
!     first check info array to make sure all elements of info
!     are either zero or one.
      do 10 i=2,11
         if(info(i).ne.0.and.info(i).ne.1)go to 701
   10    continue
!
      if(neq.le.0)go to 702
!
!     check and compute maximum order
      mxord=5
      if(info(9).eq.0)go to 20
         mxord=iwork(lmxord)
         if(mxord.lt.1.or.mxord.gt.5)go to 703
   20    iwork(lmxord)=mxord
!
!     compute mtype,lenpd,lenrw.check ml and mu.
      if(info(6).ne.0)go to 40
         lenpd=neq**2
         lenrw=40+(iwork(lmxord)+4)*neq+lenpd
         if(info(5).ne.0)go to 30
            iwork(lmtype)=2
            go to 60
   30       iwork(lmtype)=1
            go to 60
   40 if(iwork(lml).lt.0.or.iwork(lml).ge.neq)go to 717
      if(iwork(lmu).lt.0.or.iwork(lmu).ge.neq)go to 718
      lenpd=(2*iwork(lml)+iwork(lmu)+1)*neq
      if(info(5).ne.0)go to 50
         iwork(lmtype)=5
         mband=iwork(lml)+iwork(lmu)+1
         msave=(neq/mband)+1
         lenrw=40+(iwork(lmxord)+4)*neq+lenpd+2*msave
         go to 60
   50    iwork(lmtype)=4
         lenrw=40+(iwork(lmxord)+4)*neq+lenpd
!
!     check lengths of rwork and iwork
   60 leniw=20+neq
      iwork(lnpd)=lenpd
      if(lrw.lt.lenrw)go to 704
      if(liw.lt.leniw)go to 705
!
!     check to see that tout is different from t
      if(tout .eq. t)go to 719
!
!     check hmax
      if(info(7).eq.0)go to 70
         hmax=rwork(lhmax)
         if(hmax.le.0.0d0)go to 710
   70 continue
!
!     initialize counters
      iwork(lnst)=0
      iwork(lnre)=0
      iwork(lnje)=0
!
      iwork(lnstl)=0
      idid=1
      go to 200
!
!-----------------------------------------------------------------------
!     this block is for continuation calls
!     only. here we check info(1), and if the
!     last step was interrupted we check whether
!     appropriate action was taken.
!-----------------------------------------------------------------------
!
  100 continue
      if(info(1).eq.1)go to 110
      if(info(1).ne.-1)go to 701
!
!     if we are here, the last step was interrupted
!     by an error condition from ddastp, and
!     appropriate action was not taken. this
!     is a fatal error.
      write (xern1, '(i8)') idid
      call xermsg ('slatec', 'ddassl',                                  &
     &   'the last step terminated with a negative value of idid = ' // &
     &   xern1 // ' and no appropriate action was taken.  ' //          &
     &   'run terminated', -998, 2)
      return
  110 continue
      iwork(lnstl)=iwork(lnst)
!
!-----------------------------------------------------------------------
!     this block is executed on all calls.
!     the error tolerance parameters are
!     checked, and the work array pointers
!     are set.
!-----------------------------------------------------------------------
!
  200 continue
!     check rtol,atol
      nzflg=0
      rtoli=rtol(1)
      atoli=atol(1)
      do 210 i=1,neq
         if(info(2).eq.1)rtoli=rtol(i)
         if(info(2).eq.1)atoli=atol(i)
         if(rtoli.gt.0.0d0.or.atoli.gt.0.0d0)nzflg=1
         if(rtoli.lt.0.0d0)go to 706
         if(atoli.lt.0.0d0)go to 707
  210    continue
      if(nzflg.eq.0)go to 708
!
!     set up rwork storage.iwork storage is fixed
!     in data statement.
      le=ldelta+neq
      lwt=le+neq
      lphi=lwt+neq
      lpd=lphi+(iwork(lmxord)+1)*neq
      lwm=lpd
      ntemp=npd+iwork(lnpd)
      if(info(1).eq.1)go to 400
!
!-----------------------------------------------------------------------
!     this block is executed on the initial call
!     only. set the initial step size, and
!     the error weight vector, and phi.
!     compute initial yprime, if necessary.
!-----------------------------------------------------------------------
!
      tn=t
      idid=1
!
!     set error weight vector wt
      call ddawts(neq,info(2),rtol,atol,y,rwork(lwt),rpar,ipar)
      do 305 i = 1,neq
         if(rwork(lwt+i-1).le.0.0d0) go to 713
  305    continue
!
!     compute unit roundoff and hmin
      uround = d1mach(4)
      ! uround = epsilon(uround)
      rwork(lround) = uround
      hmin = 4.0d0*uround*max(abs(t),abs(tout))
!
!     check initial interval to see that it is long enough
      tdist = abs(tout - t)
      if(tdist .lt. hmin) go to 714
!
!     check ho, if this was input
      if (info(8) .eq. 0) go to 310
         ho = rwork(lh)
         if ((tout - t)*ho .lt. 0.0d0) go to 711
         if (ho .eq. 0.0d0) go to 712
         go to 320
  310  continue
!
!     compute initial stepsize, to be used by either
!     ddastp or ddaini, depending on info(11)
      ho = 0.001d0*tdist
      ypnorm = ddanrm(neq,yprime,rwork(lwt),rpar,ipar)
      if (ypnorm .gt. 0.5d0/ho) ho = 0.5d0/ypnorm
      ho = sign(ho,tout-t)
!     adjust ho if necessary to meet hmax bound
  320 if (info(7) .eq. 0) go to 330
         rh = abs(ho)/rwork(lhmax)
         if (rh .gt. 1.0d0) ho = ho/rh
!     compute tstop, if applicable
  330 if (info(4) .eq. 0) go to 340
         tstop = rwork(ltstop)
         if ((tstop - t)*ho .lt. 0.0d0) go to 715
         if ((t + ho - tstop)*ho .gt. 0.0d0) ho = tstop - t
         if ((tstop - tout)*ho .lt. 0.0d0) go to 709
!
!     compute initial derivative, updating tn and y, if applicable
  340 if (info(11) .eq. 0) go to 350
      call eqn%ddaini(tn,y,yprime,neq,                                      &
     &  ho,rwork(lwt),idid,rpar,ipar,                           &
     &  rwork(lphi),rwork(ldelta),rwork(le),                            &
     &  rwork(lwm),iwork(liwm),hmin,rwork(lround),                      &
     &  info(10),ntemp)
      if (idid .lt. 0) go to 390
!
!     load h with ho.  store h in rwork(lh)
  350 h = ho
      rwork(lh) = h
!
!     load y and h*yprime into phi(*,1) and phi(*,2)
      itemp = lphi + neq
      do 370 i = 1,neq
         rwork(lphi + i - 1) = y(i)
  370    rwork(itemp + i - 1) = h*yprime(i)
!
  390 go to 500
!
!-------------------------------------------------------
!     this block is for continuation calls only. its
!     purpose is to check stop conditions before
!     taking a step.
!     adjust h if necessary to meet hmax bound
!-------------------------------------------------------
!
  400 continue
      uround=rwork(lround)
      done = .false.
      tn=rwork(ltn)
      h=rwork(lh)
      if(info(7) .eq. 0) go to 410
         rh = abs(h)/rwork(lhmax)
         if(rh .gt. 1.0d0) h = h/rh
  410 continue
      if(t .eq. tout) go to 719
      if((t - tout)*h .gt. 0.0d0) go to 711
      if(info(4) .eq. 1) go to 430
      if(info(3) .eq. 1) go to 420
      if((tn-tout)*h.lt.0.0d0)go to 490
      call ddatrp(tn,tout,y,yprime,neq,iwork(lkold),                    &
     &  rwork(lphi),rwork(lpsi))
      t=tout
      idid = 3
      done = .true.
      go to 490
  420 if((tn-t)*h .le. 0.0d0) go to 490
      if((tn - tout)*h .ge. 0.0d0) go to 425
      call ddatrp(tn,tn,y,yprime,neq,iwork(lkold),                      &
     &  rwork(lphi),rwork(lpsi))
      t = tn
      idid = 1
      done = .true.
      go to 490
  425 continue
      call ddatrp(tn,tout,y,yprime,neq,iwork(lkold),                    &
     &  rwork(lphi),rwork(lpsi))
      t = tout
      idid = 3
      done = .true.
      go to 490
  430 if(info(3) .eq. 1) go to 440
      tstop=rwork(ltstop)
      if((tn-tstop)*h.gt.0.0d0) go to 715
      if((tstop-tout)*h.lt.0.0d0)go to 709
      if((tn-tout)*h.lt.0.0d0)go to 450
      call ddatrp(tn,tout,y,yprime,neq,iwork(lkold),                    &
     &   rwork(lphi),rwork(lpsi))
      t=tout
      idid = 3
      done = .true.
      go to 490
  440 tstop = rwork(ltstop)
      if((tn-tstop)*h .gt. 0.0d0) go to 715
      if((tstop-tout)*h .lt. 0.0d0) go to 709
      if((tn-t)*h .le. 0.0d0) go to 450
      if((tn - tout)*h .ge. 0.0d0) go to 445
      call ddatrp(tn,tn,y,yprime,neq,iwork(lkold),                      &
     &  rwork(lphi),rwork(lpsi))
      t = tn
      idid = 1
      done = .true.
      go to 490
  445 continue
      call ddatrp(tn,tout,y,yprime,neq,iwork(lkold),                    &
     &  rwork(lphi),rwork(lpsi))
      t = tout
      idid = 3
      done = .true.
      go to 490
  450 continue
!     check whether we are within roundoff of tstop
      if(abs(tn-tstop).gt.100.0d0*uround*                               &
     &   (abs(tn)+abs(h)))go to 460
      call ddatrp(tn,tstop,y,yprime,neq,iwork(lkold),                   &
     &  rwork(lphi),rwork(lpsi))
      idid=2
      t=tstop
      done = .true.
      go to 490
  460 tnext=tn+h
      if((tnext-tstop)*h.le.0.0d0)go to 490
      h=tstop-tn
      rwork(lh)=h
!
  490 if (done) go to 580
!
!-------------------------------------------------------
!     the next block contains the call to the
!     one-step integrator ddastp.
!     this is a looping point for the integration steps.
!     check for too many steps.
!     update wt.
!     check for too much accuracy requested.
!     compute minimum stepsize.
!-------------------------------------------------------
!
  500 continue
!     check for failure to compute initial yprime
      if (idid .eq. -12) go to 527
!
!     check for too many steps
      if((iwork(lnst)-iwork(lnstl)).lt.500)                             &
     &   go to 510
           idid=-1
           go to 527
!
!     update wt
  510 call ddawts(neq,info(2),rtol,atol,rwork(lphi),                    &
     &  rwork(lwt),rpar,ipar)
      do 520 i=1,neq
         if(rwork(i+lwt-1).gt.0.0d0)go to 520
           idid=-3
           go to 527
  520 end do
!
!     test for too much accuracy requested.
      r=ddanrm(neq,rwork(lphi),rwork(lwt),rpar,ipar)*                   &
     &   100.0d0*uround
      if(r.le.1.0d0)go to 525
!     multiply rtol and atol by r and return
      if(info(2).eq.1)go to 523
           rtol(1)=r*rtol(1)
           atol(1)=r*atol(1)
           idid=-2
           go to 527
  523 do 524 i=1,neq
           rtol(i)=r*rtol(i)
  524      atol(i)=r*atol(i)
      idid=-2
      go to 527
  525 continue
!
!     compute minimum stepsize
      hmin=4.0d0*uround*max(abs(tn),abs(tout))
!
!     test h vs. hmax
      if (info(7) .ne. 0) then
         rh = abs(h)/rwork(lhmax)
         if (rh .gt. 1.0d0) h = h/rh
      endif
!
      call eqn%ddastp(tn,y,yprime,neq,                                      &
     &   h,rwork(lwt),info(1),idid,rpar,ipar,                   &
     &   rwork(lphi),rwork(ldelta),rwork(le),                           &
     &   rwork(lwm),iwork(liwm),                                        &
     &   rwork(lalpha),rwork(lbeta),rwork(lgamma),                      &
     &   rwork(lpsi),rwork(lsigma),                                     &
     &   rwork(lcj),rwork(lcjold),rwork(lhold),                         &
     &   rwork(ls),hmin,rwork(lround),                                  &
     &   iwork(lphase),iwork(ljcalc),iwork(lk),                         &
     &   iwork(lkold),iwork(lns),info(10),ntemp)
  527 if(idid.lt.0)go to 600
!
!--------------------------------------------------------
!     this block handles the case of a successful return
!     from ddastp (idid=1).  test for stop conditions.
!--------------------------------------------------------
!
      if(info(4).ne.0)go to 540
           if(info(3).ne.0)go to 530
             if((tn-tout)*h.lt.0.0d0)go to 500
             call ddatrp(tn,tout,y,yprime,neq,                          &
     &         iwork(lkold),rwork(lphi),rwork(lpsi))
             idid=3
             t=tout
             go to 580
  530        if((tn-tout)*h.ge.0.0d0)go to 535
             t=tn
             idid=1
             go to 580
  535        call ddatrp(tn,tout,y,yprime,neq,                          &
     &         iwork(lkold),rwork(lphi),rwork(lpsi))
             idid=3
             t=tout
             go to 580
  540 if(info(3).ne.0)go to 550
      if((tn-tout)*h.lt.0.0d0)go to 542
         call ddatrp(tn,tout,y,yprime,neq,                              &
     &     iwork(lkold),rwork(lphi),rwork(lpsi))
         t=tout
         idid=3
         go to 580
  542 if(abs(tn-tstop).le.100.0d0*uround*                               &
     &   (abs(tn)+abs(h)))go to 545
      tnext=tn+h
      if((tnext-tstop)*h.le.0.0d0)go to 500
      h=tstop-tn
      go to 500
  545 call ddatrp(tn,tstop,y,yprime,neq,                                &
     &  iwork(lkold),rwork(lphi),rwork(lpsi))
      idid=2
      t=tstop
      go to 580
  550 if((tn-tout)*h.ge.0.0d0)go to 555
      if(abs(tn-tstop).le.100.0d0*uround*(abs(tn)+abs(h)))go to 552
      t=tn
      idid=1
      go to 580
  552 call ddatrp(tn,tstop,y,yprime,neq,                                &
     &  iwork(lkold),rwork(lphi),rwork(lpsi))
      idid=2
      t=tstop
      go to 580
  555 call ddatrp(tn,tout,y,yprime,neq,                                 &
     &   iwork(lkold),rwork(lphi),rwork(lpsi))
      t=tout
      idid=3
      go to 580
!
!--------------------------------------------------------
!     all successful returns from ddassl are made from
!     this block.
!--------------------------------------------------------
!
  580 continue
      rwork(ltn)=tn
      rwork(lh)=h
      return
!
!-----------------------------------------------------------------------
!     this block handles all unsuccessful
!     returns other than for illegal input.
!-----------------------------------------------------------------------
!
  600 continue
      itemp=-idid
      go to (610,620,630,690,690,640,650,660,670,675,                   &
     &  680,685), itemp
!
!     the maximum number of steps was taken before
!     reaching tout
  610 write (xern3, '(1p,d15.6)') tn
      call xermsg ('slatec', 'ddassl',                                  &
     &   'at current t = ' // xern3 // ' 500 steps taken on this ' //   &
     &   'call before reaching tout', idid, 1)
      go to 690
!
!     too much accuracy for machine precision
  620 write (xern3, '(1p,d15.6)') tn
      call xermsg ('slatec', 'ddassl',                                  &
     &   'at t = ' // xern3 // ' too much accuracy requested for ' //   &
     &   'precision of machine. rtol and atol were increased to ' //    &
     &   'appropriate values', idid, 1)
      go to 690
!
!     wt(i) .le. 0.0 for some i (not at start of problem)
  630 write (xern3, '(1p,d15.6)') tn
      call xermsg ('slatec', 'ddassl',                                  &
     &   'at t = ' // xern3 // ' some element of wt has become .le. ' //&
     &   '0.0', idid, 1)
      go to 690
!
!     error test failed repeatedly or with h=hmin
  640 write (xern3, '(1p,d15.6)') tn
      write (xern4, '(1p,d15.6)') h
      call xermsg ('slatec', 'ddassl',                                  &
     &   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //         &
     &   ' the error test failed repeatedly or with abs(h)=hmin',       &
     &   idid, 1)
      go to 690
!
!     corrector convergence failed repeatedly or with h=hmin
  650 write (xern3, '(1p,d15.6)') tn
      write (xern4, '(1p,d15.6)') h
      call xermsg ('slatec', 'ddassl',                                  &
     &   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //         &
     &   ' the corrector failed to converge repeatedly or with ' //     &
     &   'abs(h)=hmin', idid, 1)
      go to 690
!
!     the iteration matrix is singular
  660 write (xern3, '(1p,d15.6)') tn
      write (xern4, '(1p,d15.6)') h
      call xermsg ('slatec', 'ddassl',                                  &
     &   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //         &
     &   ' the iteration matrix is singular', idid, 1)
      go to 690
!
!     corrector failure preceded by error test failures.
  670 write (xern3, '(1p,d15.6)') tn
      write (xern4, '(1p,d15.6)') h
      call xermsg ('slatec', 'ddassl',                                  &
     &   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //         &
     &   ' the corrector could not converge.  also, the error test ' // &
     &   'failed repeatedly.', idid, 1)
      go to 690
!
!     corrector failure because ires = -1
  675 write (xern3, '(1p,d15.6)') tn
      write (xern4, '(1p,d15.6)') h
      call xermsg ('slatec', 'ddassl',                                  &
     &   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //         &
     &   ' the corrector could not converge because ires was equal ' // &
     &   'to minus one', idid, 1)
      go to 690
!
!     failure because ires = -2
  680 write (xern3, '(1p,d15.6)') tn
      write (xern4, '(1p,d15.6)') h
      call xermsg ('slatec', 'ddassl',                                  &
     &   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //         &
     &   ' ires was equal to minus two', idid, 1)
      go to 690
!
!     failed to compute initial yprime
  685 write (xern3, '(1p,d15.6)') tn
      write (xern4, '(1p,d15.6)') ho
      call xermsg ('slatec', 'ddassl',                                  &
     &   'at t = ' // xern3 // ' and stepsize h = ' // xern4 //         &
     &   ' the initial yprime could not be computed', idid, 1)
      go to 690
!
  690 continue
      info(1)=-1
      t=tn
      rwork(ltn)=tn
      rwork(lh)=h
      return
!
!-----------------------------------------------------------------------
!     this block handles all error returns due
!     to illegal input, as detected before calling
!     ddastp. first the error message routine is
!     called. if this happens twice in
!     succession, execution is terminated
!
!-----------------------------------------------------------------------
  701 call xermsg ('slatec', 'ddassl',                                  &
     &   'some element of info vector is not zero or one', 1, 1)
      go to 750
!
  702 write (xern1, '(i8)') neq
      call xermsg ('slatec', 'ddassl',                                  &
     &   'neq = ' // xern1 // ' .le. 0', 2, 1)
      go to 750
!
  703 write (xern1, '(i8)') mxord
      call xermsg ('slatec', 'ddassl',                                  &
     &   'maxord = ' // xern1 // ' not in range', 3, 1)
      go to 750
!
  704 write (xern1, '(i8)') lenrw
      write (xern2, '(i8)') lrw
      call xermsg ('slatec', 'ddassl',                                  &
     &   'rwork length needed, lenrw = ' // xern1 //                    &
     &   ', exceeds lrw = ' // xern2, 4, 1)
      go to 750
!
  705 write (xern1, '(i8)') leniw
      write (xern2, '(i8)') liw
      call xermsg ('slatec', 'ddassl',                                  &
     &   'iwork length needed, leniw = ' // xern1 //                    &
     &   ', exceeds liw = ' // xern2, 5, 1)
      go to 750
!
  706 call xermsg ('slatec', 'ddassl',                                  &
     &   'some element of rtol is .lt. 0', 6, 1)
      go to 750
!
  707 call xermsg ('slatec', 'ddassl',                                  &
     &   'some element of atol is .lt. 0', 7, 1)
      go to 750
!
  708 call xermsg ('slatec', 'ddassl',                                  &
     &   'all elements of rtol and atol are zero', 8, 1)
      go to 750
!
  709 write (xern3, '(1p,d15.6)') tstop
      write (xern4, '(1p,d15.6)') tout
      call xermsg ('slatec', 'ddassl',                                  &
     &   'info(4) = 1 and tstop = ' // xern3 // ' behind tout = ' //    &
     &   xern4, 9, 1)
      go to 750
!
  710 write (xern3, '(1p,d15.6)') hmax
      call xermsg ('slatec', 'ddassl',                                  &
     &   'hmax = ' // xern3 // ' .lt. 0.0', 10, 1)
      go to 750
!
  711 write (xern3, '(1p,d15.6)') tout
      write (xern4, '(1p,d15.6)') t
      call xermsg ('slatec', 'ddassl',                                  &
     &   'tout = ' // xern3 // ' behind t = ' // xern4, 11, 1)
      go to 750
!
  712 call xermsg ('slatec', 'ddassl',                                  &
     &   'info(8)=1 and h0=0.0', 12, 1)
      go to 750
!
  713 call xermsg ('slatec', 'ddassl',                                  &
     &   'some element of wt is .le. 0.0', 13, 1)
      go to 750
!
  714 write (xern3, '(1p,d15.6)') tout
      write (xern4, '(1p,d15.6)') t
      call xermsg ('slatec', 'ddassl',                                  &
     &   'tout = ' // xern3 // ' too close to t = ' // xern4 //         &
     &   ' to start integration', 14, 1)
      go to 750
!
  715 write (xern3, '(1p,d15.6)') tstop
      write (xern4, '(1p,d15.6)') t
      call xermsg ('slatec', 'ddassl',                                  &
     &   'info(4)=1 and tstop = ' // xern3 // ' behind t = ' // xern4,  &
     &   15, 1)
      go to 750
!
  717 write (xern1, '(i8)') iwork(lml)
      call xermsg ('slatec', 'ddassl',                                  &
     &   'ml = ' // xern1 // ' illegal.  either .lt. 0 or .gt. neq',    &
     &   17, 1)
      go to 750
!
  718 write (xern1, '(i8)') iwork(lmu)
      call xermsg ('slatec', 'ddassl',                                  &
     &   'mu = ' // xern1 // ' illegal.  either .lt. 0 or .gt. neq',    &
     &   18, 1)
      go to 750
!
  719 write (xern3, '(1p,d15.6)') tout
      call xermsg ('slatec', 'ddassl',                                  &
     &  'tout = t = ' // xern3, 19, 1)
      go to 750
!
  750 idid=-33
      if(info(1).eq.-1) then
         call xermsg ('slatec', 'ddassl',                               &
     &      'repeated occurrences of illegal input$$' //                &
     &      'run terminated. apparent infinite loop', -999, 2)
      endif
!
      info(1)=-1
      return
!-----------end of subroutine ddassl------------------------------------
      end subroutine ddassl
!deck ddaini
      subroutine ddaini (eqn, x, y, yprime, neq, h, wt, idid, rpar,&
     &   ipar, phi, delta, e, wm, iwm, hmin, uround, nonneg, ntemp)
!***begin prologue  ddaini
!***subsidiary
!***purpose  initialization routine for ddassl.
!***library   slatec (dassl)
!***type      double precision (sdaini-s, ddaini-d)
!***author  petzold, linda r., (llnl)
!***description
!-----------------------------------------------------------------
!     ddaini takes one step of size h or smaller
!     with the backward euler method, to
!     find yprime.  x and y are updated to be consistent with the
!     new step.  a modified damped newton iteration is used to
!     solve the corrector iteration.
!
!     the initial guess for yprime is used in the
!     prediction, and in forming the iteration
!     matrix, but is not involved in the
!     error test. this may have trouble
!     converging if the initial guess is no
!     good, or if g(x,y,yprime) depends
!     nonlinearly on yprime.
!
!     the parameters represent:
!     x --         independent variable
!     y --         solution vector at x
!     yprime --    derivative of solution vector
!     neq --       number of equations
!     h --         stepsize. imder may use a stepsize
!                  smaller than h.
!     wt --        vector of weights for error
!                  criterion
!     idid --      completion code with the following meanings
!                  idid= 1 -- yprime was found successfully
!                  idid=-12 -- ddaini failed to find yprime
!     rpar,ipar -- real and integer parameter arrays
!                  that are not altered by ddaini
!     phi --       work space for ddaini
!     delta,e --   work space for ddaini
!     wm,iwm --    real and integer arrays storing
!                  matrix information
!
!-----------------------------------------------------------------
!***routines called  ddajac, ddanrm, ddaslv
!***revision history  (yymmdd)
!   830315  date written
!   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
!   901019  merged changes made by c. ulrich with slatec 4.0 format.
!   901026  added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (fnf)
!   901030  minor corrections to declarations.  (fnf)
!***end prologue  ddaini
!
      integer  neq, idid, ipar(*), iwm(*), nonneg, ntemp
      double precision                                                  &
     &   x, y(*), yprime(*), h, wt(*), rpar(*), phi(neq,*), delta(*),   &
     &   e(*), wm(*), hmin, uround
      class(problem_ddassl) :: eqn

      integer  i, ier, ires, jcalc, lnje, lnre, m, maxit, mjac, ncf,    &
     &   nef, nsf
      double precision                                                  &
     &   cj, damp, delnrm, err, oldnrm, r, rate, s, xold, ynorm
      logical  convgd
!
      parameter (lnre=12)
      parameter (lnje=13)
!
      data maxit/10/,mjac/5/
      data damp/0.75d0/
!
!
!---------------------------------------------------
!     block 1.
!     initializations.
!---------------------------------------------------
!
!***first executable statement  ddaini
      idid=1
      nef=0
      ncf=0
      nsf=0
      xold=x
      ynorm=ddanrm(neq,y,wt,rpar,ipar)
!
!     save y and yprime in phi
      do 100 i=1,neq
         phi(i,1)=y(i)
  100    phi(i,2)=yprime(i)
!
!
!----------------------------------------------------
!     block 2.
!     do one backward euler step.
!----------------------------------------------------
!
!     set up for start of corrector iteration
  200 cj=1.0d0/h
      x=x+h
!
!     predict solution and derivative
      do 250 i=1,neq
  250   y(i)=y(i)+h*yprime(i)
!
      jcalc=-1
      m=0
      convgd=.true.
!
!
!     corrector loop.
  300 iwm(lnre)=iwm(lnre)+1
      ires=0
!
      call eqn%res(x,y,yprime,delta,ires,rpar,ipar)
      if (ires.lt.0) go to 430
!
!
!     evaluate the iteration matrix
      if (jcalc.ne.-1) go to 310
      iwm(lnje)=iwm(lnje)+1
      jcalc=0
      call eqn%ddajac(neq,x,y,yprime,delta,cj,h,                            &
     &   ier,wt,e,wm,iwm,ires,                                      &
     &   uround,rpar,ipar,ntemp)
!
      s=1000000.d0
      if (ires.lt.0) go to 430
      if (ier.ne.0) go to 430
      nsf=0
!
!
!
!     multiply residual by damping factor
  310 continue
      do 320 i=1,neq
  320    delta(i)=delta(i)*damp
!
!     compute a new iterate (back substitution)
!     store the correction in delta
!
      call ddaslv(neq,delta,wm,iwm)
!
!     update y and yprime
      do 330 i=1,neq
         y(i)=y(i)-delta(i)
  330    yprime(i)=yprime(i)-cj*delta(i)
!
!     test for convergence of the iteration.
!
      delnrm=ddanrm(neq,delta,wt,rpar,ipar)
      if (delnrm.le.100.d0*uround*ynorm)                                &
     &   go to 400
!
      if (m.gt.0) go to 340
         oldnrm=delnrm
         go to 350
!
  340 rate=(delnrm/oldnrm)**(1.0d0/m)
      if (rate.gt.0.90d0) go to 430
      s=rate/(1.0d0-rate)
!
  350 if (s*delnrm .le. 0.33d0) go to 400
!
!
!     the corrector has not yet converged. update
!     m and and test whether the maximum
!     number of iterations have been tried.
!     every mjac iterations, get a new
!     iteration matrix.
!
      m=m+1
      if (m.ge.maxit) go to 430
!
      if ((m/mjac)*mjac.eq.m) jcalc=-1
      go to 300
!
!
!     the iteration has converged.
!     check nonnegativity constraints
  400 if (nonneg.eq.0) go to 450
      do 410 i=1,neq
  410    delta(i)=min(y(i),0.0d0)
!
      delnrm=ddanrm(neq,delta,wt,rpar,ipar)
      if (delnrm.gt.0.33d0) go to 430
!
      do 420 i=1,neq
         y(i)=y(i)-delta(i)
  420    yprime(i)=yprime(i)-cj*delta(i)
      go to 450
!
!
!     exits from corrector loop.
  430 convgd=.false.
  450 if (.not.convgd) go to 600
!
!
!
!-----------------------------------------------------
!     block 3.
!     the corrector iteration converged.
!     do error test.
!-----------------------------------------------------
!
      do 510 i=1,neq
  510    e(i)=y(i)-phi(i,1)
      err=ddanrm(neq,e,wt,rpar,ipar)
!
      if (err.le.1.0d0) return
!
!
!
!--------------------------------------------------------
!     block 4.
!     the backward euler step failed. restore x, y
!     and yprime to their original values.
!     reduce stepsize and try again, if
!     possible.
!---------------------------------------------------------
!
  600 continue
      x = xold
      do 610 i=1,neq
         y(i)=phi(i,1)
  610    yprime(i)=phi(i,2)
!
      if (convgd) go to 640
      if (ier.eq.0) go to 620
         nsf=nsf+1
         h=h*0.25d0
         if (nsf.lt.3.and.abs(h).ge.hmin) go to 690
         idid=-12
         return
  620 if (ires.gt.-2) go to 630
         idid=-12
         return
  630 ncf=ncf+1
      h=h*0.25d0
      if (ncf.lt.10.and.abs(h).ge.hmin) go to 690
         idid=-12
         return
!
  640 nef=nef+1
      r=0.90d0/(2.0d0*err+0.0001d0)
      r=max(0.1d0,min(0.5d0,r))
      h=h*r
      if (abs(h).ge.hmin.and.nef.lt.10) go to 690
         idid=-12
         return
  690    go to 200
!
!-------------end of subroutine ddaini----------------------
      end subroutine ddaini
!deck ddajac
      subroutine ddajac (eqn, neq, x, y, yprime, delta, cj, h, ier, wt, e,   &
     &   wm, iwm, ires, uround, rpar, ipar, ntemp)
!***begin prologue  ddajac
!***subsidiary
!***purpose  compute the iteration matrix for ddassl and form the
!            lu-decomposition.
!***library   slatec (dassl)
!***type      double precision (sdajac-s, ddajac-d)
!***author  petzold, linda r., (llnl)
!***description
!-----------------------------------------------------------------------
!     this routine computes the iteration matrix
!     pd=dg/dy+cj*dg/dyprime (where g(x,y,yprime)=0).
!     here pd is computed by the user-supplied
!     routine jac if iwm(mtype) is 1 or 4, and
!     it is computed by numerical finite differencing
!     if iwm(mtype)is 2 or 5
!     the parameters have the following meanings.
!     y        = array containing predicted values
!     yprime   = array containing predicted derivatives
!     delta    = residual evaluated at (x,y,yprime)
!                (used only if iwm(mtype)=2 or 5)
!     cj       = scalar parameter defining iteration matrix
!     h        = current stepsize in integration
!     ier      = variable which is .ne. 0
!                if iteration matrix is singular,
!                and 0 otherwise.
!     wt       = vector of weights for computing norms
!     e        = work space (temporary) of length neq
!     wm       = real work space for matrices. on
!                output it contains the lu decomposition
!                of the iteration matrix.
!     iwm      = integer work space containing
!                matrix information
!     res      = name of the external user-supplied routine
!                to evaluate the residual function g(x,y,yprime)
!     ires     = flag which is equal to zero if no illegal values
!                in res, and less than zero otherwise.  (if ires
!                is less than zero, the matrix was not completed)
!                in this case (if ires .lt. 0), then ier = 0.
!     uround   = the unit roundoff error of the machine being used.
!     jac      = name of the external user-supplied routine
!                to evaluate the iteration matrix (this routine
!                is only used if iwm(mtype) is 1 or 4)
!-----------------------------------------------------------------------
!***routines called  dgbfa, dgefa
!***revision history  (yymmdd)
!   830315  date written
!   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
!   901010  modified three max calls to be all on one line.  (fnf)
!   901019  merged changes made by c. ulrich with slatec 4.0 format.
!   901026  added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (fnf)
!   901101  corrected purpose.  (fnf)
!***end prologue  ddajac
!
      integer  neq, ier, iwm(*), ires, ipar(*), ntemp
      double precision                                                  &
     &   x, y(*), yprime(*), delta(*), cj, h, wt(*), e(*), wm(*),       &
     &   uround, rpar(*)
      class(problem_ddassl) :: eqn
!

      integer  i, i1, i2, ii, ipsave, isave, j, k, l, lenpd, lipvt,     &
     &   lml, lmtype, lmu, mba, mband, meb1, meband, msave, mtype, n,   &
     &   npd, npdm1, nrow
      double precision  del, delinv, squr, ypsave, ysave
!
      parameter (npd=1)
      parameter (lml=1)
      parameter (lmu=2)
      parameter (lmtype=4)
      parameter (lipvt=21)
!
!***first executable statement  ddajac
      ier = 0
      npdm1=npd-1
      mtype=iwm(lmtype)
      go to (100,200,300,400,500),mtype
!
!
!     dense user-supplied matrix
  100 lenpd=neq*neq
      do i=1,lenpd
         wm(npdm1+i)=0.0d0
      end do
      call eqn%jac(x,y,yprime,wm(npd),cj,rpar,ipar)
      go to 230
!
!
!     dense finite-difference-generated matrix
  200 ires=0
      nrow=npdm1
      squr = sqrt(uround)
      do i=1,neq
         del=squr*max(abs(y(i)),abs(h*yprime(i)),abs(wt(i)))
         del=sign(del,h*yprime(i))
         del=(y(i)+del)-y(i)
         ysave=y(i)
         ypsave=yprime(i)
         y(i)=y(i)+del
         yprime(i)=yprime(i)+cj*del
         call eqn%res(x,y,yprime,e,ires,rpar,ipar)
         if (ires .lt. 0) return
         delinv=1.0d0/del
         do l=1,neq
            wm(nrow+l)=(e(l)-delta(l))*delinv
         end do
         nrow=nrow+neq
         y(i)=ysave
         yprime(i)=ypsave
      end do

!
!     do dense-matrix lu decomposition on pd
  230    call dgefa(wm(npd),neq,neq,iwm(lipvt),ier)
      return
!
!
!     dummy section for iwm(mtype)=3
  300 return
!
!
!     banded user-supplied matrix
  400 lenpd=(2*iwm(lml)+iwm(lmu)+1)*neq
      do 410 i=1,lenpd
  410    wm(npdm1+i)=0.0d0
      call eqn%jac(x,y,yprime,wm(npd),cj,rpar,ipar)
      meband=2*iwm(lml)+iwm(lmu)+1
      go to 550
!
!
!     banded finite-difference-generated matrix
  500 mband=iwm(lml)+iwm(lmu)+1
      mba=min(mband,neq)
      meband=mband+iwm(lml)
      meb1=meband-1
      msave=(neq/mband)+1
      isave=ntemp-1
      ipsave=isave+msave
      ires=0
      squr=sqrt(uround)
      do 540 j=1,mba
         do 510 n=j,neq,mband
          k= (n-j)/mband + 1
          wm(isave+k)=y(n)
          wm(ipsave+k)=yprime(n)
          del=squr*max(abs(y(n)),abs(h*yprime(n)),abs(wt(n)))
          del=sign(del,h*yprime(n))
          del=(y(n)+del)-y(n)
          y(n)=y(n)+del
  510     yprime(n)=yprime(n)+cj*del
      call eqn%res(x,y,yprime,e,ires,rpar,ipar)
      if (ires .lt. 0) return
      do 530 n=j,neq,mband
          k= (n-j)/mband + 1
          y(n)=wm(isave+k)
          yprime(n)=wm(ipsave+k)
          del=squr*max(abs(y(n)),abs(h*yprime(n)),abs(wt(n)))
          del=sign(del,h*yprime(n))
          del=(y(n)+del)-y(n)
          delinv=1.0d0/del
          i1=max(1,(n-iwm(lmu)))
          i2=min(neq,(n+iwm(lml)))
          ii=n*meb1-iwm(lml)+npdm1
          do 520 i=i1,i2
  520       wm(ii+i)=(e(i)-delta(i))*delinv
  530    continue
  540 end do
!
!
!     do lu decomposition of banded pd
  550 call dgbfa(wm(npd),meband,neq,                                    &
     &    iwm(lml),iwm(lmu),iwm(lipvt),ier)
      return
!------end of subroutine ddajac------
      end subroutine ddajac
!deck ddanrm
      double precision function ddanrm (neq, v, wt, rpar, ipar)
!***begin prologue  ddanrm
!***subsidiary
!***purpose  compute vector norm for ddassl.
!***library   slatec (dassl)
!***type      double precision (sdanrm-s, ddanrm-d)
!***author  petzold, linda r., (llnl)
!***description
!-----------------------------------------------------------------------
!     this function routine computes the weighted
!     root-mean-square norm of the vector of length
!     neq contained in the array v,with weights
!     contained in the array wt of length neq.
!        ddanrm=sqrt((1/neq)*sum(v(i)/wt(i))**2)
!-----------------------------------------------------------------------
!***routines called  (none)
!***revision history  (yymmdd)
!   830315  date written
!   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
!   901019  merged changes made by c. ulrich with slatec 4.0 format.
!   901026  added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (fnf)
!***end prologue  ddanrm
!
      integer  neq, ipar(*)
      double precision  v(neq), wt(neq), rpar(*)
!
      integer  i
      double precision  sum, vmax
!
!***first executable statement  ddanrm
      ddanrm = 0.0d0
      vmax = 0.0d0
      do 10 i = 1,neq
        if(abs(v(i)/wt(i)) .gt. vmax) vmax = abs(v(i)/wt(i))
   10   continue
      if(vmax .le. 0.0d0) go to 30
      sum = 0.0d0
      do 20 i = 1,neq
   20   sum = sum + ((v(i)/wt(i))/vmax)**2
      ddanrm = vmax*sqrt(sum/neq)
   30 continue
      return
!------end of function ddanrm------
      end function ddanrm
!deck ddaslv
      subroutine ddaslv (neq, delta, wm, iwm)
!***begin prologue  ddaslv
!***subsidiary
!***purpose  linear system solver for ddassl.
!***library   slatec (dassl)
!***type      double precision (sdaslv-s, ddaslv-d)
!***author  petzold, linda r., (llnl)
!***description
!-----------------------------------------------------------------------
!     this routine manages the solution of the linear
!     system arising in the newton iteration.
!     matrices and real temporary storage and
!     real information are stored in the array wm.
!     integer matrix information is stored in
!     the array iwm.
!     for a dense matrix, the linpack routine
!     dgesl is called.
!     for a banded matrix,the linpack routine
!     dgbsl is called.
!-----------------------------------------------------------------------
!***routines called  dgbsl, dgesl
!***revision history  (yymmdd)
!   830315  date written
!   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
!   901019  merged changes made by c. ulrich with slatec 4.0 format.
!   901026  added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (fnf)
!***end prologue  ddaslv
!
      integer  neq, iwm(*)
      double precision  delta(*), wm(*)
!
      external  dgbsl, dgesl
!
      integer  lipvt, lml, lmu, lmtype, meband, mtype, npd
      parameter (npd=1)
      parameter (lml=1)
      parameter (lmu=2)
      parameter (lmtype=4)
      parameter (lipvt=21)
!
!***first executable statement  ddaslv
      mtype=iwm(lmtype)
      go to(100,100,300,400,400),mtype
!
!     dense matrix
  100 call dgesl(wm(npd),neq,neq,iwm(lipvt),delta,0)
      return
!
!     dummy section for mtype=3
  300 continue
      return
!
!     banded matrix
  400 meband=2*iwm(lml)+iwm(lmu)+1
      call dgbsl(wm(npd),meband,neq,iwm(lml),                           &
     &  iwm(lmu),iwm(lipvt),delta,0)
      return
!------end of subroutine ddaslv------
      end subroutine ddaslv
!deck ddastp
      subroutine ddastp (eqn, x, y, yprime, neq, h, wt, jstart,    &
     &   idid, rpar, ipar, phi, delta, e, wm, iwm, alpha, beta, gamma,  &
     &   psi, sigma, cj, cjold, hold, s, hmin, uround, iphase, jcalc, k,&
     &   kold, ns, nonneg, ntemp)
!***begin prologue  ddastp
!***subsidiary
!***purpose  perform one step of the ddassl integration.
!***library   slatec (dassl)
!***type      double precision (sdastp-s, ddastp-d)
!***author  petzold, linda r., (llnl)
!***description
!-----------------------------------------------------------------------
!     ddastp solves a system of differential/
!     algebraic equations of the form
!     g(x,y,yprime) = 0,  for one step (normally
!     from x to x+h).
!
!     the methods used are modified divided
!     difference,fixed leading coefficient
!     forms of backward differentiation
!     formulas. the code adjusts the stepsize
!     and order to control the local error per
!     step.
!
!
!     the parameters represent
!     x  --        independent variable
!     y  --        solution vector at x
!     yprime --    derivative of solution vector
!                  after successful step
!     neq --       number of equations to be integrated
!     res --       external user-supplied subroutine
!                  to evaluate the residual.  the call is
!                  call res(x,y,yprime,delta,ires,rpar,ipar)
!                  x,y,yprime are input.  delta is output.
!                  on input, ires=0.  res should alter ires only
!                  if it encounters an illegal value of y or a
!                  stop condition.  set ires=-1 if an input value
!                  of y is illegal, and ddastp will try to solve
!                  the problem without getting ires = -1.  if
!                  ires=-2, ddastp returns control to the calling
!                  program with idid = -11.
!     jac --       external user-supplied routine to evaluate
!                  the iteration matrix (this is optional)
!                  the call is of the form
!                  call jac(x,y,yprime,pd,cj,rpar,ipar)
!                  pd is the matrix of partial derivatives,
!                  pd=dg/dy+cj*dg/dyprime
!     h --         appropriate step size for next step.
!                  normally determined by the code
!     wt --        vector of weights for error criterion.
!     jstart --    integer variable set 0 for
!                  first step, 1 otherwise.
!     idid --      completion code with the following meanings:
!                  idid= 1 -- the step was completed successfully
!                  idid=-6 -- the error test failed repeatedly
!                  idid=-7 -- the corrector could not converge
!                  idid=-8 -- the iteration matrix is singular
!                  idid=-9 -- the corrector could not converge.
!                             there were repeated error test
!                             failures on this step.
!                  idid=-10-- the corrector could not converge
!                             because ires was equal to minus one
!                  idid=-11-- ires equal to -2 was encountered,
!                             and control is being returned to
!                             the calling program
!     rpar,ipar -- real and integer parameter arrays that
!                  are used for communication between the
!                  calling program and external user routines
!                  they are not altered by ddastp
!     phi --       array of divided differences used by
!                  ddastp. the length is neq*(k+1),where
!                  k is the maximum order
!     delta,e --   work vectors for ddastp of length neq
!     wm,iwm --    real and integer arrays storing
!                  matrix information such as the matrix
!                  of partial derivatives,permutation
!                  vector, and various other information.
!
!     the other parameters are information
!     which is needed internally by ddastp to
!     continue from step to step.
!
!-----------------------------------------------------------------------
!***routines called  ddajac, ddanrm, ddaslv, ddatrp
!***revision history  (yymmdd)
!   830315  date written
!   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
!   901019  merged changes made by c. ulrich with slatec 4.0 format.
!   901026  added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (fnf)
!   951030  reset psi(1), phi(*,2) at 690. (ach)
!   000711  fixed newton convergence test below 360 (ach)
!***end prologue  ddastp
!
      integer  neq, jstart, idid, ipar(*), iwm(*), iphase, jcalc, k,    &
     &   kold, ns, nonneg, ntemp
      double precision                                                  &
     &   x, y(*), yprime(*), h, wt(*), rpar(*), phi(neq,*), delta(*),   &
     &   e(*), wm(*), alpha(*), beta(*), gamma(*), psi(*), sigma(*), cj,&
     &   cjold, hold, s, hmin, uround
      class(problem_ddassl) :: eqn

      integer  i, ier, ires, j, j1, kdiff, km1, knew, kp1, kp2, lctf,   &
     &   letf, lmxord, lnje, lnre, lnst, m, maxit, ncf, nef, nsf, nsp1
      double precision                                                  &
     &   alpha0, alphas, cjlast, ck, delnrm, enorm, erk, erkm1,         &
     &   erkm2, erkp1, err, est, hnew, oldnrm, pnorm, r, rate, temp1,   &
     &   temp2, terk, terkm1, terkm2, terkp1, xold, xrate
      logical  convgd
!
      parameter (lmxord=3)
      parameter (lnst=11)
      parameter (lnre=12)
      parameter (lnje=13)
      parameter (letf=14)
      parameter (lctf=15)
!
      data maxit/4/
      data xrate/0.25d0/
!
!
!
!
!
!-----------------------------------------------------------------------
!     block 1.
!     initialize. on the first call,set
!     the order to 1 and initialize
!     other variables.
!-----------------------------------------------------------------------
!
!     initializations for all calls
!***first executable statement  ddastp
      idid=1
      xold=x
      ncf=0
      nsf=0
      nef=0
      if(jstart .ne. 0) go to 120
!
!     if this is the first step,perform
!     other initializations
      iwm(letf) = 0
      iwm(lctf) = 0
      k=1
      kold=0
      hold=0.0d0
      jstart=1
      psi(1)=h
      cjold = 1.0d0/h
      cj = cjold
      s = 100.d0
      jcalc = -1
      delnrm=1.0d0
      iphase = 0
      ns=0
  120 continue
!
!
!
!
!
!-----------------------------------------------------------------------
!     block 2
!     compute coefficients of formulas for
!     this step.
!-----------------------------------------------------------------------
  200 continue
      kp1=k+1
      kp2=k+2
      km1=k-1
      xold=x
      if(h.ne.hold.or.k .ne. kold) ns = 0
      ns=min(ns+1,kold+2)
      nsp1=ns+1
      if(kp1 .lt. ns)go to 230
!
      beta(1)=1.0d0
      alpha(1)=1.0d0
      temp1=h
      gamma(1)=0.0d0
      sigma(1)=1.0d0
      do 210 i=2,kp1
         temp2=psi(i-1)
         psi(i-1)=temp1
         beta(i)=beta(i-1)*psi(i-1)/temp2
         temp1=temp2+h
         alpha(i)=h/temp1
         sigma(i)=(i-1)*sigma(i-1)*alpha(i)
         gamma(i)=gamma(i-1)+alpha(i-1)/h
  210    continue
      psi(kp1)=temp1
  230 continue
!
!     compute alphas, alpha0
      alphas = 0.0d0
      alpha0 = 0.0d0
      do 240 i = 1,k
        alphas = alphas - 1.0d0/i
        alpha0 = alpha0 - alpha(i)
  240   continue
!
!     compute leading coefficient cj
      cjlast = cj
      cj = -alphas/h
!
!     compute variable stepsize error coefficient ck
      ck = abs(alpha(kp1) + alphas - alpha0)
      ck = max(ck,alpha(kp1))
!
!     decide whether new jacobian is needed
      temp1 = (1.0d0 - xrate)/(1.0d0 + xrate)
      temp2 = 1.0d0/temp1
      if (cj/cjold .lt. temp1 .or. cj/cjold .gt. temp2) jcalc = -1
      if (cj .ne. cjlast) s = 100.d0
!
!     change phi to phi star
      if(kp1 .lt. nsp1) go to 280
      do 270 j=nsp1,kp1
         do 260 i=1,neq
  260       phi(i,j)=beta(j)*phi(i,j)
  270    continue
  280 continue
!
!     update time
      x=x+h
!
!
!
!
!
!-----------------------------------------------------------------------
!     block 3
!     predict the solution and derivative,
!     and solve the corrector equation
!-----------------------------------------------------------------------
!
!     first,predict the solution and derivative
  300 continue
      do 310 i=1,neq
         y(i)=phi(i,1)
  310    yprime(i)=0.0d0
      do 330 j=2,kp1
         do 320 i=1,neq
            y(i)=y(i)+phi(i,j)
  320       yprime(i)=yprime(i)+gamma(j)*phi(i,j)
  330 end do
      pnorm = ddanrm (neq,y,wt,rpar,ipar)
!
!
!
!     solve the corrector equation using a
!     modified newton scheme.
      convgd= .true.
      m=0
      iwm(lnre)=iwm(lnre)+1
      ires = 0
      call eqn%res(x,y,yprime,delta,ires,rpar,ipar)
      if (ires .lt. 0) go to 380
!
!
!     if indicated,reevaluate the
!     iteration matrix pd = dg/dy + cj*dg/dyprime
!     (where g(x,y,yprime)=0). set
!     jcalc to 0 as an indicator that
!     this has been done.
      if(jcalc .ne. -1)go to 340
      iwm(lnje)=iwm(lnje)+1
      jcalc=0
      call eqn%ddajac(neq,x,y,yprime,delta,cj,h,                            &
     & ier,wt,e,wm,iwm,ires,uround,rpar,                        &
     & ipar,ntemp)
      cjold=cj
      s = 100.d0
      if (ires .lt. 0) go to 380
      if(ier .ne. 0)go to 380
      nsf=0
!
!
!     initialize the error accumulation vector e.
  340 continue
      do 345 i=1,neq
  345    e(i)=0.0d0
!
!
!     corrector loop.
  350 continue
!
!     multiply residual by temp1 to accelerate convergence
      temp1 = 2.0d0/(1.0d0 + cj/cjold)
      do 355 i = 1,neq
  355   delta(i) = delta(i) * temp1
!
!     compute a new iterate (back-substitution).
!     store the correction in delta.
      call ddaslv(neq,delta,wm,iwm)
!
!     update y, e, and yprime
      do 360 i=1,neq
         y(i)=y(i)-delta(i)
         e(i)=e(i)-delta(i)
  360    yprime(i)=yprime(i)-cj*delta(i)
!
!     test for convergence of the iteration
      delnrm=ddanrm(neq,delta,wt,rpar,ipar)
      if (m .gt. 0) go to 365
         if (delnrm .le. 100.d0*uround*pnorm) go to 375
         oldnrm = delnrm
         go to 367
  365 rate = (delnrm/oldnrm)**(1.0d0/m)
      if (rate .gt. 0.90d0) go to 370
      s = rate/(1.0d0 - rate)
  367 if (s*delnrm .le. 0.33d0) go to 375
!
!     the corrector has not yet converged.
!     update m and test whether the
!     maximum number of iterations have
!     been tried.
      m=m+1
      if(m.ge.maxit)go to 370
!
!     evaluate the residual
!     and go back to do another iteration
      iwm(lnre)=iwm(lnre)+1
      ires = 0
      call eqn%res(x,y,yprime,delta,ires,                                   &
     &  rpar,ipar)
      if (ires .lt. 0) go to 380
      go to 350
!
!
!     the corrector failed to converge in maxit
!     iterations. if the iteration matrix
!     is not current,re-do the step with
!     a new iteration matrix.
  370 continue
      if(jcalc.eq.0)go to 380
      jcalc=-1
      go to 300
!
!
!     the iteration has converged.  if nonnegativity of solution is
!     required, set the solution nonnegative, if the perturbation
!     to do it is small enough.  if the change is too large, then
!     consider the corrector iteration to have failed.
  375 if(nonneg .eq. 0) go to 390
      do 377 i = 1,neq
  377    delta(i) = min(y(i),0.0d0)
      delnrm = ddanrm(neq,delta,wt,rpar,ipar)
      if(delnrm .gt. 0.33d0) go to 380
      do 378 i = 1,neq
  378    e(i) = e(i) - delta(i)
      go to 390
!
!
!     exits from block 3
!     no convergence with current iteration
!     matrix,or singular iteration matrix
  380 convgd= .false.
  390 jcalc = 1
      if(.not.convgd)go to 600
!
!
!
!
!
!-----------------------------------------------------------------------
!     block 4
!     estimate the errors at orders k,k-1,k-2
!     as if constant stepsize was used. estimate
!     the local error at order k and test
!     whether the current step is successful.
!-----------------------------------------------------------------------
!
!     estimate errors at orders k,k-1,k-2
      enorm = ddanrm(neq,e,wt,rpar,ipar)
      erk = sigma(k+1)*enorm
      terk = (k+1)*erk
      est = erk
      knew=k
      if(k .eq. 1)go to 430
      do 405 i = 1,neq
  405   delta(i) = phi(i,kp1) + e(i)
      erkm1=sigma(k)*ddanrm(neq,delta,wt,rpar,ipar)
      terkm1 = k*erkm1
      if(k .gt. 2)go to 410
      if(terkm1 .le. 0.5d0*terk)go to 420
      go to 430
  410 continue
      do 415 i = 1,neq
  415   delta(i) = phi(i,k) + delta(i)
      erkm2=sigma(k-1)*ddanrm(neq,delta,wt,rpar,ipar)
      terkm2 = (k-1)*erkm2
      if(max(terkm1,terkm2).gt.terk)go to 430
!     lower the order
  420 continue
      knew=k-1
      est = erkm1
!
!
!     calculate the local error for the current step
!     to see if the step was successful
  430 continue
      err = ck * enorm
      if(err .gt. 1.0d0)go to 600
!
!
!
!
!
!-----------------------------------------------------------------------
!     block 5
!     the step is successful. determine
!     the best order and stepsize for
!     the next step. update the differences
!     for the next step.
!-----------------------------------------------------------------------
      idid=1
      iwm(lnst)=iwm(lnst)+1
      kdiff=k-kold
      kold=k
      hold=h
!
!
!     estimate the error at order k+1 unless:
!        already decided to lower order, or
!        already using maximum order, or
!        stepsize not constant, or
!        order raised in previous step
      if(knew.eq.km1.or.k.eq.iwm(lmxord))iphase=1
      if(iphase .eq. 0)go to 545
      if(knew.eq.km1)go to 540
      if(k.eq.iwm(lmxord)) go to 550
      if(kp1.ge.ns.or.kdiff.eq.1)go to 550
      do 510 i=1,neq
  510    delta(i)=e(i)-phi(i,kp2)
      erkp1 = (1.0d0/(k+2))*ddanrm(neq,delta,wt,rpar,ipar)
      terkp1 = (k+2)*erkp1
      if(k.gt.1)go to 520
      if(terkp1.ge.0.5d0*terk)go to 550
      go to 530
  520 if(terkm1.le.min(terk,terkp1))go to 540
      if(terkp1.ge.terk.or.k.eq.iwm(lmxord))go to 550
!
!     raise order
  530 k=kp1
      est = erkp1
      go to 550
!
!     lower order
  540 k=km1
      est = erkm1
      go to 550
!
!     if iphase = 0, increase order by one and multiply stepsize by
!     factor two
  545 k = kp1
      hnew = h*2.0d0
      h = hnew
      go to 575
!
!
!     determine the appropriate stepsize for
!     the next step.
  550 hnew=h
      temp2=k+1
      r=(2.0d0*est+0.0001d0)**(-1.0d0/temp2)
      if(r .lt. 2.0d0) go to 555
      hnew = 2.0d0*h
      go to 560
  555 if(r .gt. 1.0d0) go to 560
      r = max(0.5d0,min(0.9d0,r))
      hnew = h*r
  560 h=hnew
!
!
!     update differences for next step
  575 continue
      if(kold.eq.iwm(lmxord))go to 585
      do 580 i=1,neq
  580    phi(i,kp2)=e(i)
  585 continue
      do 590 i=1,neq
  590    phi(i,kp1)=phi(i,kp1)+e(i)
      do 595 j1=2,kp1
         j=kp1-j1+1
         do 595 i=1,neq
  595    phi(i,j)=phi(i,j)+phi(i,j+1)
      return
!
!
!
!
!
!-----------------------------------------------------------------------
!     block 6
!     the step is unsuccessful. restore x,psi,phi
!     determine appropriate stepsize for
!     continuing the integration, or exit with
!     an error flag if there have been many
!     failures.
!-----------------------------------------------------------------------
  600 iphase = 1
!
!     restore x,phi,psi
      x=xold
      if(kp1.lt.nsp1)go to 630
      do 620 j=nsp1,kp1
         temp1=1.0d0/beta(j)
         do 610 i=1,neq
  610       phi(i,j)=temp1*phi(i,j)
  620    continue
  630 continue
      do 640 i=2,kp1
  640    psi(i-1)=psi(i)-h
!
!
!     test whether failure is due to corrector iteration
!     or error test
      if(convgd)go to 660
      iwm(lctf)=iwm(lctf)+1
!
!
!     the newton iteration failed to converge with
!     a current iteration matrix.  determine the cause
!     of the failure and take appropriate action.
      if(ier.eq.0)go to 650
!
!     the iteration matrix is singular. reduce
!     the stepsize by a factor of 4. if
!     this happens three times in a row on
!     the same step, return with an error flag
      nsf=nsf+1
      r = 0.25d0
      h=h*r
      if (nsf .lt. 3 .and. abs(h) .ge. hmin) go to 690
      idid=-8
      go to 675
!
!
!     the newton iteration failed to converge for a reason
!     other than a singular iteration matrix.  if ires = -2, then
!     return.  otherwise, reduce the stepsize and try again, unless
!     too many failures have occurred.
  650 continue
      if (ires .gt. -2) go to 655
      idid = -11
      go to 675
  655 ncf = ncf + 1
      r = 0.25d0
      h = h*r
      if (ncf .lt. 10 .and. abs(h) .ge. hmin) go to 690
      idid = -7
      if (ires .lt. 0) idid = -10
      if (nef .ge. 3) idid = -9
      go to 675
!
!
!     the newton scheme converged, and the cause
!     of the failure was the error estimate
!     exceeding the tolerance.
  660 nef=nef+1
      iwm(letf)=iwm(letf)+1
      if (nef .gt. 1) go to 665
!
!     on first error test failure, keep current order or lower
!     order by one.  compute new stepsize based on differences
!     of the solution.
      k = knew
      temp2 = k + 1
      r = 0.90d0*(2.0d0*est+0.0001d0)**(-1.0d0/temp2)
      r = max(0.25d0,min(0.9d0,r))
      h = h*r
      if (abs(h) .ge. hmin) go to 690
      idid = -6
      go to 675
!
!     on second error test failure, use the current order or
!     decrease order by one.  reduce the stepsize by a factor of
!     four.
  665 if (nef .gt. 2) go to 670
      k = knew
      r = 0.25d0
      h = r*h
      if (abs(h) .ge. hmin) go to 690
      idid = -6
      go to 675
!
!     on third and subsequent error test failures, set the order to
!     one and reduce the stepsize by a factor of four.
  670 k = 1
      r = 0.25d0
      h = r*h
      if (abs(h) .ge. hmin) go to 690
      idid = -6
      go to 675
!
!
!     for all crashes, restore y to its last value,
!     interpolate to find yprime at last x, and return
  675 continue
      call ddatrp(x,x,y,yprime,neq,k,phi,psi)
      return
!
!
!     go back and try this step again.
!     if this is the first step, reset psi(1) and rescale phi(*,2).
  690 if (kold .eq. 0) then
        psi(1) = h
        do 695 i = 1,neq
  695     phi(i,2) = r*phi(i,2)
        endif
      go to 200
!
!------end of subroutine ddastp------
      end subroutine ddastp
!deck ddatrp
      subroutine ddatrp (x, xout, yout, ypout, neq, kold, phi, psi)
!***begin prologue  ddatrp
!***subsidiary
!***purpose  interpolation routine for ddassl.
!***library   slatec (dassl)
!***type      double precision (sdatrp-s, ddatrp-d)
!***author  petzold, linda r., (llnl)
!***description
!-----------------------------------------------------------------------
!     the methods in subroutine ddastp use polynomials
!     to approximate the solution. ddatrp approximates the
!     solution and its derivative at time xout by evaluating
!     one of these polynomials, and its derivative,there.
!     information defining this polynomial is passed from
!     ddastp, so ddatrp cannot be used alone.
!
!     the parameters are:
!     x     the current time in the integration.
!     xout  the time at which the solution is desired
!     yout  the interpolated approximation to y at xout
!           (this is output)
!     ypout the interpolated approximation to yprime at xout
!           (this is output)
!     neq   number of equations
!     kold  order used on last successful step
!     phi   array of scaled divided differences of y
!     psi   array of past stepsize history
!-----------------------------------------------------------------------
!***routines called  (none)
!***revision history  (yymmdd)
!   830315  date written
!   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
!   901019  merged changes made by c. ulrich with slatec 4.0 format.
!   901026  added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (fnf)
!***end prologue  ddatrp
!
      integer  neq, kold
      double precision  x, xout, yout(*), ypout(*), phi(neq,*), psi(*)
!
      integer  i, j, koldp1
      double precision  c, d, gamma, temp1
!
!***first executable statement  ddatrp
      koldp1=kold+1
      temp1=xout-x
      do 10 i=1,neq
         yout(i)=phi(i,1)
   10    ypout(i)=0.0d0
      c=1.0d0
      d=0.0d0
      gamma=temp1/psi(1)
      do 30 j=2,koldp1
         d=d*gamma+c/psi(j-1)
         c=c*gamma
         gamma=(temp1+psi(j-1))/psi(j)
         do 20 i=1,neq
            yout(i)=yout(i)+c*phi(i,j)
   20       ypout(i)=ypout(i)+d*phi(i,j)
   30    continue
      return
!
!------end of subroutine ddatrp------
      end subroutine ddatrp
!deck ddawts
      subroutine ddawts (neq, iwt, rtol, atol, y, wt, rpar, ipar)
!***begin prologue  ddawts
!***subsidiary
!***purpose  set error weight vector for ddassl.
!***library   slatec (dassl)
!***type      double precision (sdawts-s, ddawts-d)
!***author  petzold, linda r., (llnl)
!***description
!-----------------------------------------------------------------------
!     this subroutine sets the error weight vector
!     wt according to wt(i)=rtol(i)*abs(y(i))+atol(i),
!     i=1,-,n.
!     rtol and atol are scalars if iwt = 0,
!     and vectors if iwt = 1.
!-----------------------------------------------------------------------
!***routines called  (none)
!***revision history  (yymmdd)
!   830315  date written
!   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
!   901019  merged changes made by c. ulrich with slatec 4.0 format.
!   901026  added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (fnf)
!***end prologue  ddawts
!
      integer  neq, iwt, ipar(*)
      double precision  rtol(*), atol(*), y(*), wt(*), rpar(*)
!
      integer  i
      double precision  atoli, rtoli
!
!***first executable statement  ddawts
      rtoli=rtol(1)
      atoli=atol(1)
      do 20 i=1,neq
         if (iwt .eq.0) go to 10
           rtoli=rtol(i)
           atoli=atol(i)
   10      wt(i)=rtoli*abs(y(i))+atoli
   20      continue
      return
!-----------end of subroutine ddawts------------------------------------
      end subroutine ddawts
!deck xermsg
      subroutine xermsg (librar, subrou, messg, nerr, level)
!-----------------------------------------------------------------------
! subroutines xermsg, xsetf, xsetun, and the function routine ixsav, as
! given here, constitute a simplified version of the slatec error
! handling package.  written by a. c. hindmarsh, 18 november 1992.
!
! all arguments are input arguments.
! librar = library name (character array).  prefixed to message.
! subrou = routine name (character array).  prefixed to message.
! messg  = the message (character array).
! nerr   = integer error number.  prefixed to message.
! level  = the error level..
!          0 or 1 means recoverable (control returns to caller).
!          2 means fatal (run is aborted--see note below).
!
! note..  this routine has been simplified in the following ways..
! 1. a single prefix line is printed with nerr, subrou, and librar.
! 2. the message in messg is printed, unaltered, on lines of up to 72
!    characters each using a format of (a).
! 3. if level = 2, control passes to the statement   stop
!    to abort the run.  this statement may be machine-dependent.
!
! for a different default logical unit number, change the data
! statement in function routine ixsav.
! for a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
! subroutines called by xermsg.. none
! function routines called by xermsg.. ixsav
! intrinsic function used by xermsg.. len
!-----------------------------------------------------------------------
      character*(*) librar, subrou, messg
      integer nerr, level
      integer i1, i2, il, lenmsg, llen, lunit, mesflg, nlines
      parameter (llen = 72)
!
! get message print flag and logical unit number. ----------------------
      lunit = ixsav (1, 0, .false.)
      mesflg = ixsav (2, 0, .false.)
      if (mesflg .eq. 0) go to 100
! write nerr, subrou, and librar. --------------------------------------
      i1 = len(subrou)
      i2 = len(librar)
      write (lunit, 10) nerr, subrou(1:i1), librar(1:i2)
   10 format(/,'***error number ',i6,' from ',a,' in library ',a,'***')
! write the message. ---------------------------------------------------
      lenmsg = len(messg)
      nlines = ( (lenmsg - 1)/llen ) + 1
      do 20 il = 1,nlines
        i1 = 1 + (il - 1)*llen
        i2 = min(il*llen,lenmsg)
        write (lunit,'(a)') messg(i1:i2)
   20   continue
! abort the run if level = 2. ------------------------------------------
  100 if (level .ne. 2) return
      stop
!----------------------- end of subroutine xermsg ----------------------
      end subroutine xermsg
!deck xsetun
      subroutine xsetun (lun)
!-----------------------------------------------------------------------
! this routine resets the logical unit number for messages.
!
! subroutines called by xsetun.. none
! function routine called by xsetun.. ixsav
!-----------------------------------------------------------------------
      integer lun, junk
!
      if (lun .gt. 0) junk = ixsav (1,lun,.true.)
      return
!----------------------- end of subroutine xsetun ----------------------
      end subroutine xsetun
!deck xsetf
      subroutine xsetf (mflag)
!-----------------------------------------------------------------------
! this routine resets the print control flag mflag.
!
! subroutines called by xsetf.. none
! function routine called by xsetf.. ixsav
!-----------------------------------------------------------------------
      integer mflag, junk
!
      if (mflag .eq. 0 .or. mflag .eq. 1) junk = ixsav (2,mflag,.true.)
      return
!----------------------- end of subroutine xsetf -----------------------
      end subroutine xsetf
!deck ixsav
      integer function ixsav (ipar, ivalue, iset)
      logical iset
      integer ipar, ivalue
!-----------------------------------------------------------------------
! ixsav saves and recalls one of two error message parameters:
!   lunit, the logical unit number to which messages are printed, and
!   mesflg, the message print flag.
! this is a modification of the slatec library routine j4save.
!
! saved local variables..
!  lunit  = logical unit number for messages.
!           the default is 6 (machine-dependent).
!  mesflg = print control flag..
!           1 means print all messages (the default).
!           0 means no printing.
!
! on input..
!   ipar   = parameter indicator (1 for lunit, 2 for mesflg).
!   ivalue = the value to be set for the parameter, if iset = .true.
!   iset   = logical flag to indicate whether to read or write.
!            if iset = .true., the parameter will be given
!            the value ivalue.  if iset = .false., the parameter
!            will be unchanged, and ivalue is a dummy argument.
!
! on return..
!   ixsav = the (old) value of the parameter.
!
! subroutines/functions called by ixsav.. none
!-----------------------------------------------------------------------
      integer lunit, mesflg
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this routine.
!-----------------------------------------------------------------------
      save lunit, mesflg
      data lunit/6/, mesflg/1/
!
      if (ipar .eq. 1) then
        ixsav = lunit
        if (iset) lunit = ivalue
        endif
!
      if (ipar .eq. 2) then
        ixsav = mesflg
        if (iset) mesflg = ivalue
        endif
!
      return
!----------------------- end of function ixsav -------------------------
      end function ixsav


!deck d1mach
      double precision function d1mach (i)
      implicit none
      integer :: i
      double precision :: b, x
!***begin prologue  d1mach
!***purpose  return floating point machine dependent constants.
!***library   slatec
!***category  r1
!***type      single precision (d1mach-s, d1mach-d)
!***keywords  machine constants
!***author  fox, p. a., (bell labs)
!           hall, a. d., (bell labs)
!           schryer, n. l., (bell labs)
!***description
!
!   d1mach can be used to obtain machine-dependent parameters for the
!   local machine environment.  it is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        a = d1mach(i)
!
!   where i=1,...,5.  the (output) value of a above is determined by
!   the (input) value of i.  the results for various values of i are
!   discussed below.
!
!   d1mach(1) = b**(emin-1), the smallest positive magnitude.
!   d1mach(2) = b**emax*(1 - b**(-t)), the largest magnitude.
!   d1mach(3) = b**(-t), the smallest relative spacing.
!   d1mach(4) = b**(1-t), the largest relative spacing.
!   d1mach(5) = log10(b)
!
!   assume single precision numbers are represented in the t-digit,
!   base-b form
!
!              sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
!
!   where 0 .le. x(i) .lt. b for i=1,...,t, 0 .lt. x(1), and
!   emin .le. e .le. emax.
!
!   the values of b, t, emin and emax are provided in i1mach as
!   follows:
!   i1mach(10) = b, the base.
!   i1mach(11) = t, the number of base-b digits.
!   i1mach(12) = emin, the smallest exponent e.
!   i1mach(13) = emax, the largest exponent e.
!
!
!***references  p. a. fox, a. d. hall and n. l. schryer, framework for
!                 a portable library, acm transactions on mathematical
!                 software 4, 2 (june 1978), pp. 177-188.
!***routines called  xermsg
!***revision history  (yymmdd)
!   790101  date written
!   960329  modified for fortran 90 (be after suggestions by ehg)
!***end prologue  d1mach
!
      x = 1.0d0
      b = radix(x)
      select case (i)
        case (1)
          d1mach = b**(minexponent(x)-1) ! the smallest positive magnitude.
        case (2)
          d1mach = huge(x)               ! the largest magnitude.
        case (3)
          d1mach = b**(-digits(x))       ! the smallest relative spacing.
        case (4)
          d1mach = b**(1-digits(x))      ! the largest relative spacing.
        case (5)
          d1mach = log10(b)
        case default
          write (*, fmt = 9000)
 9000     format ('1error    1 in d1mach - i out of bounds')
          stop
      end select
      return
      end function d1mach

end module ddassl_mod
