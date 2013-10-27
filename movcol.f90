module movcol_mod

  use ddassl_mod

  ! numerical parameters
  real(8), parameter ::&
       &zero  = 0.0_8  ,&
       &quart = 0.25_8 ,&
       &half  = 0.5_8  ,&
       &one   = 1.0_8  ,&
       &two   = 2.0_8  ,&
       &three = 3.0_8  ,&
       &four  = 4.0_8  ,&
       &six   = 6.0_8

  ! algorithm parameters?
  real(8), parameter ::&
       &tau1 = 1.e-4_8,&
       &s1 = 0.211324865347452d0,&
       &s2 = 0.788675134652548d0

  ! real(16), parameter ::&
  !      &tau2 = 0.0001_16

  private

  type, abstract, public, extends(problem_ddassl) :: problem_movcol
     ! work arrays
     real(8), allocatable, private :: rwork(:)
     integer, allocatable, private :: iwork(:)

     ! pointers to the physical quantities
     real(8), pointer :: x(:), u(:,:), ux(:,:), y(:,:), ydot(:,:)

     ! print times
     real(8), allocatable :: touta(:)

     ! parameters
     real(8) :: par(20)            !I don't know why par has size 20
     real(8), pointer :: job       !par(1)
     real(8), pointer :: output    !par(2)
     real(8), pointer :: mmpde     !par(3)
     real(8), pointer :: tau       !par(4)
     real(8), pointer :: gamma     !par(5)
     real(8), pointer :: ip        !par(6)
     real(8), pointer :: stpmax    !par(7)
     real(8), pointer :: left_end  !par(8)
     real(8), pointer :: right_end !par(9)

     ! system size and number of mesh points
     integer :: npts, npde

     ! output file id
     integer :: nprnt = 111
   contains
     private

     ! main procedure
     procedure, public :: solve => movcol_solve
     ! initialization function
     procedure, public :: init => movcol_init

     ! other private  procedures
     procedure :: movcol
     procedure :: slnout
     procedure :: respde
     procedure :: respd1
     procedure :: evlmnt
     procedure :: resmeq
     ! procedures defined for parent problem_ddassl
     procedure :: res => resode
     procedure :: jac => jacode
     ! these procedures must be defined in a child class
     procedure(defivs_i), deferred :: defivs
     procedure(defout_i), deferred :: defout
     procedure(defmsh_i), deferred :: defmsh
     procedure(defpde_i), deferred :: defpde
     procedure(defbcp_i), deferred :: defbcp
     procedure(defbcm_i), deferred :: defbcm
     procedure(defmnt_i), deferred :: defmnt
  end type problem_movcol

  abstract interface

     subroutine defivs_i (eqn, x, u, ux)
       import problem_movcol
       class(problem_movcol)          :: eqn
       real(8) :: x, u(eqn%npde), ux(eqn%npde)
     end subroutine defivs_i

     subroutine defout_i (eqn, npde, npts, t, xmesh, xmesht, u, ux, ut, tstep,&
          touta, ntouta, istop, index, nts)
       import problem_movcol
       class(problem_movcol)          :: eqn
       real(8), dimension(npts, npde) :: u,  ux, ut
       real(8), dimension(npts)       :: xmesh, xmesht
       real(8), dimension(ntouta)     :: touta
       real(8)                        :: t, tstep
     end subroutine defout_i

     subroutine defmsh_i (eqn, xmesh)
       import problem_movcol
       class(problem_movcol)  :: eqn
       real(8) :: xmesh(eqn%npts)
     end subroutine defmsh_i

     subroutine defpde_i (eqn, index, t, x, u, ux, ut, uxt, fg)
       import problem_movcol
       class(problem_movcol) :: eqn
       integer :: index
       real(8) :: t, x
       real(8), dimension(eqn%npde) ::  u, ux, ut, uxt, fg
     end subroutine defpde_i

     subroutine defbcp_i (eqn, index, t, x, xt, u, ux, uxx, ut, uxt, res)
       import problem_movcol
       class(problem_movcol)  :: eqn
       integer :: index
       real(8) :: t, x, xt
       real(8), dimension(eqn%npde) ::  u, ux, uxx, ut, uxt, res
     end subroutine defbcp_i

     subroutine defbcm_i (eqn, index, t, x, xt, u, ux, uxx, ut, uxt, res)
       import problem_movcol
       class(problem_movcol)  :: eqn
       integer :: index
       real(8) :: t, x, xt, res
       real(8), dimension(eqn%npde) ::  u, ux, uxx, ut, uxt
     end subroutine defbcm_i

     subroutine defmnt_i (eqn, t, x, u, ux, uxx, ut, uxt, fmntr)
       import problem_movcol
       class(problem_movcol)  :: eqn
       real(8) :: t, x, fmntr
       real(8), dimension(eqn%npde) :: u, ux, uxx, ut, uxt
     end subroutine defmnt_i

  end interface

  contains

    ! allocate the work tables
    subroutine movcol_init(eqn, npde, npts)
      class(problem_movcol), target :: eqn
      integer :: npde, npts

      ! length of work arrays
      integer :: lrw, liw

      ! ipmax?
      integer :: ipmax = 4

      ! number of dependent values per mesh point
      integer :: m

      eqn%npde  = npde
      eqn%npts  = npts

      lrw = 60 + 6 * npde+ (22 + 3 * ipmax) * npts&
           & + (55 +  12 * ipmax) * npts * npde&
           & + (24 + 12 * ipmax) * npts * npde**2
      liw = 20 + npts + 2 * npts * npde

      allocate(eqn%rwork(lrw))
      allocate(eqn%iwork(liw))

      ! initialize work arrays with zeroes
      eqn%rwork = 0.0_8
      eqn%iwork = 0


      ! if the touta table was not allocate, allocate it as an empty table
      if(.not. allocated(eqn%touta)) then
         allocate(eqn%touta(0))
      end if

      ! arrange the parameter pointers
      eqn%job       => eqn%par(1)
      eqn%output    => eqn%par(2)
      eqn%mmpde     => eqn%par(3)
      eqn%tau       => eqn%par(4)
      eqn%gamma     => eqn%par(5)
      eqn%ip        => eqn%par(6)
      eqn%stpmax    => eqn%par(7)
      eqn%left_end  => eqn%par(8)
      eqn%right_end => eqn%par(9)

      ! set the default values for parameters
      ! set the whole table to -1
      eqn%par       = -1
      ! then set the particular parameters
      eqn%output    = 6         !par(2)
      eqn%left_end  = 0         !par(8)
      eqn%right_end = 1         !par(9)

      m = 2*npde+1
      eqn%y   (1:m,1:npts) => eqn%rwork(1       :  m*npts)
      eqn%ydot(1:m,1:npts) => eqn%rwork(m*npts+1:2*npts*m)

      eqn%x  => eqn%y(m,            :)
      eqn%u  => eqn%y(1:npde,       :)
      eqn%ux => eqn%y(npde+1:2*npde,:)

    end subroutine movcol_init

    subroutine movcol_solve(eqn, atol, rtol, touta, iflag, filename)
      class(problem_movcol) :: eqn
      real(8) :: atol, rtol, touta(:)
      integer :: iflag
      character(len=*) :: filename

      real :: tcpu, timaray(2)

      open(newunit = eqn%nprnt, file = filename)

      tcpu = etime (timaray)

      call eqn%movcol(&
           &atol   = atol,&
           &rtol   = rtol,&
           &touta  = touta,&
           &iflag  = iflag)

      tcpu = etime (timaray) - tcpu

      !-----
      ! output cpu times
      !-----

      write (6, 99995) tcpu
      close (eqn%nprnt)
      !
99995 format(/' **** cpu time for mesh movement: ', f12.2/)

    end subroutine movcol_solve


!
! version: nov. 22, 1995
!
! 1/17/96  change integer phypde to logical phypde
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
! purpose :                                                           c
!                                                                     c
!          movcol is primarily intended to solve systems of           c
!          second-order parabolic pdes in one space dimension. it     c
!          is also capable of solving hyperbolic pdes with suitably   c
!          smooth solutions [hr1].                                    c
!                                                                     c
!          movcol uses a method of lines approach based upon a        c
!          moving collocation method [hr1].  the physical pdes are    c
!          discretized in space with a cubic hermite colloction-      c
!          type method, and the mmpdes (moving mesh pdes) for         c
!          computing the moving mesh points [hrr1] are discretized    c
!          with a 3-point finite difference method.  the resulting    c
!          ode system is integrated in time with the dae solver       c
!          dassl developed by l. petzold [petzold].                   c
!                                                                     c
!          the physical pdes are considered in the general divergence c
!          form                                                       c
!                                                                     c
!          f(t, x, u, ux, ut, uxt) = [ g(t, x, u, ux, ut, uxt) ]_x    c
!                                                                     c
!                                            x^l(t) < x < x^r(t)      c
!                                            t_a < t < t_b            c
!                                                                     c
!          subject to the separated boundary conditions               c
!                                                                     c
!          b^l(t, x, xt, u, ux, uxx, ut, uxt) = 0   at x = x^l(t)     c
!          b^r(t, x, xt, u, ux, uxx, ut, uxt) = 0   at x = x^r(t)     c
!                                                                     c
!          and the initial conditions                                 c
!                                                                     c
!          u(x,t_a) = u(x)                                            c
!                                                                     c
!          where f, g, b^l, b^r, u and u are vector-valued functions  c
!          of length npde. it is assumed that the boundaries x^l(t)   c
!          and x^r(t) are determined by the two conditions            c
!                                                                     c
!          b^l_m(t, x, xt, u, ux, uxx, ut, uxt) = 0   at x = x^l(t)   c
!          b^r_m(t, x, xt, u, ux, uxx, ut, uxt) = 0   at x = x^r(t)   c
!                                                                     c
!          where b^l_m and b^r_m are scalar functions.                c
!                                                                     c
!          it is the user's responsibility to ensure that the pdes    c
!          are well posed.                                            c
!                                                                     c
!          notation :                                                 c
!                                                                     c
!          t      -- time                                             c
!          x      -- physical spatial coordinate                      c
!          x^l(t) -- left boundary                                    c
!          x^r(t) -- right boundary                                   c
!          u      -- physical solution variable                       c
!          xt     -- mesh speed x_t                                   c
!          ux     -- 1st derivative of u with respect to x            c
!          uxx    -- 2nd derivative of u with respect to x            c
!          ut     -- 1st derivative of u with respect to t            c
!          uxt    -- 2nd mixed derivative                             c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
! authors:                                                            c
!                                                                     c
!          weizhang huang                                             c
!          department of mathematics                                  c
!          the university of kansas                                   c
!          lawrence, ks 66044                                         c
!          usa                                                        c
!          whuang@math.ku.edu                                      c
!                                                                     c
!          robert d. russell                                          c
!          department of mathematics and statistics                   c
!          simon fraser university                                    c
!          burnaby, b.c. v5a 1s6, canada                              c
!          rdr@cs.sfu.ca                                              c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
! references :                                                        c
!                                                                     c
! 1. [hr1] w. huang and r. d. russell,                                c
!    a moving collocation method for solving time dependent partial   c
!    differential equations,                                          c
!    appl. numer. math. 20 (1996), 101-116.                           c
! 2. [hr2] w. huang and r. d. russell,                                c
!    analysis of moving mesh partial differential equations with      c
!    spatial smoothing,                                               c
!    siam j. numer. anal. 34 (1997), 1106-1126.                       c
! 3. [hrr1] w. huang, y. ren and r. d. russell,                       c
!    moving mesh partial differential equations (mmpdes) based on the c
!    equidistribution principle,                                      c
!    siam j. numer. anal. 31 (1994), 709--730.                        c
! 4. [hrr2] w. huang, y. ren and r. d. russell,                       c
!    moving mesh methods based on moving mesh partial differential    c
!    equations,                                                       c
!    j. comput. phys. 113 (1994), 279--290.                           c
! 5. [petzold] l. r. petzold,                                         c
!    a description of dassl: a differential/algebraic system solver,  c
!    sand 82-8637, sandia labs, livermore, usa                        c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
! usage :                                                             c
!                                                                     c
!      integer          npde, npts, ntouta, iflag, lrw, liw,          c
!     +                 iwork(liw)                                    c
!      double precision atol, rtol, touta(ntouta), par(20),           c
!     +                 rwork(lrw)                                    c
!                                                                     c
!      call movcol(npde, npts, atol, rtol, touta, ntouta, par,        c
!     +            iflag, rwork, lrw, iwork, liw)                     c
!                                                                     c
!      note :                                                         c
!                                                                     c
!      lrw  >= 60 + 6*npde + (22+3*ip)*npts                           c
!            + (55+12*ip)*npts*npde + (24+12*ip)*npts*npde**2         c
!      liw  >= 20 + npts + 2*npts*npde                                c
!                                                                     c
!      where ip := par(6)                                             c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
! description of variables ( i : input, o : output)                   c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! 1. description of parameters
!
! npde:i    (integer) number of physical pdes
!
! npts:i    (integer) number of mesh nodes
!
! atol:i    (real) absolute tolerance for time integration
!
! rtol:i    (real) relative tolerance for time integration
!
! touta:i   (real array) prescribed times for outputting information
!           about performance of time integration with dassl.
!           the ode system is integrated from touta(1) := t_a to
!           touta(ntouta) := t_b
!
! ntouta:i  (integer) dimension of touta
!
! par:i     (real array) array of length 20 containing optional
!           parameters
!
!   note:   if a parameter does not have a value within its permitted
!           range, then its default value is used. therefore, if a
!           default value is desired, the parameter can be set to a
!           value outside of its permitted range, for example, to a
!           negative number (such as -1)
!
!   par  name     meaning
!
!   1   job   - parameter specifying whether or not an equidistributed
!                mesh is to be generated with respect to the initial
!                solution
!
!                job = 1: solve the physical pdes starting with a mesh
!                         provided by the user in the subroutine defmsh
!                job = 2: solve the physical pdes after first computing
!                         an equidistributed mesh
!                job = 3: generate an equidistributed mesh and stop
!
!                note: in the cases job = 2 and 3, the physical region
!                must be defined by setting par(8) = x^l(t_a) and
!                par(9) = x^r(t_a). the mesh generation uses the
!                initial value u(x) of the physical solution provided
!                by the user in the subroutine defivs. the mesh is
!                generated internally by solving the user-prescribed
!                mmpde and the pde u_t = u(x), 0 <= t <= 1, with the
!                user-provided monitor function
!
!                permitted values: 1 <= job <= 3
!
!                default value: job = 2
!
!   2   inform - parameter for selection of output information about
!                performance of the time integration code dassl
!
!                permitted values: 1 <= inform
!
!                default value: inform = 6
!
!   3   mmpde  - parameter used for specifying the basic moving mesh
!                pde
!
!                mmpde = 1: fixed mesh
!                mmpde = 2: mmpde4 in [hrr1, hrr2]
!                mmpde = 3: mmpde6 in [hrr1, hrr2]
!                mmpde = 4: mmpde(24) with spatial smoothing in [hr2],
!                           a smoothed version of mmpde4 in terms of
!                           the mesh concentration function
!
!                permitted values: 1 <= mmpde <= 4
!
!                default value: mmpde = 4
!
!   4   tau    - time smoothing parameter in mmpdes
!
!                permitted values: 0 <= tau
!
!                recommended values: 1e-6 <= tau <= 1e-3
!
!                default value: tau = 1e-4
!
!   5   gamma  - spatial smoothing parameter used when mmpde = 4
!
!                permitted values: 0 < gamma
!
!                recommended values: 1 <= gamma <= 2
!
!                default value: gamma = 1
!
!   6   ip     - spatial smoothing parameter used when mmpde = 2 or 3
!
!                permitted values: 0 <= ip
!
!                recommended values: 2 <= ip <= 4
!
!                default value: ip = 2
!
!   7   stpmax - the allowed maximal time stepsize for the time
!                integration
!
!                permitted values: 0 < stpmax
!
!                default value: no restriction on time stepsize
!
!   8   xl     - the left end of the initial phyical region, i.e. xl =
!                x_l(t_a). this value is used when job = 2 and 3
!
!   9   xr     - the right end of the initial phyical region, i.e. xr =
!                x_r(t_a). this value is used when job = 2 and 3
!
!   10 - 20    - not used
!
! iflag:o   (integer) an indicator reporting the code's performance
!
!           iflag = +1:  successful run and default values of some of
!                        parameters par(1) -- par(7) have been used
!           iflag =  0:  successful run
!           iflag = -1:  at least one of the workspace arrays, rwork
!                        or iwork, has a length less than that required
!           iflag = -2:  at least one of the parameters, npde, npts,
!                        atol, rtol or ntouta, is incorrectly input
!           iflag = -3:  time integration stops since istop is set
!                        less than zero by the user in defout
!           iflag = -4:  time integration stops since idid <= -2 in
!                        dassl (see documentation in dassl)
!
! 2. workspace arrays
!
! rwork     (real array) workspace
!
! lrw:i     (integer) dimension of rwork. it must be at least
!           60 + 6*npde + (22+3*ip)*npts + (55+12*ip)*npts*npde
!           + (24+12*ip)*npts*npde**2, where ip := par(6)
!           (if ipar(3)=mmpde.ne.2 or 3, then one can consider
!           ip to be 0 in this lower bound for lrw)
!
! iwork     (integer array) workspace
!
! liw:i     (integer) dimension of iwork. it must be at least
!           20 + npts + 2*npts*npde
!
! 3. subroutines for problem definition
!
! defpde    define the physical pdes
! defbcp    define the physical boundary conditions
! defivs    define the initial values for the physical solution and
!           its first spatial derivative
! defmsh    define the initial mesh
! defbcm    define the boundary conditions for the mesh
! defmnt    define the monitor function
!
! see section 5 below for a detailed description of these subroutines
!
! 4. subroutines for solution output
!
! defout    output the physical solution and mesh
! deffnc    interpolate the physical solutions and first spatial
!           derivative via cubic hermite interpolation. this is
!           provided for user evaluation of solutions at non-mesh
!           points
!
! see section 5 below for a detailed description of these subroutines
!
! 5. subroutines (i: included, u: user must provide)
!
!    note : no input variable values should be changed in the user
!           provided subroutines
!
! defpde:u  subroutine defpde (npde, index, t, x, u, ux, ut, uxt, fg)
!           integer          npde, index
!           double precision t, x, u(npde), ux(npde), ut(npde),
!                            uxt(npde), fg(npde)
!
!           this subroutine is used to define the physical pdes. the f
!           and g terms, given by fg, must be defined in terms of input
!           variables npde, t, x, u, ux, ut and uxt by
!
!           fg := f(t, x, u, ux, ut, uxt)      if index < 0
!           fg := g(t, x, u, ux, ut, uxt)      if index > 0
!
! defbcp:u  subroutine defbcp(npde, index, t, x, xt, u, ux, uxx, ut,
!                             uxt, res)
!           integer          npde, index
!           double precision t, x, xt, u(npde), ux(npde), uxx(npde),
!                            ut(npde), uxt(npde), res(npde)
!
!           this subroutine is used to define the boundary conditions
!           for the physical solution variables. the residual res must
!           be defined in terms of input variables npde, index, t, x,
!           xt, u, ux uxx, ut and uxt by
!
!           res(k) := b^l_k(t, x, xt, u, ux, uxx, ut, uxt) if index < 0
!           res(k) := b^r_k(t, x, xt, u, ux, uxx, ut, uxt) if index > 0
!
!           for k = 1, ..., npde
!
! defivs:u  subroutine defivs(npde, x, u, ux)
!           integer          npde
!           double precision x, u(npde), ux(npde)
!
!           this subroutine is used to define the initial values for
!           the physical solution and its first spatial derivative,
!           where u and ux must be defined as functions of x (input
!           variable)
!
! defmsh:u  subroutine defmsh(npts, xmesh)
!           integer          npts
!           double precision xmesh(npts)
!
!           this subroutine is used to define the initial mesh, where
!           the entire initial mesh xmesh(1) < xmesh(2) < ...
!           < xmesh(npts) must be defined if job = 1.
!
! defbcm:u  subroutine defbcm(npde, index, t, x, xt, u, ux, uxx, ut,
!                             uxt, res)
!           integer          npde, index
!           double precision t, x, xt, u(npde), ux(npde), uxx(npde),
!                            ut(npde), uxt(npde), res
!
!           this subroutine is used to define the boundary conditions
!           for the mesh (the coordinate transformation). the residual
!           res must be defined in terms of input variables npde, index,
!           t, x, xt, u, ux, uxx, ut and uxt by
!
!           res := b^l_m(t, x, xt, u, ux, uxx, ut, uxt)   if index < 0
!           res := b^r_m(t, x, xt, u, ux, uxx, ut, uxt)   if index > 0
!
!           note: for fixed endpoints, simply use res := x - x^l if
!           index < 0 and res := x - x^r if index > 0
!
! defmnt:u  subroutine defmnt(npde, t, x, u, ux, uxx, ut, uxt, fmntr)
!           integer          npde
!           double precision t, x, u(npde), ux(npde), uxx(npde),
!                            ut(npde), uxt(npde), fmntr
!
!           this subroutine is used to define the monitor function. the
!           value fmntr of the monitor function must be defined in
!           terms of input variables npde, t, x, u, ux, uxx, ut and uxt
!
!           note: for the arclength monitor function,
!           fmntr := sqrt( 1 + ux(1)^2 + ... + ux(npde)^2 )
!
! defout:u  subroutine defout(npde, npts, t, xmesh, xmesht, u, ux, ut,
!                             tstep, touta, ntouta, istop, index, nts)
!           integer          npde, npts, i, ntouta, istop, index, nts
!           double precision xmesh(npts), xmesht(npts), u(npts,npde),
!                            ux(npts,npde), ut(npts,npde), tstep,
!                            touta(ntouta)
!
!           this subroutine is used to output the mesh points (xmesh)
!           and mesh speeds (xmesht) and the solution values (u, ut, ux)
!           at time t, where u(i,k) := u_k(x_i(t),t) (k=1,...,mpdes;
!           i=1,...,npts), tstep is the current time stepsize and touta
!           is the array of output times prescribed by the user. the use
!           can stop the time integration when some condition is met by
!           setting istop < 0.  c           index is a variable such tha
!
!           the code generates a mesh                   when index < 0
!                               (for job = 2 or 3)
!           the code is solves the physical pde         when index > 0
!                               (for job = 1 or 2)
!
!           all variables are input variables
!
! deffnc:i  subroutine deffnc(x, k, npde, npts, xmesh, u, ux, uval,
!                             uxval)
!           integer          k, npde, npts
!           double precision x, xmesh(npts), u(npts,npde),
!                            ux(npts,npde), uval, uxval
!
!           this subroutine provides the value of the k_th component
!           of the solution and its first spatial derivative at x,
      subroutine movcol (eqn, atol, rtol, touta, iflag)
!           using a cubic hermite spline which interpolates the
!           given solution values (u and ux) at mesh points (xmesh)
!           (see defout). it is provided for convenient user
!           evaluation of the solution and its first spatial
!           derivative at non-mesh points
!
! movcl1:i  this is the key subroutine (movcol is simply a user
!           interface subroutine)
!
! resode:i  defines the ode system for dassl
!
! jacode:i  defines a dummy routine for dassl
!
! respde:i  computes the residuals for the discretizations of the
!           physical pdes and their bcs
!
! respd1:i  computes the residuals for the discretizations of the
!           (artificial) pde and bcs during the mesh generation
!
! evlmnt:i  computes the monitor function
!
! smtmnt:i  smooths the monitor function (only used for mmpde = 2
!           and 3)
!
! resmeq:i  computes the residuals of the discrete mesh equations and
!           their bcs
!
! drvtvs:i  computes values of u, u_x, u_xx, u_t and u_xt for the
!           cubic hermite interpolatory spline solution at point x
!
! fncshp:i  defines the shape function and derivative values for
!           cubic hermite polynomial interpolation
!
! zermch:i  computes the machine unit roundoff in double precision
!
! slnout:i  driver for the subroutine defout for solution output
!
! ddassl    double precision version of dassl. this and certain
!           subroutines (from linpack) used by dassl must be included.
!           for a description see [petzold]
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
! end of the description of variables                                 c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!...this subroutine is the code's user interface which checks input
!...parameters, sets the default values, simplifies the declaration
!...of workspace arrays and calls the main subroutine movcl1
!
!
      implicit real (8)(a - h, o - z)
      real(8) :: rtol, atol
      class(problem_movcol), target :: eqn
      real(8) :: touta(:)

      integer :: liw, lrw, ntouta, npts, npde, i
      real(8) :: tau
      real(8), pointer :: par(:), rwork(:)
      integer, pointer :: iwork(:)

      par   => eqn%par
      iwork => eqn%iwork
      rwork => eqn%rwork

      ntouta = size(touta)
      liw = size(iwork)
      lrw = size(rwork)
      npts = eqn%npts
      npde = eqn%npde
!
      zrmch = epsilon(dum)
      iflag = 0
!
!...check basic parameters job, inform, mmpde, tau, gamma, ip and
!...stpmax. if their values are not within the permitted ranges,
!...default values are set and iflag = 1.
!
      job = par (1) + zrmch
      if (par (1) .le. - zrmch) job = - 1
      if (job.le.0.or.job.ge.4) then
         iflag = 1
         job = 2
      endif
!
      inform = par (2) + zrmch
      if (par (2) .le. - zrmch) inform = - 1
      if (inform.le.0) then
         iflag = 1
         inform = 6
      endif
!
      mmpde = par (3) + zrmch
      if (par (3) .le. - zrmch) mmpde = - 1
      if (mmpde.le.0.or.mmpde.ge.5) then
         iflag = 1
         mmpde = 4
      endif
!
      tau = par (4) + zrmch
      if (tau.le.zero) then
         iflag = 1
         tau = tau1
      endif
!
      gamma = par (5) + zrmch
      if (gamma.le.zero) then
         iflag = 1
         gamma = one
      endif
!
      ip = par (6) + zrmch
      if (par (6) .le. - zrmch) ip = - 1
      if (ip.lt.0) then
         iflag = 1
         ip = 2
      endif
!
      stpmax = par (7)
      if (stpmax.le.zrmch) then
         iflag = 1
         stpmax = one / zrmch
      endif
!
!...check basic parameters npde, npts, atol, rtol, ntouta,
!...xl and xr
!
      if (npde.le.0) iflag = - 1
      if (npts.le.0) iflag = - 1
      if (abs (atol) + abs (rtol) .le.zrmch) iflag = - 1
      if (ntouta.le.1) iflag = - 1
      if (job.eq.2.or.job.eq.3) then
         xl = par (8)
         xr = par (9)
         if (xl.gt.xr) iflag = - 1
      endif
!
!...if iflag = - 1, at least one of the parameters npde, npts,
!...atol, rtol and ntouta has been incorrectly input, so the
!...computation is halted
!
      if (iflag.eq. - 1) then
         write (inform, 99910) iflag
         write (inform, 99990) npde, npts, job, inform, mmpde, tau,     &
         gamma, ip, atol, rtol, stpmax, ntouta
         return
      endif
!
!...check the dimensions of the workspace arrays
!
      neq = (2 * npde+1) * npts
      if (mmpde.eq.2.or.mmpde.eq.3) then
         lrw0 = 40 + 12 * neq + 3 * (ip + 2) * (2 * npde+1) * neq
      else
         lrw0 = 40 + 12 * neq + 3 * (0 + 2) * (2 * npde+1) * neq
      endif
      lrw1 = 20 + 6 * npde+ (3 * npde+2) * npts
      if (lrw.lt.2 * npde+lrw0 + lrw1) iflag = - 2
      if (liw.lt.20 + npts * (2 * npde+1) ) iflag = - 2
!
!...when iflag = - 2, at least one of the dimensions of the
!...workspace arrays rwork and iwork is incorrectly input,
!...so the computation is stopped
!
      if (iflag.eq. - 2) then
         write (inform, 99920) iflag
         write (inform, 99990) npde, npts, job, inform, mmpde, tau,     &
         gamma, ip, atol, rtol, stpmax, ntouta
         return
      endif
!
!...call movcl1
!
      open(unit=111, file="par.dat")
      do i = 1, 20
         write(111,"(i5,f10.5)") i, par(i)
      end do
      close(111)

      call movcl1 (eqn, npde, npts, job, inform, iflag, xl, xr, mmpde,&
           & tau, gamma, ip, atol, rtol, stpmax, touta, ntouta,&
           & rwork (1            : neq                      ),&      !y
           & rwork (1+neq        : neq        + neq         ),&      !ydot
           & rwork (1+2*neq      : 2*neq      + lrw0        ),&      !rwk
           & lrw0, iwork, liw,&
           & rwork (1+2*neq+lrw0 : 2*neq+lrw0 + lrw1 ), lrw1)        !rwk1
      print *, lrw, lrw0, lrw1
      !
      return
!
99910 format(                                                           &
     &/' *****************************************************'         &
     &/' iflag = ',i3                                                   &
     &/' halt......at least one of the parameters npde, npts'           &
     &/' atol, rtol, ntouta, xl, and xr is incorrectly input.'          &
     &/' they must satisfy :'                                           &
     &/' npde >= 1                    npts >= 1'                        &
     &/' |atol| + |rtol| > 0          ntouta >= 2'                      &
     &/' xl <= xr (when job = 2 and 3)'                                 &
     &/' *****************************************************'/)
99920 format(                                                           &
     &/' *****************************************************'         &
     &/' iflag = ',i3                                                   &
     &/' halt......at least one of the workspace arrays rwork'          &
     &/' and iwork has insufficient length. check the values'           &
     &/' of lrw and liw. also check the values for npde, npts'          &
     &/' and ip := par(6)'                                              &
     &/' *****************************************************'/)
99990 format(                                                           &
     &/' *****************************************************'         &
     &/'  npde =',i12,'  npts =',i12,                                   &
     &'   job =',i12,' inform=', i10                                    &
     &/' mmpde =',i12,'   tau =',e12.5,                                 &
     &' gamma =',e12.5,'    ip =',i10                                   &
     &/'  atol =',e12.5,'  rtol =', e12.5,                              &
     &' stpmax=',e12.4,' ntouta=',i10                                   &
     &/' *****************************************************'/)
      end subroutine movcol
!
!***********************************************************************
!***********************************************************************
!
      subroutine movcl1 (eqn, npde, npts, job, inform, iflag, xl, xr, mmpde, &
      tau, gamma, ip, atol, rtol, stpmax, touta, ntouta, y, ydot, rwk,  &
      lrw, iwk, liw, rwk1, lrw1)
!
!...this is the main subroutine.
!...the relationships between y, ydot, and (x, u) are
!...
!...       y(m*(i-1)+2*(k-1)+1) := u_k(x_i(t),t)
!...       y(m*(i-1)+2*(k-1)+2) := (u_k)_x(x_i(t),t)
!...       y(m*(i-1)        +m) := x_i(t)
!...
!...                               k = 1, ..., npde
!...                               i = 1, ..., npts
!...
!...       ydot := y_t
!...
!...and m = 2*npde+1
!...
!...job = 1: solve the physical pdes starting with the given mesh
!...job = 2: first generate a mesh corresponding to the given
!...         initial value of the solution and then solve
!...         the physical pdes starting with this mesh
!...job = 3: only generate a mesh corresponding to the given
!...         initial value of the solution
!
      implicit real (8)(a - h, o - z)
      class(problem_movcol) :: eqn
      real(8) :: rtol1(1), atol1(1)
      dimension iwk (liw), info (15), iwk1 (20)
      dimension touta (ntouta), y (npts * (2 * npde+1) ), rwk (lrw)
      dimension rwk1 (lrw1), ydot (npts * (2 * npde+1) )
      logical phypde
!
!...define some basic parameters
!
      m = 2 * npde+1
      neq = m * npts
      m11 = 20
      m12 = 20 + 1 * npde
      m13 = 20 + 2 * npde
      m14 = 20 + 3 * npde
      m15 = 20 + 4 * npde
      m16 = 20 + 5 * npde
      m21 = 20 + 6 * npde
      m22 = 20 + 6 * npde+1 * npts
      m31 = 20 + 6 * npde+2 * npts
      m32 = 20 + 6 * npde+2 * npts + 1 * npde * npts
      m33 = 20 + 6 * npde+2 * npts + 2 * npde * npts
      zrmch = epsilon (1.0d0)
!     phypde = 'true' for solving the physical pdes
!     phypde = 'false' for generating a mesh
      if (job.eq.1) phypde = .true.
      if (job.eq.2) phypde = .false.
      if (job.eq.3) phypde = .false.
!
!...compute the initial mesh
!
!     for solving the physical pde
      if (phypde) then
         call eqn%defmsh(eqn%x)
      endif
!
!     the mesh generation case: starting from a uniform mesh
      if (.not.phypde) then
         do i = 1, npts
            y (m * i) = xl + (xr - xl) * (i - 1) / (npts - 1)
         end do
      endif

!
!...print the headers
!
      write (inform, 99920)
!  this is the entry point for one of the gotos
   20 if (phypde) then
         write (inform, 99930)
      else
         write (inform, 99940)
      endif

!
!...values for t, y and ydot which are consistent with the ode
!...system (arising in the method of lines procedure) are
!...needed initially to ensure convergence of ddassl at the first
!...step of the numerical integration. these initial values are
!...first obtained
!
!...the initial time derivatives ydot are computed by taking two
!...time integration steps with the first-order backward
!...differentiation formula and suppressing ddassl's time
!...truncation error control during these steps
!
!...compute the initial values of the physical solution
!
      if (phypde) then
         print *, "phypde = ",phypde
!        the case of solving the physical pdes
         do i = 1, npts
            call eqn%defivs(eqn%x(i),eqn%u(:,i),eqn%ux(:,i))
         end do
      else
!        the mesh generation case
         eqn%u  = 0.0_8
         eqn%ux = 0.0_8
      endif

      close(111)
!
!...define the input parameters for subroutine resode
!
      iwk1 (1) = npde
      iwk1 (2) = npts
      iwk1 (3) = mmpde
      iwk1 (4) = ip
      if (phypde) then
         iwk1 (5) = + 1
      else
         iwk1 (5) = - 1
      endif
      rwk1 (1) = gamma
      rwk1 (2) = tau
!
!...define the input parameters for ddassl
!
      do i = 1, 15
         info (i) = 0
      end do

      info (3) = 1
      info (6) = 1
!        define the lower and upper bandwidths for the jacobian
      if (mmpde.eq.2.or.mmpde.eq.3) then
         iwk (1) = min (m * (ip + 2), neq - 1)
      else
         iwk (1) = min (m * (0 + 2), neq - 1)
      endif
      iwk (2) = iwk (1)
      info (7) = 1
!        define the maximal time stepsize
      rwk (2) = min (sqrt (zrmch), atol, rtol)
      info (8) = 1
!        define the initial time stepsize
      rwk (3) = rwk (2) * half
      info (9) = 1
!        define the order of the bdf method (1st order)
      iwk (3) = 1
      info (11) = 1
!        define the initial values for ydot
      do i = 1, neq
         ydot (i) = zero
      end do
!
!...take two time integration steps using ddassl
!
      if (phypde) then
         t = touta (1)
         tout = touta (1) + one
      else
         t = zero
         tout = one
      endif
      atol1 = one
      rtol1 = one
      do itout = 1, 2

         ! open(unit=111, file="inside_resode.dat")
         ! ! print *, npde, npts, mmpde, ip, gamma, tau, m
         ! ! print *, zero, quart, half, one, two, three, four, six
         ! do i = 1, eqn%npts
         !    write(111,"(i10,(f30.10))") i, eqn%x(i)
         ! end do
         ! close(111)
         ! stop

         call ddassl (eqn, neq, t, y, ydot, tout, info, rtol1, atol1,&
         idid, rwk, lrw, iwk, liw, rwk1, iwk1)
      end do
!
!...now begin the actual numerical integration
!
!...the physical and mesh pdes are discretized in space with
!...a moving collocation method. the resulting ode system is
!...integrated in time with ddassl
!
!...reset the initial data for the mesh and solution
!
!     for the mesh
!
!     for solving the physical pde
      if (job.eq.1) then
         call eqn%defmsh ( eqn%x )
      endif
!
!     the mesh generation case: starting from a uniform mesh
      if (.not.phypde) then
         do i = 1, npts
            eqn%x(i) = xl + (xr - xl) * (i - 1) / (npts - 1)
         end do
      endif
!
!     for the physical solution
      if (phypde) then
!        the case of solving the physical pdes
         do i = 1, npts
            call eqn%defivs (eqn%x(i), eqn%u(:,i), eqn%ux(:,i))
         end do

      else
!        the mesh generation case
         eqn%u  = zero
         eqn%ux = zero
      endif
!
!...define the input parameters for subroutine resode
!
      iwk1 (1) = npde
      iwk1 (2) = npts
      iwk1 (3) = mmpde
      iwk1 (4) = ip
      if (phypde) then
         iwk1 (5) = + 1
      else
         iwk1 (5) = - 1
      endif
      rwk1 (1) = gamma
      rwk1 (2) = tau
!
!...define the input parameters for ddassl
!
      do i = 1, 15
         info (i) = 0
      end do
      info (3) = 1
      info (6) = 1
!        define the lower and upper bandwidths for the jacobian
      if (mmpde.eq.2.or.mmpde.eq.3) then
         iwk (1) = min (m * (ip + 2), neq - 1)
      else
         iwk (1) = min (m * (0 + 2), neq - 1)
      endif
      iwk (2) = iwk (1)
      info (7) = 1
!        define the maximal time stepsize
      rwk (2) = stpmax
      info (8) = 1
!        define the initial time stepsize
      temp = abs (touta (ntouta) )
      if (.not.phypde) temp = one
      rwk (3) = min (stpmax, rtol, atol, temp, one)
      info (11) = 1
!
!...print the initial solution
!
      if (phypde) then
         t = touta (1)
         index = + 1
      else
         t = zero
         index = - 1
      endif
      tstep = zero
      istop = 0
      nts = 0
      call eqn%slnout (npde, npts, t, y, ydot, touta, ntouta, tstep, istop, &
      index, nts, rwk1 (m11 + 1), rwk1 (m12 + 1), rwk1 (m13 + 1),       &
      rwk1 (m14 + 1), rwk1 (m15 + 1), rwk1 (m21 + 1), rwk1 (m22 + 1),   &
      rwk1 (m31 + 1), rwk1 (m32 + 1), rwk1 (m33 + 1))
!     if istop < 0, stop the computation
      if (istop.lt.0) then
         iflag = - 3
         write (inform, 99950) iflag
         write (inform, 99990) npde, npts, job, inform, mmpde, tau,     &
         gamma, ip, atol, rtol, stpmax, ntouta
         return
      endif
!
!...perform the time integration with ddassl
!
      if (phypde) then
         itoutm = ntouta
      else
         itoutm = 2
      endif
      do itout = 2, itoutm
         if (phypde) then
            tout = touta (itout)
         else
            tout = one
         endif
  110    call ddassl (eqn, neq, t, y, ydot, tout, info, rtol1, atol1,  &
         idid, rwk, lrw, iwk, liw, rwk1, iwk1)
!
!...print the solution at t
!
         if (phypde) then
            index = + 1
         else
            index = - 1
         endif
         tstep = rwk (7)
         istop = 0
         nts = iwk (11)
         call eqn%slnout (npde, npts, t, y, ydot, touta, ntouta, tstep,     &
         istop, index, nts, rwk1 (m11 + 1), rwk1 (m12 + 1), rwk1 (m13 + &
         1), rwk1 (m14 + 1), rwk1 (m15 + 1), rwk1 (m21 + 1), rwk1 (m22 +&
         1), rwk1 (m31 + 1), rwk1 (m32 + 1), rwk1 (m33 + 1))
!     if istop < 0, stop the computation
         if (istop.lt.0) then
            iflag = - 3
            write (inform, 99950) iflag
            write (inform, 99990) npde, npts, job, inform, mmpde, tau,  &
            gamma, ip, atol, rtol, stpmax, ntouta
            return
         endif
!
!...check the performance of the time integrator
!
         if (idid.eq.1) goto 110
!        a step was successfully taken in the intermediate-output
!        mode. ddassl has not yet reached tout. continue the time
!        integration
         if (idid.eq. - 1) then
!        a large amount of work has been expended. reset the
!        limitation and continue the time integration
            info (1) = 1
            goto 110
         endif
         if (idid.ge. - 1.and.idid.le.3) then
!        the time integration is successful to tout. output the
!        performance information
            write (inform, 99970) t, idid, iwk (14), iwk (15), iwk (11),&
            iwk (13), iwk (12), iwk (8), rwk (7), rwk (3)
            info (1) = 1
         else
!        the time integration fails and the computation is stopped
            iflag = - 4
            write (inform, 99980) iflag, idid
            write (inform, 99990) npde, npts, job, inform, mmpde, tau,  &
            gamma, ip, atol, rtol, stpmax, ntouta
            return
         endif
      end do
!
!...if job = 2, go back to solve the physical pdes
!
      if (job.eq.2.and. (.not.phypde) ) then
         phypde = .true.
         goto 20
      endif
!
!...the computation is completed. output summary information and
!...return
!
      write (inform, 99990) npde, npts, job, inform, mmpde, tau, gamma, &
      ip, atol, rtol, stpmax, ntouta
      return
99920 format(                                                           &
     &/' *****************************************************'         &
     &/' ***                                               ***'         &
     &/' ***                  movcol                       ***'         &
     &/' ***                                               ***'         &
     &/' *****************************************************'         &
     &/' nts  : time steps taken'                                       &
     &/' jac  : jacobian evaluations'                                   &
     &/' etf  : error test failures'                                    &
     &/' cfn  : convergence test failures in newton iteration'          &
     &/' idid : indicator of the performance of dassl'                  &
     &/' nres : number of the ode residual evaluations'                 &
     &/' ord  : order of the method for the time integration'           &
     &/' steplast : time step used in last step'                        &
     &/' stepnow :  time step used in current step'                     &
     &/' *****************************************************'/)
99930 format(                                                           &
     &/' *****************************************************'         &
     &/' ***                                               ***'         &
     &/' ***           solve the physical pdes             ***'         &
     &/' ***                                               ***'         &
     &/' *****************************************************'/)
99940 format(                                                           &
     &/' *****************************************************'         &
     &/' ***                                               ***'         &
     &/' ***       generate an equidistributed mesh        ***'         &
     &/' ***                                               ***'         &
     &/' *****************************************************'/)
99950 format(                                                           &
     &/' *****************************************************'         &
     &/' iflag = ',i3                                                   &
     &/' halt......istop is set to be < 0 by the user in defout'        &
     &/' *****************************************************'/)
99970 format(                                                           &
     &/' ---------------------- information for time integration',      &
     &' ---------------------'/                                         &
     &'  time =',e12.5,'  idid =',i12,'   etf =',i12,'   cfn =',i10/    &
     &'   nts =',i12,'   jac =',i12,'  nres =',i12,'   ord =',i10/      &
     &' steplast=',e12.5,' stepnow=',e12.5,                             &
     &/' -----------------------------------------------------',        &
     &'------------------------')
99980 format(                                                           &
     &/' *****************************************************'         &
     &/' iflag = ',i3,' idid =', i3                                     &
     &/' halt......idid <= -2 in dassl (see documentation in'           &
     &/' dassl)'                                                        &
     &/' *****************************************************'/)
99990 format(                                                           &
     &/' *****************************************************',        &
     &'************************'                                        &
     &/'  npde =',i12,'  npts =',i12,                                   &
     &'   job =',i12,' inform=', i10                                    &
     &/' mmpde =',i12,'   tau =',e12.5,                                 &
     &' gamma =',e12.5,'    ip =',i10                                   &
     &/'  atol =',e12.5,'  rtol =', e12.5,                              &
     &' stpmax=',e12.4,' ntouta=',i10                                   &
     &/' *****************************************************',        &
     &'************************'/)
      end subroutine movcl1
!
!***********************************************************************
!***********************************************************************
!
      subroutine jacode (eqn, t, y, ydot, pd, cj, rwk, iwk)
!
!...this dummy routine is required by ddassl
!
      class(problem_movcol) :: eqn
      integer, dimension(*) :: iwk
      real(8) :: t, cj
      real(8), dimension(*) :: y, ydot, rwk
      real(8), dimension(10,*) :: pd
      return
      end subroutine jacode
!
!***********************************************************************
!***********************************************************************
!
      subroutine resode (eqn, t, y, ydot, res, ires, rwk, iwk)
!
!...this subroutine used by ddassl computes the residuals for
!...the ode system resulting from the discretization of the
!...physical and mesh pdes.
!...the relationships between y, ydot, res c...and (x, u) are
!...
!...       y(m*i)               := x_i(t)
!...       y(m*(i-1)+2*(k-1)+1) := u_k(x_i(t),t)
!...       y(m*(i-1)+2*(k-1)+2) := (u_k)_x(x_i(t),t)
!...
!...                               k = 1, ..., npde
!...                               i = 1, ..., npts
!...
!...       ydot := y_t
!...       res  := f(t, x, u, ux, ut, uxt)-[ g(t, x, u, ux, ut, uxt) ]_x
!...
!...and m = 2*npde+1
!
      class(problem_movcol) :: eqn
      integer :: ires
      integer, dimension(*) :: iwk
      real(8) :: t
      real(8), dimension(*) :: y, ydot, res, rwk

      ! local variables
      real(8) :: tau, gamma
      logical phypde
!
!...define some basic parameters
!

      npde = iwk (1)
      npts = iwk (2)
      mmpde = iwk (3)
      ip = iwk (4)
      if (iwk (5) .gt.0) then
!        the case for solving the physical pdes
         phypde = .true.
      else
!        the mesh generation case
         phypde = .false.
      endif
      gamma = rwk (1)
      tau = rwk (2)

      m = 2 * npde+1
      m11 = 20
      m12 = 20 + 1 * npde
      m13 = 20 + 2 * npde
      m14 = 20 + 3 * npde
      m15 = 20 + 4 * npde
      m16 = 20 + 5 * npde
      m21 = 20 + 6 * npde
      m22 = 20 + 6 * npde+npts
!
!...compute residuals of the physical pdes and their bcs
!
      if (phypde) then
!        the case of solving the physical pdes
         call eqn%respde (npde, npts, t, y, ydot, res, rwk (m11 + 1),       &
         rwk (m12 + 1), rwk (m13 + 1), rwk (m14 + 1), rwk (m15 + 1),    &
         rwk (m16 + 1))
      else
!        the mesh generation case
         call eqn%respd1 (npde, npts, t, y, ydot, res, rwk (m11 + 1),       &
         rwk (m12 + 1), rwk (m13 + 1), rwk (m14 + 1), rwk (m15 + 1) )
      endif
!
!...compute and smooth the monitor function
!
!     rwk(m21+i) := fmntr(i)
      if (mmpde.ge.2) call eqn%evlmnt (npde, npts, t, y, ydot, rwk (m21 + 1)&
      , rwk (m11 + 1), rwk (m12 + 1), rwk (m13 + 1), rwk (m14 + 1),     &
      rwk (m15 + 1) )
      if (mmpde.eq.2.or.mmpde.eq.3) call smtmnt (npts, ip, rwk (m21 + 1)&
      , rwk (m22 + 1) )
!
!...compute residuals for the discrete mesh equations and their bcs.
!
      call eqn%resmeq ( npde, npts, mmpde, tau, gamma, y, ydot, res, rwk (   &
      m21 + 1), rwk (m11 + 1), rwk (m12 + 1), rwk (m13 + 1), rwk (m14 + &
      1), rwk (m15 + 1))
      if (.not.phypde) then
!        the mesh generation case (for which boundaries are fixed)
         res (m) = ydot (m)
         res (m * npts) = ydot (m * npts)
      endif
!
      return
      end subroutine resode
!
!***********************************************************************
!***********************************************************************
!
      subroutine respde (eqn, npde, npts, t, y, ydot, res, u, ux, uxx, ut,   &
      uxt, rw)
!
!...this subroutine computes the residuals for the discretizations
!...of the physical pdes and their bcs.
!...the relationships between y, ydot, res, and (x, u) are
!...
!...       y(m*i)               := x_i(t)
!...       y(m*(i-1)+2*(k-1)+1) := u_k(x_i(t),t)
!...       y(m*(i-1)+2*(k-1)+2) := (u_k)_x(x_i(t),t)
!...
!...                               k = 1, ..., npde
!...                               i = 1, ..., npts
!...
!...       ydot := y_t
!...       res  := f(t, x, u, ux, ut, uxt)-[ g(t, x, u, ux, ut, uxt) ]_x
!...
!...and m = 2*npde+1
!
      implicit real (8)(a - h, o - z)
      class(problem_movcol) :: eqn
      dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      dimension res ( (2 * npde+1) * npts), u (npde), ux (npde),        &
      uxx (npde)
      dimension ut (npde), uxt (npde), rw (npde), w (3, 2)
!
      m = 2 * npde+1
!
!...define the weights for the cell-average collocation approximation
!...to the right-hand-side term
!
      w (1, 1) = - one-two * (s2 - s1)
      w (2, 1) = four * (s2 - s1)
      w (3, 1) = one-two * (s2 - s1)
      w (1, 2) = - one+two * (s2 - s1)
      w (2, 2) = - four * (s2 - s1)
      w (3, 2) = one+two * (s2 - s1)
!
!...compute the residuals for the discretization of the physical pdes
!...at the interior (gauss) points
!
      do 50 i = 1, npts - 1
         h = y (m * (i + 1) ) - y (m * i)
!
!...compute the left-hand-side term f(...) using collocation
!
         do 20 j = 1, 2
            if (j.eq.1) x = y (m * i) + s1 * h
            if (j.eq.2) x = y (m * i) + s2 * h
            call drvtvs (npde, npts, i, x, y, ydot, u, ux, uxx, ut, uxt)
            call eqn%defpde (- 1, t, x, u, ux, ut, uxt, rw)
            do 10 k = 1, npde
               if (j.eq.1) res (m * (i - 1) + 2 * (k - 1) + 2) = rw (k)
               if (j.eq.2) res (m * (i) + 2 * (k - 1) + 1) = rw (k)
   10       end do
   20    end do
!
!...compute the right-hand-side term [ g(...) ]_x using cell-average
!...collocation
!
         do 40 j = 1, 3
            if (j.eq.1) x = y (m * i)
            if (j.eq.2) x = y (m * i) + half * h
            if (j.eq.3) x = y (m * (i + 1) )
            call drvtvs (npde, npts, i, x, y, ydot, u, ux, uxx, ut, uxt)
            call eqn%defpde (+ 1, t, x, u, ux, ut, uxt, rw)
            do 30 k = 1, npde
               res (m * (i - 1) + 2 * (k - 1) + 2) = res (m * (i - 1)   &
               + 2 * (k - 1) + 2) - rw (k) * w (j, 1) / h
               res (m * (i) + 2 * (k - 1) + 1) = res (m * (i) + 2 *     &
               (k - 1) + 1) - rw (k) * w (j, 2) / h
   30       end do
   40    end do
   50 end do
!
!...compute the bcs at x = x^l (left end) using collocation
!
      x = y (m)
      xt = ydot (m)
      call drvtvs (npde, npts, 1, x, y, ydot, u, ux, uxx, ut, uxt)
      call eqn%defbcp (- 1, t, x, xt, u, ux, uxx, ut, uxt, rw)
      do 60 k = 1, npde
         res (2 * (k - 1) + 1) = rw (k)
   60 end do
!
!...compute the bcs at x = x^r (right end) using collocation
!
      x = y (m * npts)
      xt = ydot (m * npts)
      call drvtvs (npde, npts, npts - 1, x, y, ydot, u, ux, uxx, ut,    &
      uxt)
      call eqn%defbcp (+ 1, t, x, xt, u, ux, uxx, ut, uxt, rw)
      do 70 k = 1, npde
         res (m * (npts - 1) + 2 * (k - 1) + 2) = rw (k)
   70 end do
!
      return
      end subroutine respde
!
!***********************************************************************
!***********************************************************************
!
      subroutine respd1 (eqn, npde, npts, t, y, ydot, res, u, ux, uxx, ut,   &
      uxt)
!
!...for the mesh generation case, the (artificial) physical pde is
!...defined as
!...          u_t = u(x), 0 < t <= 1
!...
!...where u(x) is the given initial solution defined in subroutine
!...defivs
!...
!...this subroutine computes the residuals for the discretization
!...to this pde and for the bcs.
!...the relationships between y, ydot, res, and (x, u) are
!...
!...       y(m*i)               := x_i(t)
!...       y(m*(i-1)+2*(k-1)+1) := u_k(x_i(t),t)
!...       y(m*(i-1)+2*(k-1)+2) := (u_k)_x(x_i(t),t)
!...
!...                               k = 1, ..., npde
!...                               i = 1, ..., npts
!...
!...       ydot := y_t
!...       res  := u_t - u(x)
!...
!...and m = 2*npde+1
!
      implicit real (8)(a - h, o - z)
      class(problem_movcol) :: eqn
      dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      dimension res ( (2 * npde+1) * npts), u (npde), ux (npde),        &
      uxx (npde)
      dimension ut (npde), uxt (npde)
!
      m = 2 * npde+1
!
!...compute the residuals for the discretizations to the pdes at
!...the interior (gauss) points
!
      do 30 i = 1, npts - 1
         h = y (m * (i + 1) ) - y (m * i)
         ht = ydot (m * (i + 1) ) - ydot (m * i)
         do 20 j = 1, 2
            if (j.eq.1) then
               x = y (m * i) + s1 * h
               xt = ydot (m * i) + s1 * ht
            endif
            if (j.eq.2) then
               x = y (m * i) + s2 * h
               xt = ydot (m * i) + s2 * ht
            endif
            call drvtvs (npde, npts, i, x, y, ydot, u, ux, uxx, ut, uxt)
            call eqn%defivs (x, u, uxx)
            do 10 k = 1, npde
!
!...1st order upwind treatment
!
               ut (k) = ut (k) + ux (k) * xt
               temp0 = y (m * (i - 1) + 2 * (k - 1) + 1) * fncshp (1, 0,&
               zero) + y (m * (i - 1) + 2 * (k - 1) + 2) * fncshp (2, 0,&
               zero) * h + y (m * (i) + 2 * (k - 1) + 1) * fncshp (3, 0,&
               zero) + y (m * (i) + 2 * (k - 1) + 2) * fncshp (4, 0,    &
               zero) * h
               temp1 = y (m * (i - 1) + 2 * (k - 1) + 1) * fncshp (1, 0,&
               s1) + y (m * (i - 1) + 2 * (k - 1) + 2) * fncshp (2, 0,  &
               s1) * h + y (m * (i) + 2 * (k - 1) + 1) * fncshp (3, 0,  &
               s1) + y (m * (i) + 2 * (k - 1) + 2) * fncshp (4, 0, s1)  &
               * h
               temp2 = y (m * (i - 1) + 2 * (k - 1) + 1) * fncshp (1, 0,&
               s2) + y (m * (i - 1) + 2 * (k - 1) + 2) * fncshp (2, 0,  &
               s2) * h + y (m * (i) + 2 * (k - 1) + 1) * fncshp (3, 0,  &
               s2) + y (m * (i) + 2 * (k - 1) + 2) * fncshp (4, 0, s2)  &
               * h
               temp3 = y (m * (i - 1) + 2 * (k - 1) + 1) * fncshp (1, 0,&
               one) + y (m * (i - 1) + 2 * (k - 1) + 2) * fncshp (2, 0, &
               one) * h + y (m * (i) + 2 * (k - 1) + 1) * fncshp (3, 0, &
               one) + y (m * (i) + 2 * (k - 1) + 2) * fncshp (4, 0, one)&
               * h
               if (j.eq.1) then
                  if (xt.ge.zero) temp = (temp2 - temp1) / (s2 - s1)    &
                  / h
                  if (xt.le.zero) temp = (temp1 - temp0) / s1 / h
                  res (m * (i - 1) + 2 * (k - 1) + 2) = (ut (k) - temp *&
                  xt) - u (k)
               endif
               if (j.eq.2) then
                  if (xt.ge.zero) temp = (temp3 - temp2) / (one-s2)     &
                  / h
                  if (xt.le.zero) temp = (temp2 - temp1) / (s2 - s1)    &
                  / h
                  res (m * (i) + 2 * (k - 1) + 1) = (ut (k) - temp * xt)&
                  - u (k)
               endif
   10       end do
   20    end do
   30 end do
!
!...compute the residual for the collocation solution of the
!...pdes at x = x^l (left end)
!
      x = y (m)
      xt = ydot (m)
      call drvtvs (npde, npts, 1, x, y, ydot, u, ux, uxx, ut, uxt)
      call eqn%defivs (x, ux, uxx)
!     here ux represents u(x)
      do 40 k = 1, npde
         res (2 * (k - 1) + 1) = u (k) - t * ux (k)
   40 end do
!
!...compute the residual for the collocation solution of the
!...pdes at x = x^r (right end)
!
!
      x = y (m * npts)
      xt = ydot (m * npts)
      call drvtvs (npde, npts, npts - 1, x, y, ydot, u, ux, uxx, ut,    &
      uxt)
      call eqn%defivs (x, ux, uxx)
!     here ux represents u(x)
      do 50 k = 1, npde
         res (m * (npts - 1) + 2 * (k - 1) + 2) = u (k) - t * ux (k)
   50 end do
!
      return
      end subroutine respd1
!
!***********************************************************************
!***********************************************************************
!
      subroutine resmeq (eqn, npde, npts, mmpde, tau, gamma, y, ydot, res,   &
      fmntr, u, ux, uxx, ut, uxt)
!
!...this subroutine computes the residuals for the finite difference
!...approximations of the mesh equations c...and their bcs.
!...the relationships between y, ydot, res, and (x, u) are
!...
!...       y(m*i)               := x_i(t)
!...       y(m*(i-1)+2*(k-1)+1) := u_k(x_i(t),t)
!...       y(m*(i-1)+2*(k-1)+2) := (u_k)_x(x_i(t),t)
!...
!...                               k = 1, ..., npde
!...                               i = 1, ..., npts
!...
!...       ydot := y_t
!...
!...and m = 2*npde+1
!
      implicit real (8)(a - h, o - z)
      class(problem_movcol) :: eqn
      real(8) :: tau
      dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      dimension res ( (2 * npde+1) * npts), fmntr (npts), u (npde)
      dimension ux (npde), uxx (npde), ut (npde), uxt (npde)
!
      m = 2 * npde+1
!
!...mmpde = 1: fixed mesh
!
      if (mmpde.eq.1) then
         do 10 i = 1, npts
            res (m * i) = ydot (m * i)
   10    end do
         return
      endif
!
!...compute bcs for the coordinate transformation (the mesh)
!...for mpde = 2, 3, and 4
!
      x = y (m)
      xt = ydot (m)
      call drvtvs (npde, npts, 1, x, y, ydot, u, ux, uxx, ut, uxt)
      call eqn%defbcm ( - 1, t, x, xt, u, ux, uxx, ut, uxt, temp)
      res (m) = temp
      x = y (m * npts)
      xt = ydot (m * npts)
      call drvtvs (npde, npts, npts - 1, x, y, ydot, u, ux, uxx, ut,    &
      uxt)
      call eqn%defbcm ( + 1, t, x, xt, u, ux, uxx, ut, uxt, temp)
      res (m * npts) = temp
!
!...mmpde = 2: mmpde4 in [hrr1, hrr2] at the interior mesh points
!
      if (mmpde.eq.2) then
         do 20 i = 2, npts - 1
            temp1 = fmntr (i - 1)
            temp2 = fmntr (i)
            res (m * i) = temp2 * (y (m * (i + 1) ) - y (m * i) )       &
            - temp1 * (y (m * i) - y (m * (i - 1) ) ) + (temp2 *        &
            (ydot (m * (i + 1) ) - ydot (m * i) ) - temp1 * (ydot (m *  &
            i) - ydot (m * (i - 1) ) ) ) * tau
   20    end do
         return
      endif
!
!...mmpde = 3: mmpde6 in [hrr1, hrr2] at the interior mesh points
!
      if (mmpde.eq.3) then
         do 30 i = 2, npts - 1
            temp1 = fmntr (i - 1)
            temp2 = fmntr (i)
            res (m * i) = temp2 * (y (m * (i + 1) ) - y (m * i) )       &
            - temp1 * (y (m * i) - y (m * (i - 1) ) ) + ( (ydot (m *    &
            (i + 1) ) - ydot (m * i) ) - (ydot (m * i) - ydot (m *      &
            (i - 1) ) ) ) * tau
   30    end do
         return
      endif
!
!...mmpde = 4: mmpde(24) in [hr2], which is a smoothed version
!...of mmpde4 in terms of the mesh concentration function at the
!...interior mesh points
!
      xlambd = gamma * (gamma + one)
      if (mmpde.eq.4) then
         do 40 i = 2, npts - 1
            temp1 = fmntr (i - 1)
            temp2 = fmntr (i)
            ym1 = - (ydot (m * i) - ydot (m * (i - 1) ) ) / (y (m * i)  &
            - y (m * (i - 1) ) ) **2 * tau + one / (y (m * i) - y (m *  &
            (i - 1) ) )
            yp1 = - (ydot (m * (i + 1) ) - ydot (m * i) ) / (y (m *     &
            (i + 1) ) - y (m * i) ) **2 * tau + one / (y (m * (i + 1) ) &
            - y (m * i) )
            if (i - 2.ge.1) then
               ym3 = - (ydot (m * (i - 1) ) - ydot (m * (i - 2) ) )     &
               / (y (m * (i - 1) ) - y (m * (i - 2) ) ) **2 * tau + one &
               / (y (m * (i - 1) ) - y (m * (i - 2) ) )
            else
               ym3 = ym1
            endif
            if (i + 2.le.npts) then
               yp3 = - (ydot (m * (i + 2) ) - ydot (m * (i + 1) ) )     &
               / (y (m * (i + 2) ) - y (m * (i + 1) ) ) **2 * tau + one &
               / (y (m * (i + 2) ) - y (m * (i + 1) ) )
            else
               yp3 = yp1
            endif
            res (m * i) = (yp1 - xlambd * (yp3 - two * yp1 + ym1) )     &
            / temp2 - (ym1 - xlambd * (yp1 - two * ym1 + ym3) ) / temp1
   40    end do
         return
      endif
!
      return
      end subroutine resmeq
!
!***********************************************************************
!***********************************************************************
!
      subroutine drvtvs (npde, npts, i, x, y, ydot, u, ux, uxx, ut, uxt)
!
!...this subroutine computes the value (u) and partial derivatives (ux,
!...uxx, ut, uxt) of the cubic hermite interpolatory spline at any
!...point (x) in the i_th interval [x_i, x_(i+1)], i = 1, ..., npts-1,
!...by using nodal values y and ydot which are defined as
!...
!...       y(m*i)               := x_i(t)
!...       y(m*(i-1)+2*(k-1)+1) := u_k(x_i(t),t)
!...       y(m*(i-1)+2*(k-1)+2) := (u_k)_x(x_i(t),t)
!...
!...                               k = 1, ..., npde
!...                               i = 1, ..., npts
!...
!...       ydot := y_t
!...
!...and m = 2*npde+1
!
      implicit real (8)(a - h, o - z)
      dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      dimension u (npde), ux (npde), uxx (npde), ut (npde), uxt (npde)
      dimension temp (4, 0:2)
!
      m = 2 * npde+1
      h = y (m * (i + 1) ) - y (m * i)
      ht = ydot (m * (i + 1) ) - ydot (m * i)
      s = (x - y (m * i) ) / h
!
      do j = 1, 4
         do k = 0, 2
            temp (j, k) = fncshp (j, k, s)
         end do
      end do

!
      do 10 k = 1, npde
         u (k) = y (m * (i - 1) + 2 * (k - 1) + 1) * temp (1, 0)        &
         + y (m * (i - 1) + 2 * (k - 1) + 2) * temp (2, 0) * h + y (m * &
         (i) + 2 * (k - 1) + 1) * temp (3, 0) + y (m * (i) + 2 *        &
         (k - 1) + 2) * temp (4, 0) * h
         ux (k) = (y (m * (i - 1) + 2 * (k - 1) + 1) * temp (1, 1)      &
         + y (m * (i - 1) + 2 * (k - 1) + 2) * temp (2, 1) * h + y (m * &
         (i) + 2 * (k - 1) + 1) * temp (3, 1) + y (m * (i) + 2 *        &
         (k - 1) + 2) * temp (4, 1) * h) / h
         ut (k) = ydot (m * (i - 1) + 2 * (k - 1) + 1) * temp (1, 0)    &
         + ydot (m * (i - 1) + 2 * (k - 1) + 2) * temp (2, 0) * h +     &
         ydot (m * (i) + 2 * (k - 1) + 1) * temp (3, 0) + ydot (m *     &
         (i) + 2 * (k - 1) + 2) * temp (4, 0) * h + y (m * (i - 1)      &
         + 2 * (k - 1) + 2) * temp (2, 0) * ht + y (m * (i) + 2 *       &
         (k - 1) + 2) * temp (4, 0) * ht - ux (k) * (ydot (m * (i - 1)  &
         + m) + s * ht)
         uxx (k) = (y (m * (i - 1) + 2 * (k - 1) + 1) * temp (1, 2)     &
         + y (m * (i - 1) + 2 * (k - 1) + 2) * temp (2, 2) * h + y (m * &
         (i) + 2 * (k - 1) + 1) * temp (3, 2) + y (m * (i) + 2 *        &
         (k - 1) + 2) * temp (4, 2) * h) / h / h
         uxt (k) = (ydot (m * (i - 1) + 2 * (k - 1) + 1) * temp (1, 1)  &
         + ydot (m * (i - 1) + 2 * (k - 1) + 2) * temp (2, 1) * h +     &
         ydot (m * (i) + 2 * (k - 1) + 1) * temp (3, 1) + ydot (m *     &
         (i) + 2 * (k - 1) + 2) * temp (4, 1) * h + y (m * (i - 1)      &
         + 2 * (k - 1) + 2) * temp (2, 1) * ht + y (m * (i) + 2 *       &
         (k - 1) + 2) * temp (4, 1) * ht) / h - ux (k) * ht / h - uxx ( &
         k) * (ydot (m * (i - 1) + m) + s * ht)
   10 end do
!
      return
      end subroutine drvtvs
!
!***********************************************************************
!***********************************************************************
!
      real(8) function fncshp (j, k, s)
!
!...this function subroutine computes, at a point s in [0, 1], the
!...value for the k_th derivative (k =  0, 1 or 2) of the j_th shape
!...function (j = 1, 2, 3 or 4) for the cubic hermite interpolate
!
      implicit real (8)(a - h, o - z)
!
      if (j.eq.1) then
         if (k.eq.0) fncshp = (one+two * s) * (one-s) * (one-s)
         if (k.eq.1) fncshp = - six * s * (one-s)
         if (k.eq.2) fncshp = - six * (one-two * s)
      endif
      if (j.eq.2) then
         if (k.eq.0) fncshp = s * (one-s) * (one-s)
         if (k.eq.1) fncshp = (one-s) * (one-three * s)
         if (k.eq.2) fncshp = six * s - four
      endif
      if (j.eq.3) then
         if (k.eq.0) fncshp = (three-two * s) * s * s
         if (k.eq.1) fncshp = six * s * (one-s)
         if (k.eq.2) fncshp = six * (one-two * s)
      endif
      if (j.eq.4) then
         if (k.eq.0) fncshp = (s - one) * s * s
         if (k.eq.1) fncshp = s * (three * s - two)
         if (k.eq.2) fncshp = six * s - two
      endif
!
      return
      end function fncshp
!
!***********************************************************************
!***********************************************************************
!
      subroutine evlmnt (eqn, npde, npts, t, y, ydot, fmntr, u, ux, uxx, ut, &
      uxt)
!
!...this subroutine computes the monitor function value
!...fmntr(i) := c...m_(i+1/2).
!...the relationships between y, ydot, and (x, u) are
!...
!...       y(m*i)               := x_i(t)
!...       y(m*(i-1)+2*(k-1)+1) := u_k(x_i(t),t)
!...       y(m*(i-1)+2*(k-1)+2) := (u_k)_x(x_i(t),t)
!...
!...                               k = 1, ..., npde
!...                               i = 1, ..., npts
!...
!...       ydot := y_t
!...
!...and m = 2*npde+1
!
      implicit real (8)(a - h, o - z)
      class(problem_movcol) :: eqn
      dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      dimension u (npde), ux (npde), uxx (npde), ut (npde), uxt (npde)
      dimension fmntr (npts)
!
      m = 2 * npde+1
!
      do 20 i = 1, npts - 1
         fmntr (i) = zero
         h = y (m * (i + 1) ) - y (m * i)
         do 10 j = 1, 2
            if (j.eq.1) x = y (m * i) + s1 * h
            if (j.eq.2) x = y (m * i) + s2 * h
            call drvtvs (npde, npts, i, x, y, ydot, u, ux, uxx, ut, uxt)
            call eqn%defmnt (t, x, u, ux, uxx, ut, uxt, temp)
            fmntr (i) = fmntr (i) + half * temp
   10    end do
   20 end do
!
!...approximate the monitor function more accurately in the end
!...subintervals (i.e., for i = 1 and npts-1)
!
      do 40 i = 1, npts - 1, npts - 2
         fmntr (i) = half * fmntr (i)
         do 30 j = i, i + 1
            x = y (m * j)
            call drvtvs (npde, npts, i, x, y, ydot, u, ux, uxx, ut, uxt)
            call eqn%defmnt (t, x, u, ux, uxx, ut, uxt, temp)
            fmntr (i) = fmntr (i) + quart * temp
   30    end do
   40 end do
!
      return
      end subroutine evlmnt
!
!***********************************************************************
!***********************************************************************
!
      subroutine smtmnt (npts, ip, fmntr, rw)
!
!...this subroutine smooths the monitor function. it is called
!...only when mmpde = 2 or mmpde = 3
!
      implicit real (8)(a - h, o - z)
      dimension fmntr (npts), rw (npts)

      real(8) :: gamma = 2.0
!
      do 10 i = 1, npts - 1
         rw (i) = fmntr (i)
   10 end do
      do 30 i = 1, npts - 1
         temp = zero
         fmntr (i) = zero
         do 20 j = max (1, i - ip), min (npts - 1, i + ip)
            temp1 = (gamma / (gamma + one) ) **iabs (j - i)
            temp = temp + temp1
            fmntr (i) = fmntr (i) + rw (j) * temp1
   20    end do
         fmntr (i) = fmntr (i) / temp
   30 end do
!
      return
      end subroutine smtmnt
!
!***********************************************************************
!***********************************************************************
!
      subroutine deffnc (x, k, npde, npts, xmesh, u, ux, uval, uxval)
!
!...this subroutine computes the interpolatory value (uval) and first
!...spatial derivative (uxval) at any point (x) for the k_th component
!...of the solution by using the mesh points (xmesh) and corresponding
!...mesh values for u and ux
!
      implicit real (8)(a - h, o - z)
      dimension xmesh (npts), u (npts, npde), ux (npts, npde)
!
      do 10 i1 = 1, npts - 1
         if (xmesh (i1) .le.x.and.x.le.xmesh (i1 + 1) ) then
            i = i1
            goto 20
         endif
   10 end do
      if (x.le.xmesh (2) ) i = 1
      if (x.ge.xmesh (npts - 1) ) i = npts - 1
   20 h = xmesh (i + 1) - xmesh (i)
      s = (x - xmesh (i) ) / h
      uval = u (i, k) * fncshp (1, 0, s) + ux (i, k) * fncshp (2, 0, s) &
      * h + u (i + 1, k) * fncshp (3, 0, s) + ux (i + 1, k) * fncshp (4,&
      0, s) * h
      uxval = (u (i, k) * fncshp (1, 1, s) + ux (i, k) * fncshp (2, 1,  &
      s) * h + u (i + 1, k) * fncshp (3, 1, s) + ux (i + 1, k) * fncshp &
      (4, 1, s) * h) / h
      return
      end subroutine deffnc
!
!***********************************************************************
!***********************************************************************
!
      subroutine slnout (eqn, npde, npts, t, y, ydot, touta, ntouta, tstep,  &
      istop, index, nts, u, ux, uxx, ut, uxt, xmesh, xmesht, uu, uux,   &
      uut)
!
!...this subroutine outputs the solution values x_t, u, u_x,
!...and u_t at point x and time t in the format prescribed by
!...subroutine defout
!
      implicit real (8)(a - h, o - z)
      class(problem_movcol) :: eqn
      dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      dimension touta (ntouta), u (npde), ux (npde), uxx (npde)
      dimension ut (npde), uxt (npde), xmesh (npts), xmesht (npts)
      dimension uu (npts, npde), uux (npts, npde), uut (npts, npde)
!
      m = 2 * npde+1
      do 20 i = 1, npts
!
!...compute u, ux, ut at time t and x = x_i
!
         x = y (m * i)
         i1 = min (i, npts - 1)
         call drvtvs (npde, npts, i1, x, y, ydot, u, ux, uxx, ut, uxt)
!
!...store x_i, xdot_i, u(k), ux(k), ut(k) in xmesh(i), xmesht(i),
!...uu(i,k), uux(i,k) and uut(i,k), respectively
!
         xmesh (i) = x
         xmesht (i) = ydot (m * i)
         do 10 k = 1, npde
            uu (i, k) = u (k)
            uux (i, k) = ux (k)
            uut (i, k) = ut (k)
   10    end do
   20 end do
!
!...output the solution values and current time stepsize tstep
!
      call eqn%defout (npde, npts, t, xmesh, xmesht, uu, uux, uut, tstep,   &
      touta, ntouta, istop, index, nts)
!
      return
      end subroutine slnout
!
!***********************************************************************
!***********************************************************************
!
end module movcol_mod
