!> @file   movcol.f90
!! @author Pawel Biernat <pawel.biernat@gmail.com>
!! @date   Fri Nov  8 16:23:57 2013
!!
!! @brief
!!
!! @todo add y_current i ydot_current?
!!
module movcol_mod

  use ddassl_mod

  ! algorithm parameters?
  real(8), parameter ::&
       &tau1 = 1.e-4,&
       &s1 = 0.211324865347452,&
       &s2 = 0.788675134652548

  private

  type :: temp_storage
     ! of size npde
     real(8), allocatable, dimension(:)   :: u, ux, uxx, ut, uxt, rw
     ! of size npts
     real(8), allocatable, dimension(:)   :: xmesh, xmesht
     ! of size npts x npde
     real(8), allocatable, dimension(:,:) :: uu, uux, uut
     ! work array for ddassl
     real(8), allocatable, dimension(:)   :: ddassl_rwork(:)
     integer, allocatable, dimension(:)   :: ddassl_iwork(:)
     ! work arrays for movcol
     real(8), allocatable, dimension(:) :: movcol_rwork(:)
     ! integer, allocatable, dimension(:) :: movcol_iwork(:)
     ! old phypde
     logical :: physpde
  end type temp_storage


  type, abstract, public, extends(problem_ddassl) :: problem_movcol
     ! work arrays
     type(temp_storage), private :: tmp

     ! pointers to the physical quantities
     real(8), pointer     :: x(:), u(:,:), ux(:,:)
     real(8), pointer     :: xt(:), ut(:,:), uxt(:,:)
     real(8), pointer     :: resx(:), resu(:,:), resux(:,:)
     real(8), allocatable :: y(:,:), ydot(:,:)
     ! flat counterparts to y and ydot
     real(8), pointer     :: y_flat(:), ydot_flat(:)

     ! print times
     real(8), allocatable :: touta(:)

     ! parameters
     integer :: job       = 2     !par(1)
     integer :: output    = -1    !par(2)
     integer :: mmpde     = 4     !par(3)
     real(8) :: tau       = 1e-4!par(4)
     real(8) :: gamma     = 1.0 !par(5)
     integer :: ip        = 2.0 !par(6)
     real(8) :: stpmax    = 1.0/epsilon(0.0) !par(7)
     real(8) :: left_end  = 0.0 !par(8)
     real(8) :: right_end = 1.0 !par(9)

     ! relative and absolute error tolerances
     real(8) :: rtol = -1.0, atol = -1.0

     ! parameters not originally in movcol
     real(8) :: relaxation_time = 1.0

     ! system size and number of mesh points
     integer :: npts, npde

     ! number of steps taken with ddassl
     integer :: nsteps = 0

     ! error flags
     integer :: iflag = 0
     integer :: idid

     ! output file id
     integer :: nprnt = 111
     ! information is written to this unit
     integer :: inform = 6
   contains
     private

     ! main procedure
     procedure, public :: solve => movcol_solve
     ! initialization function
     procedure, public :: init => movcol_init

     ! other private  procedures
     procedure :: movcl1
     procedure :: solution_out
     procedure :: respde
     procedure :: respd1
     procedure :: evaluate_monitor
     procedure :: resmeq
     ! procedures defined for parent problem_ddassl
     procedure, public :: res => resode
     procedure, public :: jac => jacode
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

     subroutine defout_i (eqn, t, xmesh, xmesht, u, ux, ut, tstep,&
          touta, ntouta, istop, index, nts)
       import problem_movcol
       class(problem_movcol)          :: eqn
       integer                        :: ntouta, istop, nts, index
       real(8), dimension(eqn%npts, eqn%npde) :: u,  ux, ut
       real(8), dimension(eqn%npts)       :: xmesh, xmesht
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
    subroutine movcol_init(eqn)
      class(problem_movcol), target :: eqn

      integer :: npde, npts

      ! length of work arrays
      integer :: liw

      ! number of dependent values per mesh point
      integer :: m

      ! total number of equations for ddassl
      integer :: neq

      ! size of ddassl workarray and movcol workarrays
      integer :: lrw0, lrw1

      !
      !...check basic parameters job, inform, mmpde, tau, gamma, ip and
      !...stpmax. if their values are not within the permitted ranges,
      !...default values are set and iflag = 1.
      !
      eqn%iflag = 0

      if(all(eqn%job /= [1,2,3])) then
         eqn%iflag = 1
         print *, "job has to be [1,2,3]"
         return
      end if
      !
      if(all(eqn%mmpde /= [1,2,3,4])) then
         eqn%iflag = 1
         print *, "mmpde has to be [1,2,3,4]"
         return
      end if
      !
      if(eqn%ip < 0) then
         eqn%iflag = 1
         print *, "ip < 0"
         return
      endif
      !
      if(eqn%stpmax <= 0) then
         eqn%iflag = 1
         print *, "stpmax <= 0"
         return
      endif
      !
      !...check basic parameters npde, npts, atol, rtol, ntouta,
      !...xl and xr
      !
      if (eqn%npde <= 0) then
         eqn%iflag = - 1
         print *, "npde <= 0"
         return
      end if

      if (eqn%npts.le.0) then
         eqn%iflag = - 1
         print *, "npts <= 0"
         return
      end if

      if (abs (eqn%atol) + abs (eqn%rtol) .le. epsilon(0.0)) then
         eqn%iflag = - 1
         print *, "|atol| + |rtol| <= epsilon"
         return
      end if

      if (eqn%job.eq.2.or.eqn%job.eq.3) then
         if (eqn%left_end .gt. eqn%right_end) then
            eqn%iflag = - 1
            print *, "job = 2 or 3 and left_end > right_end"
            return
         end if
      endif
      !
      !...if iflag = - 1, at least one of the parameters npde, npts,
      !...atol, rtol and ntouta has been incorrectly input, so the
      !...computation is halted
      !
      if (eqn%iflag.eq. - 1) then
         write (eqn%inform, 99910) eqn%iflag
         write (eqn%inform, 99990) eqn%npde, eqn%npts, eqn%job, eqn%output, eqn%mmpde, eqn%tau,     &
         eqn%gamma, eqn%ip, eqn%atol, eqn%rtol, eqn%stpmax
         return
      endif

      npde  = eqn%npde
      npts  = eqn%npts

      m = 2*npde+1

      liw = 20 + npts*m


      ! allocate the ddassl work array
      neq = (2 * npde+1) * npts
      if (eqn%mmpde == 2 .or. eqn%mmpde == 3) then
         lrw0 = 40 + 12 * neq + 3 * (eqn%ip + 2) * (2 * npde+1) * neq
      else
         lrw0 = 40 + 12 * neq + 3 * (0 + 2) * (2 * npde+1) * neq
      endif
      allocate(eqn%tmp%ddassl_rwork(lrw0))
      allocate(eqn%tmp%ddassl_iwork(liw))
      eqn%tmp%ddassl_rwork = 0.0
      eqn%tmp%ddassl_iwork = 0

      ! allocate the work array with a size decreased by double the
      ! size of eqn%y and eqn%ydot
      lrw1 = 20 + 6 * npde+ (3 * npde+2) * npts
      allocate(eqn%tmp%movcol_rwork(lrw1))
      ! allocate(eqn%tmp%movcol_iwork(liw))

      ! initialize work arrays with zeroes
      eqn%tmp%movcol_rwork = 0.0
      ! eqn%tmp%movcol_iwork = 0

      ! if the touta table was not allocate, allocate it as an empty table
      if(.not. allocated(eqn%touta)) then
         allocate(eqn%touta(1))
         eqn%touta=1.0
      end if

      allocate(eqn%y   (m,npts))
      allocate(eqn%ydot(m,npts))

      ! flat counterparts to y and ydot
      eqn%y_flat(1:m*npts) => eqn%y
      eqn%ydot_flat(1:m*npts) => eqn%ydot

      eqn%x  => eqn%y(m,            :)
      eqn%u  => eqn%y(1:npde,       :)
      eqn%ux => eqn%y(npde+1:2*npde,:)

      eqn%xt => eqn%ydot(m,            :)
      eqn%ut => eqn%ydot(1:npde,       :)
      eqn%uxt=> eqn%ydot(npde+1:2*npde,:)

      ! temporary storage, should replace the work arrays
      allocate(eqn%tmp%u(npde))         !rwk1(m11+1)
      allocate(eqn%tmp%ux(npde))        !rwk1(m12+1)
      allocate(eqn%tmp%uxx(npde))       !rwk1(m13+1)
      allocate(eqn%tmp%ut(npde))        !rwk1(m14+1)
      allocate(eqn%tmp%uxt(npde))       !rwk1(m15+1)
      allocate(eqn%tmp%rw(npde))        !rwk1(m16+1)
      allocate(eqn%tmp%xmesh(npts))     !rwk1(m21+1)
      allocate(eqn%tmp%xmesht(npts))    !rwk1(m22+1)
      allocate(eqn%tmp%uu(npts,npde))   !rwk1(m31+1)
      allocate(eqn%tmp%uux(npts,npde))  !rwk1(m32+1)
      allocate(eqn%tmp%uut(npts,npde))  !rwk1(m33+1)

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

    end subroutine movcol_init

    subroutine movcol_solve(eqn, filename)
      class(problem_movcol) :: eqn
      character(len=*) :: filename

      real :: tcpu, tcpu1, tcpu2, timaray(2)

      if( eqn%iflag /= 0 ) then
         print *, "ERROR: cannot start solver, check iflag"
         return
      end if

      open(newunit = eqn%nprnt, file = filename)

      call etime (timaray, tcpu1)

      !
      !...call movcl1
      !
      call eqn%movcl1()

      call etime(timaray, tcpu2)
      tcpu = tcpu2 - tcpu1

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
      ! subroutine movcol (eqn, atol, rtol, touta, iflag)
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
! evaluate_monitor:i  computes the monitor function
!
! smoothen_monitor:i  smooths the monitor function (only used for mmpde = 2
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
! solution_out:i  driver for the subroutine defout for solution output
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
      ! real(8) :: rtol, atol
      ! class(problem_movcol), target :: eqn
      ! real(8) :: touta(:)
      ! ! this subroutine doesn't do anything now, do not call it
      ! print *, "Do not call subroutine movcol"
      ! return
      ! end subroutine movcol
!
!***********************************************************************
!***********************************************************************
!
      subroutine movcl1 (eqn)
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

      class(problem_movcol), target :: eqn


      ! arrays allocated for ddassl
      integer :: info (15), iwk1 (0)
      real(8) :: rtol1(1), atol1(1)

      ! this variable controls the execution flow for the whole movcl1
      ! subroutine.  When phypd = TRUE we are solving the physical
      ! system, otherwise we are relaxing the mesh

      integer :: npde, npts, liw, ip, mmpde, job, lrw, ntouta, index

      integer :: m, neq, i, itout, istop, nts, itoutm, inform
      real(8) :: xl, xr, gamma, stpmax, zrmch, t, tau, temp, tstep
      real(8) :: tout, rtol, atol, rwk1(0)

      real(8), pointer :: rwk(:)
      real(8), pointer :: touta(:), y(:), ydot(:)
      integer, pointer :: iwk(:)

      ! bind the local variables to the variables stored in eqn
      job = eqn%job
      xl  = eqn%left_end
      xr  = eqn%right_end
      mmpde  = eqn%mmpde
      gamma  = eqn%gamma
      ip     = eqn%ip
      stpmax = eqn%stpmax
      tau    = eqn%tau
      inform = eqn%inform
      npde   = eqn%npde
      npts   = eqn%npts
      rtol   = eqn%rtol
      atol   = eqn%atol

      rwk    => eqn%tmp%ddassl_rwork
      iwk    => eqn%tmp%ddassl_iwork
      ! iwk    => eqn%tmp%movcol_iwork
      y      => eqn%y_flat
      ydot   => eqn%ydot_flat
      touta  => eqn%touta

      lrw    = size(rwk)
      liw    = size(iwk)
      ntouta = size(touta)

      ! so far so good, so set the error flag to OK
      eqn%iflag = 0
!
!...define some basic parameters
!
      m = 2 * npde+1
      neq = m * npts
      zrmch = epsilon (1.0)
!     phypde = 'true' for solving the physical pdes
!     phypde = 'false' for generating a mesh
      if (job.eq.1) eqn%tmp%physpde = .true.
      if (job.eq.2) eqn%tmp%physpde = .false.
      if (job.eq.3) eqn%tmp%physpde = .false.
!
!...compute the initial mesh
!
!     for solving the physical pde
      if (eqn%tmp%physpde) then
         call eqn%defmsh(eqn%x)
      else
!     the mesh generation case: starting from a uniform mesh
         do i = 1, npts
            eqn%x (i) = xl + (xr - xl) * (i - 1) / (npts - 1)
         end do
      endif

!
!...print the headers
!
      write (inform, 99920)
!  this is the entry point for one of the gotos
   20 if (eqn%tmp%physpde) then
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
      if (eqn%tmp%physpde) then
!        the case of solving the physical pdes
         do i = 1, npts
            call eqn%defivs(eqn%x(i), eqn%u(:,i), eqn%ux(:,i))
         end do
      else
!        the mesh generation case
         eqn%u  = 0.0
         eqn%ux = 0.0
      endif
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
      rwk (3) = rwk (2) * 0.5
      info (9) = 1
!        define the order of the bdf method (1st order)
      iwk (3) = 1
      info (11) = 1
!        define the initial values for ydot
      do i = 1, neq
         ydot (i) = 0.0
      end do
!
!...take two time integration steps using ddassl
!
      if (eqn%tmp%physpde) then
         t = touta (1)
         tout = touta (1) + 1.0
      else
         t = 0.0
         tout = 1.0
      endif

      atol1 = 1.0
      rtol1 = 1.0
      do itout = 1, 2
         call ddassl (eqn, neq, t, y, ydot, tout, info, rtol1, atol1,&
         eqn%idid, rwk, lrw, iwk, liw, rwk1, iwk1)
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
      if (.not.eqn%tmp%physpde) then
         do i = 1, npts
            eqn%x(i) = xl + (xr - xl) * (i - 1) / (npts - 1)
         end do
         eqn%u  = 0.0
         eqn%ux = 0.0
      endif
!
!     for the physical solution
      if (eqn%tmp%physpde) then
!        the case of solving the physical pdes
         do i = 1, npts
            call eqn%defivs(eqn%x(i), eqn%u(:,i), eqn%ux(:,i))
         end do
      endif
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
      if (.not.eqn%tmp%physpde) temp = 1.0
      rwk (3) = min (stpmax, rtol, atol, temp, 1.0)
      info (11) = 1
!
!...print the initial solution
!
      if (eqn%tmp%physpde) then
         t = touta (1)
         index = + 1
      else
         t = 0.0
         index = - 1
      endif
      tstep = 0.0
      istop = 0
      nts = 0
      call eqn%solution_out (t, tstep, istop, index, nts)
!     if istop < 0, stop the computation
      if (istop.lt.0) then
         eqn%iflag = - 3
         write (inform, 99950) eqn%iflag
         write (inform, 99990) npde, npts, job, inform, mmpde, tau,     &
         gamma, ip, atol, rtol, stpmax, ntouta
         return
      endif
!
!...perform the time integration with ddassl
!
      if (eqn%tmp%physpde) then
         itoutm = ntouta
      else
         itoutm = 2
      endif
      do itout = 2, itoutm
         if (eqn%tmp%physpde) then
            tout = touta (itout)
         else
            tout = eqn%relaxation_time
         endif

110      continue

         rtol1 = eqn%rtol
         atol1 = eqn%atol
         if( eqn%tmp%physpde ) eqn%nsteps = eqn%nsteps+1

         call ddassl (eqn, neq, t, y, ydot, tout, info, rtol1, atol1,&
              eqn%idid, rwk, lrw, iwk, liw, rwk1, iwk1)

!
!...print the solution at t
!
         if (eqn%tmp%physpde) then
            index = + 1
         else
            index = - 1
         endif
         tstep = rwk (7)
         istop = 0
         nts = iwk (11)
         call eqn%solution_out (t, tstep, istop, index, nts)
!        call eqn%solution_out (npde, npts, t, y, ydot, touta, ntouta,
!         tstep, istop, index, nts)

!     if istop < 0, stop the computation
         if (istop.lt.0) then
            eqn%iflag = - 3
            write (inform, 99950) eqn%iflag
            write (inform, 99990) npde, npts, job, inform, mmpde, tau,  &
            gamma, ip, atol, rtol, stpmax, ntouta
            return
         endif
!
!...check the performance of the time integrator
!
         if (eqn%idid.eq.1) goto 110
!        a step was successfully taken in the intermediate-output
!        mode. ddassl has not yet reached tout. continue the time
!        integration
         if (eqn%idid.eq. - 1) then
!        a large amount of work has been expended. reset the
!        limitation and continue the time integration
            info (1) = 1
            goto 110
         endif
         if (eqn%idid.ge. - 1.and.eqn%idid.le.3) then
!        the time integration is successful to tout. output the
!        performance information
            write (inform, 99970) t, eqn%idid, iwk (14), iwk (15), iwk (11),&
            iwk (13), iwk (12), iwk (8), rwk (7), rwk (3)
            info (1) = 1
         else
!        the time integration fails and the computation is stopped
            eqn%iflag = - 4
            write (inform, 99980) eqn%iflag, eqn%idid
            write (inform, 99990) npde, npts, job, inform, mmpde, tau,  &
            gamma, ip, atol, rtol, stpmax, ntouta
            return
         endif
      end do
!
!...if job = 2, go back to solve the physical pdes
!
      if (job.eq.2.and. (.not.eqn%tmp%physpde) ) then
         eqn%tmp%physpde = .true.
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
      real(8), dimension(*), target :: y, ydot, res, rwk

      ! local variables
      integer :: npde, npts, mmpde, ip, m
      real(8) :: tau, gamma
      real(8), pointer :: y2d(:,:), ydot2d(:,:), res2d(:,:)
!
!...define some basic parameters
!

      npde = eqn%npde
      npts = eqn%npts
      mmpde = eqn%mmpde
      ip = eqn%ip

      gamma = eqn%gamma
      tau = eqn%tau

      m = 2 * npde+1

      y2d   (1:m,1:npts) => y(1:m*npts)
      ydot2d(1:m,1:npts) => ydot(1:m*npts)
      res2d (1:m,1:npts) => res(1:m*npts)

      ! update the pointers before calling subroutines defined by user
      eqn%x    => y2d   (m,            :)
      eqn%u    => y2d   (1:npde,       :)
      eqn%ux   => y2d   (npde+1:2*npde,:)

      eqn%xt   => ydot2d(m,            :)
      eqn%ut   => ydot2d(1:npde,       :)
      eqn%uxt  => ydot2d(npde+1:2*npde,:)

      eqn%resx => res2d (m,            :)
      eqn%resu => res2d (1:npde,       :)
      eqn%resux=> res2d (npde+1:2*npde,:)

!
!...compute residuals of the physical pdes and their bcs
!
      if (eqn%tmp%physpde) then
!        the case of solving the physical pdes
         call eqn%respde (npde, npts, t, y, ydot, res, eqn%tmp%u,       &
         eqn%tmp%ux, eqn%tmp%uxx, eqn%tmp%ut, eqn%tmp%uxt, eqn%tmp%rw)
      else
!        the mesh generation case
         call eqn%respd1 (npde, npts, t, y, ydot, res, eqn%tmp%u,       &
         eqn%tmp%ux, eqn%tmp%uxx, eqn%tmp%ut, eqn%tmp%uxt )
      endif
!
!...compute and smooth the monitor function
!
!     rwk(m21+i) := fmntr(i)
      if (mmpde.ge.2) call eqn%evaluate_monitor (npde, npts, t, y, ydot,&
           & eqn%tmp%xmesh, eqn%tmp%u, eqn%tmp%ux, eqn%tmp%uxx, eqn%tmp%ut,     &
           & eqn%tmp%uxt )
      if (mmpde.eq.2.or.mmpde.eq.3) call smoothen_monitor (npts, ip, eqn%tmp%xmesh&
      , eqn%tmp%xmesht )
!
!...compute residuals for the discrete mesh equations and their bcs.
!
      call eqn%resmeq ( npde, npts, mmpde, tau, gamma, y, ydot, res, eqn%tmp%xmesh,&
           & eqn%tmp%u, eqn%tmp%ux, eqn%tmp%uxx, eqn%tmp%ut, eqn%tmp%uxt)
      if (.not.eqn%tmp%physpde) then
!        the mesh generation case (for which boundaries are fixed)
         eqn%resx(   1) = eqn%xt(   1)
         eqn%resx(npts) = eqn%xt(npts)
      endif
!

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
      class(problem_movcol) :: eqn
      real(8) :: t, y, ydot, res, u, ux, uxx, ut, uxt, rw
      integer :: npde, npts

      dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      dimension res ( (2 * npde+1) * npts), u (npde), ux (npde),        &
      uxx (npde)
      dimension ut (npde), uxt (npde), rw (npde)

      ! local variables
      real(8) :: w(3,2), h, x, xt
      integer :: i, j, k, m

!
      m = 2 * npde+1
!
!...define the weights for the cell-average collocation approximation
!...to the right-hand-side term
!
      w (1, 1) =  -1.0-2.0 * (s2 - s1)
      w (2, 1) =       4.0 * (s2 - s1)
      w (3, 1) =   1.0-2.0 * (s2 - s1)
      w (1, 2) =  -1.0+2.0 * (s2 - s1)
      w (2, 2) =     - 4.0 * (s2 - s1)
      w (3, 2) =   1.0+2.0 * (s2 - s1)
!
!...compute the residuals for the discretization of the physical pdes
!...at the interior (gauss) points
!
      do 50 i = 1, npts - 1
         h = eqn%x(i+1)-eqn%x(i)
!
!...compute the left-hand-side term f(...) using collocation
!
         do 20 j = 1, 2
            if (j.eq.1) x = eqn%x(i) + s1 * h
            if (j.eq.2) x = eqn%x(i) + s2 * h
            call drvtvs (npde, npts, i, x, y, ydot, u, ux, uxx, ut, uxt)
            call eqn%defpde (- 1, t, x, u, ux, ut, uxt, rw)
            do 10 k = 1, npde
               if (j.eq.1) res(m*(i-1) + 2*(k-1) + 2) = rw (k)
               if (j.eq.2) res(m*(i)   + 2*(k-1) + 1) = rw (k)
   10       end do
   20    end do
!
!...compute the right-hand-side term [ g(...) ]_x using cell-average
!...collocation
!
         do 40 j = 1, 3
            if (j.eq.1) x = eqn%x(i)
            if (j.eq.2) x = eqn%x(i) + h/2.0
            if (j.eq.3) x = eqn%x(i+1)
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
      class(problem_movcol) :: eqn
      real(8) :: t, y, ydot, res, u, ux, uxx, ut, uxt
      integer :: npde, npts
      dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      dimension res ( (2 * npde+1) * npts), u (npde), ux (npde),        &
      uxx (npde)
      dimension ut (npde), uxt (npde)

      ! local variables
      integer :: m, i, j, k
      real(8) :: h, ht, x, xt, temp, temp0, temp1, temp2, temp3

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
               0.0) + y (m * (i - 1) + 2 * (k - 1) + 2) * fncshp (2, 0,&
               0.0) * h + y (m * (i) + 2 * (k - 1) + 1) * fncshp (3, 0,&
               0.0) + y (m * (i) + 2 * (k - 1) + 2) * fncshp (4, 0,    &
               0.0) * h
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
               1.0) + y (m * (i - 1) + 2 * (k - 1) + 2) * fncshp (2, 0, &
               1.0) * h + y (m * (i) + 2 * (k - 1) + 1) * fncshp (3, 0, &
               1.0) + y (m * (i) + 2 * (k - 1) + 2) * fncshp (4, 0, 1.0)&
               * h
               if (j.eq.1) then
                  if (xt.ge.0.0) temp = (temp2 - temp1) / (s2 - s1)    &
                  / h
                  if (xt.le.0.0) temp = (temp1 - temp0) / s1 / h
                  res (m * (i - 1) + 2 * (k - 1) + 2) = (ut (k) - temp *&
                  xt) - u (k)
               endif
               if (j.eq.2) then
                  if (xt.ge.0.0) temp = (temp3 - temp2) / (1.0-s2)     &
                  / h
                  if (xt.le.0.0) temp = (temp2 - temp1) / (s2 - s1)    &
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
      class(problem_movcol) :: eqn
      real(8) :: tau, y, ydot, res, fmntr, u, ux, uxx, ut, uxt, gamma
      integer :: npde, npts, mmpde
      dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      dimension res ( (2 * npde+1) * npts), fmntr (npts), u (npde)
      dimension ux (npde), uxx (npde), ut (npde), uxt (npde)

      ! local variables
      integer :: m, i
      real(8) :: x, xt, t, temp, temp1, temp2, ym1, yp1, ym3, yp3, xlambd

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
      xlambd = gamma * (gamma + 1.0)
      if (mmpde.eq.4) then
         do 40 i = 2, npts - 1
            temp1 = fmntr (i - 1)
            temp2 = fmntr (i)
            ym1 = - (ydot (m * i) - ydot (m * (i - 1) ) ) / (y (m * i)  &
            - y (m * (i - 1) ) ) **2 * tau + 1.0 / (y (m * i) - y (m *  &
            (i - 1) ) )
            yp1 = - (ydot (m * (i + 1) ) - ydot (m * i) ) / (y (m *     &
            (i + 1) ) - y (m * i) ) **2 * tau + 1.0 / (y (m * (i + 1) ) &
            - y (m * i) )
            if (i - 2.ge.1) then
               ym3 = - (ydot (m * (i - 1) ) - ydot (m * (i - 2) ) )     &
               / (y (m * (i - 1) ) - y (m * (i - 2) ) ) **2 * tau + 1.0 &
               / (y (m * (i - 1) ) - y (m * (i - 2) ) )
            else
               ym3 = ym1
            endif
            if (i + 2.le.npts) then
               yp3 = - (ydot (m * (i + 2) ) - ydot (m * (i + 1) ) )     &
               / (y (m * (i + 2) ) - y (m * (i + 1) ) ) **2 * tau + 1.0 &
               / (y (m * (i + 2) ) - y (m * (i + 1) ) )
            else
               yp3 = yp1
            endif
            res (m * i) = (yp1 - xlambd * (yp3 - 2.0 * yp1 + ym1) )     &
            / temp2 - (ym1 - xlambd * (yp1 - 2.0 * ym1 + ym3) ) / temp1
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
      integer :: npde, npts, i, j, k
      real(8) :: y, ydot, u, ux, uxx, ut, uxt, x
      dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      dimension u (npde), ux (npde), uxx (npde), ut (npde), uxt (npde)

      ! local variables
      integer :: m
      real(8) :: h, ht, temp(4,0:2), s
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
      real(8) :: s
      integer :: j, k
!
      if (j.eq.1) then
         if (k.eq.0) fncshp = (1.0+2.0 * s) * (1.0-s) * (1.0-s)
         if (k.eq.1) fncshp = - 6.0 * s * (1.0-s)
         if (k.eq.2) fncshp = - 6.0 * (1.0-2.0 * s)
      endif
      if (j.eq.2) then
         if (k.eq.0) fncshp = s * (1.0-s) * (1.0-s)
         if (k.eq.1) fncshp = (1.0-s) * (1.0-3.0 * s)
         if (k.eq.2) fncshp = 6.0 * s - 4.0
      endif
      if (j.eq.3) then
         if (k.eq.0) fncshp = (3.0-2.0 * s) * s * s
         if (k.eq.1) fncshp = 6.0 * s * (1.0-s)
         if (k.eq.2) fncshp = 6.0 * (1.0-2.0 * s)
      endif
      if (j.eq.4) then
         if (k.eq.0) fncshp = (s - 1.0) * s * s
         if (k.eq.1) fncshp = s * (3.0 * s - 2.0)
         if (k.eq.2) fncshp = 6.0 * s - 2.0
      endif
!
      return
      end function fncshp
!
!***********************************************************************
!***********************************************************************
!
      subroutine evaluate_monitor (eqn, npde, npts, t, y, ydot, fmntr, u, ux, uxx, ut, &
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
      class(problem_movcol) :: eqn
      integer :: npde, npts
      real(8) ::t, y, ydot, fmntr, u, ux, uxx, ut, uxt
      dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      dimension u (npde), ux (npde), uxx (npde), ut (npde), uxt (npde)
      dimension fmntr (npts)

      ! local variables
      integer :: m, i, j
      real(8) :: h, x, temp
!
      m = 2 * npde+1
      fmntr = 0.0
!
      do i = 1, npts - 1
         h = eqn%x(i+1)-eqn%x(i)
         do j = 1, 2
            if (j.eq.1) x = eqn%x(i) + s1 * h
            if (j.eq.2) x = eqn%x(i) + s2 * h
            call drvtvs (npde, npts, i, x, y, ydot, u, ux, uxx, ut, uxt)
            call eqn%defmnt (t, x, u, ux, uxx, ut, uxt, temp)
            fmntr (i) = fmntr (i) + 0.5 * temp
         end do
      end do
!
!...approximate the monitor function more accurately in the end
!...subintervals (i.e., for i = 1 and npts-1)
!
      do i = 1, npts - 1, npts - 2
         fmntr (i) = 0.5 * fmntr (i)
         do j = i, i + 1
            x = y (m * j)
            call drvtvs (npde, npts, i, x, y, ydot, u, ux, uxx, ut, uxt)
            call eqn%defmnt (t, x, u, ux, uxx, ut, uxt, temp)
            fmntr (i) = fmntr (i) + 0.25 * temp
         end do
      end do
!
      return
      end subroutine evaluate_monitor
!
!***********************************************************************
!***********************************************************************
!
      subroutine smoothen_monitor (npts, ip, fmntr, rw)
!
!...this subroutine smooths the monitor function. it is called
!...only when mmpde = 2 or mmpde = 3
!
      real(8) :: fmntr, rw
      integer :: ip
      integer :: npts
      dimension fmntr (npts), rw (npts)

      ! local variables
      integer :: i, j
      real(8) :: gamma = 2.0, temp, temp1

      rw = fmntr
      temp  = 0.0
      fmntr = 0.0
      do i = 1, npts - 1
         do j = max(1, i - ip), min(npts - 1, i + ip)
            temp1 = (gamma / (gamma + 1.0) )**abs(j - i)
            temp = temp + temp1
            fmntr (i) = fmntr (i) + rw (j) * temp1
         end do
         fmntr (i) = fmntr (i) / temp
      end do
!
      return
      end subroutine smoothen_monitor
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
      integer :: npts, npde, k
      real(8) :: xmesh, u, ux, x, uval, uxval
      dimension xmesh (npts), u (npts, npde), ux (npts, npde)

      ! local variables
      integer :: i1, i
      real(8) :: h, s
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
      subroutine solution_out (eqn, t, tstep, istop, index, nts)
      ! subroutine solution_out (eqn, npde, npts, t, y, ydot, touta, ntouta, tstep,  &
      ! istop, index, nts, u, ux, uxx, ut, uxt, xmesh, xmesht, uu, uux,   &
      ! uut)
!
!...this subroutine outputs the solution values x_t, u, u_x,
!...and u_t at point x and time t in the format prescribed by
!...subroutine defout
!
      class(problem_movcol), target :: eqn
      integer :: index, istop, nts
      real(8) :: t, tstep
      ! integer :: npde, npts, ntouta, istop, index, nts
      ! real(8) :: y, ydot, touta, u, ux, uxx, ut, uxt, xmesh, xmesht, uu, uux, uut
      ! real(8) :: t, tstep
      ! dimension y ( (2 * npde+1) * npts), ydot ( (2 * npde+1) * npts)
      ! dimension touta (ntouta), u (npde), ux (npde), uxx (npde)
      ! dimension ut (npde), uxt (npde), xmesh (npts), xmesht (npts)
      ! dimension uu (npts, npde), uux (npts, npde), uut (npts, npde)

      ! local variables
      integer :: npde, npts, m, i, i1
      real(8) :: x
      type(temp_storage), pointer :: tmp
      tmp => eqn%tmp
      npde = eqn%npde
      npts = eqn%npts
!
      m = 2 * npde+1

      do i = 1, npts
!
!...compute u, ux, ut at time t and x = x_i
!
         x = eqn%x(i)
         i1 = min (i, npts - 1)
         call drvtvs (npde, npts, i1, x, eqn%y_flat, eqn%ydot_flat,&
              & tmp%u, tmp%ux, tmp%uxx, tmp%ut, tmp%uxt)
!
!...store x_i, xdot_i, u(k), ux(k), ut(k) in xmesh(i), xmesht(i),
!...uu(i,k), uux(i,k) and uut(i,k), respectively
!
         tmp%xmesh (i) = x
         tmp%xmesht (i) = eqn%ydot (m, i)
         tmp%uu (i,:) = tmp%u
         tmp%uux(i,:) = tmp%ux
         tmp%uut(i,:) = tmp%ut
      end do
!
!...output the solution values and current time stepsize tstep
!
      call eqn%defout (t, tmp%xmesh, tmp%xmesht, tmp%uu, tmp%uux, tmp%uut, tstep,   &
      eqn%touta, size(eqn%touta), istop, index, nts)
!
      return
      end subroutine solution_out
!
!***********************************************************************
!***********************************************************************
!
end module movcol_mod
