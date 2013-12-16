module my_problem_mod

  use movcol_mod

  type, extends(problem_movcol) :: my_problem
     real    :: dim = 7
     real    :: amplitude = 0.0
   contains
     procedure :: defivs
     procedure :: defout
     procedure :: defmsh
     procedure :: defpde
     procedure :: defbcp
     procedure :: defbcm
     procedure :: defmnt
  end type my_problem

  ! this variable is used to communicate with the defout function
  ! where it controls weather the result will be printed or not

contains
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                     c
  !                 define the physical problem                         c
  !            subroutines: defpde, defbcp, defivs, defmsh              c
  !                                                                     c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !-----
  ! define the physical pdes, where
  ! fg := f(t, x, u, ux, ut, uxt)      if index < 0
  ! fg := g(t, x, u, ux, ut, uxt)      if index > 0
  !-----
  !
  ! f=(g)_x
  subroutine defpde (eqn, index, t, x, u, ux, ut, uxt, fg)
    class(my_problem) :: eqn
    integer :: index
    real(8) :: t, x
    real(8), dimension(eqn%npde) ::  u, ux, ut, uxt, fg

    associate(d => eqn%dim)
      ! if( x < 1.e-3 ) then
      !    print *, x, d-3, x**(d-3), (x**2*ut+(d-1)/2.0*sin(2*u)), ux
      ! end if
      ! if( x < epsilon(1.0) ) then
      !    fg = 0.0
      ! else
      !    if (index.lt.0) fg = x**(d-3)*(x**2*ut+(d-1)/2.0*sin(2*u))
      !    if (index.gt.0) fg = x**(d-1)*ux
      ! end if

      if (index.lt.0) fg = x**(d-3)*(x**2*ut+(d-1)/2.0*sin(2*u))
      if (index.gt.0) fg = x**(d-1)*ux
      ! if (index.lt.0) fg = (d-1)/x*ux-(d-1)/2.0/x**2*sin(2*u)
      ! if (index.gt.0) fg = -ux
    end associate

  end subroutine defpde
  !
  !-----
  ! define the physical boundary conditions, where
  ! res(k) := b^l_k(t, x, xt, u, ux, uxx, ut, uxt)  if index < 0
  ! res(k) := b^r_k(t, x, xt, u, ux, uxx, ut, uxt)  if index > 0
  ! for k = 1, ..., npde
  !-----
  !
  subroutine defbcp (eqn, index, t, x, xt, u, ux, uxx, ut, uxt,res)
    class(my_problem) :: eqn
    integer :: index
    real(8) :: t, x, xt
    real(8), dimension(eqn%npde) ::  u, ux, uxx, ut, uxt, res
    if (index.lt.0) res = u
    if (index.gt.0) res = u - acos(-1.0)
  end subroutine defbcp


  !
  !-----
  ! define the initial values
  !-----
  !
  subroutine defivs (eqn, x, u, ux)
    class(my_problem) :: eqn
    real(8) :: x, u(eqn%npde), ux(eqn%npde)

    associate(a => eqn%amplitude)
      u  = x + a*sin(x)
      ux = 1.0 + a*cos(x)
    end associate

  end subroutine defivs


  !
  !-----
  ! define the initial mesh. this mesh is used only when job = 1
  !-----
  !
  subroutine defmsh (eqn, xmesh)
    class(my_problem) :: eqn
    real(8) :: xmesh(eqn%npts)
    integer :: i

    associate( npts => eqn%npts )
      xmesh = [(acos(-1.0)*(i-1.0)/(npts-1.0), i=1,npts)]
    end associate

  end subroutine defmsh
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                     c
  !          define bcs and monitor function for moving mesh            c
  !               subroutines: defbcm, defmnt                           c
  !                                                                     c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !-----
  ! define the boundary conditions for the coordinate transformation,
  ! where
  ! res := b^l_m(t, x, xt, u, ux, uxx, ut, uxt)  if index < 0
  ! res := b^r_m(t, x, xt, u, ux, uxx, ut, uxt)  if index > 0
  !-----
  !
  subroutine defbcm (eqn, index, t, x, xt, u, ux, uxx, ut, uxt, res)
    class(my_problem) :: eqn
    integer :: index
    real(8) :: t, x, xt, res
    real(8), dimension(eqn%npde) ::  u, ux, uxx, ut, uxt
    if (index.lt.0) res = x
    if (index.gt.0) res = x - acos(-1.0)
  end subroutine defbcm

  !
  !-----
  ! define the monitor function
  !-----
  !
  subroutine defmnt (eqn, t, x, u, ux, uxx, ut, uxt, fmntr)
    class(my_problem) :: eqn
    real(8) :: t, x, fmntr
    real(8), dimension(eqn%npde) :: u, ux, uxx, ut, uxt
    !     define the arclength monitor function
    fmntr = dsqrt (1. + ux (1) **2)
    fmntr = sqrt(10.0+ux(1)**2)
    fmntr = abs(ux(1))
    fmntr = 1.0
    ! fmntr = 10.0+abs(ux(1))+sqrt(abs(ux(1)))
  end subroutine defmnt

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                     c
  !                        output solutions                             c
  !                       subroutine: defout                            c
  !                                                                     c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !-----
  ! output solutions. the time integration may be stopped if istop is set
  ! less than zero.
  ! the code is solving the physical pdes   when index > 0
  ! the code is generating a mesh           when index < 0
  !-----
  !
  subroutine defout (eqn, t, xmesh, xmesht, u, ux, ut, tstep,&
       touta, ntouta, istop, index, nts)
    class(my_problem) :: eqn
    integer :: ntouta, istop, index, nts
    real(8), dimension(eqn%npts, eqn%npde) :: u,  ux, ut
    real(8), dimension(eqn%npts)       :: xmesh, xmesht
    real(8), dimension(ntouta)     :: touta
    real(8)                        :: t, tstep

    integer :: i, npts
    npts = eqn%npts

    !
    !-----
    ! do not output solutions when generating a mesh with respect to
    ! the given initial values of the physical solution
    !-----
    !
    if (index.lt.0) return
    !
    !-----
    ! output solutions and errors at t = touta(1), ..., touta(ntouta)
    !-----
    !
    if( mod(eqn%nsteps,50) == 0 )  then

       write(eqn%nprnt, *)
       write(eqn%nprnt, *)
       write(eqn%nprnt, '("# t = ", g0)') t
       do i = 1, npts
          write(eqn%nprnt, *) xmesh(i), u(i,1), ux(i,1)
       end do

    end if

    if( u(2,1) > 1.0 ) istop = -1

  end subroutine defout

end module my_problem_mod


program ex1
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!                           example 1                                 c
!                                                                     c
! this is an example of using movcol to solve burgers' equation with  c
! an exact solution :                                                 c
!                                                                     c
!     u_t = [ 10^(-4) u_x - u^2/2 ]_x, 0 < x < 1, 0 < t <= 1          c
!     u(0, t) = exact(0, t), u(1, t) = exact(1,t)                     c
!     u(x, 0) = exact(x, 0)                                           c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! version: nov. 22, 1995
!
!-----
! set values for some basic parameters and declare arrays
!-----
!
  use movcol_mod
  use my_problem_mod

  real(8), allocatable :: touta(:)
  real(8) :: atol, rtol
  integer :: iflag, npde, npts

  type(my_problem) :: my_eqn

  ! size of mesh and equation number
  my_eqn%npde  =  1
  my_eqn%npts  = 201

  ! set the parameters
  my_eqn%left_end  = 0.0
  my_eqn%right_end = acos(-1.0)

  ! error tolerances
  my_eqn%atol = 1.d-5
  my_eqn%rtol = 1.d-6

  ! tau
  my_eqn%tau  = 1.e-6
  my_eqn%mmpde= 3
  my_eqn%job  = 2
  ! my_eqn%ip   = 0

  ! output times
  allocate(my_eqn%touta(2))
  my_eqn%touta = [0.00, 10.00]

  ! initialize solver
  call my_eqn%init()

  ! solve the equations
  call my_eqn%solve("soln.dat")

end program ex1
