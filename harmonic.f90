module my_problem_mod

  use movcol_mod

  type, extends(problem_movcol) :: my_problem
     integer    :: dim = 7
     integer    :: k   = 1
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
    real :: t, x
    real, dimension(eqn%npde) ::  u, ux, ut, uxt, fg

    associate(d => eqn%dim, k => eqn%k)
      if (index.lt.0) fg = ut+k*(d+k-2)/2.0/x**2*sin(2*u)-(d-1)/x*ux
      if (index.gt.0) fg = ux
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
    real :: t, x, xt
    real, dimension(eqn%npde) ::  u, ux, uxx, ut, uxt, res
    if (index.lt.0) res = ut
    if (index.gt.0) res = u-x

  end subroutine defbcp


  !
  !-----
  ! define the initial values
  !-----
  !
  subroutine defivs (eqn, x, u, ux)
    class(my_problem) :: eqn
    real :: x, u(eqn%npde), ux(eqn%npde)

    u  = x
    ux = 1.0

  end subroutine defivs

  !
  !-----
  ! define the initial mesh. this mesh is used only when job = 1
  !-----
  !
  subroutine defmsh (eqn, xmesh)
    class(my_problem) :: eqn
    real :: xmesh(eqn%npts)
    integer :: i

    associate( npts => eqn%npts )
      xmesh = [(eqn%left_end+eqn%right_end*(i-1.0)/(npts-1.0), i=1,npts)]
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
    real :: t, x, xt, res
    real, dimension(eqn%npde) ::  u, ux, uxx, ut, uxt
    if (index.lt.0) res = xt
    if (index.gt.0) res = xt
  end subroutine defbcm

  !
  !-----
  ! define the monitor function
  !-----
  !
  subroutine defmnt (eqn, t, x, u, ux, uxx, ut, uxt, fmntr)
    class(my_problem) :: eqn
    real :: t, x, fmntr
    real, dimension(eqn%npde) :: u, ux, uxx, ut, uxt
    !     define the arclength monitor function
    fmntr = dsqrt (1. + ux (1) **2)
    ! fmntr = sqrt(10.0+ux(1)**2)
    ! fmntr = abs(ux(1))
    ! fmntr = 1.0
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
    real, dimension(eqn%npts, eqn%npde) :: u,  ux, ut
    real, dimension(eqn%npts)       :: xmesh, xmesht
    real, dimension(ntouta)     :: touta
    real                        :: t, tstep

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

       write(eqn%nprnt(1), *)
       write(eqn%nprnt(1), *)
       write(eqn%nprnt(1), '("# t = ", g0)') t
       do i = 1, npts
          write(eqn%nprnt(1), *) xmesh(i), u(i,1), ux(i,1)
       end do

    end if


    if( mod(eqn%nsteps,1) == 0 )  then
       write(eqn%nprnt(2), *) t, ux(1,1)
    end if

    if( ux(2,1) > 1.0e7 ) istop = -1
    ! if( u(2,1) > .05 ) istop = -1

  end subroutine defout

end module my_problem_mod


program ex1
!
!-----
! set values for some basic parameters and declare arrays
!-----
!
  use movcol_mod
  use my_problem_mod

  real, allocatable :: touta(:)
  real :: atol, rtol
  integer :: iflag, npde, npts

  type(my_problem), allocatable :: my_eqn
  character(len=200) :: dirname

  do npts = 400, 400, 100
     allocate(my_eqn)

     ! size of mesh and equation number
     my_eqn%npde  =  1
     my_eqn%npts  = npts

     ! set the parameters
     my_eqn%left_end  = 0.0
     my_eqn%right_end = acos(-1.0)

     ! error tolerances
     my_eqn%atol = 1.d-5
     my_eqn%rtol = 1.d-6

     ! tau
     my_eqn%tau  = 1.e-10
     my_eqn%mmpde= 4
     my_eqn%job  = 1
     ! my_eqn%ip   = 0

     ! output times
     allocate(my_eqn%touta(2))
     my_eqn%touta = [0.00, 10.00]

     ! output files
     allocate(character(len=200) :: my_eqn%filenames(2))

     write(dirname, '("data/d",i0.3,"/k",i0.3,"/n",i0.4,"/")') my_eqn%dim, my_eqn%k, my_eqn%npts
     ! make sure the directory exists
     call system('rm -rf '//trim(dirname))
     call system('mkdir -p '//trim(dirname))
     ! append the particular file names
     my_eqn%filenames = trim(dirname) // ["soln.dat","at0.dat"]

     ! initialize solver
     call my_eqn%init()

     ! solve the equations
     call my_eqn%solve()

     deallocate(my_eqn)

  end do

end program ex1
