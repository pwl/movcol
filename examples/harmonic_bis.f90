module my_problem_mod

  use movcol_mod

  type, extends(problem_movcol) :: my_problem
     ! integer    :: dim = 27
     ! integer    :: k   = 3
     integer    :: dim = 8
     integer    :: k   = 1
     real       :: amplitude
     integer    :: bis = 0
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
    if (index.gt.0) res = u

  end subroutine defbcp


  !
  !-----
  ! define the initial values
  !-----
  !
  subroutine defivs (eqn, x, u, ux)
    class(my_problem) :: eqn
    real :: x, u(eqn%npde), ux(eqn%npde)

    associate( a => eqn%amplitude )
      u  = a*sin(x)
      ux = a*cos(x)
    end associate

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
    fmntr = sqrt (1. + ux (1) **2)
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
  subroutine defout (eqn, t, xmesh, xmesht, u, ux, ut, uxt, uxx, tstep,&
       touta, ntouta, istop, index, nts)
    class(my_problem) :: eqn
    integer :: ntouta, istop, index, nts
    real, dimension(eqn%npts, eqn%npde) :: u,  ux, ut, uxt, uxx
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
    if( mod(eqn%nsteps,5) == 0 )  then

       write(eqn%nprnt(1), *)
       write(eqn%nprnt(1), *)
       write(eqn%nprnt(1), '("# t = ", ES50.32)') t
       do i = 1, npts
          write(eqn%nprnt(1), '(3(ES50.32," "))') xmesh(i), u(i,1), ux(i,1), uxx(i,1)
       end do

    end if


    if( mod(eqn%nsteps,1) == 0 )  then
       write(eqn%nprnt(2), '(2(ES50.32," "))') t, ux(1,1)
    end if

    if( any(u(:,1)>1.75) ) then
       istop = -1
       eqn%bis = +1
    end if

    if( all(u(:,1)<acos(-1.0)/2.0) ) then
       istop = -1
       eqn%bis = -1
    end if

  end subroutine defout


  function bis(amplitude)
    real, intent(in) :: amplitude
    real :: bis

    type(my_problem), allocatable :: my_eqn
    character(len=200) :: dirname

    allocate(my_eqn)

    my_eqn%amplitude = amplitude

    ! size of mesh and equation number
    my_eqn%npde  =  1
    my_eqn%npts  = 400

    ! set the parameters
    my_eqn%left_end  = 0.0
    my_eqn%right_end = acos(-1.0)

    ! error tolerances
    my_eqn%atol = 1.d-6
    my_eqn%rtol = 1.d-5

    ! tau
    my_eqn%tau  = 1.e-10
    my_eqn%mmpde= 4
    my_eqn%job  = 2

    ! output times
    allocate(my_eqn%touta(2))
    my_eqn%touta = [0.00, 10.00]

    ! output files
    allocate(character(len=200) :: my_eqn%filenames(2))

    write(dirname, '("data_bis/d",i0.3,"/k",i0.3,"/n",i0.4,"/amp",f0.32,"/")') my_eqn%dim, my_eqn%k, my_eqn%npts, my_eqn%amplitude
    ! make sure the directory exists
    call system('mkdir -p '//trim(dirname))
    ! append the particular file names
    my_eqn%filenames = trim(dirname) // ["soln.dat","at0.dat"]

    ! initialize solver
    call my_eqn%init()

    ! solve the equations
    call my_eqn%solve()

    if( my_eqn%idid <= -2 ) then
       bis = 0.0
    else
       bis = my_eqn%bis
    end if

    deallocate(my_eqn)

  end function bis

end module my_problem_mod

program harmonic_bis
  use my_problem_mod

  real :: amp0 = 0.1
  real :: amp1 = 3.0
  real :: ampn, valn, val0, val1
  real :: dist
  integer :: i

  dist = amp1-amp0
  i = 0
  val0 = bis(amp0)
  val1 = bis(amp1)


  if( val0*val1 > 0 ) then
     print *, "No zero fund in the initial interval"
  else
     do while( dist > 1.e-30 )
        print '(i0.4,", ",3(ES50.32,", "))', i, amp1, amp0, dist
        ampn = (amp0+amp1)/2.0
        valn = bis(ampn)

        if( valn == 0 ) then
           print *, "Error occured when solving pde."
           exit
        end if

        if( val0*valn > 0 ) then
           amp0 = ampn
           val0 = valn
        else
           amp1 = ampn
        end if

        dist = amp1-amp0
        i = i+1
     end do
  end if

  print '(ES50.32)', (amp1+amp0)/2

end program harmonic_bis
