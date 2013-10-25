module my_problem_mod

  use movcol_mod

  type, extends(problem_movcol) :: my_problem
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
  subroutine defpde (eqn, index, t, x, u, ux, ut, uxt, fg)
    class(my_problem) :: eqn
    integer :: index
    real(8) :: t, x
    real(8), dimension(eqn%npde) ::  u, ux, ut, uxt, fg

    real(8) :: eps
    eps = 1.e-4
    if (index.lt.0) fg (1) = ut (1)
    if (index.gt.0) fg (1) = eps * ux (1) - u (1) **2 / 2.d0
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
    if (index.lt.0) res (1) = u (1) - exact (0.d0, t)
    if (index.gt.0) res (1) = u (1) - exact (1.d0, t)
  end subroutine defbcp


  !
  !-----
  ! define the initial values
  !-----
  !
  subroutine defivs (eqn, x, u, ux)
    class(my_problem) :: eqn
    real(8) :: x, u(eqn%npde), ux(eqn%npde)
    u (1) = exact (x, 0.d0)
    ux (1) = dexact (x, 0.d0)
  end subroutine defivs


  !
  !-----
  ! define the initial mesh. this mesh is used only when job = 1
  !-----
  !
  subroutine defmsh (eqn, xmesh)
    class(my_problem) :: eqn
    real(8) :: xmesh(eqn%npts)
    integer :: i, npts
    npts = eqn%npts
    do i = 1, npts
       !     define a uniform mesh
       xmesh (i) = (i - 1) / (npts - 1)
    end do
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
    if (index.lt.0) res = x - 0.d0
    if (index.gt.0) res = x - 1.d0
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
    fmntr = dsqrt (1.d0 + ux (1) **2)
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
  subroutine defout (eqn, npde, npts, t, xmesh, xmesht, u, ux, ut, tstep,&
       touta, ntouta, istop, index, nts)
    class(my_problem) :: eqn
    integer :: npts, npde, ntouta, istop, index, nts
    real(8), dimension(npts, npde) :: u,  ux, ut
    real(8), dimension(npts)       :: xmesh, xmesht
    real(8), dimension(ntouta)     :: touta
    real(8)                        :: t, tstep

    real(8) :: x, uval, uexac
    integer :: n, i

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
    do n = 1, ntouta
       if (dabs (t - touta (n) ) .le.1.d-15) then
          write (eqn%nprnt, 99997) t
          do i = 1, npts
             x = xmesh (i)
             !         call deffnc(x, 1, npde, npts, xmesh, u, ux, uval, uxval)
             uval = u (i, 1)
             uexac = exact (x, t)
             write (eqn%nprnt, 99999) x, uval, uexac, dabs (uval - uexac)
          end do
       endif
    end do
    return
99997 format(/'# t = ', f10.4/)
99999 format(6(1x,e12.5))
  end subroutine defout

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                     c
  !                      define the exact solution                      c
  !                                                                     c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !-----
  ! define the exact solution for the current problem
  !-----
  !
  real(8) function exact (x, t)
    real(8) :: x,t
    real(8) :: eps, a, b, c, temp, ea, eb, ec
    eps = 1.d-4
    a = - 0.05d0 / eps * (x - 0.5d0 + 4.95d0 * t)
    b = - 0.25d0 / eps * (x - 0.5d0 + 0.75d0 * t)
    c = - 0.50d0 / eps * (x - 0.375d0)
    temp = dmax1 (a, b, c)
    ea = 0.d0
    eb = 0.d0
    ec = 0.d0
    if (a - temp.ge. - 35.d0) ea = dexp (a - temp)
    if (b - temp.ge. - 35.d0) eb = dexp (b - temp)
    if (c - temp.ge. - 35.d0) ec = dexp (c - temp)
    exact = (0.1d0 * ea + 0.5d0 * eb + ec) / (ea + eb + ec)
  end function exact
  !
  !-----
  ! define the derivative of the exact solution for the current
  ! problem
  !-----
  !
  real(8) function dexact (x, t)
    real(8) :: x, t
    real(8) :: eps, a, b, c, temp, r1, r2, r3, r1x, r2x, r3x
    eps = 1.d-4
    a = - 0.05d0 / eps * (x - 0.5d0 + 4.95d0 * t)
    b = - 0.25d0 / eps * (x - 0.5d0 + 0.75d0 * t)
    c = - 0.50d0 / eps * (x - 0.375d0)
    temp = dmax1 (a, b, c)
    r1 = 0.d0
    r2 = 0.d0
    r3 = 0.d0
    if (a - temp.ge. - 35.d0) r1 = dexp (a - temp)
    if (b - temp.ge. - 35.d0) r2 = dexp (b - temp)
    if (c - temp.ge. - 35.d0) r3 = dexp (c - temp)
    r1x = - 0.05d0 * r1 / eps
    r2x = - 0.25d0 * r2 / eps
    r3x = - 0.50d0 * r3 / eps
    dexact = ( ( - 0.4d0 * r2 - 0.9d0 * r3) * r1x + (0.4d0 * r1 -     &
         0.5d0 * r3) * r2x + (0.9d0 * r1 + 0.5d0 * r2) * r3x) / (r1 + r2 + &
         r3) **2
  end function dexact
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
  npde  =  1
  npts  = 41

  ! initialize solver
  call my_eqn%init(npde,npts)

  ! set the parameters
  my_eqn%left_end  = 0.0
  my_eqn%right_end = 1.0

  ! error tolerances
  atol = 1.d-3
  rtol = 1.d-5

  ! output times
  allocate(touta(4))
  touta = [0.00, 0.25, 0.55, 1.00]

  ! solve the equations
  call my_eqn%solve(atol, rtol, touta, iflag, "soln.dat")

end program ex1
