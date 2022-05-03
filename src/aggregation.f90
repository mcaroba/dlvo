! Copyright by Miguel A. Caro (2022)

program aggregation

  use potentials
  use integrators

  implicit none

!**************************************************************************************
! Variable declarations
! Time and temperature stuff
  real*8 :: dt, tau, T
  integer :: Nt

! Particle arrays
  real*8, allocatable :: pos(:,:), vel(:,:), mass(:), force(:,:), Z(:), sigma(:), epsilon(:), &
                         xi_prev(:,:)
  real*8 :: Ep, Ek, xi_array(1:2, 1:3)

! Constants
  real*8, parameter :: kB = 8.6173303d-5, eps0 = 5.526349406d-3, pi = 3.141592653589793d0

! Colloidal particle
  real*8 :: epsilon_col, sigma_col, Z_col, m_col
  integer :: n_col

! Solvent particles
  real*8 :: epsilon_sol, sigma_sol, Z_sol, Z1_sol, Z2_sol, m_sol
  integer :: n_sol, n1_sol, n2_sol

! Box
  real*8 :: Lx, L(1:3)

! Miscellaneous
  real*8 :: econst, fi(1:3), this_Ep, new_pos(1:3), new_vel(1:3), d, dist(1:3)
  integer :: i, j, Np, k, n_write, i2, step, iostatus
  logical :: crazy_forces, too_close, PBC(1:3)
  character*32 :: keyword, junk
!**************************************************************************************

  econst = 1.d0 / 4.d0 / pi / eps0
! Define the time parameters
!  Nt = 10000 ! No. of time steps
  dt = 0.5d0 ! Time step

! Thermostat
!  T = 1000.d0 ! Target temperature
  tau = 100.d0 ! Time constant

! Write every these many steps
  n_write = 10

! We hardcode sol parameters
  sigma_sol = 3.d0 ! Angstrom
  epsilon_sol = 0.01d0 ! eV
  Z_sol = 0.1d0 ! electron charges
!  n_sol = 100 ! number of particles with each charge
  m_sol = 10.d0 ! amu

! We hardcode some col parameters
  sigma_col = 10.d0 ! Angstrom
  epsilon_col = 1.d0 ! eV
  Z_col = 0.5d0 ! electron charges
  m_col = 100.d0 ! amu

! Read in number of colloidal particles
  open(unit=10, file="input", status="old")
  iostatus = 0
  do while( iostatus == 0 )
    read(10,*,iostat=iostatus) keyword
    if( keyword == "Nt" )then
      backspace(10)
      read(10,*,iostat=iostatus) junk, junk, Nt
    else if( keyword == "T" )then
      backspace(10)
      read(10,*,iostat=iostatus) junk, junk, T
    else if( keyword == "n_sol" )then
      backspace(10)
      read(10,*,iostat=iostatus) junk, junk, n_sol
    else if( keyword == "n_col" )then
      backspace(10)
      read(10,*,iostat=iostatus) junk, junk, n_col
    end if      
  end do
  close(10)

! We need to make sure that charge neutrality is preserved
  n1_sol = n_sol
  if( Z_sol /= 0.d0 )then
    n2_sol = int((Z_sol*dfloat(n_sol) + Z_col*dfloat(n_col)) / Z_sol)
  else
    n2_sol = n1_sol
  end if
  Z_sol = (Z_sol*dfloat(n1_sol+n2_sol)) / dfloat(n1_sol+n2_sol)
  Z1_sol = Z_sol
  Z2_sol = -Z_sol

  Np = n_col + n1_sol + n2_sol


! Estimate box size from vdW radii (properly, we should use a barostat)
  Lx = (5.d0 * (sigma_sol**3*(n1_sol+n2_sol) + sigma_col**3*n_col))**(1.d0/3.d0)
  L = [Lx, Lx, Lx]


! Allocate particle arrays
  allocate( pos(1:3, 1:Np) )
  allocate( xi_prev(1:3, 1:Np) )
  allocate( vel(1:3, 1:Np) )
  allocate( force(1:3, 1:Np) )
  allocate( mass(1:Np) )
  allocate( Z(1:Np) )
  allocate( sigma(1:Np) )
  allocate( epsilon(1:Np) )

! Initialize the arrays
  vel = 0.d0
  mass(1:n_col) = m_col
  mass(n_col+1:Np) = m_sol
  Z(1:n_col) = Z_col
  Z(n_col+1:n_col+n1_sol) = Z1_sol
  Z(n_col+n1_sol+1:Np) = Z2_sol
  sigma(1:n_col) = sigma_col
  sigma(n_col+1:Np) = sigma_sol
  epsilon(1:n_col) = epsilon_col
  epsilon(n_col+1:Np) = epsilon_sol
  pos = 0.d0
  do i = 1, Np
    too_close = .true.
    do while( too_close )
      too_close = .false.
      call random_number( pos(1:3, i) )
      pos(1:3, i) = pos(1:3, i) * L(1:3)
      do j = 1, i-1
        call get_distance(pos(1:3, i), pos(1:3, j), L, [.true., .true., .true.], new_pos, d)
        if( d < (sigma(i)+sigma(j) / 2.d0) )then
          too_close = .true.
          exit
        end if
      end do
    end do
  end do

  open(unit=10, file="out.xyz", status="unknown")
  open(unit=20, file="thermo.log", status="unknown")

  write(*,*) "Progress:"
  write(*,*) "|----------------------------------------|"
  write(*,'(A)', advance="no") "  "
  k = 1
  do step = 0, Nt
!   Write to file every n_write steps
    if( mod(step, n_write) == 0 )then
      write(10,*) Np
      write(10,*) 'Lattice="', L(1), 0.d0, 0.d0, 0.d0, L(2), 0.d0, 0.d0, 0.d0, L(3), '"'
      do i = 1, Np
        call get_distance(0.5d0*L, pos(1:3, i), L, [.true., .true., .true.], new_pos, d)
        if( i <= n_col )then
          write(10,*) "Col", new_pos, force(1:3, i)
        else if( i <= n_col + n1_sol )then
          write(10,*) "Sol1", new_pos, force(1:3, i)
        else
          write(10,*) "Sol2", new_pos, force(1:3, i)
        end if
      end do
    end if

!   Print progress bar
    if( step/(dfloat(Nt)/40.d0) >= k )then
      write(*,'(A)', advance="no") "*"
      k = k + 1
    end if

!   Compute energy and force
    force = 0.d0
    Ep = 0.d0
    Ek = 0.d0
    crazy_forces = .false.
    do i = 1, Np
      do j = i+1, Np
        call pairwise_electrostatic_potential(pos(1:3, i), pos(1:3, j), Z(i), Z(j), econst, &
                                              L, [.true., .true., .true.], this_Ep, fi)
        Ep = Ep + this_Ep
        force(1:3, i) = force(1:3, i) + fi
        force(1:3, j) = force(1:3, j) - fi
        call lj_potential(pos(1:3, i), pos(1:3, j), sigma(i), sigma(j), epsilon(i), epsilon(j), &
                          L, [.true., .true., .true.], this_Ep, fi)
        Ep = Ep + this_Ep
        force(1:3, i) = force(1:3, i) + fi
        force(1:3, j) = force(1:3, j) - fi
      end do
    end do


!   Integrate the equations of motion
    do i = 1, Np
      if( step == 0 )then
        xi_prev(1:3, i) = pos(1:3, i)
        call lazy_man(pos(1:3, i), vel(1:3, i), force(1:3, i), mass(i), dt, new_pos, new_vel)
      else
        xi_array(1, 1:3) = xi_prev(1:3, i)
        xi_array(2, 1:3) = pos(1:3, i)
        call verlet(xi_array, force(1:3, i), mass(i), dt, new_pos)
      end if
      xi_prev(1:3, i) = pos(1:3, i)
      pos(1:3, i) = new_pos(1:3)
!     Interpolate the velocity to estimate the temperature
      vel(1:3, i) = (pos(1:3, i) - xi_prev(1:3, i)) / dt
      Ek = Ek + 0.5d0 * mass(i) * dot_product(vel(1:3, i), vel(1:3, i))
    end do
!   Remove CM velocity (should be close to zero since we initialized properly)
    call remove_cm_vel(vel, mass)
!   Rescale the velocities to control the temperature
    call berendsen_thermostat(vel, T, 2.d0/3.d0/dfloat(Np-1)/kB*Ek, tau, dt)
!   Readjust the positions according to new velocities
    pos = xi_prev + vel*dt

!   Log the thermodynamic properties
    write(20,*) step, Ep, Ek
  end do
  write(*,*)

  close(10)
  close(20)

end program
