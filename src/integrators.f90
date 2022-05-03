! Copyright by Miguel A. Caro (2022)

module integrators

  implicit none

  contains

! Lazy Man's approach
  subroutine lazy_man(x_in, v_in, F, m, dt, x_out, v_out)

    implicit none

    real*8, intent(in) :: x_in(1:3), v_in(1:3)
    real*8, intent(out) :: x_out(1:3), v_out(1:3)
    real*8, intent(in) :: F(1:3), m, dt

    x_out(1:3) = x_in(1:3) + v_in(1:3)*dt + F(1:3)/m/2.d0*dt**2
    v_out(1:3) = v_in(1:3) + F(1:3)/m*dt

  end subroutine


! Regular Verlet
  subroutine verlet(x_in, F, m, dt, x_out)

    implicit none

    real*8, intent(in) :: x_in(1:2,1:3)
    real*8, intent(out) :: x_out(1:3)
    real*8, intent(in) :: F(1:3), m, dt

    x_out(1:3) = 2.d0*x_in(2,1:3) - x_in(1,1:3) + F(1:3)/m*dt**2

  end subroutine




  subroutine berendsen_thermostat(vel, T0, T, tau, dt)

    implicit none

    real*8, intent(inout) :: vel(:,:)
    real*8, intent(in) :: T0, T, tau, dt

    vel = vel * dsqrt(1.d0 + dt/tau * (T0/T - 1.d0))

  end subroutine




  subroutine remove_cm_vel(vel, M)

!   I should adapt this code to mixed boundary conditions, where
!   the CM velocity can be removed per Cartesian dimension independently

    implicit none

    real*8, intent(inout) :: vel(:,:)
    real*8, intent(in) :: M(:)
    real*8 :: cm_pos(1:3), cm_vel(1:3), total_mass
    integer :: Np, i

    Np = size(vel, 2)

    cm_vel = 0.d0
    total_mass = 0.d0
    do i = 1, Np
      cm_vel(1:3) = cm_vel(1:3) + M(i)*vel(1:3,i)
      total_mass = total_mass + M(i)
    end do
    cm_vel = cm_vel / total_mass
    do i = 1, Np
      vel(1:3,i) = vel(1:3,i) - cm_vel(1:3)
    end do

  end subroutine

end module
