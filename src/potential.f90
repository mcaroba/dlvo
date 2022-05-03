! Copyright by Miguel A. Caro (2022).

module potentials

  implicit none

  contains

!**********************************************************************************************
! This subroutine returns the distance between ri and rj under
! certain boundary conditions
  subroutine get_distance(posi, posj, L, PBC, dist, d)

    implicit none

    real*8, intent(in) :: posi(1:3), posj(1:3), L(1:3)
    logical, intent(in) :: PBC(1:3)
    real*8, intent(out) :: d
    real*8, intent(out) :: dist(1:3)
    real*8 :: d2
    integer :: i

    d2 = 0.d0
    do i = 1, 3
      if( PBC(i) )then
        dist(i) = modulo(posj(i) - posi(i), L(i))
        if( dist(i) > L(i)/2.d0 )then
          dist(i) = dist(i) - L(i)
        end if
      else
        dist(i) = posj(i) - posi(i)
      end if
      d2 = d2 + dist(i)**2
    end do
    d = dsqrt(d2)

  end subroutine get_distance
!**********************************************************************************************




!**********************************************************************************************
! This returns potential energy and force for a 1/r type potential. G is a
! constant prefactor; it could be the gravity constant or e^2/(4 pi eps0), etc.
! The interaction is smoothed out to zero at a distance equivalent to half the
! shortest dimension of the simulation box. The attenuation length for the cutoff
! function is the whole length to avoid artifacts in the force close to the cutoff

  subroutine pairwise_electrostatic_potential(posi, posj, Zi, Zj, G, L, PBC, &
                                              Epot, fi)

    implicit none

    real*8, intent(in) :: posi(1:3), posj(1:3), Zi, Zj, G, L(1:3)
    logical, intent(in) :: PBC(1:3)
    real*8, intent(out) :: Epot, fi(1:3)
    real*8 :: d, dist(1:3), fcut, dfcut, pi, Lmin, start

    pi = dacos(-1.d0)

    call get_distance(posi, posj, L, PBC, dist, d)

!   We enable the cutoff function only if under periodic boundary conditions
    if( any(PBC) )then
      Lmin = minval(L) / 2.d0
      start = 0.d0*Lmin
      if( d >= Lmin )then
        fcut = 0.d0
        dfcut = 0.d0
      else if( d <= start )then
        fcut = 1.d0
        dfcut = 0.d0
      else
        fcut = 0.5d0 * (1.d0 + dcos((d-start)*pi/(Lmin-start)))
        dfcut = -0.5d0 * sin((d-start)*pi/(Lmin-start)) * pi / (Lmin-start)
      end if
    else
      fcut = 1.d0
      dfcut = 0.d0
    end if

    Epot = G * Zi*Zj / d * fcut

!   The force on i is calculated assuming the convention that dist(1) = xj - xi
    fi(1:3) = G * Zi*Zj * (-1.d0/d**2*fcut + 1.d0/d*dfcut) * dist(1:3)/d

  end subroutine
!**********************************************************************************************




!**********************************************************************************************
! This subroutine returns the interaction energy between two particles according to the
! Lennard-Jones potential. Rcut is chosen as half the minimum dimension in L
  subroutine lj_potential(posi, posj, sigmai, sigmaj, epsiloni, epsilonj, &
                          L, PBC, Epot, fi)

    implicit none

    real*8, intent(in) :: posi(1:3), posj(1:3), sigmai, sigmaj, &
                          epsiloni, epsilonj, L(1:3)
    logical, intent(in) :: PBC(1:3)
    real*8, intent(out) :: Epot, fi(1:3)
    real*8 :: d, epsilon, sigma, pi, dist(1:3), E0, Rcut

    pi = dacos(-1.d0)

    Rcut = minval(L)/2.d0

    sigma = 0.5d0 * (sigmai + sigmaj)
    epsilon = dsqrt(epsiloni*epsilonj)

    E0 = 4.d0*epsilon * ( (sigma/Rcut)**12 - (sigma/Rcut)**6 )
    
    call get_distance(posi, posj, L, PBC, dist, d)

    if( d < Rcut )then
      Epot = 4.d0*epsilon * ( (sigma/d)**12 - (sigma/d)**6 ) - E0
!     The force on i is calculated assuming the convention that dist(1) = xj - xi
      fi(1:3) = 4.d0*epsilon/d**2 * (-12.d0*(sigma/d)**12 + &
                6.d0*(sigma/d)**6) * dist(1:3)
    else
!     There is no discontinuity of the potential at Rcut, but there
!     is a discontinuity of the force
      Epot = 0.d0
      fi(1:3) = 0.d0
    end if

  end subroutine lj_potential
!**********************************************************************************************



end module
