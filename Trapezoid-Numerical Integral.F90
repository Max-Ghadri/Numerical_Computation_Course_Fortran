program numerical_integration
  implicit none
  ! Define constants
  real(8) :: gamma, omega, a, l
  integer :: num_steps_r, num_steps_z
  real(8) :: dr, dz, r, z, inner_integral, outer_integral
  real(8) :: f_r, f_r_prev, f_z, f_z_prev

  ! Assign values to constants
  gamma = 1.0d0
  omega = 50.0d-6
  a = 0.002d0
  l = 0.02d0

  ! Set number of steps for numerical integration
  num_steps_r = 1000
  num_steps_z = 1000

  ! Calculate step sizes
  dr = a / num_steps_r
  dz = l / num_steps_z

  ! Initialize integrals
  outer_integral = 0.0d0

  ! Outer integral over z using the trapezoidal rule
  f_z_prev = exp(-gamma * 0.0d0)  ! Value at z = 0
  do z = dz, l, dz
    f_z = exp(-gamma * z)

    ! Inner integral over r using the trapezoidal rule
    inner_integral = 0.0d0
    f_r_prev = exp(-2.0d0 * (0.0d0)**2 / omega**2) * 0.0d0  ! Value at r = 0
    do r = dr, a, dr
      f_r = exp(-2.0d0 * r**2 / omega**2) * r
      inner_integral = inner_integral + 0.5d0 * (f_r + f_r_prev) * dr
      f_r_prev = f_r  ! Update previous value
    end do

    ! Add contribution from inner integral to the outer integral
    outer_integral = outer_integral + 0.5d0 * (f_z + f_z_prev) * inner_integral * dz
    f_z_prev = f_z  ! Update previous value
  end do

  ! Print the final result
  print *, "The value of the integral G is: ", outer_integral

end program numerical_integration
