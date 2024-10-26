program simpson_integration
  implicit none
  ! Define constants
  real(8) :: gamma, omega, a, l
  integer :: num_steps_r, num_steps_z
  real(8) :: dr, dz, r, z, inner_integral, outer_integral
  integer :: i

  ! Assign values to constants
  gamma = 1.0d0
  omega = 50.0d-6
  a = 0.002d0
  l = 0.02d0

  ! Set number of steps for numerical integration (must be even for Simpson's rule)
  num_steps_r = 1000
  num_steps_z = 1000

  ! Calculate step sizes
  dr = a / num_steps_r
  dz = l / num_steps_z

  ! Initialize integrals
  outer_integral = 0.0d0

  ! Outer integral using Simpson's rule
  do i = 0, num_steps_z, 2
    z = i * dz
    if (i == 0 .or. i == num_steps_z) then
      outer_integral = outer_integral + simpson_inner_integral(z)  ! f(a) or f(b)
    elseif (mod(i, 2) == 0) then
      outer_integral = outer_integral + 2.0d0 * simpson_inner_integral(z)  ! Even terms
    else
      outer_integral = outer_integral + 4.0d0 * simpson_inner_integral(z)  ! Odd terms
    end if
  end do
  outer_integral = outer_integral * dz / 3.0d0

  ! Print the final result
  print *, "The value of the integral G is: ", outer_integral

contains

  ! Function to perform the inner integral using Simpson's rule
  function simpson_inner_integral(z) result(inner_result)
    real(8), intent(in) :: z
    real(8) :: inner_result, r
    integer :: j
    inner_result = 0.0d0

    do j = 0, num_steps_r, 2
      r = j * dr
      if (j == 0 .or. j == num_steps_r) then
        inner_result = inner_result + exp(-2.0d0 * r**2 / omega**2) * r  ! f(a) or f(b)
      elseif (mod(j, 2) == 0) then
        inner_result = inner_result + 2.0d0 * exp(-2.0d0 * r**2 / omega**2) * r  ! Even terms
      else
        inner_result = inner_result + 4.0d0 * exp(-2.0d0 * r**2 / omega**2) * r  ! Odd terms
      end if
    end do
    inner_result = inner_result * dr / 3.0d0
  end function simpson_inner_integral

end program simpson_integration
