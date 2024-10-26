program numerical_integration
  implicit none
  ! Define constants
  real(8) :: gama, omegaf, radius, length
  integer :: num_steps_r, num_steps_z
  real(8) :: deltar, deltaz, r, z, r_integral, z_integral

  nr = 1000
  nz = 1000
  
  gama       = 1.
  radius     = 0.002
  length     = 0.02
  omegaf     = 50.e-6
  z_integral = 0.

  deltar = radius / nr
  deltaz = length / nz

  do k = 0, nz
	 z = k * deltaz 
	 
    r_integral = 0.
	do j = 0, nr
	   r = j * deltar
       r_integral = r_integral + exp(-2 * r**2 / omegaf**2) * r * deltar
    end do

    z_integral= z_integral+ exp(-gama * z) * r_integral * deltaz
  end do
  G = z_integral
  
  wirte(*,*) "The value of the Normalization Constant (G): ", G

end program numerical_integration
