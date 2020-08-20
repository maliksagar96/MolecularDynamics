program main

    use class_module

    implicit none
    integer :: i,j,k,ii,jj,kk,monteCarloSteps = 1000
    real*8 :: kineticEnergy, pe, TE
    !Initializing everyhing
    !call posInit()
    x(1) = 1.0
    y(1) = 1.0
    z(1) = 1.0
    x(2) = 14.8
    y(2) = 1.0
    z(2) = 1.0

    open(14, file = 'velocities.dat')

    do i = 1, particles
      velocity_x(i) = 0
      velocity_y(i) = 0
      velocity_z(i) = 0
    end do
    !call initVelocities()
    call initForce()

    do i = 1, 1000
      call verletAlgorithm()
      print *,i,x(1)
    end do

    close(14)

end program main
