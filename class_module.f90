module class_module
  implicit none

    real*8 :: length = 15.0                                           !Box is cubical box of Boxsize = boxsize*boxsize*boxsize
    integer, parameter :: particles = 2
    real, parameter :: dt = 0.001
    real*8, parameter :: sigma = 1.0d0, eps = 1.0d0, mass = 1.0d0, temp = 1.0d0
    real*8, parameter :: rc = 2.5*sigma
    real*8 :: velocity_x(particles), velocity_y(particles), velocity_z(particles)
    real*8 :: x(particles), y(particles), z(particles)
    integer, parameter :: seed_x = 2550, seed_y = 456, seed_z = 7982  !To check random number calling for differernt seeds
    real :: dt2by2 = (dt**2)/(mass*2.0), dtby2 = dt/(2.0*mass)
    real*8 :: potential = 0
    real*8 :: force_x(particles), force_y(particles), force_z(particles)
    real*8:: f_rc = eps*((12*(sigma**12)/(rc**13)) -(6*(sigma**6)/rc**7))
    real*8 :: u_rc = eps*((sigma**12)/(rc**12) - (sigma**6)/(rc**6))
    real :: force_flag



contains
  !Defining all the funtions first then defining all the subroutines

    function velInit(arrayDim, seed) result(v)
      implicit none
      integer, intent(in) :: arrayDim, seed
      real*8, dimension(arrayDim) :: v
      integer :: i, j = 0
      real :: vc = 10   !velocity constant
      call srand(seed)
      do i = 1, size(v)
        v(i) = vc*(rand() - 0.05d0)
      end do
    end function velInit

    !Minimum image convention
    function minDist(x1, x2) result(dist)
      implicit none
      real*8, intent(in) :: x1, x2
      real*8 :: dist
      
        if((x1 - x2) > length /2) then
          dist = x1 -x2 - length

        else if ((x1 - x2)< (-length /2)) then
          dist = x1 - x2 + length

        else
          dist = x1 - x2
        end if
    end function minDist

    ! Always provide square root of r^2
    function calcForce(xi, xj, r) result(force)
      implicit none
      real*8, intent(in) :: xi, xj,r
      real*8::force, sig6, xDiff
      xDiff = minDist(xj, xi)
      sig6 = (sigma**6)/(r**6)
      force = (6*sig6*eps*(2*sig6-1)*(xDiff)/(r**2) -f_rc*(xDiff)/r)
    end function calcForce

    function boundaryBounce(x1) result(xBox)
      implicit none
      real*8 :: x1
      real*8 :: xBox
      xBox = x1
      if(x1 > length) then
        xBox = xBox - length * float(int(xBox/length))
      else if (x1 < 0) then
        xBox = xBox + length *  float(int(-xBox/length) + 1)
      else
        xBox = x1
      end if

    end function


    function calcKE(vx, vy, vz) result(KE)
      implicit none
      real*8, intent(in):: vx, vy, vz
      real*8 ::KE
      KE = 0.5 * mass * (vx**2 + vy**2 + vz**2)
    end function calcKE

    function calcPE(r) result(PE)
      implicit none
      real*8, intent(in) :: r
      real*8 :: PE, sig6
      sig6 = (sigma**6)/(r**6)
      PE = eps*sig6*(sig6 - 1) - u_rc + r*f_rc
    end function calcPE

    ! Initializing position with a equally spaced particles with a constant separation between 2 adjacent molecules in any direction
    subroutine posInit()

      implicit none
      integer :: i, j, k, pcounter, dummyCounter, counter
      real :: sepFactor = 1.2
      integer :: maxCount
      sepFactor = sigma*sepFactor
      maxCount = floor(length/sepFactor)
      x(1) = sepFactor
      y(1) = sepFactor
      z(1) = sepFactor
      counter = 0
      do i = 1, particles-1
        counter = counter + 1
        if(counter == maxCount) then
          counter = 0
          x(i+1) = sepFactor
        else
          x(i+1) = x(i) + sepFactor
        end if
        !print *,mod(i, (maxCount)), x(i)
      end do

      counter = 0
      do i = 1, particles-1
        counter = counter + 1
        !print *,y(i),i
        if(counter == maxCount) then
          counter = 0
          y(i+1) = y(i) + sepFactor
        else
          y(i+1) = y(i)
        end if

        if(y(i+1) > sepFactor*maxCount) then
          y(i+1) = sepFactor
        end if
      end do

      counter = 0
      do i = 1, particles-1
        counter = counter + 1
        !print *,y(i),i
        if(counter == (maxCount**2)) then
          counter = 0
          z(i+1) = z(i) + sepFactor
        else
          z(i+1) = z(i)
        end if
      end do
    end subroutine posInit

    subroutine updatePos()
      implicit none
      integer :: i,j
      real*8::r

        do i = 1, particles
            !print*, boundaryBounce(x(i) + velocity_x(i) * dt + force_x(i) * dt2by2)
            x(i) = boundaryBounce(x(i) + velocity_x(i) * dt + force_x(i) * dt2by2)
            y(i) = boundaryBounce(y(i) + velocity_y(i) * dt + force_y(i) * dt2by2)
            z(i) = boundaryBounce(z(i) + velocity_z(i) * dt + force_z(i) * dt2by2)

        end do
        !! Impose boundary cross coundition
    end subroutine updatePos

    subroutine initVelocities()
      implicit none
      integer :: i
      real*8 :: avgVx, avgVy, avgVz, total_x
      !initializing without normalizing and Making Velocity of COM zero
      velocity_x = velInit(size(velocity_x), seed_x)
      velocity_y = velInit(size(velocity_y), seed_y)
      velocity_z = velInit(size(velocity_z), seed_z)
      !averaging velocities
      do i = 1, size(velocity_x)
        avgVx = velocity_x(i)+avgVx
        avgVy = velocity_y(i)+avgVy
        avgVz = velocity_z(i)+avgVz
      end do

      avgVx = avgVx/float(particles)
      avgVy = avgVy/float(particles)
      avgVz = avgVz/float(particles)

      !Making Velocity of COM of zero
      do i = 1, size(velocity_x)
        velocity_x(i) = avgVx - velocity_x(i)
        velocity_y(i) = avgVy - velocity_y(i)
        velocity_z(i) = avgVz - velocity_z(i)
      end do
    end subroutine initVelocities

    subroutine updateVelocity()
      implicit none
      integer :: i,j
      real*8::r
      do i = 1, particles
        velocity_x(i) = velocity_x(i) + force_x(i) * dtby2         !Heart of the subroutine. Mass is incorporated in dt2by2
        velocity_y(i) = velocity_y(i) + force_y(i) * dtby2
        velocity_z(i) = velocity_z(i) + force_z(i) * dtby2
      end do
    end subroutine updateVelocity

    subroutine initForce()
      implicit none
      !Calculating force on ith element
      integer :: i,j
      real*8 :: r
      do i = 1, size(force_x)
        force_x(i) = 0
        force_y(i) = 0
        force_z(i) = 0
      end do
      call updateForce()
    end subroutine initForce

    subroutine updateForce()
      implicit none
      integer :: i,j
      real*8 :: r, fx, fy, fz
      integer :: counter = 0
      do i = 1, (size(force_x)-1)
        do j = (i+1), size(force_x)
          r = sqrt(minDist(x(i),x(j))**2 + minDist(y(i), y(j))**2 + minDist(z(i), z(j))**2)
          !print *, y(i), y(j)
          if(r<rc) then
            !force on ith particle due to all other particles

            counter = counter + 1
            fx = calcForce(x(i), x(j), r)
            fy = calcForce(y(i), y(j), r)
            fz = calcForce(z(i), z(j), r)

            force_x(i) = force_x(i) + fx
            force_y(i) = force_y(i) + fy
            force_z(i) = force_z(i) + fz

            !force on jth particle due to ith
            force_x(j) = force_x(j) - fx
            force_y(j) = force_y(j) - fy
            force_z(j) = force_z(j) - fz
            !print *, counter,force_x(i), r
          end if
        end do
      end do
    end subroutine

    subroutine verletAlgorithm()
      call updatePos()
      call updateVelocity()
      call updateForce()
    end subroutine verletAlgorithm


end module
