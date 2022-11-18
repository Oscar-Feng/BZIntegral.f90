program temp
  use BZInt3D
  implicit none

  ! Testing Quad3DRuleΘ and Quad3DRuleδ. In order to compare between Fortran and
  ! Julia, one may need to add a space between "!" and "$OMP..." in the
  ! corresponding subroutine in bzint3d.f90 to disable parallelism
  integer :: i,j,k,iter
  real(wp), allocatable, dimension(:,:,:) :: Emesh
  real(wp), allocatable, dimension(:,:,:) :: Wmesh
  real(wp),dimension(9) :: arr
  real :: t1,t2
  real(wp), parameter :: PI=4*atan(1.0_wp)

  allocate(Emesh(9,9,9),Wmesh(9,9,9))

  iter=3
  arr=(/-1.0_wp, -0.75_wp, -0.5_wp, -0.25_wp, 0.0_wp, 0.25_wp, 0.5_wp, 0.75_wp, 1.0_wp/)

  do i=1,9
    do j=1,9
      do k=1,9
        Emesh(i,j,k)=(arr(i)*arr(i)+arr(j)*arr(j)+arr(k)*arr(k))*0.5_wp
      end do
    end do
  end do

  call CPU_TIME(t1)
  call Quad3DRuleTheta(Emesh, 0.5_wp, Wmesh, iter)
  call CPU_TIME(t2)
  print*, "calculated volume of radius 1 sphere: ", sum(Wmesh)*8
  print*, "relative error: ", (sum(Wmesh)*8-4*PI/3)/(4*PI/3)*100,"%"
  print*, "time used", t2-t1

  call CPU_TIME(t1)
  call Quad3DRuleDelta(Emesh, 0.5_wp, Wmesh, iter)
  call CPU_TIME(t2)
  print*, "calculated surface area of radius 1 sphere: ", sum(Wmesh)*8
  print*, "relative error: ", (sum(Wmesh)*8-4*PI)/(4*PI)*100,"%"
  print*, "time used", t2-t1

  deallocate(Emesh)
  deallocate(Wmesh)

! An alternative way to test time:
! integer :: t1,t2,COUNT_RATE
! call system_clock(t1,COUNT_RATE)
! call ...(subroutine to be tested)
! call system_clock(t2,COUNT_RATE)
! print*, real(t2-t1)/COUNT_RATE
end program temp