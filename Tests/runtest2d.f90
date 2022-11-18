program temp
  use BZInt2D
  implicit none

! Testing Quad2DRuleΘ and Quad2DRuleδ. In order to compare between Fortran and
! Julia, one may need to add a space between "!" and "$OMP..." in the
! corresponding subroutine in bzint2d.f90 to disable parallelism
  integer :: i,j,iter
  real(wp), allocatable, dimension(:,:) :: Emesh
  real(wp), allocatable, dimension(:,:) :: Wmesh
  real(wp),dimension(9) :: arr
  real :: t1, t2
  real(wp), parameter :: PI=4*atan(1.0_wp)

  allocate(Emesh(9,9),Wmesh(9,9))

  iter=3
  arr=(/-1.0_wp, -0.75_wp, -0.5_wp, -0.25_wp, 0.0_wp, 0.25_wp, 0.5_wp, 0.75_wp, 1.0_wp/)

  do i=1,9
    do j=1,9
      Emesh(i,j)=arr(i)*arr(i)+arr(j)*arr(j)
    end do
  end do

  call CPU_TIME(t1)
  call Quad2DRuleTheta(Emesh, 1.0_wp, Wmesh, iter)
  call CPU_TIME(t2)
  print*, "calculated area of radius 1 circle:", sum(Wmesh)*4
  print*, "relative error: ", (sum(Wmesh)*4-pi)/pi*100, "%"
  print*, "time used:", t2-t1, "s"

  call CPU_TIME(t1)
  call Quad2DRuleDelta(Emesh, 1.0_wp, Wmesh, iter)
  call CPU_TIME(t2)
  print*, "calculated perimeter of radius 1 circle:", sum(Wmesh)*4
  print*, "relative error: ", (sum(Wmesh)*4-pi)/pi*100, "%"
  print*, "time used:", t2-t1 , "s"

  deallocate(Emesh)
  deallocate(Wmesh)

! Another way to test time:
! integer :: t1,t2,COUNT_RATE
! call system_clock(t1,COUNT_RATE)
! call ...(subroutine to be tested)
! call system_clock(t2,COUNT_RATE)
! print*, real(t2-t1)/COUNT_RATE
end program temp