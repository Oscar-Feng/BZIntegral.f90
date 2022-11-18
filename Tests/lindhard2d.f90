program temp
  use BZInt2D
  implicit none

! Testing Quad2DRuleΘ𝔇 and Quad2DRuleΘδ.
!                  fₖ - fₖ₊ₚ
! χ(p,ω) = ∫d²k --------------
!                ϵₖ-ϵₖ₊ₚ+ω+iη

  real(wp), allocatable, dimension(:,:) :: Emesh,Eplusmesh,Dmesh,Wmesh_1,Wmesh_2,Wmesh
  real(wp) :: q, eF, v, vol
  real(wp), allocatable, dimension(:) :: arr_x
  real(wp), dimension(9) :: arr_y
  integer :: Nx, i, j, k, iter, qq
  real(wp), parameter :: PI=4*atan(1.0_wp)
  integer :: t1,t2

  eF=0.5_wp
  v=0.0_wp
  iter=4
  q=0.2_wp

  Nx=ceiling(8+4*q)
  if (mod(Nx,2)==0) then
    Nx=Nx+1
  end if

  allocate(Emesh(9,Nx),Eplusmesh(9,Nx),Dmesh(9,Nx),Wmesh_1(9,Nx),&
           Wmesh_2(9,Nx),Wmesh(9,Nx))
  allocate(arr_x(Nx))

  vol = (1.120_wp+1.125_wp+q)*(1.120_wp+1.125_wp)
  arr_y=(/-1.125_wp,-0.844375_wp,-0.56375_wp,-0.283125_wp,-0.0025_wp,&
       0.278125_wp,0.55875_wp,0.839375_wp,1.12_wp/)
  do i=1,Nx
    arr_x(i)=(-1.125_wp-q)+(i-1)*(2.245_wp+q)/(Nx-1)
  end do

  do i=1,9
    do k=1,Nx
      Emesh(i,k)=(arr_y(i)*arr_y(i)+arr_x(k)*arr_x(k))/2
    end do
  end do

  do i=1,9
    do k=1,Nx
      Eplusmesh(i,k)=(arr_y(i)*arr_y(i)+(arr_x(k)+q)*(arr_x(k)+q))/2
    end do
  end do

  do i=1,9
    do k=1,Nx
      Dmesh(i,k)=(v+Emesh(i,k)-Eplusmesh(i,k))*2*PI
    end do
  end do

  call Quad2DRuleThetaFrakD(Emesh,eF,Dmesh, Wmesh_1,iter)
  call Quad2DRuleThetaFrakD(Eplusmesh,eF,Dmesh, Wmesh_2,iter)

  Wmesh=Wmesh_1-Wmesh_2

  print*, "q=", q, "v=", v, "iter=", iter
  print*, "Lindhard3D_Re, result=", sum(Wmesh)*vol

  call Quad3DRuleThetaDelta(Emesh,eF,Dmesh, Wmesh_1,iter)
  call Quad3DRuleThetaDelta(Eplusmesh,eF,Dmesh, Wmesh_2,iter)

  Wmesh=Wmesh_1-Wmesh_2

  print*, "q=", q, "v=", v, "iter=", iter
  print*, "Lindhard3D_Im, result=", -sum(Wmesh)*vol*PI

  deallocate(Emesh,Wmesh,Dmesh,Eplusmesh,Wmesh_1,Wmesh_2)
  deallocate(arr_x)

end program temp
