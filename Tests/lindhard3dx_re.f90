program temp
  use BZInt3D
  implicit none

! Testing Quad3DRuleΘ𝔇 and Quad3DRuleΘΘ𝔇.
!                fₖ (1- fₖ₊ₚ)      fₖ (1- fₖ₋ₚ)
! χ(p,ω) = ∫d³k ------------- +  -------------
!                ω-ϵₖ+ϵₖ₊ₚ+iη      -ω-ϵₖ+ϵₖ₋ₚ-iη

  real(wp), allocatable, dimension(:,:,:) :: Emesh,mnsEmesh,Eplusmesh,Emnusmesh,&
                            D1mesh,D2mesh,Wmesh_1,Wmesh_2,Wmesh_3,Wmesh_4,Wmesh
  real(wp) :: q, eF, v, vol
  real(wp), dimension(9) :: arr_y
  integer :: i, j, k, iter
  real(wp), parameter :: PI=4*atan(1.0_wp)
  real :: t1,t2

  eF=0.5_wp
  v=0.5_wp
  q=0.5_wp
  iter=3

  allocate(Emesh(9,9,9),mnsEmesh(9,9,9),Eplusmesh(9,9,9),Emnusmesh(9,9,9),&
           D1mesh(9,9,9),D2mesh(9,9,9),Wmesh_1(9,9,9),Wmesh_2(9,9,9),&
          Wmesh_3(9,9,9),Wmesh_4(9,9,9),Wmesh(9,9,9))

  vol = (1.120_wp+1.125_wp)*(1.120_wp+1.125_wp)*(1.120_wp+1.125_wp)
  arr_y=(/-1.125_wp,-0.844375_wp,-0.56375_wp,-0.283125_wp,-0.0025_wp,&
       0.278125_wp,0.55875_wp,0.839375_wp,1.12_wp/)

  do i=1,9
    do j=1,9
      do k=1,9
        Emesh(i,j,k)=(arr_y(i)*arr_y(i)+arr_y(j)*arr_y(j)+arr_y(k)*arr_y(k))/2
      end do
    end do
  end do

  mnsEmesh=eF-Emesh

  do i=1,9
    do j=1,9
      do k=1,9
        Eplusmesh(i,j,k)=(arr_y(i)*arr_y(i)+arr_y(j)*arr_y(j)+(arr_y(k)+q)*(arr_y(k)+q))/2
      end do
    end do
  end do

  do i=1,9
    do j=1,9
      do k=1,9
        Emnusmesh(i,j,k)=(arr_y(i)*arr_y(i)+arr_y(j)*arr_y(j)+(arr_y(k)-q)*(arr_y(k)-q))/2
      end do
    end do
  end do

  D1mesh=(v-Emesh+Eplusmesh)*4*PI
  D2mesh=(-v-Emesh+Emnusmesh)*4*PI
  Eplusmesh=eF-Eplusmesh
  Emnusmesh=eF-Emnusmesh

  call CPU_TIME(t1)
  call Quad3DRuleThetaFrakD(Emesh,eF,D1mesh,Wmesh_1,iter)
  call Quad3DRuleThetaThetaFrakD(mnsEmesh,Eplusmesh,D1mesh, Wmesh_2,iter)
  call Quad3DRuleThetaFrakD(Emesh,eF,D2mesh,Wmesh_3,iter)
  call Quad3DRuleThetaThetaFrakD(mnsEmesh,Emnusmesh,D2mesh, Wmesh_4,iter)
  call CPU_TIME(t2)

  Wmesh=Wmesh_1-Wmesh_2+Wmesh_3-Wmesh_4

  print*, "q=", q, "v=", v, "iter=", iter
  print*, "Lindhard3D_ReX, result=", -sum(Wmesh)*vol
  print*, "time used:", t2-t1 , "s"

  deallocate(Emesh,mnsEmesh,Wmesh,D1mesh,D2mesh,Eplusmesh,Emnusmesh,Wmesh_1,Wmesh_2,Wmesh_3,Wmesh_4)
end program temp
