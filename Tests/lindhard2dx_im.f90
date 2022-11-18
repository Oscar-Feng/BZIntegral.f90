program temp
  use BZInt2D
  implicit none

! Testing Quad2DRuleΘδ and Quad2DRuleΘΘδ.
!                fₖ (1- fₖ₊ₚ)      fₖ (1- fₖ₋ₚ)
! χ(p,ω) = ∫d²k ------------- +  -------------
!                ω-ϵₖ+ϵₖ₊ₚ+iη      -ω-ϵₖ+ϵₖ₋ₚ-iη

  real(wp) :: q,v,eF,vol
  integer :: i,j,iter
  real(wp), dimension(9) :: arr_y
  real(wp), allocatable, dimension(:,:) :: Emesh,mnsEmesh,Eplusmesh,Emnsmesh,Dmesh_1,Dmesh_2,Wmesh_1,Wmesh_2,Wmesh_3,Wmesh_4,Wmesh
  real(wp), parameter :: PI=4*atan(1.0_wp)
  real :: t1, t2

  q=0.5_wp
  v=0.5_wp
  eF=0.5_wp
  iter=2

  arr_y=(/-1.125_wp,-0.844375_wp,-0.56375_wp,-0.283125_wp,-0.0025_wp,&
       0.278125_wp,0.55875_wp,0.839375_wp,1.12_wp/)

  vol = (1.120_wp+1.125_wp)*(1.120_wp+1.125_wp)
  allocate(Emesh(9,9),mnsEmesh(9,9),Eplusmesh(9,9),Emnsmesh(9,9),Dmesh_1(9,9),Dmesh_2(9,9),&
         Wmesh_1(9,9),Wmesh_2(9,9),Wmesh_3(9,9),Wmesh_4(9,9),Wmesh(9,9))

  do i=1,9
    do j=1,9
      Emesh(i,j)=(arr_y(i)*arr_y(i)+arr_y(j)*arr_y(j))/2
    end do
  end do

  do i=1,9
    do j=1,9
      Eplusmesh(i,j)=(arr_y(i)*arr_y(i)+(arr_y(j)+q)*(arr_y(j)+q))/2
    end do
  end do

  do i=1,9
    do j=1,9
      Emnsmesh(i,j)=(arr_y(i)*arr_y(i)+(arr_y(j)-q)*(arr_y(j)-q))/2
    end do
  end do

  do i=1,9
    do j=1,9
      Dmesh_1(i,j)=(v-Emesh(i,j)+Eplusmesh(i,j))*2*PI
    end do
  end do

  do i=1,9
    do j=1,9
      Dmesh_2(i,j)=(-v-Emesh(i,j)+Emnsmesh(i,j))*2*PI
    end do
  end do

  do i=1,9
    do j=1,9
      Eplusmesh(i,j)=eF-Eplusmesh(i,j)
    end do
  end do

  do i=1,9
    do j=1,9
      Emnsmesh(i,j)=eF-Emnsmesh(i,j)
    end do
  end do

  do i=1,9
    do j=1,9
      mnsEmesh(i,j)=eF-Emesh(i,j)
    end do
  end do

  call CPU_TIME(t1)
  call Quad2DRuleThetaDelta(Emesh, eF, Dmesh_1, Wmesh_1, iter)
  call Quad2DRuleThetaDelta(Emesh, eF, Dmesh_2, Wmesh_3, iter)
  call Quad2DRuleThetaThetaDelta(mnsEmesh,Eplusmesh,Dmesh_1,Wmesh_2,iter)
  call Quad2DRuleThetaThetaDelta(mnsEmesh,Emnsmesh,Dmesh_2,Wmesh_4,iter)
  call CPU_TIME(t2)

  do i=1,9
    do j=1,9
      Wmesh(i,j)=Wmesh_1(i,j)-Wmesh_2(i,j)-Wmesh_3(i,j)+Wmesh_4(i,j)
    end do
  end do

  print*, "q=", q, "v=", v
  print*, "result=", sum(Wmesh)*vol*PI
  print*, "time used:", t2-t1 , "s"

  deallocate(Emesh,mnsEmesh,Eplusmesh,Emnsmesh,Dmesh_1,Dmesh_2,Wmesh_1,Wmesh_2,Wmesh_3,Wmesh_4,Wmesh)
end program temp
