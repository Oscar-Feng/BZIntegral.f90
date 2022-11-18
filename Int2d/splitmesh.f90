module SplitMesh
  implicit none
  integer, dimension(2,6,2), parameter :: box2qtrig=reshape((/1,1, 3,1, 1,3,&
    2,1, 1,2, 2,2,  3,3, 1,3, 3,1, 2,3, 3,2, 2,2/),(/2,6,2/))
  integer, dimension(2,3,2), parameter :: rect2trig=reshape((/1,1, 2,1, 1,2,&
    2,2, 1,2, 2,1/),(/2,3,2/))
contains
  ! In order for this splitting to be possible, mesh should have
  ! even numbers of space in each dimension. (So the input integers
  ! should be odd numbers)

  ! In the integer arrays QTrigs(Trigs) and boxes, the elements in the first
  ! dimension stand for x,y separately
  subroutine Mesh2Trig(Nxplus1, Nyplus1, Trigs)
    integer, intent(in) :: Nxplus1, Nyplus1
    integer, allocatable, dimension(:,:,:),intent(out):: Trigs
    integer, allocatable, dimension(:,:) :: boxes
    integer :: leng, i,j,count

    leng=(Nxplus1-1)*(Nyplus1-1)
    allocate(boxes(2,leng))

    count=1
    do i=0, Nyplus1-2
      do j=0, Nxplus1-2
        boxes(1,count)=j
        boxes(2,count)=i
        count=count+1
      end do
    end do

    allocate(Trigs(2,3,leng*2))
    do i=1,leng
      do j=1,2
        count=j+(i-1)*2
        Trigs(1,:,count)=boxes(1,i)+rect2trig(1,:,j)
        Trigs(2,:,count)=boxes(2,i)+rect2trig(2,:,j)
      end do
    end do
    deallocate(boxes)
  end subroutine Mesh2Trig

  subroutine Mesh2QTrig(Nxplus1, Nyplus1, QTrigs)
    integer, intent(in) :: Nxplus1, Nyplus1
    integer, allocatable, dimension(:,:,:),intent(out):: QTrigs
    integer, allocatable, dimension(:,:) :: boxes
    integer :: leng, i,j,count

    leng=(Nxplus1-1)*(Nyplus1-1)/4
    allocate(boxes(2,leng))

    count=1
    do i=0, Nyplus1-2, 2
      do j=0, Nxplus1-2, 2
         boxes(1,count)=j
         boxes(2,count)=i
         count=count+1
       end do
     end do

     allocate(QTrigs(2,6,leng*2))
     do i=1,leng
       do j=1,2
         count=j+(i-1)*2
         QTrigs(1,:,count)=boxes(1,i)+box2qtrig(1,:,j)
         QTrigs(2,:,count)=boxes(2,i)+box2qtrig(2,:,j)
       end do
     end do
     deallocate(boxes)
  end subroutine Mesh2QTrig
end module SplitMesh
