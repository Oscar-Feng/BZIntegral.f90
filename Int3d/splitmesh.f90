module SplitMesh
  implicit none
  integer, dimension(3,4,6), parameter :: box2tetra=reshape((/0,0,0,  1,0,0,&
      0,1,0,  1,0,1,    1,1,0,  1,0,0,  0,1,0,  1,0,1,    0,0,0,  0,0,1,&
      0,1,0,  1,0,1,    0,1,1,  0,0,1,  0,1,0,  1,0,1,    1,1,0,  1,1,1,&
      0,1,0,  1,0,1,    0,1,1,  1,1,1,  0,1,0,  1,0,1/),(/3,4,6/))
  integer, dimension(3,10,6), parameter :: box2qtetra=reshape((/0,0,0,  2,0,0,&
      0,2,0,  2,0,2,  1,0,0,  0,1,0,  1,0,1,  1,1,0,  1,1,1,  2,0,1,&
      2,2,0,  2,0,0,  0,2,0,  2,0,2,  2,1,0,  1,2,0,  2,1,1,  1,1,0,  1,1,1,&
      2,0,1,      0,0,0,  0,0,2,  0,2,0,  2,0,2,  0,0,1,  0,1,0,  1,0,1,&
      0,1,1,  1,1,1,  1,0,2,      0,2,2,  0,0,2,  0,2,0,  2,0,2,  0,1,2,&
      0,2,1,  1,1,2,  0,1,1,  1,1,1,  1,0,2,      2,2,0,  2,2,2,  0,2,0,&
      2,0,2,  2,2,1,  1,2,0,  2,1,1,  1,2,1,  1,1,1,  2,1,2,      0,2,2,&
      2,2,2,  0,2,0,  2,0,2,  1,2,2,  0,2,1,  1,1,2,  1,2,1,  1,1,1,&
      2,1,2/),(/3,10,6/))
contains
! Split uniform mesh with open boundary into tetrahedrons.
! Usually we get periodic mesh at start, then for this function to
!   work, we need add 1 to size of mesh in each dimensions.
! In the integer arrays Qtetras and boxes, the elements in the first
!   dimension stand for x,y,z separately
   subroutine Mesh2Tetra(Nxplus1, Nyplus1, Nzplus1, Tetras)
     integer, intent(in) :: Nxplus1, Nyplus1, Nzplus1
     integer, allocatable, dimension(:,:,:), intent(out) :: Tetras
     integer, allocatable, dimension(:,:) :: boxes
     integer :: leng, i, j, k, count

     leng=(Nxplus1-1)*(Nyplus1-1)*(Nzplus1-1)
     allocate(boxes(3,leng))

     count=1
     do i=1, Nzplus1-1
       do j=1, Nyplus1-1
         do k=1, Nxplus1-1
           boxes(1,count)=k
           boxes(2,count)=j
           boxes(3,count)=i
           count=count+1
         end do
       end do
     end do

     allocate(Tetras(3,4,leng*6))
     do i=1,leng
       do j=1,6
         count=j+(i-1)*6
         Tetras(1,:,count)=boxes(1,i)+box2tetra(1,:,j)
         Tetras(2,:,count)=boxes(2,i)+box2tetra(2,:,j)
         Tetras(3,:,count)=boxes(3,i)+box2tetra(3,:,j)
       end do
     end do
     deallocate(boxes)
   end subroutine Mesh2Tetra

! In order for this splitting to be possible, mesh should have
!   even numbers of space in each dimension. (So the input integers
!   should be odd numbers)
   subroutine Mesh2QTetra(Nxplus1, Nyplus1, Nzplus1, Qtetras)
     integer, intent(in) :: Nxplus1, Nyplus1, Nzplus1
     integer, allocatable, dimension(:,:,:), intent(out) :: QTetras
     integer, allocatable, dimension(:,:) :: boxes
     integer :: leng, i, j, k, count

     leng=(Nxplus1-1)*(Nyplus1-1)*(Nzplus1-1)/8
     allocate(boxes(3,leng))

     count=1
     do i=1, Nzplus1-1, 2
       do j=1, Nyplus1-1, 2
         do k=1, Nxplus1-1, 2
           boxes(1,count)=k
           boxes(2,count)=j
           boxes(3,count)=i
           count=count+1
         end do
       end do
     end do

     allocate(QTetras(3,10,leng*6))
     do i=1,leng
       do j=1,6
         count=j+(i-1)*6
         QTetras(1,:,count)=boxes(1,i)+box2qtetra(1,:,j)
         QTetras(2,:,count)=boxes(2,i)+box2qtetra(2,:,j)
         QTetras(3,:,count)=boxes(3,i)+box2qtetra(3,:,j)
       end do
     end do
     deallocate(boxes)
   end subroutine Mesh2QTetra
end module SplitMesh
