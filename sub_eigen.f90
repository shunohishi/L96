subroutine eigen(n,A,evalue,evector)

  implicit none

  integer i,j
  integer,intent(in) :: n
  double precision,intent(in) :: A(n,n)        !Symmetric Matrix
  double precision,intent(out) :: evalue(n)    !Eigen value
  double precision,intent(out) :: evector(n,n) !Eigen vector
  integer inan
  
  integer lwork,info
  double precision,allocatable :: work(:)
  
  lwork=3*n-1
  allocate(work(lwork))
  evector(:,:)=A(:,:)
  call dsyev('V','U',n,evector,n,evalue,work,lwork,info)

  if(info /= 0)then
     write(*,*) "***Error: Eigen",info
     inan=-99
     stop
  end if
  deallocate(work)
  
end subroutine eigen
