!--------------------------------------------------------------
! Calculate inverce matrix |
!---------------------------
!
! N: Number of column &  A
! A(N,N) : Matrix
!
!--------------------------------------------------------------

subroutine inverce_matrix(N,A)
  
  implicit none

  integer,parameter :: iswitch=1 !1:Use eigenvalue decomposition, 2: lapack
  
  integer,intent(in) :: N
  double precision,intent(inout) :: A(N,N)
  
  integer i,inan
  double precision evalue(N),evector(N,N)

  integer info,ipiv(N),lwork
  double precision :: work(64*N)
  
  if(iswitch == 1)then
     
     !---Use eigen vector & value
     call eigen(N,A,evalue,evector)
     A(:,:)=0.d0
     do i=1,N
        A(i,i)=evalue(i)**(-1.d0)
     end do
     A(:,:)=matmul(A(:,:),transpose(evector(:,:)))
     A(:,:)=matmul(evector(:,:),A(:,:))
        
  else if(iswitch == 2)then
     
     !---Use lappack     
     call dgetrf(N,N,A,N,ipiv,info)
     if(info .ne. 0)then
        A(:,:)=0.d0
        inan=-99
        write(6,*)"***Error at inverse matrix dgetrf: info=",info
!        stop
     endif
     
     lwork=int(64*N)
     call dgetri(N,A,N,ipiv,work,lwork,info)
     if(info .ne. 0)then
        A(:,:)=0.d0
        inan=-99
        write(6,*)"***Error at inverse matrix dgetri: info=",info
!        stop
     endif

  end if
  
end subroutine inverce_matrix

