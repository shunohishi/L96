!------------------------------------------------------------------------------------------
! Make random number |
!------------------------------------------------------------------------------------------

subroutine make_original_random(idt,n,out)

  implicit none

  integer,intent(in) :: idt
  integer,intent(in) :: n
  integer seedsize
  integer,allocatable :: seed(:)
  
  double precision,intent(out) :: out(n)
  
  !Random number
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  call random_seed(get=seed)
  seed(:)=seed(:)+idt
  call random_seed(put=seed)
  call random_number(out(:))
  deallocate(seed)
  
end subroutine make_original_random

!----------------------------------------------

subroutine make_random(idt,mean,sigma,n,out)

  implicit none

  integer,intent(in) :: idt
  integer i
  integer,intent(in) :: n
  integer seedsize
  integer,allocatable :: seed(:)
  
  double precision ave,sd
  double precision,intent(in) :: mean,sigma
  double precision,intent(out) :: out(n)
  
  !Random number
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  call random_seed(get=seed)
  seed(:)=seed(:)+idt
  call random_seed(put=seed)
  call random_number(out)
  deallocate(seed)

  !Average
  ave=0.d0
  do i=1,n
     ave=ave+out(i)/dble(n)
  end do

  !Standard deviation
  sd=0.d0
  do i=1,n
     sd=sd+(out(i)-ave)**2.d0/dble(n)
  end do
  sd=sqrt(sd)
  
  !Normalization
  do i=1,n
     out(i)=(out(i)-ave)/sd
  end do

  !Add/Substract mean, inflate/deflate standard deviation
  do i=1,n
     out(i)=out(i)*sigma+mean
  end do
  
!  write(*,*) "Mean:",ave,"--->",mean
!  write(*,*) "SD:",sd,"--->",sigma
  
end subroutine make_random

!------------------------------------------------------------------------------------------
! Gaussian-shaped random number |
!--------------------------------
!
! Use of Box Muller's method
! Mean:0, Sd:1
!
!------------------------------------------------------------------------------------------

subroutine make_gaussian_random(idt,mean,sigma,n,out)

  use parameter  
  implicit none
    
  integer i
  integer,intent(in) :: idt !For different random number
  integer,intent(in) :: n
  integer seedsize
  integer,allocatable :: seed(:)
  
  double precision ave,sd,ske
  double precision,intent(in) :: mean,sigma
  double precision noise1(n),noise2(n)
  double precision,intent(out) :: out(n)
  
  !Random number
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  !Noise1
  call random_seed(get=seed)
  seed(:)=seed(:)+idt
  call random_seed(put=seed)
  call random_number(noise1)
  !Noise2
  call random_seed(get=seed)
  seed(:)=seed(:)+idt+1
  call random_seed(put=seed)
  call random_number(noise2)
  deallocate(seed)

  
  do i=1,n
     out(i)=sqrt(-2.d0*dlog(1.d0-noise1(i)))*dcos(2.d0*pi*noise2(i))
  end do
  
  !Average
  ave=0.d0
  do i=1,n
     ave=ave+out(i)/dble(n)
  end do

  !Standard deviation
  sd=0.d0
  do i=1,n
     sd=sd+(out(i)-ave)**2.d0/dble(n)
  end do
  sd=sqrt(sd)

  !Skewness
  ske=0.d0
  do i=1,n
     ske=ske+((out(i)-ave)/sd)**3.d0
  end do
  ske=dble(n)/(dble(n-1.d0)*dble(n-2.d0))*ske
  
  !Normalization
  do i=1,n
     out(i)=(out(i)-ave)/sd
  end do

  !Add/Substract mean, inflate/deflate standard deviation
  do i=1,n
     out(i)=out(i)*sigma+mean
  end do  
  
!  write(*,*) "Mean:",ave,"--->",mean
!  write(*,*) "SD:",sd,"--->",sigma
!  write(*,*) "Skewness:",ske
  
end subroutine make_gaussian_random

!-----------------------------------------------------------------------------------------------------
! Integer random number without duplication |
!-----------------------------------------------------------------------------------------------------

subroutine make_int_random(idt,snum,enum,nnum,out)

  implicit none

  integer inum,jnum
  
  integer,intent(in) :: idt
  integer,intent(in) :: snum,enum !Start,End edge
  integer,intent(in) :: nnum      !Number of random number
  integer,intent(out) :: out(nnum)
  
  integer count
  integer jdt

  integer seedsize
  integer,allocatable :: seed(:)
  
  double precision random(nnum)
  double precision dtmp(1)
  
  !Initial random number
  !Random number
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  call random_seed(get=seed)
  seed(:)=seed(:)+idt
  call random_seed(put=seed)
  call random_number(random)
  !Make integer random number
  out(:)=int(random(:)*(enum-snum+1)+snum)
  deallocate(seed)
  
  !Check duplication
  do inum=1,nnum

     count=0
     jdt=0
     
     do while(count < nnum)   
        do jnum=1,nnum
           
           count=count+1
           
           if(inum == jnum)cycle           
           if(out(inum) == out(jnum))then
              call make_original_random(jdt,1,dtmp)
              out(inum)=int(dtmp(1)*(enum-snum+1)+snum)
              count=0
              jdt=jdt+1
              exit
           end if           
           
        end do !jnum
        
     end do !while
  end do !inum
  
end subroutine make_int_random
