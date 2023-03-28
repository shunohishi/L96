subroutine int_sort(n,x)

  implicit none

  integer i,j
  integer itmp

  integer,intent(in) :: n
  integer x(n)


  !Sort
  do i=1,n

     j=i+1
     
     do while(j <= n)

        if(x(j) < x(i))then
           itmp=x(i)
           x(i)=x(j)
           x(j)=itmp
        else
           j=j+1
        end if        

     end do
  end do
  
end subroutine int_sort
