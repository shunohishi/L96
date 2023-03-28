subroutine make_pco(npco,pco)

  implicit none

  integer ipco
  
  integer,intent(in) :: npco
  double precision,intent(out) :: pco(npco)

  do ipco=1,npco
     pco(ipco)=-1.d0+0.1*dble(ipco-1)
  end do
  
end subroutine make_pco

!------------------------------------------------------------------

subroutine make_xferr(nens,xt,xf,xferr)

  use parameter
  implicit none

  integer ix,iens
  double precision xfmean(nx)
  
  integer,intent(in) :: nens

  double precision,intent(in) :: xt(nx),xf(nx,nens)
  double precision,intent(out) :: xferr(nx)


  xfmean(:)=0.d0
  do iens=1,nens
     do ix=1,nx
        xfmean(ix)=xfmean(ix)+xf(ix,iens)
     end do
  end do

  xfmean(:)=xfmean(:)/dble(nens)

  xferr(:)=xfmean(:)-xt(:)
  
end subroutine make_xferr

!-----------------------------------------------------------------

subroutine make_inf_nc(inf_parm)

  use parameter
  implicit none

  integer iinf
  double precision inf_parm(ninf)

  do iinf=1,ninf
     if(1 <= iinf .and. iinf <= 11)then     !1.00, 1.01, ..., 1.10
        inf_parm(iinf)=1.d0+0.01d0*dble(iinf-1)
     else if(12 <= iinf .and. iinf <= 20)then !1.20, 1.30, ..., 2.00
        inf_parm(iinf)=1.2d0+0.1d0*dble(iinf-12)
     else if(21 <= iinf .and. iinf <= 23)then   !3.00, ..., 5.00
        inf_parm(iinf)=3.0d0+1.0d0*dble(iinf-21)
     end if
  end do  
  
end subroutine make_inf_nc

!-----------------------------------------------------------------

subroutine make_R(iR,nobs,R)

  implicit none

  integer iobs
  
  integer,intent(in) :: iR,nobs
  double precision,intent(out) :: R(nobs,nobs)

  R(:,:)=0.d0
  do iobs=1,nobs
     if(1 <= iR .and. iR <= 11)then !R = 1.0, 1.1, ..., 2.0
        R(iobs,iobs)=1.d0+0.1d0*dble(iR-1)
     else if(12 <= iR .and. iR <= 18)then !R=3.0, ..., 20.0
        R(iobs,iobs)=3.d0+1.0d0*dble(iR-12)
     end if
  end do
  
end subroutine make_R

!-----------------------------------------------------------------

subroutine est_id_inf(rmse,id_inf)

  use parameter
  implicit none

  integer ipco,iinf
  
  double precision rmse_min(npco)
  
  double precision,intent(in) :: rmse(ninf,npco)
  integer,intent(out) :: id_inf(npco)

  id_inf(:)=0
  
  do ipco=1,npco
     rmse_min(ipco)=minval(rmse(:,ipco))   
  end do

  do ipco=1,npco
     do iinf=1,ninf
        if(rmse(iinf,ipco) == rmse_min(ipco))then
           id_inf(ipco)=iinf
        end if
     end do
  end do

end subroutine est_id_inf

!---------------------------------------------------------------

subroutine est_id_inf_R(rmse,id_inf,id_R)

  use parameter
  implicit none

  integer ipco,iR,iinf
  
  double precision rmse_min(npco)
  
  double precision,intent(in) :: rmse(ninf,nR,npco)
  integer,intent(out) :: id_inf(npco),id_R(npco)

  id_inf(:)=0
  id_R(:)=0
  
  do ipco=1,npco
     rmse_min(ipco)=minval(rmse(:,:,ipco))   
  end do

  do ipco=1,npco
     do iR=1,nR
        do iinf=1,ninf
           if(rmse(iinf,iR,ipco) == rmse_min(ipco))then
              id_inf(ipco)=iinf
              id_R(ipco)=iR
           end if
        end do
     end do
  end do
  
end subroutine est_id_inf_R

!-------------------------------------------------------------

subroutine est_id2_R(rmse,id2_R)

  use parameter
  implicit none

  integer ipco,iR,iinf
  double precision rmse_min(ninf,npco)
  
  double precision,intent(in) :: rmse(ninf,nR,npco)
  integer,intent(out) :: id2_R(ninf,npco)

  id2_R(:,:)=0

  do ipco=1,npco
     do iinf=1,ninf
        rmse_min(iinf,ipco)=minval(rmse(iinf,:,ipco))
     end do
  end do

  do ipco=1,npco
     do iR=1,nR
        do iinf=1,ninf
           if(rmse(iinf,iR,ipco) == rmse_min(iinf,ipco))then
              id2_R(iinf,ipco)=iR
           end if
        end do
     end do
  end do
  
end subroutine est_id2_R
