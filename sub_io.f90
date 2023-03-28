!------------------------------------------------------------------------------
! KF/3DVAR |
!------------------------------------------------------------------------------

subroutine write_data(nobs,iinf,it_io,it,xt,x,xf,xa,obs,inf_parm)

  use parameter
  implicit none

  integer,intent(in) :: nobs,iinf
  integer,intent(in) :: it_io,it
  integer ix
  
  double precision,intent(in) :: xt(nx),x(nx),xf(nx),xa(nx),obs(nx)
  double precision,intent(in) :: inf_parm
    
  if(nobs == 1 .and. iinf == 1 .and. it_io == 0)then
     open(1,file=trim(dir(execute_da))//"/x.dat",status="replace")
  else
     open(1,file=trim(dir(execute_da))//"/x.dat",access="append")
  end if
  do ix=1,nx
     write(1,'(i6,E15.4,i6,6f12.5)') nobs,inf_parm,ix,dt*dble(it)/int_day,xt(ix),x(ix),xf(ix),xa(ix),obs(ix)
  end do
  close(1)

end subroutine write_data

!---------------------------

subroutine write_P(nobs,iinf,it_io,it,Pf,Pa,inf_parm)

  use parameter
  implicit none

  integer,intent(in) :: nobs,iinf
  integer,intent(in) :: it_io,it
  integer ix
  
  double precision,intent(in) :: Pf(nx,nx),Pa(nx,nx)
  double precision,intent(in) :: inf_parm
    
  if(nobs == 1 .and. iinf == 1 .and. it_io == 0)then
     open(11,file=trim(dir(execute_da))//"/P.dat",status="replace")
  else
     open(11,file=trim(dir(execute_da))//"/P.dat",access="append")
  end if
  do ix=1,nx
     write(11,'(i6,E15.4,i6,3f12.5)') nobs,inf_parm,ix,dt*dble(it)/int_day,Pf(ix,ix),Pa(ix,ix)
  end do
  close(11)
  
end subroutine write_P

!-----------------------------

subroutine write_bias_rmse(nobs,iinf,st_obs,it,bias,biasf,biasa,biasobs,rmse,rmsef,rmsea,rmseobs,inf_parm)

  use parameter
  implicit none

  integer,intent(in) :: nobs,iinf
  integer,intent(in) :: st_obs,it
  double precision,intent(in) :: bias,biasf,biasa,biasobs
  double precision,intent(in) :: rmse,rmsef,rmsea,rmseobs
  double precision,intent(in) :: inf_parm
  
  if(nobs == 1 .and. iinf == 1 .and. st_obs == it)then
     open(31,file=trim(dir(execute_da))//"/bias_rmse.dat",status="replace")
  else
     open(31,file=trim(dir(execute_da))//"/bias_rmse.dat",access="append")
  end if
  write(31,'(i6,E15.4,9f12.5)') nobs,inf_parm,dt*dble(it)/int_day,bias,biasf,biasa,biasobs,rmse,rmsef,rmsea,rmseobs
  close(31)

end subroutine write_bias_rmse

!--------------------------------

subroutine write_ave_bias_rmse(nobs,iinf,bias,biasf,biasa,biaso, &
     & rmse,rmsef,rmsea,rmseo,inf_parm)

  use parameter
  implicit none

  integer,intent(in) :: nobs,iinf

  double precision,intent(in) :: bias,biasf,biasa,biaso
  double precision,intent(in) :: rmse,rmsef,rmsea,rmseo
  double precision,intent(in) :: inf_parm

  if(nobs == 1 .and. iinf == 1)then
     open(32,file=trim(dir(execute_da))//"/ave_bias_rmse.dat",status="replace")
  else
     open(32,file=trim(dir(execute_da))//"/ave_bias_rmse.dat",access="append")
  end if
  write(32,'(i6,E15.4,8f12.5)') nobs,inf_parm,bias,biasf,biasa,biaso,rmse,rmsef,rmsea,rmseo
  close(32)

end subroutine write_ave_bias_rmse

!--------------------------------------------

subroutine write_ave_Pf_Pa(nobs,iinf,ifo,sfo,inf_parm, &
                   & Pfdave,Pfdstd,Padave,Padstd, &
                   & Pfcdave,Pfcdstd,Pacdave,Pacdstd, &
                   & Pfncdave,Pfncdstd,Pancdave,Pancdstd)

  use parameter
  implicit none
  
  integer,intent(in) :: nobs,iinf,ifo,sfo
  
  double precision,intent(in) :: inf_parm
  double precision,intent(in) :: Pfdave,Pfdstd,Padave,Padstd
  double precision,intent(in) :: Pfcdave,Pfcdstd,Pacdave,Pacdstd
  double precision,intent(in) :: Pfncdave,Pfncdstd,Pancdave,Pancdstd

  if(nobs == nobs .and. iinf == 1 .and. ifo == sfo)then
     open(61,file=trim(dir(execute_da))//"/ave_Pf.dat",status="replace")
     open(62,file=trim(dir(execute_da))//"/ave_Pa.dat",status="replace")
  else
     open(61,file=trim(dir(execute_da))//"/ave_Pf.dat",access="append")
     open(62,file=trim(dir(execute_da))//"/ave_Pa.dat",access="append")
  endif
  write(61,'(i6,E15.4,7f12.5)') nobs,inf_parm,0.1d0*ifo, &
       & Pfdave,Pfdstd,Pfcdave,Pfcdstd,Pfncdave,Pfncdstd
  write(62,'(i6,E15.4,7f12.5)') nobs,inf_parm,0.1d0*ifo, &
       & Padave,Padstd,Pacdave,Pacdstd,Pancdave,Pancdstd
  close(61)
  close(62)
  
end subroutine write_ave_Pf_Pa

!------------------------------------------------------------------------------------------
! EnKF |
!------------------------------------------------------------------------------------------

subroutine write_ensmean(nens,nobs,iinf,sigma_loc,it_io,it, &
     & xt,x,xf,xa,obserr,inf_parm)

  use parameter
  implicit none

  integer,intent(in) :: nens,nobs,iinf,sigma_loc
  integer,intent(in) :: it_io,it
  integer ix
  
  double precision,intent(in) :: xt(nx),x(nx)
  double precision,intent(in) :: xf(nx),xa(nx)
  double precision,intent(in) :: obserr(nx)
  double precision,intent(in) :: inf_parm
    
!  if(nens == 16 .and. nobs == nx .and. iinf == 1 .and. sigma_loc == 1 .and. cor_set == -1.d0 .and. it_io == 0)then
!     open(1,file=trim(dir(execute_da))//"/x.dat",status="replace")
!  else
  open(1,file=trim(dir(execute_da))//"/x.dat",access="append")
!  end if
  
  do ix=1,nx
     write(1,'(2i6,f12.5,2i6,6f12.5)') nens,nobs,inf_parm,sigma_loc,ix,dt*dble(it)/int_day, &
             & xt(ix),x(ix),xf(ix),xa(ix),obserr(ix)
  end do
  
  close(1)

end subroutine write_ensmean

!------------------------------

subroutine write_enssprd(nens,nobs,iinf,sigma_loc,it_io,it, &
     & xfsprd,xasprd,inf_parm)

  use parameter
  implicit none

  integer,intent(in) :: nens,nobs,iinf,sigma_loc
  integer,intent(in) :: it_io,it
  integer ix

  double precision,intent(in) :: xfsprd(nx),xasprd(nx)
  double precision xfsprdmean,xasprdmean
  double precision,intent(in) :: inf_parm
  
  xfsprdmean=0.d0
  xasprdmean=0.d0

  do ix=1,nx
     xfsprdmean=xfsprdmean+xfsprd(ix)/dble(nx)
     xasprdmean=xasprdmean+xasprd(ix)/dble(nx)
  end do
  
!  if(nens == 16 .and. nobs == nx .and. iinf == 1 .and. sigma_loc == 1 .and. it_io == 0)then
!     open(21,file=trim(dir(execute_da))//"/sprd.dat",status="replace")
!     open(22,file=trim(dir(execute_da))//"/sprdmean.dat",status="replace")
!  else
  open(21,file=trim(dir(execute_da))//"/sprd.dat",access="append")
  open(22,file=trim(dir(execute_da))//"/sprdmean.dat",access="append")
!  end if
  
  do ix=1,nx
     write(21,'(2i6,E15.4,2i6,3f12.5)') nens,nobs,inf_parm,sigma_loc,ix,dt*dble(it)/int_day,xfsprd(ix),xasprd(ix)
  end do
  write(22,'(2i6,E15.4,i6,3f12.5)') nens,nobs,inf_parm,sigma_loc,dt*dble(it)/int_day,xfsprdmean,xasprdmean
  
  close(21)
  close(22)
  
end subroutine write_enssprd

!------------------------------

subroutine write_ens_bias_rmse(nens,nobs,iinf,sigma_loc,st_obs,it, &
     & bias,biasf,biasa,biaso, &
     & rmse,rmsef,rmsea,rmseo, &
     & inf_parm)

  use parameter
  implicit none

  integer,intent(in) :: nens,nobs,iinf,sigma_loc
  integer,intent(in) :: st_obs,it
  double precision,intent(in) :: bias,biasf,biasa,biaso
  double precision,intent(in) :: rmse,rmsef,rmsea,rmseo
  double precision,intent(in) :: inf_parm
  
  !  if(nens == 16 .and. nobs == nx .and. iinf == 1 .and. sigma_loc == 1 .and. st_obs == it)then
  !     open(31,file=trim(dir(execute_da))//"/bias_rmse.dat",status="replace")
  !  else
  open(30,file=trim(dir(execute_da))//"/bias.dat",access="append")
  open(31,file=trim(dir(execute_da))//"/rmse.dat",access="append")
  !  end if
  write(30,'(2i6,E15.4,i6,5f12.5)') nens,nobs,inf_parm,sigma_loc,dt*dble(it)/int_day, &
       & bias,biasf,biasa,biaso
  write(31,'(2i6,E15.4,i6,5f12.5)') nens,nobs,inf_parm,sigma_loc,dt*dble(it)/int_day, &
       & rmse,rmsef,rmsea,rmseo
  close(30)
  close(31)

end subroutine write_ens_bias_rmse

!--------------------------------

subroutine write_ens_ave_bias_rmse_sprd(nens,nobs,iinf,sigma_loc, &
     & bias,biasf,biasa,biaso, &
     & rmse,rmsef,rmsea,rmseo, &
     & sprdf,sprda, &
     & inf_parm)

  use parameter
  implicit none

  integer,intent(in) :: nens,nobs,iinf,sigma_loc
  
  double precision,intent(in) :: bias,biasf,biasa,biaso
  double precision,intent(in) :: rmse,rmsef,rmsea,rmseo
  double precision,intent(in) :: sprdf,sprda
  double precision,intent(in) :: inf_parm

!  if(nens == 16 .and. nobs == nx .and. iinf == 1 .and. sigma_loc == 1)then
!     open(32,file=trim(dir(execute_da))//"/ave_bias_rmse.dat",status="replace")
!  else
  open(32,file=trim(dir(execute_da))//"/ave_bias.dat",access="append")
  open(33,file=trim(dir(execute_da))//"/ave_rmse.dat",access="append")
  open(34,file=trim(dir(execute_da))//"/ave_sprd.dat",access="append")
!  end if
  write(32,'(2i6,E15.4,i6,4f12.5)') nens,nobs,inf_parm,sigma_loc, &
       & bias,biasf,biasa,biaso
  write(33,'(2i6,E15.4,i6,4f12.5)') nens,nobs,inf_parm,sigma_loc, &
       & rmse,rmsef,rmsea,rmseo
  write(34,'(2i6,E15.4,i6,2f12.5)') nens,nobs,inf_parm,sigma_loc, &
       & sprdf,sprda
  close(32)
  close(33)
  close(34)

end subroutine write_ens_ave_bias_rmse_sprd

!----------------------------------

subroutine write_Neff(nens,nobs,iinf,sigma_loc,ave_Neff,std_Neff,inf_parm)

  use parameter
  implicit none

  integer,intent(in) :: nens,nobs,iinf,sigma_loc

  double precision,intent(in) :: ave_Neff,std_Neff
  double precision,intent(in) :: inf_parm

  if(nens == 16 .and. nobs == nx .and. iinf == 1 .and. sigma_loc == 1)then
     open(35,file=trim(dir(execute_da))//"/Neff.dat",status="replace")
  else
     open(35,file=trim(dir(execute_da))//"/Neff.dat",access="append")
  end if
  write(35,'(2i6,E15.4,i6,2f12.5)') nens,nobs,inf_parm,sigma_loc,ave_Neff,std_Neff
  close(35)
  
end subroutine write_Neff

!---------------------------------------------------------------------------------------
! Growth rate |
!---------------------------------------------------------------------------------------

subroutine write_growth_rate(it,ave_bias,std_bias,ave_rmse,std_rmse)

  use parameter
  implicit none
  
  integer,intent(in) :: it
  integer status,system
  
  double precision,intent(in) :: ave_bias,std_bias,ave_rmse,std_rmse

  status=system("mkdir -p GR")
  
  if(it == 1)then
     open(41,file="GR/growth_rate.dat",status="replace")
  else
     open(41,file="GR/growth_rate.dat",access="append")
  end if
  write(41,'(5f12.5)') dble(it)*dt/int_day,ave_bias,std_bias,ave_rmse,std_rmse
  close(41)
    
end subroutine write_growth_rate
