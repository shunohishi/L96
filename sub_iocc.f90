subroutine write_ensmean(nens,nobs,inf_parm,pco, &
     & it_io,it, &
     & xt,x,xf,xa,obserr)

  use parameter
  implicit none

  integer ix
  
  integer,intent(in) :: nens,nobs
  integer,intent(in) :: it_io,it
  
  double precision,intent(in) :: xt(nx),x(nx)
  double precision,intent(in) :: xf(nx),xa(nx)
  double precision,intent(in) :: obserr(nx)
  double precision,intent(in) :: inf_parm,pco
    
  open(1,file=trim(dir(execute_da))//"/x.dat",access="append")
  do ix=1,nx
     write(1,'(2i6,2f12.5,i6,6f12.5)') nens,nobs,inf_parm,pco,ix,dt*dble(it)/int_day, &
             & xt(ix),x(ix),xf(ix),xa(ix),obserr(ix)
  end do
  close(1)

end subroutine write_ensmean

!------------------------------

subroutine write_enssprd(nens,nobs,inf_parm,pco, &
     & it_io,it, &
     & xfsprd,xasprd)

  use parameter
  implicit none

  integer,intent(in) :: nens,nobs
  integer,intent(in) :: it_io,it
  integer ix

  double precision,intent(in) :: xfsprd(nx),xasprd(nx)
  double precision xfsprdmean,xasprdmean
  double precision,intent(in) :: inf_parm,pco
  
  xfsprdmean=0.d0
  xasprdmean=0.d0

  do ix=1,nx
     xfsprdmean=xfsprdmean+xfsprd(ix)/dble(nx)
     xasprdmean=xasprdmean+xasprd(ix)/dble(nx)
  end do
  
  open(1,file=trim(dir(execute_da))//"/sprd.dat",access="append")
  open(2,file=trim(dir(execute_da))//"/sprdmean.dat",access="append")
  
  do ix=1,nx
     write(1,'(2i6,2f12.5,i6,3f12.5)') nens,nobs,inf_parm,pco,ix,dt*dble(it)/int_day,xfsprd(ix),xasprd(ix)
  end do
  write(2,'(2i6,2f12.5,3f12.5)') nens,nobs,inf_parm,pco,dt*dble(it)/int_day,xfsprdmean,xasprdmean
  
  close(1)
  close(2)
  
end subroutine write_enssprd

!----------------------------------
  
subroutine write_static_ave(inf_parm,pco,id_inf, &
     & bias,biasf,biasa,biaso, &
     & rmse,rmsef,rmsea,rmseo, &
     & sprdf,sprda,cor)

  use parameter
  implicit none

  integer ipco,iinf
  integer id
  
  integer,intent(in) :: id_inf(npco)
  double precision,intent(in) :: inf_parm(ninf),pco(npco)
  double precision,intent(in) :: bias(ninf,npco),biasf(ninf,npco),biasa(ninf,npco),biaso(ninf,npco)
  double precision,intent(in) :: rmse(ninf,npco),rmsef(ninf,npco),rmsea(ninf,npco),rmseo(ninf,npco)
  double precision,intent(in) :: sprdf(ninf,npco),sprda(ninf,npco)
  double precision,intent(in) :: cor(ninf,npco)
  
  open(11,file=trim(dir(execute_da))//"/bias_ave.dat",status="replace")
  open(12,file=trim(dir(execute_da))//"/rmse_ave.dat",status="replace")
  open(13,file=trim(dir(execute_da))//"/sprd_ave.dat",status="replace")
  open(14,file=trim(dir(execute_da))//"/cor_ave.dat",status="replace")
  open(21,file=trim(dir(execute_da))//"/bias_ave_min.dat",status="replace")
  open(22,file=trim(dir(execute_da))//"/rmse_ave_min.dat",status="replace")
  open(23,file=trim(dir(execute_da))//"/sprd_ave_min.dat",status="replace")
  open(24,file=trim(dir(execute_da))//"/cor_ave_min.dat",status="replace")
  do ipco=1,npco

     do iinf=1,ninf
        write(11,'(2f12.5,4f12.5)') inf_parm(iinf),pco(ipco), &
             & bias(iinf,ipco),biasf(iinf,ipco),biasa(iinf,ipco),biaso(iinf,ipco)
        write(12,'(2f12.5,4f12.5)') inf_parm(iinf),pco(ipco), &
             & rmse(iinf,ipco),rmsef(iinf,ipco),rmsea(iinf,ipco),rmseo(iinf,ipco)
        write(13,'(2f12.5,2f12.5)') inf_parm(iinf),pco(ipco), &
             & sprdf(iinf,ipco),sprda(iinf,ipco)
        write(14,'(2f12.5,f12.5)') inf_parm(iinf),pco(ipco),cor(iinf,ipco)
     end do

     id=id_inf(ipco)
     write(21,'(2f12.5,4f12.5)') inf_parm(id),pco(ipco), &
          & bias(id,ipco),biasf(id,ipco),biasa(id,ipco),biaso(id,ipco)
     write(22,'(2f12.5,4f12.5)') inf_parm(id),pco(ipco), &
          & rmse(id,ipco),rmsef(id,ipco),rmsea(id,ipco),rmseo(id,ipco)
     write(23,'(2f12.5,2f12.5)') inf_parm(id),pco(ipco), &
          & sprdf(id,ipco),sprda(id,ipco)
     write(24,'(2f12.5,f12.5)') inf_parm(id),pco(ipco),cor(id,ipco)

  end do
  close(11)
  close(12)
  close(13)
  close(14)
  close(21)
  close(22)
  close(23)
  close(24)

end subroutine write_static_ave

!-------------------------------------------------------------------------

subroutine write_static_ave_nc(inf_parm,R,pco,id_inf,id_R,id2_R, &
     & bias,biasf,biasa,biaso, &
     & rmse,rmsef,rmsea,rmseo, &
     & sprdf,sprda,cor)

  use parameter
  implicit none

  integer ipco,iR,iinf
  
  integer,intent(in) :: id_inf(npco),id_R(npco),id2_R(ninf,npco)
  double precision,intent(in) :: inf_parm(ninf),R(nR),pco(npco)
  double precision,intent(in) :: bias(ninf,nR,npco),biasf(ninf,nR,npco),biasa(ninf,nR,npco),biaso(ninf,nR,npco)
  double precision,intent(in) :: rmse(ninf,nR,npco),rmsef(ninf,nR,npco),rmsea(ninf,nR,npco),rmseo(ninf,nR,npco)
  double precision,intent(in) :: sprdf(ninf,nR,npco),sprda(ninf,nR,npco)
  double precision,intent(in) :: cor(ninf,nR,npco)
  
  open(11,file=trim(dir(execute_da))//"/bias_ave.dat",status="replace")
  open(12,file=trim(dir(execute_da))//"/rmse_ave.dat",status="replace")
  open(13,file=trim(dir(execute_da))//"/sprd_ave.dat",status="replace")
  open(14,file=trim(dir(execute_da))//"/cor_ave.dat",status="replace")
  open(21,file=trim(dir(execute_da))//"/bias_ave_min.dat",status="replace")
  open(22,file=trim(dir(execute_da))//"/rmse_ave_min.dat",status="replace")
  open(23,file=trim(dir(execute_da))//"/sprd_ave_min.dat",status="replace")
  open(24,file=trim(dir(execute_da))//"/cor_ave_min.dat",status="replace")
  open(31,file=trim(dir(execute_da))//"/bias_ave_Rmin.dat",status="replace")
  open(32,file=trim(dir(execute_da))//"/rmse_ave_Rmin.dat",status="replace")
  open(33,file=trim(dir(execute_da))//"/sprd_ave_Rmin.dat",status="replace")
  open(34,file=trim(dir(execute_da))//"/cor_ave_Rmin.dat",status="replace")

  do ipco=1,npco

     do iR=1,nR
        do iinf=1,ninf
           write(11,'(3f12.5,4f12.5)') inf_parm(iinf),R(iR),pco(ipco), &
                & bias(iinf,iR,ipco),biasf(iinf,iR,ipco),biasa(iinf,iR,ipco),biaso(iinf,iR,ipco)
           write(12,'(3f12.5,4f12.5)') inf_parm(iinf),R(iR),pco(ipco), &
                & rmse(iinf,iR,ipco),rmsef(iinf,iR,ipco),rmsea(iinf,iR,ipco),rmseo(iinf,iR,ipco)
           write(13,'(3f12.5,2f12.5)') inf_parm(iinf),R(iR),pco(ipco), &
                & sprdf(iinf,iR,ipco),sprda(iinf,iR,ipco)
           write(14,'(3f12.5,f12.5)') inf_parm(iinf),R(iR),pco(ipco),cor(iinf,iR,ipco)
        end do
     end do
     
     iinf=id_inf(ipco)
     iR=id_R(ipco)
     
     write(21,'(3f12.5,4f12.5)') inf_parm(iinf),R(iR),pco(ipco), &
          & bias(iinf,iR,ipco),biasf(iinf,iR,ipco),biasa(iinf,iR,ipco),biaso(iinf,iR,ipco)
     write(22,'(3f12.5,4f12.5)') inf_parm(iinf),R(iR),pco(ipco), &
          & rmse(iinf,iR,ipco),rmsef(iinf,iR,ipco),rmsea(iinf,iR,ipco),rmseo(iinf,iR,ipco)
     write(23,'(3f12.5,2f12.5)') inf_parm(iinf),R(iR),pco(ipco), &
          & sprdf(iinf,iR,ipco),sprda(iinf,iR,ipco)
     write(24,'(3f12.5,f12.5)') inf_parm(iinf),R(iR),pco(ipco),cor(iinf,iR,ipco)

     do iinf=1,ninf

        iR=id2_R(iinf,ipco)
        write(31,'(3f12.5,4f12.5)') inf_parm(iinf),R(iR),pco(ipco), &
             & bias(iinf,iR,ipco),biasf(iinf,iR,ipco),biasa(iinf,iR,ipco),biaso(iinf,iR,ipco)
        write(32,'(3f12.5,4f12.5)') inf_parm(iinf),R(iR),pco(ipco), &
             & rmse(iinf,iR,ipco),rmsef(iinf,iR,ipco),rmsea(iinf,iR,ipco),rmseo(iinf,iR,ipco)
        write(33,'(3f12.5,2f12.5)') inf_parm(iinf),R(iR),pco(ipco), &
             & sprdf(iinf,iR,ipco),sprda(iinf,iR,ipco)
        write(34,'(3f12.5,f12.5)') inf_parm(iinf),R(iR),pco(ipco),cor(iinf,iR,ipco)

     end do

     
  end do
  close(11)
  close(12)
  close(13)
  close(14)
  close(21)
  close(22)
  close(23)
  close(24)
  close(31)
  close(32)
  close(33)
  close(34)

end subroutine write_static_ave_nc
