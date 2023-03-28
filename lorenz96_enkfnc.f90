program main

  use parameter
  implicit none

  !---Common
  integer it,nt
  integer it_spin,nt_spin
  integer ix,jx
  integer iens,nens
  integer iinf
  integer iR
  integer status,system  

  !Obs.
  integer iobs,nobs
  integer it_obs,nt_obs
  integer st_obs,et_obs
  
  !IO
  integer it_io
  
  !---Variable
  !True
  double precision xtin(nx),xin(nx)
  double precision xtout(nx),xout(nx)

  !Forecast & Analysis
  double precision,allocatable :: xfin(:,:),xfout(:,:),xaout(:,:) !EnKF
  double precision xfinmean(nx),xfoutmean(nx),xaoutmean(nx)       !Ensemble mean
  double precision xfinsprd(nx),xfoutsprd(nx),xaoutsprd(nx)       !Ensemble sprd
  
  double precision,allocatable :: noise_ens(:,:)

  !Error
  double precision,allocatable :: xferr(:,:)    !Forecast error
  double precision,allocatable :: Hxferr(:,:)   !               in obs. space 
    
  !Obs.
  double precision,allocatable :: obs(:)
  double precision,allocatable :: gnoise(:,:)
  double precision,allocatable :: obserr(:,:)
  
  !---Error covariance matrix
  double precision,allocatable :: R(:,:)         !Observation error covariance matrix: R
  double precision Rout(nR)
  double precision,allocatable :: H(:,:),HT(:,:) !Observation operator
  double precision inf_parm(ninf)                !Multiplicative inflation

  integer id_inf(npco) !Optimal inflation parameter for each parameter a
  integer id_R(npco)   !        R
  integer id2_R(ninf,npco)
  
  !---Static variable
  double precision,allocatable :: bias(:),biasf(:),biasa(:),biaso(:)             !Bias
  double precision,allocatable :: rmse(:),rmsef(:),rmsea(:),rmseo(:)             !RMSE
  double precision,allocatable :: sprdf(:),sprda(:)                              !Spread
  
  double precision ave_bias(ninf,nR,npco),ave_biasf(ninf,nR,npco),ave_biasa(ninf,nR,npco),ave_biaso(ninf,nR,npco) !Averaged Bias
  double precision std_bias(ninf,nR,npco),std_biasf(ninf,nR,npco),std_biasa(ninf,nR,npco),std_biaso(ninf,nR,npco) !STD Bias 
  
  double precision ave_rmse(ninf,nR,npco),ave_rmsef(ninf,nR,npco),ave_rmsea(ninf,nR,npco),ave_rmseo(ninf,nR,npco) !Averaged RMSE
  double precision std_rmse(ninf,nR,npco),std_rmsef(ninf,nR,npco),std_rmsea(ninf,nR,npco),std_rmseo(ninf,nR,npco) !STD RMSE

  double precision ave_sprdf(ninf,nR,npco),ave_sprda(ninf,nR,npco) !Averaged SPRD
  double precision std_sprdf(ninf,nR,npco),std_sprda(ninf,nR,npco) !STD SPRD
  
  !---Noise
  double precision noise(nx)
  double precision,allocatable :: noise_pf(:,:)

  !---ETKFCC
  integer ipco
  double precision pco(npco)                              !Parameter for Correlated Observation
  double precision,allocatable :: cor(:)              !Correlation between Hxferr vs. obserr
  double precision ave_cor(ninf,nR,npco),std_cor(ninf,nR,npco)  !Average & STD Correlation
  
  !---Null
  double precision null
  double precision xnull(nx)
  double precision Pnull(nx,nx)

  !---Tmp
  double precision tmp
  double precision Ptmp(nx,nx)

  !--NaN
  integer inan
  
  !---check nx,nobs,execute_da
  if(nx <= 3)then
     write(6,'(a)') "***Error: Choose nx > 3"
     stop
  end if

  if(execute_da == 10)then
     write(6,'(a)') "Execute ETKFNC"
  else
     write(6,'(a)') "***Error: Choose execute_da=10"
     stop
  end if
  
  !---Number of time
  nt=nint(dble(nday)*int_day/dt)
  nt_spin=nint(dble(nday_spin)*int_day/dt)
  write(6,'(a,i6)') "Number of time:",nt

  !--Make directory
  status=system("mkdir -p "//trim(dir(execute_da)))

  !---Null
  xnull(:)=0.d0
  Pnull(:,:)=0.d0
  
  !---DA run: Kalman Filter or 3D-VAR  
  write(6,'(a)') "===Start DA run==="
  
  !---Observation
  !Time
  st_obs=nint((dble(yyyys_obs)*365.d0+dble(dds_obs)+dble(hhs_obs)/24.d0)*int_day/dt)
  et_obs=nint((dble(yyyye_obs)*365.d0+dble(dde_obs)+dble(hhe_obs)/24.d0)*int_day/dt)
  nt_obs=(et_obs-st_obs)/nint(int_obs/dt)+1
  write(6,'(a,i6,a,i6,a)') "Obs time:",st_obs,"-",et_obs
  write(6,'(a,i6)') "Number of obs time:",nt_obs
  allocate(bias(nt_obs),biasf(nt_obs),biasa(nt_obs),biaso(nt_obs))
  allocate(rmse(nt_obs),rmsef(nt_obs),rmsea(nt_obs),rmseo(nt_obs))
  allocate(sprdf(nt_obs),sprda(nt_obs))
  
  !---Random number
  call make_gaussian_random(0,noise_mean_da,noise_sd_da,nx,noise)
  allocate(gnoise(nx,nt_obs))
  allocate(obserr(nx,nt_obs),xferr(nx,nt_obs))
  do ix=1,nx
     call make_gaussian_random(ix,0.d0,err_obs,nt_obs,gnoise(ix,:))
  end do
  
  !---Covariance inflation
  call make_inf_nc(inf_parm)
  call make_pco(npco,pco)
  
  !---Ensemble size
  nens=40
  write(*,'(a,i6)') "Ensemble size:",nens

  allocate(xfin(nx,nens),xfout(nx,nens),xaout(nx,nens))
  allocate(noise_ens(nx,nens))
  do iens=1,nens
     call make_gaussian_random(iens,noise_mean_da,noise_sd_da,nx,noise_ens(:,iens))
  end do
  
  !_____________________________________________________________________________________________________
  !---Start do loop
  !!!--- Nobs loop ---!!!
  do nobs=nx,nx !Number of obs.
     
  allocate(obs(nobs),R(nobs,nobs),H(nobs,nx),HT(nx,nobs))
  allocate(Hxferr(nobs,nt_obs))
  allocate(cor(nobs))
  
  !!!--- Inflation loop ---!!!
  do iinf=1,ninf !Covariance inflation

  write(*,'(a,f12.5)') "Inflation:",inf_parm(iinf)
     
  do iR=1,nR

  call make_R(iR,nobs,R)
  Rout(iR)=R(1,1)
  write(*,'(a,f12.5)') "R:",Rout(iR)
  
  !!!--- Parameter for Correlated Observation loop --- !!!
  do ipco=1,npco
        
  write(*,'(a,f12.5)') "Parameter a:",pco(ipco)
  
  !!!--- Initialization ---!!!
  it_obs=0
  it_io=0
  inan=0
  !True run
  call initialization(amp_t,deg_t,xnull,xtin) !noise=xnull=0.d0
  !Free run
  call initialization(amp_da,deg_da,noise,xin)
  !Ensemble run
  do iens=1,nens
     call initialization(amp_da,deg_da,noise_ens(:,iens),xfin(:,iens))
  end do
  !Ens mean/sprd
  do ix=1,nx
     call ave_std(nens,xfin(ix,:),xfinmean(ix),xfinsprd(ix))
  end do

  !Write enseman/sprd***
  !call write_ensmean(nens,nobs,inf_parm(iinf),pco(ipco),it_io,0, &
  !     & xtin,xin,xfinmean,xnull,xnull)
  !call write_enssprd(nens,nobs,inf_parm(iinf),pco(ipco),it_io,0, &
  !     & xfinsprd,xnull)
  
  !Matrix
  !call initialization_ECM(nobs,R,Ptmp)

  !---------------------------------------------------------------------------------------------
  !   Main loop |
  !---------------------------------------------------------------------------------------------

  
  do it=1,nt

     if(mod(it,1000) == 0 .or. mod(it,nt) == 0)then
        write(6,'(5(a5,i6,a1,i6,x))') &
             & "obs: ",nobs,"/",nx,&
             & "inf: ",iinf,"/",ninf, &
             & "R:   ",iR,"/",nR, &
             & "pco: ",ipco,"/",npco, &
             & "time:",it,"/",nt
     end if
        
     !---Forwading Lorenz model
     call forward_lorenz96(xtin,xtout)
     call forward_lorenz96(xin,xout)
     do iens=1,nens
        call forward_lorenz96(xfin(:,iens),xfout(:,iens))
     end do

     !---Spin-up
     if(it == st_obs)then
        do it_spin=1,nt_spin
           
           call forward_lorenz96(xtin,xtout)
           call forward_lorenz96(xin,xout)
           do iens=1,nens
              call forward_lorenz96(xfin(:,iens),xfout(:,iens))
           end do
           
           if(mod(it_spin,nint(int_obs/dt)) == 0)then
              
              call obs_operator(it,nobs,H,HT)        
              obs(:)=matmul(H(:,:),xtout(:)+gnoise(:,it_obs))
              call etkfcc(nens,nobs,xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                   & obs,R,H,0.d0,inf_parm(iinf))
              
              xfout(:,:)=xaout(:,:)
              
           end if
           
        end do !it_spin
     end if
        
     !---Ensemble Kalman Filter
     if(mod(it,nint(int_obs/dt)) == 0)then
        
        call obs_operator(it,nobs,H,HT)
        
        if(st_obs <= it .and. it <= et_obs)then !Assimilation period
           
           it_obs=it_obs+1
           
           !Make correlated obs.
           call make_xferr(nens,xtout,xfout,xferr(:,it_obs))
           obserr(:,it_obs)=gnoise(:,it_obs)+pco(ipco)*xferr(:,it_obs)
           
           obs(:)=matmul(H(:,:),xtout(:)+obserr(:,it_obs))
           
           if(execute_da == 10)then
              call etkfcc(nens,nobs,xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                   & obs,R,H,0.d0,inf_parm(iinf))
           end if
           
           !---Bias & RMSD
           call bias_rmse(nx,xout,xtout,bias(it_obs),rmse(it_obs))
           call bias_rmse(nx,xfoutmean,xtout,biasf(it_obs),rmsef(it_obs))
           call bias_rmse(nx,xaoutmean,xtout,biasa(it_obs),rmsea(it_obs))
           call bias_rmse(nobs,obs,matmul(H(:,:),xtout(:)),biaso(it_obs),rmseo(it_obs))
           call ave_std(nx,xfoutsprd,sprdf(it_obs),tmp)
           call ave_std(nx,xaoutsprd,sprda(it_obs),tmp)
           
           !call write_ens_bias_rmse(nens,nobs,iinf,ipco,it, &
           !& bias(it_obs),biasf(it_obs),biasa(it_obs),biaso(it_obs), &
           !& rmse(it_obs),rmsef(it_obs),rmsea(it_obs),rmseo(it_obs), &
           !& inf_parm(iinf))
           
           !xfout <-- xaout
           xfout(:,:)=xaout(:,:)
           
        else !Free run period
           
           do ix=1,nx
              call ave_std(nens,xfout(ix,:),xfoutmean(ix),xfoutsprd(ix))
           end do
           
           obs(:)=0.d0
           xaoutmean(:)=0.d0
           xaoutsprd(:)=0.d0
              
        end if !st_obs <= it <= et_obs
        
        !---IO
        it_io=it_io+1
        
        !call write_ensmean(nens,nobs,inf_parm(iinf),pco(ipco),it_io,it, &
        !     & xtout,xout,xfoutmean,xaoutmean,obserr)
        !call write_enssprd(nens,nobs,inf_parm(iinf),pco(ipco),it_io,it, &
        !     & xfoutsprd,xaoutsprd)
        
     end if !it_obs
     
     !xfin <-- xfout
     xtin(:)=xtout(:)
     xin(:)=xout(:)
     xfin(:,:)=xfout(:,:)
     
  end do !it

  !---Bias & RMSD
  call ave_std(nt_obs,bias,ave_bias(iinf,iR,ipco),std_bias(iinf,iR,ipco))
  call ave_std(nt_obs,biasf,ave_biasf(iinf,iR,ipco),std_biasf(iinf,iR,ipco))
  call ave_std(nt_obs,biasa,ave_biasa(iinf,iR,ipco),std_biasa(iinf,iR,ipco))
  call ave_std(nt_obs,biaso,ave_biaso(iinf,iR,ipco),std_biaso(iinf,iR,ipco))     
  call ave_std_rmse(nt_obs,rmse,ave_rmse(iinf,iR,ipco),std_rmse(iinf,iR,ipco))
  call ave_std_rmse(nt_obs,rmsef,ave_rmsef(iinf,iR,ipco),std_rmsef(iinf,iR,ipco))
  call ave_std_rmse(nt_obs,rmsea,ave_rmsea(iinf,iR,ipco),std_rmsea(iinf,iR,ipco))
  call ave_std_rmse(nt_obs,rmseo,ave_rmseo(iinf,iR,ipco),std_rmseo(iinf,iR,ipco))
  call ave_std_rmse(nt_obs,sprdf,ave_sprdf(iinf,iR,ipco),std_sprdf(iinf,iR,ipco))
  call ave_std_rmse(nt_obs,sprda,ave_sprda(iinf,iR,ipco),std_sprda(iinf,iR,ipco))
     
  write(6,'(a,2(f12.5,a))') "Free run: Bias(",ave_bias(iinf,iR,ipco),") RMSE(",ave_rmse(iinf,iR,ipco),")"
  write(6,'(a,3(f12.5,a))') "Forecast: Bias(",ave_biasf(iinf,iR,ipco),") RMSE(",ave_rmsef(iinf,iR,ipco), &
       & ") SPREAD(",ave_sprdf(iinf,iR,ipco),")"
  write(6,'(a,3(f12.5,a))') "Analysis: Bias(",ave_biasa(iinf,iR,ipco),") RMSE(",ave_rmsea(iinf,iR,ipco), &
       & ") SPREAD(",ave_sprda(iinf,iR,ipco),")"
  write(6,'(a,2(f12.5,a))') "Observation: Bias(",ave_biaso(iinf,iR,ipco),") RMSE(",ave_rmseo(iinf,iR,ipco),")"

  !---Correlation
  Hxferr(:,:)=matmul(H(:,:),xferr(:,:))
  do iobs=1,nobs
     call correlation(nt_obs,Hxferr(iobs,:),obserr(iobs,:),cor(iobs))
  end do
  call ave_std(nobs,cor(:),ave_cor(iinf,iR,ipco),std_cor(iinf,iR,ipco))
     
  end do !ipco
  end do !iR
  end do !iinf

  !---Minimum forecast RMSE and statistic values for each a
  call est_id_inf_R(ave_rmsef,id_inf,id_R)
  call est_id2_R(ave_rmsef,id2_R)

  !---IO
  call write_static_ave_nc(inf_parm,Rout,pco,id_inf,id_R,id2_R, &
       & ave_bias,ave_biasf,ave_biasa,ave_biaso, &
       & ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo, &
       & ave_sprdf,ave_sprda, &
       & ave_cor)
  
  deallocate(obs,R,H,HT)
  deallocate(Hxferr)
  deallocate(cor)
  
  end do !nobs

  deallocate(xfin,xfout,xaout)
  deallocate(noise_ens)
  
  write(6,'(a)') "===End DA run==="
    
  deallocate(bias,biasf,biasa,biaso)
  deallocate(rmse,rmsef,rmsea,rmseo)
  deallocate(sprdf,sprda)
  deallocate(gnoise)
  deallocate(obserr)

  write(6,'(a)') "=== All End === "

end program main
