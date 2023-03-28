program main

  use parameter
  implicit none

  !---Common
  integer it,nt
  integer ix
  integer iens,nens
  integer iens_exp
  integer iinf
  integer status,system  
  !Obs.
  integer iobs,nobs
  integer it_obs,nt_obs
  integer st_obs,et_obs

  !IO
  integer it_io

  integer sigma_loc !Localization scale [grid]
  integer nsigma_loc

  !Ensmeble size
  integer ens_size(nens_exp)
  
  !---Variable
  double precision xtin(nx),xin(nx)
  double precision xtout(nx),xout(nx)
  double precision,allocatable :: xfin(:,:),xfout(:,:),xaout(:,:)
  double precision xfinmean(nx),xfoutmean(nx),xaoutmean(nx) !Ensemble mean
  double precision xfinsprd(nx),xfoutsprd(nx),xaoutsprd(nx) !Ensemble sprd
  double precision,allocatable :: noise_ens(:,:)
  
  double precision,allocatable :: obs(:)
  double precision,allocatable :: obserr(:,:)
  double precision,allocatable :: obserr_ens(:,:,:)
  
  !---Error covariance matrix
  double precision,allocatable :: R(:,:)         !Observation error covarinance matrix: R
  double precision,allocatable :: H(:,:),HT(:,:) !Observation operator
  double precision inf_parm(ninf)              !Multiplicative inflation
  
  !---Static variable
  double precision,allocatable :: bias(:),biasf(:),biasa(:),biaso(:) !Bias
  double precision,allocatable :: rmse(:),rmsef(:),rmsea(:),rmseo(:) !RMSE
  double precision ave_bias,ave_biasf,ave_biasa,ave_biaso
  double precision std_bias,std_biasf,std_biasa,std_biaso
  double precision ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo
  double precision std_rmse,std_rmsef,std_rmsea,std_rmseo
  
  !---Noise
  double precision noise(nx)
  double precision,allocatable :: noise_pf(:,:)

  !---Particle filter
  double precision,allocatable :: Neff(:) !Effective ensemble size
  double precision ave_Neff,std_Neff
  
  !---Null
  double precision xnull(nx)
  double precision Pnull(nx,nx)

  !---Tmp
  double precision Ptmp(nx,nx)
  
  !---check nx,nobs,execute_da
  if(nx <= 3)then
     write(6,'(a)') "***Error: Choose nx > 3"
     stop
  end if
  if(execute_da == 3)then
     write(6,'(a)') "Execute ETKF"
  else if(execute_da == 4)then
     write(6,'(a)') "Execute LETKF"
  else if(execute_da == 5)then
     write(6,'(a)') "Execute EnSRF"
  else if(execute_da == 6)then
     write(6,'(a)') "Execute PO EnKF"
  else if(execute_da == 7)then
     write(6,'(a)') "Execute PF"
  else if(execute_da == 8)then
     write(6,'(a)') "Execute LPF"     
  else
     write(6,'(a)') "***Error: Choose execute_da=3-8"
     stop
  end if
  
  !---Number of time
  nt=nint(dble(nday)*int_day/dt)
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
  
  !---Random number
  call make_gaussian_random(0,noise_mean_da,noise_sd_da,nx,noise)
  allocate(obserr(nx,nt_obs))
  do ix=1,nx
     call make_gaussian_random(ix,0.d0,err_obs,nt_obs,obserr(ix,:))
  end do
  
  !---Covariance inflation
  call make_inf(inf_parm)

  !---Ensemble size
  call make_ensemble_size(nens_exp,ens_size)

  !_____________________________________________________________________________________________________
  
  !---Start do loop
  !!!Ensemble size loop!!!
!  do iens_exp=1,nens_exp !Number of ensemble
  do iens_exp=nens_exp,nens_exp !Number of ensemble

  nens=ens_size(iens_exp)
  write(*,'(a,i6)') "Ensemble size:",nens
     
  allocate(xfin(nx,nens),xfout(nx,nens),xaout(nx,nens))
  allocate(noise_ens(nx,nens))
  do iens=1,nens
     call make_gaussian_random(iens,noise_mean_da,noise_sd_da,nx,noise_ens(:,iens))
  end do
     
  if(execute_da == 6)then
     allocate(obserr_ens(nx,nens,nt_obs))
     do ix=1,nx
        do iens=1,nens
           call make_gaussian_random(nx+nens+ix+iens,0.d0,err_obs,nt_obs,obserr_ens(ix,iens,:))
        end do
     end do
  end if
        
  !Particle filter
  if(execute_da == 7 .or. execute_da == 8)then
     allocate(noise_pf(nx,nens))
     allocate(Neff(nt_obs))
  end if

  !!!Nobs loop!!!
  do nobs=nx,nx !Number of obs.
     
  allocate(obs(nobs),R(nobs,nobs),H(nobs,nx),HT(nx,nobs))

  !!!Inflation loop!!!
  do iinf=1,ninf !Covariance inflation
  !           do iinf=ninf,ninf !Covariance inflation
  !           inf_parm(iinf)=10.d-3
     
  write(*,'(a,E15.4)') "Inflation:",inf_parm(iinf)
     
  if(execute_da == 3 .or. execute_da == 7)then
     nsigma_loc=1
  else
     nsigma_loc=10
  end if

  !!! Localization loop"""
  !do sigma_loc=1,nsigma_loc !Localization scale
  do sigma_loc=nsigma_loc,nsigma_loc !Localization scale
     
  !---------------------------------------------------------------------------------------------
  !---Main loop
  !----------------------------------------------------------------------------------------------

  !---Initialization
  it_obs=0
  it_io=0
  !True run
  call initialization(amp_t,deg_t,xnull,xtin) !noise=xnull=0.d0
  !Free run
  call initialization(amp_da,deg_da,noise,xin)
  !Ensemble run
  do iens=1,nens
     call initialization(amp_da,deg_da,noise_ens(:,iens),xfin(:,iens))
  end do
  do ix=1,nx
     call ave_std(nens,xfin(ix,:),xfinmean(ix),xfinsprd(ix))
  end do
  
  call write_ensmean(nens,nobs,iinf,sigma_loc, &
       & it_io,0,xtin,xin,xfinmean,xnull,xnull,inf_parm(iinf))
  call write_enssprd(nens,nobs,iinf,sigma_loc, &
       & it_io,0,xfinsprd,xnull,inf_parm(iinf))
  !Matrix
  call initialization_ECM(nobs,R,Ptmp)

  do it=1,nt

     if(mod(it,1000) == 0 .or. mod(it,nt) == 0) &
          & write(6,*) &
          & "ens:",iens_exp,"/",nens_exp,"obs:",nobs,"/",nx,&
          & "inf:",iinf,"/",ninf,"loc:",sigma_loc,"/",nsigma_loc,"time:",it,"/",nt
        
        !---Forwading Lorenz model
        call forward_lorenz96(xtin,xtout)
        call forward_lorenz96(xin,xout)
        do iens=1,nens
           call forward_lorenz96(xfin(:,iens),xfout(:,iens))
        end do
        
        !---Ensemble Kalman Filter
        if(mod(it,nint(int_obs/dt)) == 0)then
           
           call obs_operator(it,nobs,H,HT)
           
           if(st_obs <= it .and. it <= et_obs)then
              
              it_obs=it_obs+1
              obs(:)=matmul(H(:,:),xtout(:))+matmul(H(:,:),obserr(:,it_obs))
              
              if(execute_da == 3)then
                 call etkf(nens,nobs,xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                      & obs,R,H,inf_parm(iinf))
              else if(execute_da == 4)then
                 call letkf(nens,nobs,sigma_loc, &
                      & xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                      & obs,R,H,inf_parm(iinf))
              else if(execute_da == 5)then
                 call ensrf(nens,nobs,sigma_loc, &
                      & xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd,obs,R,H,inf_parm(iinf))
              else if(execute_da == 6)then
                 call poenkf(nens,nobs,sigma_loc, &
                      & xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                      & obs,obserr_ens(:,:,it_obs),R,H,inf_parm(iinf))
              else if(execute_da == 7)then
                 call pf(it,nens,nobs, &
                      & xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                      & obs,R,H,Neff(it_obs))
              else if(execute_da == 8)then
                 call lpf(it,nens,nobs,sigma_loc, &
                      & xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                      & obs,R,H,Neff(it_obs))
              end if

              !---Bias & RMSD
              call bias_rmse(nx,xout,xtout,bias(it_obs),rmse(it_obs))
              call bias_rmse(nx,xfoutmean,xtout,biasf(it_obs),rmsef(it_obs))
              call bias_rmse(nx,xaoutmean,xtout,biasa(it_obs),rmsea(it_obs))
              call bias_rmse(nobs,obs,matmul(H(:,:),xtout(:)),biaso(it_obs),rmseo(it_obs))
              
              call write_ens_bias_rmse(nens,nobs,iinf,sigma_loc,st_obs,it, &
                   & bias(it_obs),biasf(it_obs),biasa(it_obs),biaso(it_obs), &
                   & rmse(it_obs),rmsef(it_obs),rmsea(it_obs),rmseo(it_obs),inf_parm(iinf))
              
              xfout(:,:)=xaout(:,:)
              
              !---Post-regulation
              if((execute_da == 7 .or. execute_da == 8) .and. add_inf)then
                 do ix=1,nx
                    call make_gaussian_random(ix+nens,0.d0,sqrt(inf_parm(iinf)),nens,noise_pf(ix,:))
                 end do
                 xfout(:,:)=xfout(:,:)+noise_pf(:,:)
              end if
              
           else
              
              do ix=1,nx
                 call ave_std(nens,xfin(ix,:),xfoutmean(ix),xfoutsprd(ix))
              end do
              
              obs(:)=0.d0
              xaoutmean(:)=0.d0
              xaoutsprd(:)=0.d0
              
           end if !st_obs <= it <= et_obs

           !---IO
           it_io=it_io+1
           call write_ensmean(nens,nobs,iinf,sigma_loc, &
                & it_io,it,xtout,xout,xfoutmean,xaoutmean,matmul(HT(:,:),obs(:)),inf_parm(iinf))
           call write_enssprd(nens,nobs,iinf,sigma_loc, &
                & it_io,it,xfoutsprd,xaoutsprd,inf_parm(iinf))
           
        end if !it_obs
        
        xtin(:)=xtout(:)
        xin(:)=xout(:)
        xfin(:,:)=xfout(:,:)
        
     end do !it

     !---Bias & RMSD
     call ave_std(nt_obs,bias,ave_bias,std_bias)
     call ave_std(nt_obs,biasf,ave_biasf,std_biasf)
     call ave_std(nt_obs,biasa,ave_biasa,std_biasa)
     call ave_std(nt_obs,biaso,ave_biaso,std_biaso)
     call ave_std(nt_obs,rmse,ave_rmse,std_rmse)
     call ave_std(nt_obs,rmsef,ave_rmsef,std_rmsef)
     call ave_std(nt_obs,rmsea,ave_rmsea,std_rmsea)
     call ave_std(nt_obs,rmseo,ave_rmseo,std_rmseo)

     write(6,'(a)') "Bias(Ave.,Std.) RMSE(Ave.,Std.)"
     write(6,'(a,4(f12.5,a))') "Free run: Bias(",ave_bias,",",std_bias,") RMSE(",ave_rmse,",",std_rmse,")"
     write(6,'(a,4(f12.5,a))') "Forecast: Bias(",ave_biasf,",",std_biasf,") RMSE(",ave_rmsef,",",std_rmsef,")"
     write(6,'(a,4(f12.5,a))') "Analysis: Bias(",ave_biasa,",",std_biasa,") RMSE(",ave_rmsea,",",std_rmsea,")"
     write(6,'(a,4(f12.5,a))') "Observation: Bias(",ave_biaso,",",std_biaso,") RMSE(",ave_rmseo,",",std_rmseo,")"


     !---IO
     call write_ens_ave_bias_rmse(nens,nobs,iinf,sigma_loc,ave_bias,ave_biasf,ave_biasa,ave_biaso, &
          & ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo,inf_parm(iinf))
     
     !--- Effective size for PF
     if(execute_da == 7 .or. execute_da == 8)then
        call ave_std(nt_obs,Neff,ave_Neff,std_Neff)
        call write_Neff(nens,nobs,iinf,sigma_loc,ave_Neff,std_Neff,inf_parm(iinf))
     end if
     
  end do !sigma_loc
  end do !iinf

  deallocate(obs,R,H,HT)

  end do !nobs

  deallocate(xfin,xfout,xaout)
  deallocate(noise_ens)

  if(execute_da == 6)then
     deallocate(obserr_ens)
  end if

  if(execute_da == 7 .or. execute_da == 8)then
     deallocate(noise_pf)
     deallocate(Neff)
  end if
     
  end do !nens
  
  write(6,'(a)') "===End DA run==="
  
  !  status=system("csh fig.csh "//trim(dir(execute_da)))
  
  deallocate(bias,biasf,biasa,biaso)
  deallocate(rmse,rmsef,rmsea,rmseo)
  deallocate(obserr)

  write(6,'(a)') "=== All End === "

end program main
