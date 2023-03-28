program main

  use parameter
  implicit none

  !---Common
  integer it,nt
  integer ix,jx
  integer iens,nens
  integer iens_exp
  integer iinf
  integer status,system  
  !Obs.
  integer iobs,nobs
  integer it_obs,nt_obs
  integer st_obs,et_obs
  integer ifo,nfo
  integer sfo,efo
  
  !IO
  integer it_io

  integer sigma_loc !Localization scale [grid]
  integer nsigma_loc

  !Ensmeble size
  integer ens_size(nens_exp)
  
  !---Variable
  !True
  double precision xtin(nx),xin(nx)
  double precision xtout(nx),xout(nx)

  !Forecast & Analysis
  double precision,allocatable :: xfin(:,:),xfout(:,:),xaout(:,:) !EnKF
  double precision xfinmean(nx),xfoutmean(nx),xaoutmean(nx)       !Ensemble mean
  double precision xfinsprd(nx),xfoutsprd(nx),xaoutsprd(nx)       !Ensemble sprd
  
  double precision,allocatable :: xfin_c(:,:),xfout_c(:,:),xaout_c(:,:)    !EnKFCC
  double precision xfinmean_c(nx),xfoutmean_c(nx),xaoutmean_c(nx)
  double precision xfinsprd_c(nx),xfoutsprd_c(nx),xaoutsprd_c(nx)

  double precision,allocatable :: xfin_nc(:,:),xfout_nc(:,:),xaout_nc(:,:) !EnKF
  double precision xfinmean_nc(nx),xfoutmean_nc(nx),xaoutmean_nc(nx)
  double precision xfinsprd_nc(nx),xfoutsprd_nc(nx),xaoutsprd_nc(nx)

  double precision,allocatable :: noise_ens(:,:)

  !Error
  double precision,allocatable :: xferr(:,:)    !Forecast error
  double precision,allocatable :: xferr_c(:,:) 
  double precision,allocatable :: xferr_nc(:,:)
    
  !Obs.
  double precision,allocatable :: obs(:)
  double precision,allocatable :: obserr(:,:)
  double precision,allocatable :: obserr_ens(:,:,:) !Perturbed Observation

  double precision,allocatable :: obs_c(:)      !Correlated obs. for EnKFCC
  double precision,allocatable :: obserr_c(:,:) !Correlated obs. error

  double precision,allocatable :: obs_nc(:)      !Correlated obs. for EnKF
  double precision,allocatable :: obserr_nc(:,:) !Correlated obs. error

  !Correlation coefficient
  double precision cor_set
  double precision cor(nx,nx)
  double precision cor_c(nx,nx)
  double precision cor_nc(nx,nx)
  
  !---Error covariance matrix
  double precision,allocatable :: R(:,:)         !Observation error covariance matrix: R
  double precision,allocatable :: H(:,:),HT(:,:) !Observation operator
  double precision inf_parm(ninf)              !Multiplicative inflation
  
  !---Static variable
  double precision,allocatable :: bias(:),biasf(:),biasa(:),biaso(:)             !Bias
  double precision,allocatable :: bias_c(:),biasf_c(:),biasa_c(:),biaso_c(:)     !      for EnKFCC
  double precision,allocatable :: bias_nc(:),biasf_nc(:),biasa_nc(:),biaso_nc(:) !      for EnKF
  double precision,allocatable :: rmse(:),rmsef(:),rmsea(:),rmseo(:)             !RMSE
  double precision,allocatable :: rmse_c(:),rmsef_c(:),rmsea_c(:),rmseo_c(:)     !     for EnKFCC
  double precision,allocatable :: rmse_nc(:),rmsef_nc(:),rmsea_nc(:),rmseo_nc(:) !     for EnKF
  double precision,allocatable :: sprdf(:),sprda(:)                              !Spread
  double precision,allocatable :: sprdf_c(:),sprda_c(:)                          !     for EnKFCC
  double precision,allocatable :: sprdf_nc(:),sprda_nc(:)                        !     for EnKF
  
  double precision ave_bias,ave_biasf,ave_biasa,ave_biaso             !Averaged Bias
  double precision std_bias,std_biasf,std_biasa,std_biaso             !STD Bias 
  double precision ave_bias_c,ave_biasf_c,ave_biasa_c,ave_biaso_c     !              for EnKFCC
  double precision std_bias_c,std_biasf_c,std_biasa_c,std_biaso_c
  double precision ave_bias_nc,ave_biasf_nc,ave_biasa_nc,ave_biaso_nc !              for EnKF
  double precision std_bias_nc,std_biasf_nc,std_biasa_nc,std_biaso_nc

  double precision ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo             !Averaged RMSE
  double precision std_rmse,std_rmsef,std_rmsea,std_rmseo             !STD RMSE
  double precision ave_rmse_c,ave_rmsef_c,ave_rmsea_c,ave_rmseo_c     !              for EnKFCC
  double precision std_rmse_c,std_rmsef_c,std_rmsea_c,std_rmseo_c
  double precision ave_rmse_nc,ave_rmsef_nc,ave_rmsea_nc,ave_rmseo_nc !              for EnKF
  double precision std_rmse_nc,std_rmsef_nc,std_rmsea_nc,std_rmseo_nc

  double precision ave_sprdf,ave_sprda       !Averaged SPRD
  double precision std_sprdf,std_sprda       !STD SPRD
  double precision ave_sprdf_c,ave_sprda_c   !              for EnKFCC
  double precision std_sprdf_c,std_sprda_c
  double precision ave_sprdf_nc,ave_sprda_nc !              for EnKF
  double precision std_sprdf_nc,std_sprda_nc
  
  !---Noise
  double precision noise(nx)
  double precision,allocatable :: noise_pf(:,:)

  !---Particle filter
  double precision,allocatable :: Neff(:) !Effective ensemble size
  double precision ave_Neff,std_Neff
  
  !---Null
  double precision null
  double precision xnull(nx)
  double precision Pnull(nx,nx)

  !---Tmp
  double precision tmp
  double precision Ptmp(nx,nx)

  !--NaN
  integer inan,inan_c,inan_nc
  
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
  allocate(sprdf(nt_obs),sprda(nt_obs))
  if(fo_cor)then
     allocate(biasf_c(nt_obs),biasa_c(nt_obs),biaso_c(nt_obs))
     allocate(rmsef_c(nt_obs),rmsea_c(nt_obs),rmseo_c(nt_obs))
     allocate(sprdf_c(nt_obs),sprda_c(nt_obs))
     allocate(biasf_nc(nt_obs),biasa_nc(nt_obs),biaso_nc(nt_obs))
     allocate(rmsef_nc(nt_obs),rmsea_nc(nt_obs),rmseo_nc(nt_obs))
     allocate(sprdf_nc(nt_obs),sprda_nc(nt_obs))
  end if
  
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
  do iens_exp=1,nens_exp !Number of ensemble
!  do iens_exp=nens_exp,nens_exp !Number of ensemble

  nens=ens_size(iens_exp)
  write(*,'(a,i6)') "Ensemble size:",nens
     
  allocate(xfin(nx,nens),xfout(nx,nens),xaout(nx,nens))
  if(fo_cor)then
     allocate(xfin_c(nx,nens),xfout_c(nx,nens),xaout_c(nx,nens))
     allocate(xfin_nc(nx,nens),xfout_nc(nx,nens),xaout_nc(nx,nens))
  end if
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
  if(execute_da == 4 .or. execute_da == 7 .or. execute_da == 8)then
     allocate(noise_pf(nx,nens))
     allocate(Neff(nt_obs))
  end if

  !!!Nobs loop!!!
  do nobs=nx,nx !Number of obs.
     
  allocate(obs(nobs),R(nobs,nobs),H(nobs,nx),HT(nx,nobs))
  if(fo_cor)then
     allocate(xferr(nx,nt_obs),xferr_c(nx,nt_obs),xferr_nc(nx,nt_obs))
     allocate(obs_c(nobs),obserr_c(nx,nt_obs))
     allocate(obs_nc(nobs),obserr_nc(nx,nt_obs))
  end if
     
  !!!Inflation loop!!!
  do iinf=1,ninf !Covariance inflation
!  do iinf=5,5 !Covariance inflation
     
  write(*,'(a,f12.5)') "Inflation:",inf_parm(iinf)
  
  !!! Localization loop!!!
    if(execute_da == 3 .or. execute_da == 7)then !ETKF or PF
     nsigma_loc=1
  else
     nsigma_loc=10
  end if
!  do sigma_loc=1,nsigma_loc !Localization scale
  do sigma_loc=1,1 !Localization scale

  !!! Correlation forecast vs. obs. loop"""
  if(fo_cor)then
     sfo=-9
     efo=9
  else
     sfo=0
     efo=0
  endif
  nfo=(efo-sfo)+1
  do ifo=sfo,efo

  cor_set=0.1d0*dble(ifo)
  write(*,'(a,f12.5)') "Correlation:",cor_set

  !---------------------------------------------------------------------------------------------
  !---Main loop
  !----------------------------------------------------------------------------------------------

  !---Initialization
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
  do ix=1,nx
     call ave_std(nens,xfin(ix,:),xfinmean(ix),xfinsprd(ix))
  end do

  if(fo_cor)then
     xfin_c(:,:)=xfin(:,:)
     xfinmean_c(:)=xfinmean(:)
     xfinsprd_c(:)=xfinsprd(:)
     xfin_nc(:,:)=xfin(:,:)
     xfinmean_nc(:)=xfinmean(:)
     xfinsprd_nc(:)=xfinsprd(:)
     inan_c=0
     inan_nc=0
  end if
  
  call write_ensmean(nens,nobs,iinf,sigma_loc,cor_set,it_io,0, &
       & xtin,xin,xfinmean,xnull,xfinmean_c,xnull,xfinmean_nc,xnull,xnull,xnull,xnull,inf_parm(iinf))
  call write_enssprd(nens,nobs,iinf,sigma_loc,cor_set,it_io,0, &
       & xfinsprd,xnull,xfinsprd_c,xnull,xfinsprd_nc,xnull,inf_parm(iinf))
  
  !Matrix
  call initialization_ECM(nobs,R,Ptmp)

  do it=1,nt

     if(mod(it,1000) == 0 .or. mod(it,nt) == 0) &
          & write(6,'(6(a5,i6,a1,i6))') &
          & "ens:",iens_exp,"/",nens_exp, &
          & "obs:",nobs,"/",nx,&
          & "inf:",iinf,"/",ninf, &
          & "loc:",sigma_loc,"/",nsigma_loc, &
          & "cor:",ifo-sfo+1,"/",nfo, &
          & "time:",it,"/",nt
     
        !---Forwading Lorenz model
        call forward_lorenz96(xtin,xtout)
        call forward_lorenz96(xin,xout)
        do iens=1,nens
           call forward_lorenz96(xfin(:,iens),xfout(:,iens))
        end do
        
        if(fo_cor)then
           do iens=1,nens
              call forward_lorenz96(xfin_c(:,iens),xfout_c(:,iens))
              call forward_lorenz96(xfin_nc(:,iens),xfout_nc(:,iens))
           end do
        end if

!        write(*,*) inan,inan_c,inan_nc
        
        !---Ensemble Kalman Filter
        if(mod(it,nint(int_obs/dt)) == 0)then
           
           call obs_operator(it,nobs,H,HT)
           
           if(st_obs <= it .and. it <= et_obs)then

              it_obs=it_obs+1

              !Make correlated obs.
              if(fo_cor)then
                 call make_cor_err(nens,xtout,xfout,xferr(:,it_obs),obserr(:,it_obs),xnull,0.d0)
                 if(inan_c == 0)then
                    call make_cor_err(nens,xtout,xfout_c,xferr_c(:,it_obs),obserr(:,it_obs),obserr_c(:,it_obs),cor_set)
                 else
                    xferr_c(:,it_obs)=0.d0
                    obserr_c(:,it_obs)=0.d0
                 end if
                 call make_cor_err(nens,xtout,xfout_nc,xferr_nc(:,it_obs),obserr(:,it_obs),obserr_nc(:,it_obs),cor_set)
              end if

              !Make obs.
              obs(:)=matmul(H(:,:),xtout(:)+obserr(:,it_obs))
              if(fo_cor)then
                 if(inan_c == 0)then
                    obs_c(:)=matmul(H(:,:),xtout(:)+obserr_c(:,it_obs))
                 else
                    obs_c(:)=0.d0
                 end if
                 obs_nc(:)=matmul(H(:,:),xtout(:)+obserr_nc(:,it_obs))
              end if
              
              if(execute_da == 3)then
                 if(fo_cor)then
                    !LETKF + orig. obs.
                    call etkfcc(nens,nobs,xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                         & obs,R,H,inf_parm(iinf),0.d0,inan)
                    !LETKFCC + cor obs.
                    call etkfcc(nens,nobs,xfout_c,xfoutmean_c,xfoutsprd_c,xaout_c,xaoutmean_c,xaoutsprd_c, &
                         & obs_c,R,H,inf_parm(iinf),cor_set,inan_c)
                    !LETKF + cor obs.
                    call etkfcc(nens,nobs,xfout_nc,xfoutmean_nc,xfoutsprd_nc,xaout_nc,xaoutmean_nc,xaoutsprd_nc, &
                         & obs_nc,R,H,inf_parm(iinf),0.d0,inan_nc)
                 else
                    call etkf(nens,nobs,xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                         & obs,R,H,inf_parm(iinf))
                 endif
              else if(execute_da == 4)then
                 if(fo_cor)then
                    !LETKF + orig obs.
                    call letkfcc(nens,nobs,sigma_loc, &
                         & xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                         & obs,R,H,inf_parm(iinf),0.d0,inan)
                    !LETKFCC + cor obs.
                    call letkfcc(nens,nobs,sigma_loc, &
                         & xfout_c,xfoutmean_c,xfoutsprd_c,xaout_c,xaoutmean_c,xaoutsprd_c, &
                         & obs_c,R,H,inf_parm(iinf),cor_set,inan_c)
                    !LETKF + cor obs.
                    call letkfcc(nens,nobs,sigma_loc, &
                         & xfout_nc,xfoutmean_nc,xfoutsprd_nc,xaout_nc,xaoutmean_nc,xaoutsprd_nc, &
                         & obs_nc,R,H,inf_parm(iinf),0.d0,inan_nc)
                 else
                    call letkf(nens,nobs,sigma_loc, &
                         & xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                         & obs,R,H,inf_parm(iinf))
                 end if
              else if(execute_da == 5)then
                 call ensrf(nens,nobs,sigma_loc, &
                      & xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd,obs,R,H,inf_parm(iinf))
              else if(execute_da == 6)then

                 if(fo_cor)then
                    
                    call poenkf(nens,nobs,sigma_loc, &
                         & xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                         & obs,obserr_ens(:,:,it_obs),R,H,inf_parm(iinf))
                    
                    call poenkf(nens,nobs,sigma_loc, &
                         & xfout_c,xfoutmean_c,xfoutsprd_c,xaout_c,xaoutmean_c,xaoutsprd_c, &
                         & obs_c,obserr_ens(:,:,it_obs),R,H,inf_parm(iinf))

                    call poenkf(nens,nobs,sigma_loc, &
                         & xfout_nc,xfoutmean_nc,xfoutsprd_nc,xaout_nc,xaoutmean_nc,xaoutsprd_nc, &
                         & obs_nc,obserr_ens(:,:,it_obs),R,H,inf_parm(iinf))
                    
                 else
                    
                    call poenkf(nens,nobs,sigma_loc, &
                      & xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd, &
                      & obs,obserr_ens(:,:,it_obs),R,H,inf_parm(iinf))
                    
                 end if
                 
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
              call ave_std(nx,xfoutsprd,sprdf(it_obs),tmp)
              call ave_std(nx,xaoutsprd,sprda(it_obs),tmp)
              
              if(fo_cor)then
                 if(inan_c == 0)then
                    call bias_rmse(nx,xfoutmean_c,xtout,biasf_c(it_obs),rmsef_c(it_obs))
                    call bias_rmse(nx,xaoutmean_c,xtout,biasa_c(it_obs),rmsea_c(it_obs))
                    call bias_rmse(nobs,obs_c,matmul(H(:,:),xtout(:)),biaso_c(it_obs),rmseo_c(it_obs))
                    call ave_std(nx,xfoutsprd_c,sprdf_c(it_obs),tmp)
                    call ave_std(nx,xaoutsprd_c,sprda_c(it_obs),tmp)
                 else
                    biasf_c(it_obs)=0.d0
                    rmsef_c(it_obs)=0.d0
                    biasa_c(it_obs)=0.d0
                    rmsea_c(it_obs)=0.d0
                    biaso_c(it_obs)=0.d0
                    rmseo_c(it_obs)=0.d0
                    sprdf_c(it_obs)=0.d0
                    sprda_c(it_obs)=0.d0                    
                 end if
                 call bias_rmse(nx,xfoutmean_nc,xtout,biasf_nc(it_obs),rmsef_nc(it_obs))
                 call bias_rmse(nx,xaoutmean_nc,xtout,biasa_nc(it_obs),rmsea_nc(it_obs))
                 call bias_rmse(nobs,obs_nc,matmul(H(:,:),xtout(:)),biaso_nc(it_obs),rmseo_nc(it_obs))
                 call ave_std(nx,xfoutsprd_nc,sprdf_nc(it_obs),tmp)
                 call ave_std(nx,xaoutsprd_nc,sprda_nc(it_obs),tmp)
              end if

              if(fo_cor)then
                 call write_ens_bias_rmse(nens,nobs,iinf,sigma_loc,cor_set,st_obs,it, &
                      & bias(it_obs),biasf(it_obs),biasa(it_obs),biaso(it_obs), &
                      & biasf_c(it_obs),biasa_c(it_obs),biaso_c(it_obs), &
                      & biasf_nc(it_obs),biasa_nc(it_obs),biaso_nc(it_obs), &
                      & rmse(it_obs),rmsef(it_obs),rmsea(it_obs),rmseo(it_obs), &
                      & rmsef_c(it_obs),rmsea_c(it_obs),rmseo_c(it_obs), &
                      & rmsef_nc(it_obs),rmsea_nc(it_obs),rmseo_nc(it_obs), &
                      & inf_parm(iinf))
              else
                 call write_ens_bias_rmse(nens,nobs,iinf,sigma_loc,cor_set,st_obs,it, &
                      & bias(it_obs),biasf(it_obs),biasa(it_obs),biaso(it_obs), &
                      & 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                      & rmse(it_obs),rmsef(it_obs),rmsea(it_obs),rmseo(it_obs), &
                      & 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                      & inf_parm(iinf))
              end if

              !xfout <-- xaout
              xfout(:,:)=xaout(:,:)
              if(fo_cor)then
                 xfout_c(:,:)=xaout_c(:,:)
                 xfout_nc(:,:)=xaout_nc(:,:)
              end if
              
              !---Post-regulation
              if((execute_da == 7 .or. execute_da == 8) .and. add_inf)then
                 do ix=1,nx
                    call make_gaussian_random(ix+nens,0.d0,sqrt(inf_parm(iinf)),nens,noise_pf(ix,:))
                 end do
                 xfout(:,:)=xfout(:,:)+noise_pf(:,:)
              end if
              
           else
              
              do ix=1,nx
                 call ave_std(nens,xfout(ix,:),xfoutmean(ix),xfoutsprd(ix))
              end do
              
              obs(:)=0.d0
              xaoutmean(:)=0.d0
              xaoutsprd(:)=0.d0

              if(fo_cor)then
                 do ix=1,nx
                    call ave_std(nens,xfout_c(ix,:),xfoutmean_c(ix),xfoutsprd_c(ix))
                    call ave_std(nens,xfout_nc(ix,:),xfoutmean_nc(ix),xfoutsprd_nc(ix))
                 end do
                 obs_c(:)=0.d0
                 obs_nc(:)=0.d0
                 xaoutmean_c(:)=0.d0
                 xaoutmean_nc(:)=0.d0
                 xaoutsprd_c(:)=0.d0
                 xaoutsprd_nc(:)=0.d0                 
              end if
              
           end if !st_obs <= it <= et_obs

           !---IO
           it_io=it_io+1

           if(fo_cor)then
              call write_ensmean(nens,nobs,iinf,sigma_loc,cor_set,it_io,it, &
                   & xtout,xout,xfoutmean,xaoutmean, &
                   & xfoutmean_c,xaoutmean_c,xfoutmean_nc,xaoutmean_nc,obserr,obserr_c,obserr_nc,inf_parm(iinf))
              call write_enssprd(nens,nobs,iinf,sigma_loc,cor_set,it_io,it, &
                   & xfoutsprd,xaoutsprd,xfoutsprd_c,xaoutsprd_c,xfoutsprd_nc,xaoutsprd_nc,inf_parm(iinf))
           else
              call write_ensmean(nens,nobs,iinf,sigma_loc,cor_set,it_io,it, &
                   & xtout,xout,xfoutmean,xaoutmean, &
                   & xnull,xnull,xnull,xnull,obserr,xnull,xnull,inf_parm(iinf))
              call write_enssprd(nens,nobs,iinf,sigma_loc,cor_set,it_io,it, &
                   & xfoutsprd,xaoutsprd,xnull,xnull,xnull,xnull,inf_parm(iinf))
           end if
                         
        end if !it_obs

        !xfin <-- xfout
        xtin(:)=xtout(:)
        xin(:)=xout(:)
        xfin(:,:)=xfout(:,:)
        if(fo_cor)then
           xfin_c(:,:)=xfout_c(:,:)
           xfin_c(:,:)=xfout_c(:,:)
           xfin_nc(:,:)=xfout_nc(:,:)
           xfin_nc(:,:)=xfout_nc(:,:)
        end if
        
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
     call ave_std(nt_obs,sprdf,ave_sprdf,std_sprdf)
     call ave_std(nt_obs,sprda,ave_sprda,std_sprda)
     
     write(6,'(a)') "Bias(Ave.,Std.) RMSE(Ave.,Std.)"
     write(6,'(a,4(f12.5,a))') "Free run: Bias(",ave_bias,",",std_bias,") RMSE(",ave_rmse,",",std_rmse,")"
     write(6,'(a,5(f12.5,a))') "Forecast: Bias(",ave_biasf,",",std_biasf,") RMSE(",ave_rmsef,",",std_rmsef,") SPREAD(",ave_sprdf,")"
     write(6,'(a,5(f12.5,a))') "Analysis: Bias(",ave_biasa,",",std_biasa,") RMSE(",ave_rmsea,",",std_rmsea,") SPREAD(",ave_sprda,")"
     write(6,'(a,4(f12.5,a))') "Observation: Bias(",ave_biaso,",",std_biaso,") RMSE(",ave_rmseo,",",std_rmseo,")"

     if(fo_cor)then

        if(inan_c == 0)then
           call ave_std(nt_obs,biasf_c,ave_biasf_c,std_biasf_c)
           call ave_std(nt_obs,biasa_c,ave_biasa_c,std_biasa_c)
           call ave_std(nt_obs,biaso_c,ave_biaso_c,std_biaso_c)
           call ave_std(nt_obs,rmsef_c,ave_rmsef_c,std_rmsef_c)
           call ave_std(nt_obs,rmsea_c,ave_rmsea_c,std_rmsea_c)
           call ave_std(nt_obs,rmseo_c,ave_rmseo_c,std_rmseo_c)
           call ave_std(nt_obs,sprdf_c,ave_sprdf_c,std_sprdf_c)
           call ave_std(nt_obs,sprda_c,ave_sprda_c,std_sprda_c)
        else
           ave_biasf_c=0.d0
           std_biasf_c=0.d0
           ave_biasa_c=0.d0
           std_biasa_c=0.d0
           ave_biaso_c=0.d0
           std_biaso_c=0.d0
           ave_rmsef_c=0.d0
           std_rmsef_c=0.d0
           ave_rmsea_c=0.d0
           std_rmsea_c=0.d0
           ave_rmseo_c=0.d0
           std_rmseo_c=0.d0
           ave_sprdf_c=0.d0
           std_sprdf_c=0.d0
           ave_sprda_c=0.d0
           std_sprda_c=0.d0
        end if
           
        write(6,'(a)') "Bias(Ave.,Std.) RMSE(Ave.,Std.) in LETKFCC with correlated obs."
        write(6,'(a,5(f12.5,a))') &
             & "Forecast: Bias(",ave_biasf_c,",",std_biasf_c,") RMSE(",ave_rmsef_c,",",std_rmsef_c,") SPREAD(",ave_sprdf_c,")"
        write(6,'(a,5(f12.5,a))') &
             & "Analysis: Bias(",ave_biasa_c,",",std_biasa_c,") RMSE(",ave_rmsea_c,",",std_rmsea_c,") SPREAD(",ave_sprda_c,")"
        write(6,'(a,4(f12.5,a))') &
             & "Observation: Bias(",ave_biaso_c,",",std_biaso_c,") RMSE(",ave_rmseo_c,",",std_rmseo_c,")"
        
        call ave_std(nt_obs,biasf_nc,ave_biasf_nc,std_biasf_nc)
        call ave_std(nt_obs,biasa_nc,ave_biasa_nc,std_biasa_nc)
        call ave_std(nt_obs,biaso_nc,ave_biaso_nc,std_biaso_nc)
        call ave_std(nt_obs,rmsef_nc,ave_rmsef_nc,std_rmsef_nc)
        call ave_std(nt_obs,rmsea_nc,ave_rmsea_nc,std_rmsea_nc)
        call ave_std(nt_obs,rmseo_nc,ave_rmseo_nc,std_rmseo_nc)
        call ave_std(nt_obs,sprdf_nc,ave_sprdf_nc,std_sprdf_nc)
        call ave_std(nt_obs,sprda_nc,ave_sprda_nc,std_sprda_nc)
        
        write(6,'(a)') "Bias(Ave.,Std.) RMSE(Ave.,Std.) in LETKF with correlated obs."
        write(6,'(a,5(f12.5,a))') &
             & "Forecast: Bias(",ave_biasf_nc,",",std_biasf_nc,") RMSE(",ave_rmsef_nc,",",std_rmsef_nc,") SPREAD(",ave_sprdf_nc,")"
        write(6,'(a,5(f12.5,a))') &
             & "Analysis: Bias(",ave_biasa_nc,",",std_biasa_nc,") RMSE(",ave_rmsea_nc,",",std_rmsea_nc,") SPREAD(",ave_sprda_nc,")"
        write(6,'(a,4(f12.5,a))') &
             & "Observation: Bias(",ave_biaso_nc,",",std_biaso_nc,") RMSE(",ave_rmseo_nc,",",std_rmseo_nc,")"
        
     end if

     !Correlation btw forecast vs. obs. error
     if(fo_cor)then
        do jx=1,nx
           do ix=1,nx
              call correlation(nt_obs,xferr(ix,:),obserr(jx,:),cor(ix,jx))
              if(inan_c == 0)then
                 call correlation(nt_obs,xferr_c(ix,:),obserr_c(jx,:),cor_c(ix,jx))
              else
                 cor_c(ix,jx)=0.d0
              end if
              call correlation(nt_obs,xferr_nc(ix,:),obserr_nc(jx,:),cor_nc(ix,jx))
           end do
        end do
     end if
     
     !---IO
     if(fo_cor)then
        call write_ens_ave_bias_rmse_sprd(nens,nobs,iinf,sigma_loc,cor_set, &
             & ave_bias,ave_biasf,ave_biasa,ave_biaso, &
             & ave_biasf_c,ave_biasa_c,ave_biaso_c, &
             & ave_biasf_nc,ave_biasa_nc,ave_biaso_nc, &
             & ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo, &
             & ave_rmsef_c,ave_rmsea_c,ave_rmseo_c, &
             & ave_rmsef_nc,ave_rmsea_nc,ave_rmseo_nc, &
             & ave_sprdf,ave_sprda, &
             & ave_sprdf_c,ave_sprda_c, &
             & ave_sprdf_nc,ave_sprda_nc, & 
             & inf_parm(iinf))
        call write_ens_cor(nens,nobs,iinf,sigma_loc,cor_set, &
             & cor,cor_c,cor_nc)

     else
        call write_ens_ave_bias_rmse_sprd(nens,nobs,iinf,sigma_loc,cor_set, &
             & ave_bias,ave_biasf,ave_biasa,ave_biaso, &
             & 0.d0,0.d0,0.d0, &
             & 0.d0,0.d0,0.d0, &
             & ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo, &
             & 0.d0,0.d0,0.d0, &
             & 0.d0,0.d0,0.d0, &
             & 0.d0,0.d0, &
             & 0.d0,0.d0, &
             & 0.d0,0.d0, &
             & inf_parm(iinf))
     end if
     
     !--- Effective size for PF
     if(execute_da == 7 .or. execute_da == 8)then
        call ave_std(nt_obs,Neff,ave_Neff,std_Neff)
        call write_Neff(nens,nobs,iinf,sigma_loc,ave_Neff,std_Neff,inf_parm(iinf))
     end if
     
  end do !ifo
  end do !sigma_loc
  end do !iinf

  deallocate(obs,R,H,HT)

  if(fo_cor)then
     deallocate(xferr,xferr_c,xferr_nc)
     deallocate(obs_c,obserr_c)
     deallocate(obs_nc,obserr_nc)
  end if
  
  end do !nobs

  deallocate(xfin,xfout,xaout)
  deallocate(noise_ens)
  if(fo_cor)then
     deallocate(xfin_c,xfout_c,xaout_c)
     deallocate(xfin_nc,xfout_nc,xaout_nc)
  end if
  
  if(execute_da == 6)then
     deallocate(obserr_ens)
  end if

  if(execute_da == 4 .or. execute_da == 7 .or. execute_da == 8)then
     deallocate(noise_pf)
     deallocate(Neff)
  end if
     
  end do !nens
  
  write(6,'(a)') "===End DA run==="
  
  !  status=system("csh fig.csh "//trim(dir(execute_da)))
  
  deallocate(bias,biasf,biasa,biaso)
  deallocate(rmse,rmsef,rmsea,rmseo)
  deallocate(sprdf,sprda)
  deallocate(obserr)
  if(fo_cor)then
     deallocate(biasf_c,biasa_c,biaso_c)
     deallocate(rmsef_c,rmsea_c,rmseo_c)
     deallocate(sprdf_c,sprda_c)
     deallocate(biasf_nc,biasa_nc,biaso_nc)
     deallocate(rmsef_nc,rmsea_nc,rmseo_nc)
     deallocate(sprdf_nc,sprda_nc)
  end if

  write(6,'(a)') "=== All End === "

end program main
