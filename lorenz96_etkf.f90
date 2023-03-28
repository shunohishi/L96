program main

  use parameter
  implicit none

    !---Common
  integer it,nt
  integer ix
  integer iens
  integer status,system
  
  !Obs.
  integer iobs
  integer it_obs,nt_obs
  integer st_obs,et_obs

  !IO
  integer it_io
  
  !---Variable
  double precision xtin(nx),xin(nx),xfin(nx,nens)
  double precision xtout(nx),xout(nx),xfout(nx,nens),xaout(nx,nens)
  double precision xfinmean(nx)                !Ensemble mean
  double precision xfoutmean(nx),xaoutmean(nx)
  double precision xfinsprd(nx)                !Ensemble sprd
  double precision xfoutsprd(nx),xaoutsprd(nx)
  double precision obs(nobs)
  double precision,allocatable :: random_obs(:,:)
  
  !---Error covariance matrix
  double precision R(nobs,nobs)         !Observation error covarinance matrix: R
  double precision H(nobs,nx),HT(nx,nobs) !Observation operator
  
  !---Static variable
  double precision,allocatable :: bias(:),biasf(:),biasa(:),biaso(:) !Bias
  double precision,allocatable :: rmse(:),rmsef(:),rmsea(:),rmseo(:) !RMSE
  double precision ave_bias,ave_biasf,ave_biasa,ave_biaso
  double precision std_bias,std_biasf,std_biasa,std_biaso
  double precision ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo
  double precision std_rmse,std_rmsef,std_rmsea,std_rmseo
  
  !---Noise
  double precision noise(nx)
  
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
  if(nobs > nx)then
     write(6,'(a)') "***Error: nobs > nx"
     stop
  end if
  if(execute_da == 3)then
     write(6,'(a)') "Execute ETKF"
  else if(execute_da == 4)then
     write(6,'(a)') "Execute LETKF"
  else
     write(6,'(a)') "***Error: Choose execute_da=3,4"
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
  !Random number
  allocate(random_obs(nobs,nt_obs))
  do iobs=1,nobs
     call make_gaussian_random(iobs,0.d0,err_obs,nt_obs,random_obs(iobs,:))
  end do
 
  do it=1,nt
     
     if(mod(it,1000) == 0 .or. mod(it,nt) == 0) write(6,'(i6,a,i6)') it,"/",nt

     !Initialization
     if(it == 1)then
        it_obs=0
        it_io=0
        !True run
        noise(:)=0.d0
        call initialization(amp_t,deg_t,noise,xtin)
        !Free run
        call make_gaussian_random(0,noise_mean_da,noise_sd_da,nx,noise)
        call initialization(amp_da,deg_da,noise,xin)
        !Ensemble run
        do iens=1,nens
           call make_gaussian_random(iens,noise_mean_da,noise_sd_da,nx,noise)
           call initialization(amp_da,deg_da,noise,xfin(:,iens))
        end do
        do ix=1,nx
           call ave_std(nens,xfin(ix,:),xfinmean(ix),xfinsprd(ix))
        end do
        call write_data(it_io,0,xtin,xin,xfinmean,xnull,xnull)
        call write_sprd(it_io,0,xfinsprd,xnull)
        !Matrix
        call initialization_ECM(R,Ptmp)
     end if
     
     !Forwading Lorenz model
     call forward_lorenz96(xtin,xtout)
     call forward_lorenz96(xin,xout)
     do iens=1,nens
        call forward_lorenz96(xfin(:,iens),xfout(:,iens))
     end do
 
     !Kalman Filter
     if(mod(it,nint(int_obs/dt)) == 0)then
        
        call obs_operator(it,H,HT)
        
        if(st_obs <= it .and. it <= et_obs)then
           
           it_obs=it_obs+1
           obs(:)=matmul(H(:,:),xtout(:))+random_obs(:,it_obs)

           if(execute_da == 3)then
              call etkf(xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd,obs,R,H)
           else if(execute_da == 4)then
              call letkf(xfout,xfoutmean,xfoutsprd,xaout,xaoutmean,xaoutsprd,obs,R,H)
           end if
           
           call bias_rmse(nx,xout,xtout,bias(it_obs),rmse(it_obs))
           call bias_rmse(nx,xfoutmean,xtout,biasf(it_obs),rmsef(it_obs))
           call bias_rmse(nx,xaoutmean,xtout,biasa(it_obs),rmsea(it_obs))
           call bias_rmse(nobs,obs,matmul(H(:,:),xtout(:)),biaso(it_obs),rmseo(it_obs))

           call write_bias_rmse(st_obs,it,bias(it_obs),biasf(it_obs),biasa(it_obs),biaso(it_obs), &
                & rmse(it_obs),rmsef(it_obs),rmsea(it_obs),rmseo(it_obs))

           xfout(:,:)=xaout(:,:)
           
        else

           do ix=1,nx
              call ave_std(nens,xfin(ix,:),xfoutmean(ix),xfoutsprd(ix))
           end do           

           obs(:)=0.d0
           xaoutmean(:)=0.d0
           xaoutsprd(:)=0.d0
           
        end if

        it_io=it_io+1
        call write_data(it_io,it,xtout,xout,xfoutmean,xaoutmean,matmul(HT(:,:),obs(:)))
        call write_sprd(it_io,it,xfoutsprd,xaoutsprd)
        
     end if
     
     xtin(:)=xtout(:)
     xin(:)=xout(:)
     xfin(:,:)=xfout(:,:)

  end do

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
  
  write(6,'(a)') "===End DA run==="
  
  !  status=system("csh fig.csh "//trim(dir(execute_da)))
  
  deallocate(bias,biasf,biasa,biaso)
  deallocate(rmse,rmsef,rmsea,rmseo)
  deallocate(random_obs)

  write(6,'(a)') "=== All End === "

end program main

