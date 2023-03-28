program main

  use parameter
  implicit none

  !---Common
  integer it,st,et,nt
  integer ix,jx
  integer iinf
  integer icor
  integer iday
  integer status,system
  integer flag,flagc,flagnc
  integer Pflag
  
  !Obs.
  integer iobs,nobs
  integer it_obs,nt_obs
  integer st_obs,et_obs
  
  !IO
  integer it_io
  
  !---Variable
  double precision xtin(nx),xin(nx),xfin(nx),xain(nx)
  double precision xtout(nx),xout(nx),xfout(nx),xaout(nx)
  double precision,allocatable :: obs(:)              !Obs. in obs. space
  double precision,allocatable :: obserr(:,:)         !Obs. error in model space
  double precision,allocatable :: Hobserr(:,:)        !           in obs space
  
  !---Error covariance matrix
  double precision,allocatable :: xferr(:,:)           !Forecast error  
  double precision Pf(nx,nx),Pftmp(nx,nx),Pfinf(nx,nx) !Forecast error covariance matrix: Pf
  double precision Pa(nx,nx)                           !Analysis error covariance matrix: Pa
  double precision,allocatable :: Pfd(:),Pad(:)        !Digonal averaged Pf,Pa
  double precision Pfdave,Pfdstd,Padave,Padstd
  double precision,allocatable :: R(:,:)               !Observation error covarinance matrix: R
  double precision,allocatable :: H(:,:),HT(:,:)       !Observation operator
  double precision,allocatable :: C(:,:)
  double precision tmp(nx,nx)
  double precision inf_parm(ninf)                    !Inflation parameter
  
  !---Static variable
  double precision,allocatable :: bias(:),biasf(:),biasa(:),biaso(:) !Bias
  double precision,allocatable :: rmse(:),rmsef(:),rmsea(:),rmseo(:) !RMSE
  double precision ave_bias,ave_biasf,ave_biasa,ave_biaso
  double precision std_bias,std_biasf,std_biasa,std_biaso
  double precision ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo
  double precision std_rmse,std_rmsef,std_rmsea,std_rmseo

  !---Forecaset and observation error correlation
  !FO/SO
  integer ifo,sfo,efo,nfo
  integer it_fo,nt_fo
  integer st_fo,et_fo
  integer,allocatable :: t_fo_obs(:)                                 

  !Include C calculation
  double precision fo_cc                              !Correlation between forecast and obs. error
  double precision xfcin(nx),xfcout(nx),xacout(nx)    !Forecast/Analysis

  double precision,allocatable :: xfcerr(:,:)         !Forecast error
  double precision Pfc(nx,nx),Pfctmp(nx,nx),Pfcinf(nx,nx)
  double precision Pac(nx,nx)
  double precision,allocatable :: Pfcd(:),Pacd(:)        !Digonal averaged Pf,Pa
  double precision Pfcdave,Pfcdstd,Pacdave,Pacdstd
  double precision,allocatable :: obsc(:)          !Correlated obs.
  double precision,allocatable :: obserrc(:,:)     !Correlated obs. error in model space
  double precision,allocatable :: Hobserrc(:,:)    !                      in obs space
  
  double precision,allocatable :: Cc(:,:)           !Forecaset and obs. error correlation
  double precision,allocatable :: biasfc(:),biasac(:),biasoc(:)
  double precision,allocatable :: rmsefc(:),rmseac(:),rmseoc(:)
  double precision ave_biasfc,ave_biasac,ave_biasoc
  double precision std_biasfc,std_biasac,std_biasoc
  double precision ave_rmsefc,ave_rmseac,ave_rmseoc
  double precision std_rmsefc,std_rmseac,std_rmseoc

  !Exclude C calculation
  double precision xfncin(nx),xfncout(nx),xancout(nx)
  double precision,allocatable :: xfncerr(:,:)
  double precision Pfnc(nx,nx),Pfnctmp(nx,nx),Pfncinf(nx,nx)
  double precision Panc(nx,nx)
  double precision,allocatable :: Pfncd(:),Pancd(:)        !Digonal averaged Pf,Pa
  double precision Pfncdave,Pfncdstd,Pancdave,Pancdstd
  double precision,allocatable :: obsnc(:)
  double precision,allocatable :: obserrnc(:,:)
  double precision,allocatable :: Hobserrnc(:,:)
  double precision,allocatable :: Cnc(:,:)
  
  double precision,allocatable :: biasfnc(:),biasanc(:),biasonc(:)
  double precision,allocatable :: rmsefnc(:),rmseanc(:),rmseonc(:)
  double precision ave_biasfnc,ave_biasanc,ave_biasonc
  double precision std_biasfnc,std_biasanc,std_biasonc
  double precision ave_rmsefnc,ave_rmseanc,ave_rmseonc
  double precision std_rmsefnc,std_rmseanc,std_rmseonc  

  !Error
  double precision avefc,stdfc !Forecast
  double precision avefnc,stdfnc !Forecast
  double precision aveoc,stdoc !Observation
  double precision aveonc,stdonc !Observation
  double precision,allocatable :: corf(:,:),corfc(:,:),corfnc(:,:)
  
  !Q
  double precision,allocatable :: q(:,:),qc(:,:)
  double precision corq(nx,nx),corqc(nx,nx)
  
  !---Noise
  double precision noise(nx)
  
  !---Null
  double precision xnull(nx)
  double precision Pnull(nx,nx)
  
  !---check nx,nobs,execute_da
  if(nx <= 3)then
     write(*,*) "***Error: Choose nx > 3"
     stop
  end if
  if(execute_da == 1)then
     write(6,'(a)') "Execute KF"
  else if(execute_da == 2)then
     write(6,'(a)') "Execute 3DVAR"
  else
     write(6,'(a)') "***Error: Choose execute_da=1,2"
     stop
  end if
  
  !---Number of time
  nt=nint(dble(nday)*int_day/dt)
  write(*,*) "Number of time:",nt

  !--Make directory
  status=system("mkdir -p "//trim(dir(execute_da)))

  !---Null
  xnull(:)=0.d0
  Pnull(:,:)=0.d0
  
  !______________________________________________________________________________________________________________
  
  !---DA run: Kalman Filter or 3D-VAR
  write(*,*) "===Start DA run==="

  !---Observation
  !Time
  st_obs=nint((dble(yyyys_obs)*365.d0+dble(dds_obs)+dble(hhs_obs)/24.d0)*int_day/dt)
  et_obs=nint((dble(yyyye_obs)*365.d0+dble(dde_obs)+dble(hhe_obs)/24.d0)*int_day/dt)
  nt_obs=(et_obs-st_obs)/nint(int_obs/dt)+1
  write(*,*) "Obs time:",st_obs,"-",et_obs
  write(*,*) "Number of obs time:",nt_obs
  allocate(bias(nt_obs),biasf(nt_obs),biasa(nt_obs),biaso(nt_obs))
  allocate(rmse(nt_obs),rmsef(nt_obs),rmsea(nt_obs),rmseo(nt_obs))
  allocate(obserr(nx,nt_obs))
  
  !---FO
  if(fo_cor)then
     st_fo=nint((dble(yyyys_fo)*365.d0+dble(dds_fo)+dble(hhs_fo)/24.d0)*int_day/dt)
     et_fo=nint((dble(yyyye_fo)*365.d0+dble(dde_fo)+dble(hhe_fo)/24.d0)*int_day/dt)
  end if

  if(so_cor)then
     st_fo=nint((dble(yyyys_so)*365.d0+dble(dds_so)+dble(hhs_so)/24.d0)*int_day/dt)
     et_fo=nint((dble(yyyye_so)*365.d0+dble(dde_so)+dble(hhe_so)/24.d0)*int_day/dt)
  end if
  
  !---FO or SO
  if(fo_cor .or. so_cor)then

     write(*,*) "FO cor. time:",st_fo,"-",et_fo  
     nt_fo=(et_fo-st_fo)/nint(int_obs/dt)+1
     write(*,*) "Number of Correlated time:",nt_fo     
     allocate(t_fo_obs(nt_fo))
     allocate(xferr(nx,nt_obs),xfcerr(nx,nt_obs),xfncerr(nx,nt_obs))
     allocate(q(nx,nt_fo),qc(nx,nt_fo))
     allocate(obserrc(nx,nt_obs),obserrnc(nx,nt_obs))
     allocate(biasfc(nt_obs),biasac(nt_obs),biasoc(nt_obs))
     allocate(biasfnc(nt_obs),biasanc(nt_obs),biasonc(nt_obs))
     allocate(rmsefc(nt_obs),rmseac(nt_obs),rmseoc(nt_obs))
     allocate(rmsefnc(nt_obs),rmseanc(nt_obs),rmseonc(nt_obs))

     allocate(Pfd(nt_obs),Pad(nt_obs))
     allocate(Pfcd(nt_obs),Pacd(nt_obs))
     allocate(Pfncd(nt_obs),Pancd(nt_obs))

     
     it_obs=0
     it_fo=0
     do it=1,nt
        if(mod(it,nint(int_obs/dt)) == 0)then
           if(st_obs <= it .and. it <= et_obs)then
              it_obs=it_obs+1
           end if
           if(st_fo <= it .and. it <= et_fo)then
              it_fo=it_fo+1
              t_fo_obs(it_fo)=it_obs
           end if
        end if
     end do
     
  end if
  
  !---Random number
  !Initial condition
  call make_gaussian_random(0,noise_mean_da,noise_sd_da,nx,noise)
  !Obs.
  do ix=1,nx
     call make_gaussian_random(ix,0.d0,err_obs,nt_obs,obserr(ix,:))
  end do
  
  !---Make covariance inflation
  call make_inf(inf_parm)
  
  !________________________________________________________________________________________________________________________________
  !---Start do loop
!  do nobs=1,nx
  do nobs=nx,nx

     allocate(obs(nobs),R(nobs,nobs),H(nobs,nx),HT(nx,nobs),C(nx,nobs))
     allocate(Hobserr(nobs,nt_obs))
     if(fo_cor .or. so_cor)then
        allocate(obsc(nobs),obsnc(nobs),Cc(nx,nobs),Cnc(nx,nobs))
        allocate(Hobserrc(nobs,nt_obs),Hobserrnc(nobs,nt_obs))
        allocate(corf(nx,nobs),corfc(nx,nobs),corfnc(nx,nobs))
     end if
     
     do iinf=1,ninf
     !do iinf=2,2
        
        !SO
        if(so_cor)then
           do ix=1,nx
              call make_gaussian_random(nx+ix,0.d0,sqrt(inf_parm(iinf)),nt_fo,q(ix,:))
              call make_cor_random(nt_fo,err_obs,sqrt(inf_parm(iinf)),sqrt(inf_parm(iinf)),fo_cc, &
                   & obserr(ix,t_fo_obs(1):t_fo_obs(nt_fo)),q(ix,1:nt_fo),qc(ix,1:nt_fo))
           end do
           do ix=1,nx
              do jx=1,nx
                 call correlation(nt_fo,q(ix,1:nt_fo),obserr(jx,t_fo_obs(1):t_fo_obs(nt_fo)),corq(ix,jx))
                 call correlation(nt_fo,qc(ix,1:nt_fo),obserr(jx,t_fo_obs(1):t_fo_obs(nt_fo)),corqc(ix,jx))
              end do
           end do
           tmp(:,:)=matmul(q(:,:),transpose(q(:,:)))/dble(nt_fo)
        end if

        if(fo_cor .or. so_cor)then
!           sfo=0
!           efo=10
           sfo=10
           efo=10
!           sfo=-10
!           efo=10
        else
           sfo=1
           efo=1
        end if
        nfo=efo-sfo+1
        
        do ifo=sfo,efo,1

           fo_cc=0.1d0*dble(ifo)
           
           do it=1,nt

              if(mod(it,1000) == 0 .or. mod(it,nt) == 0) write(*,*) nobs,"/",nx,iinf,"/",ninf,it,"/",nt

              !---Initialization
              if(it == 1)then
                 it_obs=0
                 it_io=0
                 flag=0
                 !True run
                 call initialization(amp_t,deg_t,xnull,xtin) !noise=xnull=0.d0
                 !Free run
                 call initialization(amp_da,deg_da,noise,xin)
                 call initialization(amp_da,deg_da,noise,xfin)
                 !Correlated run
                 if(fo_cor .or. so_cor)then
                    it_fo=0
                    flagc=0
                    flagnc=0
                    call initialization(amp_da,deg_da,noise,xfcin)
                    call initialization(amp_da,deg_da,noise,xfncin)
                 end if

                 !Write data
                 call write_data(nobs,iinf,it_io,0,xtin,xin,xfin,xnull,xnull,inf_parm(iinf))
                 call write_P(nobs,iinf,it_io,0,Pnull,Pnull,inf_parm(iinf))

              end if

              !---Forwading Lorenz model
              call forward_lorenz96(xtin,xtout)
              call forward_lorenz96(xin,xout)
              call forward_lorenz96(xfin,xfout)
              if(fo_cor .or. so_cor)then
                 call forward_lorenz96(xfcin,xfcout)
                 call forward_lorenz96(xfncin,xfncout)
              end if

              !---Check value
              call check_nan(xtout,xout,xfout,flag)
              if(fo_cor .or. so_cor)then
                 call check_nan(xtout,xout,xfcout,flagc)
                 call check_nan(xtout,xout,xfncout,flagnc)
              end if

              !---Initialization for covariance matrix
              if(it == st_obs) call initialization_ECM(nobs,R,Pf)
              if(it == st_obs .and. (fo_cor .or. so_cor))then
                 call initialization_ECM(nobs,R,Pfc)
                 call initialization_ECM(nobs,R,Pfnc)
              end if

              !---Lyapunov equation:Pf(t)=MPa(t-1)MT 
              if(st_obs <= it .and. it <= et_obs)then
                 if(flag == 0) call lyapunov_equation(xfout,Pf)
                 if((fo_cor .or. so_cor) .and. flagc == 0) call lyapunov_equation(xfcout,Pfc)
                 if((fo_cor .or. so_cor) .and. flagnc == 0) call lyapunov_equation(xfncout,Pfnc)
              end if

              !---Kalman Filter
              if(mod(it,nint(int_obs/dt)) == 0)then

                 call obs_operator(it,nobs,H,HT)

                 if(st_obs <= it .and. it <= et_obs)then

                    it_obs=it_obs+1
                    if((fo_cor .or. so_cor) .and. st_fo <= it .and. it <= et_fo)then
                       it_fo=it_fo+1
                    end if

                    !---Forecast error
                    if(fo_cor)then
                       xferr(:,it_obs)=xfout(:)-xtout(:)
                       xfcerr(:,it_obs)=xfcout(:)-xtout(:)
                       xfncerr(:,it_obs)=xfncout(:)-xtout(:)
                    end if

                    !---Observation
                    !Indepdendent obs.
                    obs(:)=matmul(H(:,:),xtout(:)+obserr(:,it_obs))
                    if(it_obs == 1) Hobserr(:,:)=matmul(H(:,:),obserr(:,:))
                    if(fo_cor .or. so_cor)then
                       obsc(:)=obs(:)
                       obsnc(:)=obs(:)
                       Hobserrc(:,it_obs)=matmul(H(:,:),obserrc(:,it_obs))
                       Hobserrnc(:,it_obs)=matmul(H(:,:),obserrnc(:,it_obs))
                    end if

                    !FO correlated run
                    if(fo_cor .and. st_fo <= it .and. it <= et_fo)then

                       !correlated obs. error
                       call ave_std(nx,xfcerr(:,it_obs),avefc,stdfc)
                       call ave_std(nx,obserr(:,it_obs),aveoc,stdoc)
                       call make_cor_random(nx,stdfc,stdoc,err_obs,fo_cc, &
                            & xfcerr(:,it_obs),obserr(:,it_obs),obserrc(:,it_obs))
                       call ave_std(nx,xfncerr(:,it_obs),avefnc,stdfnc)
                       call ave_std(nx,obserr(:,it_obs),aveonc,stdonc)
                       call make_cor_random(nx,stdfnc,stdonc,err_obs,fo_cc, &
                            & xfncerr(:,it_obs),obserr(:,it_obs),obserrnc(:,it_obs))

                       !Make correlated obs.
                       obsc(:)=matmul(H(:,:),xtout(:)+obserrc(:,it_obs))
                       obsnc(:)=matmul(H(:,:),xtout(:)+obserrnc(:,it_obs))
                       Hobserrc(:,it_obs)=matmul(H(:,:),obserrc(:,it_obs))
                       Hobserrnc(:,it_obs)=matmul(H(:,:),obserrnc(:,it_obs))
                    end if !fo_cor

                    !---Forecast
                    if(so_cor .and. st_fo <= it .and. it <= et_fo)then
                       xfout(:)=xfout(:)+q(:,it_fo)
                       xfcout(:)=xfcout(:)+qc(:,it_fo)
                       xfncout(:)=xfncout(:)+qc(:,it_fo)
                    end if
                    
                    call diag_ave_std(nx,Pf,Pfd(it_obs))
                    if(fo_cor .or. so_cor)then
                       call diag_ave_std(nx,Pfc,Pfcd(it_obs))
                       call diag_ave_std(nx,Pfnc,Pfncd(it_obs))
                    end if
                       
                    !---Covariance inflation
!                    call covariance_inflation(Pf,inf_parm(iinf))
!                    if(fo_cor .or. so_cor)then
!                    call covariance_inflation(Pfc,inf_parm(iinf))
!                    call covariance_inflation(Pfnc,inf_parm(iinf))
!                    end if

                    !---Check Pf
                    call check_P(it_fo,1,Pf,Pflag)
                    if(fo_cor .or. so_cor)then
                       call check_P(it_fo,1,Pfc,Pflag)
                       call check_P(it_fo,1,Pfnc,Pflag)
                    end if
                    
                    !---Kalman filter
                    !Original KF
                    if(flag == 0)then
                       C(:,:)=0.d0
                       call kalman_filter(nobs,xfout,xaout,Pf,Pa,obs,R,H,HT,C)
                    else
                       xaout(:)=0.d0
                       Pa(:,:)=0.d0
                    end if

                    !Correlated KF with C
                    if(fo_cor .or. so_cor)then
                       if(flagc == 0)then
                          if(st_fo <= it .and. it <= et_fo)then
                             if(fo_cor)then
                                call make_C_matrix_fo(nobs,fo_cc,stdfc,stdoc,Pfc,Cc)
                             end if
                             if(so_cor) call make_C_matrix_so(nobs,nt_fo,qc(:,1:nt_fo),Hobserr(:,t_fo_obs(1):t_fo_obs(nt_fo)),Cc)
                          else
                             Cc(:,:)=0.d0
                          end if

!                          if(ifo <= -2) Cc(:,:)=0.d0 !for avoid divergence
                          call kalman_filter(nobs,xfcout,xacout,Pfc,Pac,obsc,R,H,HT,Cc)
                          call check_P(it_fo,2,Pac,Pflag)
                          if(Pflag == 1)then
                             Cc(:,:)=0.d0
                             call kalman_filter(nobs,xfcout,xacout,Pfc,Pac,obsc,R,H,HT,Cc)
                             call check_P(it_fo,2,Pac,Pflag)
                             if(Pflag == 1)then
                                write(*,*) "***Error: Negative P"
                                stop
                             end if
                          end if
                       else
                          xacout(:)=0.d0
                          Pac(:,:)=0.d0
                       endif
                    end if
                    !Correlated KF without C
                    if(fo_cor .or. so_cor)then
                       if(flagnc == 0)then
                          Cnc(:,:)=0.d0
                          call kalman_filter(nobs,xfncout,xancout,Pfnc,Panc,obsnc,R,H,HT,Cnc)
                       else
                          xancout(:)=0.d0
                          Panc(:,:)=0.d0
                       end if
                    end if

                    call diag_ave_std(nx,Pa,Pad(it_obs))
                    if(fo_cor .or. so_cor)then
                       call diag_ave_std(nx,Pac,Pacd(it_obs))
                       call diag_ave_std(nx,Panc,Pancd(it_obs))
                    end if
                    
                    !---Covariance inflation
                    call covariance_inflation(Pa,inf_parm(iinf))
                    if(fo_cor .or. so_cor)then
                       call covariance_inflation(Pac,inf_parm(iinf))
                       call covariance_inflation(Panc,inf_parm(iinf))
                    end if

                    !---Check Pa
                    call check_P(it_fo,2,Pa,Pflag)
                    if(fo_cor .or. so_cor)then
                       call check_P(it_fo,2,Pac,Pflag)
                       call check_P(it_fo,2,Panc,Pflag)
                    end if
                    
                    !---Bias & RMSE at it_obs
                    call bias_rmse(nx,xout,xtout,bias(it_obs),rmse(it_obs))
                    call bias_rmse(nx,xfout,xtout,biasf(it_obs),rmsef(it_obs))
                    call bias_rmse(nx,xaout,xtout,biasa(it_obs),rmsea(it_obs))
                    call bias_rmse(nobs,obs,matmul(H(:,:),xtout(:)),biaso(it_obs),rmseo(it_obs))

                    if(fo_cor .or. so_cor)then
                       call bias_rmse(nx,xfcout,xtout,biasfc(it_obs),rmsefc(it_obs))
                       call bias_rmse(nx,xacout,xtout,biasac(it_obs),rmseac(it_obs))
                       call bias_rmse(nx,xfncout,xtout,biasfnc(it_obs),rmsefnc(it_obs))
                       call bias_rmse(nx,xancout,xtout,biasanc(it_obs),rmseanc(it_obs))
                       call bias_rmse(nobs,obsc,matmul(H(:,:),xtout(:)),biasoc(it_obs),rmseoc(it_obs))
                       call bias_rmse(nobs,obsnc,matmul(H(:,:),xtout(:)),biasonc(it_obs),rmseonc(it_obs))
                    end if
                                           
                    !---Write data
                    call write_bias_rmse(nobs,iinf,st_obs,it,bias(it_obs),biasf(it_obs),biasa(it_obs),biaso(it_obs), &
                         & rmse(it_obs),rmsef(it_obs),rmsea(it_obs),rmseo(it_obs),inf_parm(iinf))
                    
                    !---Replace data
                    xfout(:)=xaout(:)
                    Pf(:,:)=Pa(:,:)
                    if(fo_cor .or. so_cor)then
                       xfcout(:)=xacout(:)
                       Pfc(:,:)=Pac(:,:)
                       xfncout(:)=xancout(:)
                       Pfnc(:,:)=Panc(:,:)
                    end if

                 end if !st_obs <= it .and. it <= et_obs

                 it_io=it_io+1           
                 call write_data(nobs,iinf,it_io,it,xtout,xout,xfout,xaout,matmul(HT(:,:),obs(:)),inf_parm(iinf))
                 call write_P(nobs,iinf,it_io,it,Pf,Pa,inf_parm(iinf))
                 
              end if !mod(it,nint(int_obs/dt)) == 0

              xtin(:)=xtout(:)
              xin(:)=xout(:)
              xfin(:)=xfout(:)
              if(fo_cor .or. so_cor)then
                 xfcin(:)=xfcout(:)
                 xfncin(:)=xfncout(:)
              end if
              
           end do !it
           
           !---Correlation forecast and obs. error
           if(fo_cor)then
              st=t_fo_obs(1)
              et=t_fo_obs(nt_fo)
              call correlation_ef_eo(nt_fo,nobs,xferr(:,st:et),obserr(:,st:et),corf)
              call correlation_ef_eo(nt_fo,nobs,xfcerr(:,st:et),obserrc(:,st:et),corfc)
              call correlation_ef_eo(nt_fo,nobs,xfncerr(:,st:et),obserrnc(:,st:et),corfnc)
           end if
           
           !---Time averaged Bias & RMSE
           if(fo_cor .or. so_cor)then
              st=t_fo_obs(1)
              et=t_fo_obs(nt_fo)
              !Bias
              call ave_std(nt_fo,bias(st:et),ave_bias,std_bias)
              call ave_std(nt_fo,biasf(st:et),ave_biasf,std_biasf)
              call ave_std(nt_fo,biasa(st:et),ave_biasa,std_biasa)
              call ave_std(nt_fo,biaso(st:et),ave_biaso,std_biaso)
              call ave_std(nt_fo,biasfc(st:et),ave_biasfc,std_biasfc)
              call ave_std(nt_fo,biasac(st:et),ave_biasac,std_biasac)
              call ave_std(nt_fo,biasoc(st:et),ave_biasoc,std_biasoc)
              call ave_std(nt_fo,biasfnc(st:et),ave_biasfnc,std_biasfnc)
              call ave_std(nt_fo,biasanc(st:et),ave_biasanc,std_biasanc)
              call ave_std(nt_fo,biasonc(st:et),ave_biasonc,std_biasonc)
              !RMSE
              call ave_std(nt_fo,rmse(st:et),ave_rmse,std_rmse)
              call ave_std(nt_fo,rmsef(st:et),ave_rmsef,std_rmsef)
              call ave_std(nt_fo,rmsea(st:et),ave_rmsea,std_rmsea)
              call ave_std(nt_fo,rmseo(st:et),ave_rmseo,std_rmseo)
              call ave_std(nt_fo,rmsefc(st:et),ave_rmsefc,std_rmsefc)
              call ave_std(nt_fo,rmseac(st:et),ave_rmseac,std_rmseac)
              call ave_std(nt_fo,rmseoc(st:et),ave_rmseoc,std_rmseoc)
              call ave_std(nt_fo,rmsefnc(st:et),ave_rmsefnc,std_rmsefnc)
              call ave_std(nt_fo,rmseanc(st:et),ave_rmseanc,std_rmseanc)
              call ave_std(nt_fo,rmseonc(st:et),ave_rmseonc,std_rmseonc)
              !Pf,Pa
              call ave_std(nt_fo,Pfd(st:et),Pfdave,Pfdstd)
              call ave_std(nt_fo,Pfcd(st:et),Pfcdave,Pfcdstd)
              call ave_std(nt_fo,Pfncd(st:et),Pfncdave,Pfncdstd)              
              call ave_std(nt_fo,Pad(st:et),Padave,Padstd)
              call ave_std(nt_fo,Pacd(st:et),Pacdave,Pacdstd)
              call ave_std(nt_fo,Pancd(st:et),Pancdave,Pancdstd)              
           else
              call ave_std(nt_obs,bias,ave_bias,std_bias)
              call ave_std(nt_obs,biasf,ave_biasf,std_biasf)
              call ave_std(nt_obs,biasa,ave_biasa,std_biasa)
              call ave_std(nt_obs,biaso,ave_biaso,std_biaso)
              call ave_std(nt_obs,rmse,ave_rmse,std_rmse)
              call ave_std(nt_obs,rmsef,ave_rmsef,std_rmsef)
              call ave_std(nt_obs,rmsea,ave_rmsea,std_rmsea)
              call ave_std(nt_obs,rmseo,ave_rmseo,std_rmseo)
           end if

           write(*,'(a,f12.5)') "Inflation parameter:",inf_parm(iinf)
           write(*,'(a)') "Bias(Ave.,Std.) RMSE(Ave.,Std.)"
           write(*,'(a,4(f12.5,a))') "Free run: Bias(",ave_bias,",",std_bias,") RMSE(",ave_rmse,",",std_rmse,")"
           write(*,'(a,4(f12.5,a))') "Forecast: Bias(",ave_biasf,",",std_biasf,") RMSE(",ave_rmsef,",",std_rmsef,")"
           write(*,'(a,4(f12.5,a))') "Analysis: Bias(",ave_biasa,",",std_biasa,") RMSE(",ave_rmsea,",",std_rmsea,")"
           write(*,'(a,4(f12.5,a))') "Observation: Bias(",ave_biaso,",",std_biaso,") RMSE(",ave_rmseo,",",std_rmseo,")"
           if(fo_cor .or. so_cor)then
              write(*,'(a,4(f12.5,a))') "Cor. Forecast: Bias(",ave_biasfc,",",std_biasfc,") RMSE(",ave_rmsefc,",",std_rmsefc,")"
              write(*,'(a,4(f12.5,a))') "Cor. Analysis: Bias(",ave_biasac,",",std_biasac,") RMSE(",ave_rmseac,",",std_rmseac,")"
              write(*,'(a,4(f12.5,a))') "Observation: Bias(",ave_biasoc,",",std_biasoc,") RMSE(",ave_rmseoc,",",std_rmseoc,")"
              write(*,'(a,4(f12.5,a))') "Cor. Forecast: Bias(",ave_biasfnc,",",std_biasfnc,") RMSE(",ave_rmsefnc,",",std_rmsefnc,")"
              write(*,'(a,4(f12.5,a))') "Cor. Analysis: Bias(",ave_biasanc,",",std_biasanc,") RMSE(",ave_rmseanc,",",std_rmseanc,")"
              write(*,'(a,4(f12.5,a))') "Observation: Bias(",ave_biasonc,",",std_biasonc,") RMSE(",ave_rmseonc,",",std_rmseonc,")"
           endif

           if(fo_cor)then
              call write_ave_bias_rmse_fo(nobs,iinf,ifo,sfo,inf_parm(iinf), &
                   & ave_biasf,ave_biasa,ave_biaso,ave_rmsef,ave_rmsea,ave_rmseo, &
                   & ave_biasfc,ave_biasac,ave_biasoc,ave_rmsefc,ave_rmseac,ave_rmseoc, &
                   & ave_biasfnc,ave_biasanc,ave_biasonc,ave_rmsefnc,ave_rmseanc,ave_rmseonc)
              call write_cor_matrix(nobs,iinf,ifo,sfo,inf_parm(iinf),corf,corfc,corfnc)
              call write_ave_Pf_Pa(nobs,iinf,ifo,sfo,inf_parm(iinf), &
                   & Pfdave,Pfdstd,Padave,Padstd, &
                   & Pfcdave,Pfcdstd,Pacdave,Pacdstd, &
                   & Pfncdave,Pfncdstd,Pancdave,Pancdstd)                   
           else
              call write_ave_bias_rmse(nobs,iinf,ave_bias,ave_biasf,ave_biasa,ave_biaso, &
                   & ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo,inf_parm(iinf))
           end if
           
        end do !ifo
     end do !iinf

     deallocate(obs,R,H,HT,C)
     deallocate(Hobserr)
     if(fo_cor .or. so_cor)then
        deallocate(obsc,obsnc,Cc,Cnc)
        deallocate(Hobserrc,Hobserrnc)
        deallocate(corf,corfc,corfnc)
     end if
     
  end do !iobs
        
  write(*,*) "===End DA run==="
  
  !status=system("csh fig.csh "//trim(dir(execute_da)))
  
  deallocate(bias,biasf,biasa,biaso)
  deallocate(rmse,rmsef,rmsea,rmseo)
  deallocate(obserr)

  if(fo_cor .or. so_cor)then
     deallocate(t_fo_obs)
     deallocate(xferr,xfcerr,xfncerr)
     deallocate(obserrc,obserrnc)
     deallocate(q,qc)
     deallocate(biasfc,biasac,biasoc)
     deallocate(biasfnc,biasanc,biasonc)
     deallocate(rmsefc,rmseac,rmseoc)
     deallocate(rmsefnc,rmseanc,rmseonc)
     deallocate(Pfd,Pad)
     deallocate(Pfcd,Pacd)
     deallocate(Pfncd,Pancd)
  end if
  
end program main
