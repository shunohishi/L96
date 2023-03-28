program main

  use parameter
  implicit none

  !---Common
  integer it,nt
  integer ix
  integer imult
  integer iday
  integer status,system
  integer flag
  
  !Obs.
  integer iobs,nobs
  integer it_obs,nt_obs
  integer st_obs,et_obs

  !FO
  integer it_fo,nt_fo
  integer st_fo,et_fo
  integer ifo
  
  !IO
  integer it_io
  
  !---Variable
  double precision xtin(nx),xin(nx),xfin(nx),xain(nx)
  double precision xtout(nx),xout(nx),xfout(nx),xaout(nx)
  double precision,allocatable :: obs(:)              !Obs. in obs. space
  double precision,allocatable :: obserr(:,:)         !Obs. error in model space
  
  !---Error covariance matrix
  double precision Pf(nx,nx) !Forecast error covariance matrix: Pf
  double precision Pa(nx,nx) !Analysis error covariance matrix: Pa
  double precision,allocatable :: R(:,:)         !Observation error covarinance matrix: R
  double precision,allocatable :: H(:,:),HT(:,:) !Observation operator
  double precision mult_parm(nmult)  
  
  !---Static variable
  double precision,allocatable :: bias(:),biasf(:),biasa(:),biaso(:) !Bias
  double precision,allocatable :: rmse(:),rmsef(:),rmsea(:),rmseo(:) !RMSE
  double precision ave_bias,ave_biasf,ave_biasa,ave_biaso
  double precision std_bias,std_biasf,std_biasa,std_biaso
  double precision ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo
  double precision std_rmse,std_rmsef,std_rmsea,std_rmseo

  !---Forecaset and observation error correlation
  double precision fo_cc                              !Correlation between forecast and obs. error
  double precision xacout(nx),xacdout(nx),xac0out(nx)
  double precision,allocatable :: xferr(:,:)          !Forecast error
  double precision,allocatable :: obserr_random(:,:)  !Independent observation error in model space
  double precision,allocatable :: Hobserr_random(:,:) !                              in obs. space 
  double precision,allocatable :: obs_cor(:)         !Correlated obs.
  double precision,allocatable :: obserr_cor(:,:)    !Correlated obs. error
  double precision,allocatable :: Hobserr_cor(:,:)   !Obs. error in obs. space
  double precision avef,aveo,stdf,stdo
  double precision,allocatable :: cor(:,:)           !Forecaset and obs. error correlation
  double precision Pftmp(nx,nx)
  logical lC                                         !Switch C=<ef(eo)^T>
  logical C_diag                                     !diagonal/off-diagonal component 
  double precision,allocatable :: biasac(:),biasacd(:),biasac0(:),biasoc(:)
  double precision,allocatable :: rmseac(:),rmseacd(:),rmseac0(:),rmseoc(:)
  double precision ave_biasac,ave_biasacd,ave_biasac0,ave_biasoc
  double precision std_biasac,std_biasacd,std_biasac0,std_biasoc
  double precision ave_rmseac,ave_rmseacd,ave_rmseac0,ave_rmseoc
  double precision std_rmseac,std_rmseacd,std_rmseac0,std_rmseoc
  
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
     nt_fo=(et_fo-st_fo)/nint(int_obs/dt)+1
     write(*,*) "FO cor. time:",st_fo,"-",et_fo
     write(*,*) "Number of FO time:",nt_fo
     allocate(xferr(nx,nt_fo))
     allocate(obserr_random(nx,nt_fo),obserr_cor(nx,nt_fo))
     allocate(biasacd(nt_fo),biasac(nt_fo),biasac0(nt_fo),biasoc(nt_fo))
     allocate(rmseacd(nt_fo),rmseac(nt_fo),rmseac0(nt_fo),rmseoc(nt_fo))
  end if
  
  !---Random number
  !Initial condition
  call make_gaussian_random(0,noise_mean_da,noise_sd_da,nx,noise)
  !Obs.
  do ix=1,nx
     call make_gaussian_random(ix,0.d0,err_obs,nt_obs,obserr(ix,:))
  end do
  
  !---Multiplicative inflation
  call make_mult(mult_parm)

  !---Start do loop
!  do nobs=1,nx
  do nobs=nx,nx

     allocate(obs(nobs),R(nobs,nobs),H(nobs,nx),HT(nx,nobs))
     allocate(cor(nx,nobs)) !for fo
     cor(:,:)=0.d0          !original KF
     
!     do imult=2,nmult
     do imult=2,2
        
        do it=1,nt

           if(mod(it,1000) == 0 .or. mod(it,nt) == 0) write(*,*) nobs,"/",nx,imult,"/",nmult,it,"/",nt
           
           !Initialization
           if(it == 1)then
              it_obs=0
              it_io=0
              flag=0
              if(fo_cor) it_fo=0
              !True run
              call initialization(amp_t,deg_t,xnull,xtin) !noise=xnull=0.d0
              !Free run
              call initialization(amp_da,deg_da,noise,xin)
              call initialization(amp_da,deg_da,noise,xfin)
              call write_data(nobs,imult,it_io,0,xtin,xin,xfin,xnull,xnull,mult_parm(imult))
              call write_P(nobs,imult,it_io,0,Pnull,Pnull,mult_parm(imult))
           end if

           !Forwading Lorenz model
           call forward_lorenz96(xtin,xtout)
           call forward_lorenz96(xin,xout)
           call forward_lorenz96(xfin,xfout)

           call check_nan(xtout,xout,xfout,flag)
           
           !Lyapunov equation:Pf(t)=MPa(t-1)MT 
           if(it == st_obs) call initialization_ECM(nobs,R,Pf)
           if(st_obs <= it .and. it <= et_obs .and. flag == 0)then
              call lyapunov_equation(xfout,Pf)
           end if

           !Kalman Filter
           if(mod(it,nint(int_obs/dt)) == 0)then
              
              call obs_operator(it,nobs,H,HT)
              
              if(st_obs <= it .and. it <= et_obs .and. flag == 0)then

                 it_obs=it_obs+1

                 !Make obs.
                 obs(:)=matmul(H(:,:),xtout(:)+obserr(:,it_obs))

                 !Kalman filter
                 lC=.false.
                 C_diag=.false.
                 Pftmp(:,:)=Pf(:,:)
                 call kalman_filter(nobs,xfout,xaout,Pftmp,Pa,obs,R,H,HT,mult_parm(imult),cor,lC,C_diag)
                 
                 !Bias RMSE
                 call bias_rmse(nx,xout,xtout,bias(it_obs),rmse(it_obs))
                 call bias_rmse(nx,xfout,xtout,biasf(it_obs),rmsef(it_obs))
                 call bias_rmse(nx,xaout,xtout,biasa(it_obs),rmsea(it_obs))
                 call bias_rmse(nobs,obs,matmul(H(:,:),xtout(:)),biaso(it_obs),rmseo(it_obs))

                 !Write data
                 call write_bias_rmse(nobs,imult,st_obs,it,bias(it_obs),biasf(it_obs),biasa(it_obs),biaso(it_obs), &
                      & rmse(it_obs),rmsef(it_obs),rmsea(it_obs),rmseo(it_obs),mult_parm(imult))

                 !Conserve forecast error
                 if(fo_cor .and. st_fo <= it .and. it <= et_fo)then
                    it_fo=it_fo+1
                    xferr(:,it_fo)=xfout(:)-xtout(:)
                 end if

                 
                 !Replace data
                 xfout(:)=xaout(:)
                 Pf(:,:)=Pa(:,:)

              else
                 
                 obs(:)=0.d0
                 xaout(:)=0.d0
                 Pa(:,:)=0.d0
                 
              end if

              it_io=it_io+1           
              call write_data(nobs,imult,it_io,it,xtout,xout,xfout,xaout,matmul(HT(:,:),obs(:)),mult_parm(imult))
              call write_P(nobs,imult,it_io,it,Pf,Pa,mult_parm(imult))
        
           end if
           
           xtin(:)=xtout(:)
           xin(:)=xout(:)
           xfin(:)=xfout(:)
           
        end do !it

        call ave_std(nt_obs,bias,ave_bias,std_bias)
        call ave_std(nt_obs,biasf,ave_biasf,std_biasf)
        call ave_std(nt_obs,biasa,ave_biasa,std_biasa)
        call ave_std(nt_obs,biaso,ave_biaso,std_biaso)
        call ave_std(nt_obs,rmse,ave_rmse,std_rmse)
        call ave_std(nt_obs,rmsef,ave_rmsef,std_rmsef)
        call ave_std(nt_obs,rmsea,ave_rmsea,std_rmsea)
        call ave_std(nt_obs,rmseo,ave_rmseo,std_rmseo)

        write(*,'(a)') "Bias(Ave.,Std.) RMSE(Ave.,Std.)"
        write(*,'(a,4(f12.5,a))') "Free run: Bias(",ave_bias,",",std_bias,") RMSE(",ave_rmse,",",std_rmse,")"
        write(*,'(a,4(f12.5,a))') "Forecast: Bias(",ave_biasf,",",std_biasf,") RMSE(",ave_rmsef,",",std_rmsef,")"
        write(*,'(a,4(f12.5,a))') "Analysis: Bias(",ave_biasa,",",std_biasa,") RMSE(",ave_rmsea,",",std_rmsea,")"
        write(*,'(a,4(f12.5,a))') "Observation: Bias(",ave_biaso,",",std_biaso,") RMSE(",ave_rmseo,",",std_rmseo,")"

        call write_ave_bias_rmse(nobs,imult,ave_bias,ave_biasf,ave_biasa,ave_biaso, &
             & ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo,mult_parm(imult))
        
        !____________________________________________________________________________________________________________________________________
        !Correlated forecast and obs. error exp.
        if(fo_cor)then

           do ifo=0,10

              !Set correlation value
              fo_cc=0.1d0*dble(ifo)

              allocate(obs_cor(nobs))
              allocate(Hobserr_random(nobs,nt_fo),Hobserr_cor(nobs,nt_fo))

              !Prepare correlated obs. error
              do ix=1,nx
                 obserr_random(ix,1:nt_fo)=obserr(ix,(st_fo-st_obs)/nint(int_obs/dt)+1:(et_fo-st_obs)/nint(int_obs/dt)+1)
              end do

              do ix=1,nx
                 call ave_std(nt_fo,xferr(ix,:),avef,stdf)
                 call ave_std(nt_fo,obserr_random(ix,:),aveo,stdo)
                 call make_cor_random(nt_fo,stdf,stdo,err_obs,fo_cc,xferr(ix,:),obserr_random(ix,:),obserr_cor(ix,:))
              end do

              do it=1,nt

                 if(mod(it,1000) == 0 .or. mod(it,nt) == 0) &
                      & write(*,*) nobs,"/",nx,imult,"/",nmult,ifo,"/",10,it,"/",nt

                 !Initialization
                 if(it == 1)then
                    it_obs=0
                    it_fo=0
                    it_io=0
                    flag=0
                    !True run
                    call initialization(amp_t,deg_t,xnull,xtin) !noise=xnull=0.d0
                    !Free run
                    call initialization(amp_da,deg_da,noise,xin)
                    call initialization(amp_da,deg_da,noise,xfin)
                 end if

                 !Forwading Lorenz model
                 call forward_lorenz96(xtin,xtout)
                 call forward_lorenz96(xin,xout)
                 call forward_lorenz96(xfin,xfout)

                 call check_nan(xtout,xout,xfout,flag)

                 !Lyapunov equation:Pf(t)=MPa(t-1)MT 
                 if(it == st_obs) call initialization_ECM(nobs,R,Pf)
                 if(st_obs <= it .and. it <= et_obs .and. flag == 0)then
                    call lyapunov_equation(xfout,Pf)
                 end if

                 !Kalman Filter
                 if(mod(it,nint(int_obs/dt)) == 0)then

                    call obs_operator(it,nobs,H,HT)

                    if(st_obs <= it .and. it <= et_obs .and. flag == 0)then

                       it_obs=it_obs+1

                       !---Original KF using independent obs. error
                       lC=.false.
                       C_diag=.false.
                       obs(:)=matmul(H(:,:),xtout(:)+obserr(:,it_obs))
                       Pftmp(:,:)=Pf(:,:)
                       call kalman_filter(nobs,xfout,xaout,Pftmp,Pa,obs,R,H,HT,mult_parm(imult),cor,lC,C_diag)

                       !Bias RMSE
                       call bias_rmse(nx,xout,xtout,bias(it_obs),rmse(it_obs))
                       call bias_rmse(nx,xfout,xtout,biasf(it_obs),rmsef(it_obs))
                       call bias_rmse(nx,xaout,xtout,biasa(it_obs),rmsea(it_obs))                    
                       call bias_rmse(nobs,obs,matmul(H(:,:),xtout(:)),biaso(it_obs),rmseo(it_obs))

                       !Write data
                       call write_bias_rmse(nobs,imult,st_obs,it,bias(it_obs),biasf(it_obs),biasa(it_obs),biaso(it_obs), &
                            & rmse(it_obs),rmsef(it_obs),rmsea(it_obs),rmseo(it_obs),mult_parm(imult))

                       !Replace data
                       xfout(:)=xaout(:)
                       Pf(:,:)=Pa(:,:)

                    else

                       obs(:)=0.d0
                       xaout(:)=0.d0
                       Pa(:,:)=0.d0

                    end if !obs

                    if(st_fo <= it .and. it <= et_fo .and. flag == 0)then

                       it_fo=it_fo+1

                       !Correlation
                       if(nobs == nx)then
                          if(it_fo == 1)then
                             Hobserr_cor(:,:)=obserr_cor(:,:)
                             do iobs=1,nobs
                                do ix=1,nx
                                   call correlation(nt_fo,xferr(ix,:),Hobserr_cor(iobs,:),cor(ix,iobs))
                                end do
                             end do
                          else
                             cor(:,:)=cor(:,:)
                          end if
                       else
                          Hobserr_cor(:,:)=matmul(H(:,:),obserr_cor(:,:)) !H: nobs*nx, obserr_cor:nx_nt_obs
                          do iobs=1,nobs
                             do ix=1,nx
                                call correlation(nt_fo,xferr(ix,:),Hobserr_cor(iobs,:),cor(ix,iobs))
                             end do
                          end do
                       end if

                       !---Correlated KF
                       lC=.true.
                       obs_cor(:)=matmul(H(:,:),xtout(:)+obserr_cor(:,it_fo))
                       C_diag=.true.
                       Pftmp(:,:)=Pf(:,:)
                       call kalman_filter(nobs,xfout,xacdout,Pftmp,Pa,obs_cor,R,H,HT,mult_parm(imult),cor,lC,C_diag)
                       C_diag=.false.
                       Pftmp(:,:)=Pf(:,:)
                       call kalman_filter(nobs,xfout,xacout,Pftmp,Pa,obs_cor,R,H,HT,mult_parm(imult),cor,lC,C_diag)

                       !---Original KF using correlated obs. error
                       lC=.false.
                       c_diag=.true.
                       Pftmp(:,:)=Pf(:,:)
                       call kalman_filter(nobs,xfout,xac0out,Pftmp,Pa,obs_cor,R,H,HT,mult_parm(imult),cor,lC,C_diag)

                       call bias_rmse(nx,xacdout,xtout,biasacd(it_fo),rmseacd(it_fo))
                       call bias_rmse(nx,xacout,xtout,biasac(it_fo),rmseac(it_fo))
                       call bias_rmse(nx,xac0out,xtout,biasac0(it_fo),rmseac0(it_fo))                    
                       call bias_rmse(nobs,obs_cor,matmul(H(:,:),xtout(:)),biasoc(it_fo),rmseoc(it_fo))

                    end if !fo

                    it_io=it_io+1           
                    call write_data(nobs,imult,it_io,it,xtout,xout,xfout,xaout,matmul(HT(:,:),obs(:)),mult_parm(imult))
                    call write_P(nobs,imult,it_io,it,Pf,Pa,mult_parm(imult))

                 end if

                 xtin(:)=xtout(:)
                 xin(:)=xout(:)
                 xfin(:)=xfout(:)

              end do !it

              deallocate(obs_cor)
              deallocate(Hobserr_random,Hobserr_cor)

              call ave_std(nt_obs,bias,ave_bias,std_bias)
              call ave_std(nt_obs,biasf,ave_biasf,std_biasf)
              call ave_std(nt_obs,biasa,ave_biasa,std_biasa)
              call ave_std(nt_obs,biaso,ave_biaso,std_biaso)
              call ave_std(nt_obs,rmse,ave_rmse,std_rmse)
              call ave_std(nt_obs,rmsef,ave_rmsef,std_rmsef)
              call ave_std(nt_obs,rmsea,ave_rmsea,std_rmsea)
              call ave_std(nt_obs,rmseo,ave_rmseo,std_rmseo)
              
              !Monitor
              write(*,'(a)') "Bias(Ave.,Std.) RMSE(Ave.,Std.)"
              write(*,'(a,4(f12.5,a))') "Free run: Bias(",ave_bias,",",std_bias,") RMSE(",ave_rmse,",",std_rmse,")"
              write(*,'(a,4(f12.5,a))') "Forecast: Bias(",ave_biasf,",",std_biasf,") RMSE(",ave_rmsef,",",std_rmsef,")"
              write(*,'(a,4(f12.5,a))') "Analysis (KF): Bias(",ave_biasa,",",std_biasa,") RMSE(",ave_rmsea,",",std_rmsea,")"
              write(*,'(a,4(f12.5,a))') "Observation: Bias(",ave_biaso,",",std_biaso,") RMSE(",ave_rmseo,",",std_rmseo,")"
              
              call ave_std(nt_fo,biasacd,ave_biasacd,std_biasacd)
              call ave_std(nt_fo,biasac,ave_biasac,std_biasac)
              call ave_std(nt_fo,biasac0,ave_biasac0,std_biasac0)
              call ave_std(nt_fo,biasoc,ave_biasoc,std_biasoc)
              call ave_std(nt_fo,rmseacd,ave_rmseacd,std_rmseacd)
              call ave_std(nt_fo,rmseac,ave_rmseac,std_rmseac)
              call ave_std(nt_fo,rmseac0,ave_rmseac0,std_rmseac0)
              call ave_std(nt_fo,rmseoc,ave_rmseoc,std_rmseoc)

              write(*,'(a,4(f12.5,a))') &
                   & "Analysis:(CKF *C=0) Bias(",ave_biasac0,",",std_biasac0,") RMSE(",ave_rmseac0,",",std_rmseac0,")"
              write(*,'(a,4(f12.5,a))') &
                   & "Analysis:(CKF *C=diag)  Bias(",ave_biasacd,",",std_biasacd,") RMSE(",ave_rmseacd,",",std_rmseacd,")"
              write(*,'(a,4(f12.5,a))') &
                   & "Analysis:(CKF *C=full)  Bias(",ave_biasac,",",std_biasac,") RMSE(",ave_rmseac,",",std_rmseac,")"
              write(*,'(a,4(f12.5,a))') &
                   & "Correlated Obs.: Bias(",ave_biasoc,",",std_biasoc,") RMSE(",ave_rmseoc,",",std_rmseoc,")"

              write(*,'(a,5f12.5)') "Cor:",cor(1,1:5)
              write(*,'(a,5f12.5)') "Cor:",cor(2,1:5)
              write(*,'(a,5f12.5)') "Cor:",cor(3,1:5)
              write(*,'(a,5f12.5)') "Cor:",cor(4,1:5)
              write(*,'(a,5f12.5)') "Cor:",cor(5,1:5)

              call write_ave_bias_rmse_fo(nobs,imult,ifo, &
                   & ave_bias,ave_biasf,ave_biasa,ave_biaso, &
                   & ave_biasacd,ave_biasac,ave_biasac0,ave_biasoc, &
                   & ave_rmse,ave_rmsef,ave_rmsea,ave_rmseo, &
                   & ave_rmseacd,ave_rmseac,ave_rmseac0,ave_rmseoc, &
                   & mult_parm(imult),fo_cc)
              
           end do !ifo
              
        end if !fo_cor
                                
     end do !imult
     
     deallocate(obs,R,H,HT)
     deallocate(cor)
     
  end do !nobs
  !---End do loop
  
  write(*,*) "===End DA run==="
  
  !status=system("csh fig.csh "//trim(dir(execute_da)))
  
  deallocate(bias,biasf,biasa,biaso)
  deallocate(rmse,rmsef,rmsea,rmseo)
  deallocate(biasacd,biasac,biasac0,biasoc)
  deallocate(rmseacd,rmseac,rmseac0,rmseoc)
  deallocate(xferr)
  deallocate(obserr_random,obserr_cor)
  
end program main
