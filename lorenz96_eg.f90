program main

  use parameter
  implicit none

  !---Common
  integer it,nt
  integer ix
  integer iday
  integer status,system
    
  !---Variable
  double precision xtin(nx)
  double precision xtout(nx)
  
  !---Noise
  double precision noise(nx)

  !Error Growth
  integer irand_eg,itime_eg
  integer st0_eg,et0_eg,nt_eg
  integer it_eg(ntime_eg),st_eg(ntime_eg),et_eg(ntime_eg)
  
  double precision noise_eg(ntime_eg)
  double precision,allocatable :: xegin(:,:,:),xegout(:,:,:)
  double precision,allocatable :: biaseg(:,:,:),rmseeg(:,:,:)
  double precision ave_biaseg,ave_rmseeg
  double precision std_biaseg,std_rmseeg
    
  !---check nx
  if(nx <= 3)then
     write(*,*) "***Error: Choose nx > 3"
     stop
  end if
  
  !---Number of time
  nt=nint(dble(nday)*int_day/dt)
  write(*,*) "Number of time:",nt

  !--Make directory
  status=system("mkdir -p "//trim(dir(execute_da)))
  
  !Error Growth
  !Time
  st0_eg=nint((dble(yyyys_eg)*365.d0+dble(dds_eg)+dble(hhs_eg)/24.d0)*int_day/dt)
  et0_eg=nint((dble(yyyye_eg)*365.d0+dble(dde_eg)+dble(hhe_eg)/24.d0)*int_day/dt)
  nt_eg=nint(nperiod_eg*int_day/dt)
  write(*,*) "Start time for Growth rate:",st0_eg
  write(*,*) "End time for Growth rate:",et0_eg
  write(*,*) "Period for Growth rate:",nt_eg
  allocate(xegin(nx,nrand_eg,ntime_eg),xegout(nx,nrand_eg,ntime_eg))
  allocate(biaseg(nrand_eg,ntime_eg,nt_eg),rmseeg(nrand_eg,ntime_eg,nt_eg))
  call make_original_random(0,ntime_eg,noise_eg)
  do itime_eg=1,ntime_eg
     st_eg(itime_eg)=nint(st0_eg+(et0_eg-nt_eg-st0_eg)*noise_eg(itime_eg))
     et_eg(itime_eg)=st_eg(itime_eg)+nt_eg-1
  end do
  
  !---Nature run
  write(*,*) "===Start Growth rate run==="

  do it=1,nt
     
     !---Initialization
     if(it == 1)then
        noise(:)=0.d0
        call initialization(amp_t,deg_t,noise,xtin)
     end if

     do itime_eg=1,ntime_eg
        if(it == st_eg(itime_eg))then
           it_eg(itime_eg)=0
           do irand_eg=1,nrand_eg
              call make_random(irand_eg+itime_eg,noise_mean_eg,noise_sd_eg,nx,noise)              
              xegin(:,irand_eg,itime_eg)=xtin(:)+noise(:)
           end do
        end if
     end do
     
     !---Forwarding Lorenz model
     call forward_lorenz96(xtin,xtout)
     
     do itime_eg=1,ntime_eg
        
        if(st_eg(itime_eg) <= it .and. it <= et_eg(itime_eg))then
           it_eg(itime_eg)=it_eg(itime_eg)+1
           do irand_eg=1,nrand_eg
              !---Forwarding Lorenz model
              call forward_lorenz96(xegin(:,irand_eg,itime_eg),xegout(:,irand_eg,itime_eg))

              !---Bias & RMSE
              call bias_rmse(nx,xegout(:,irand_eg,itime_eg),xtout(:),&
                   & biaseg(irand_eg,itime_eg,it_eg(itime_eg)),rmseeg(irand_eg,itime_eg,it_eg(itime_eg)))
           end do
        end if
        
     end do
     
     !In <-- Out        
     xtin(:)=xtout(:)
     do itime_eg=1,ntime_eg
        if(st_eg(itime_eg) <= it .and. it <= et_eg(itime_eg))then
           xegin(:,:,itime_eg)=xegout(:,:,itime_eg)
        end if
     end do
     
  end do
  
  !Growth rate
  do it=1,nt_eg  
     call ave_std2(nrand_eg,ntime_eg,biaseg(:,:,it),ave_biaseg,std_biaseg)
     call ave_std2(nrand_eg,ntime_eg,rmseeg(:,:,it),ave_rmseeg,std_rmseeg)
     call write_growth_rate(it,ave_biaseg,std_biaseg,ave_rmseeg,std_rmseeg)
  end do
  
  deallocate(xegin,xegout)
  deallocate(biaseg,rmseeg)

  status=system("csh gr.csh")
  
  write(*,*) "===End Growth rate run==="

end program main
