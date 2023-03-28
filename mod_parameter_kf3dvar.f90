module parameter

  !---Pi
  double precision,parameter :: pi=4.d0*atan(1.d0)

  !---Initialization
  integer,parameter :: nwave = 4       !Wave number
  double precision,parameter :: amp_t=1.d0      !Infinitesimal amplitude in true run
  double precision,parameter :: deg_t=0.d0      !              degree in true run
  double precision,parameter :: amp_da=0.d0     !              amplitude in DA run
  double precision,parameter :: deg_da=0.d0     !              degree in DA run
  double precision,parameter :: noise_mean_da=0.d0 !Noise mean in DA run
  double precision,parameter :: noise_sd_da=1.d0   !      standard deviation in DA run
  
  !---Setting
  integer,parameter :: nx=40          !Number of variable(x)
  integer,parameter :: nday=365*12     !Number of time [day]
  double precision,parameter :: dt=0.01d0      !Time interval
  double precision,parameter :: int_day=0.2d0  !1day time interval
  double precision,parameter :: force=8.0d0    !External forcing
  
  !---Error Growth
  integer,parameter :: yyyys_eg=1,dds_eg=0,hhs_eg=0  !Start date
  integer,parameter :: yyyye_eg=2,dde_eg=0,hhe_eg=0  !End date
  integer,parameter :: ntime_eg=30                   !Number of time
  integer,parameter :: nrand_eg=100                  !Number of random for error growth
  double precision,parameter :: nperiod_eg=20.d0     !Integration period [day]
  double precision,parameter :: noise_mean_eg=0.d0   !Noise mean
  double precision,parameter :: noise_sd_eg=1.d-2    !Noise standard deviation
  
  !---Data assimilation
  ![1:Kalman Filter, 2:3D-Var, 3:ETKF, 4:LETKF, 5:EnSRF, 6:PO, 7:PF, 8:LPF]
  integer,parameter :: execute_da=1
  integer,parameter :: execute_pf=1                             !1: SU, 2: MN, 3: Residual
  integer,parameter :: nens_exp=1                              !Ensemble experiment
  character(5),parameter :: dir(8)=(/"KF   ","3DVAR","ETKF ","LETKF","EnSRF","PO   ","PF   ","LPF  "/)
  double precision,parameter :: err_fct=5.0d0                   !Initial Forecast error[-]
  double precision,parameter :: dx_fct=1.d-3                    !For tangent linear model
  !Inflation
  integer,parameter :: ninf=19                                 !Number of Inflation (***See make_inf***)
!  integer,parameter :: ninf=100                                 !Number of Inflation (***See make_inf***)
  logical,parameter :: mult_inf=.false.                         !Multiplicative inflation
  double precision,parameter :: mult_amp=1.d-4                  !Mult. amplitude
  logical,parameter :: add_inf=.true.                           !Additive inflation
  double precision,parameter :: add_amp=10.d-5                   !Q amplitude
!  double precision,parameter :: add_amp=10.d-4                   !Q amplitude
!  double precision,parameter :: add_amp=10.d-3                   !Q amplitude
!  double precision,parameter :: add_amp=10.d-2                   !Q amplitude
!  double precision,parameter :: add_amp=10.d-1                   !Q amplitude
  !Localization
  logical,parameter :: Pf_loc=.false.                           !Pf Localization for EnSRF/PO
  logical,parameter :: K_loc=.false.                            !K localization for EnSRF
  !EnSRF
  logical,parameter :: SRF_serial=.true.                        !Serial EnSRF or not
  
  !---Observation
  integer,parameter :: yyyys_obs=1,dds_obs=1,hhs_obs=0    !Start date
  integer,parameter :: yyyye_obs=12,dde_obs=0,hhe_obs=0    !End date
  double precision,parameter :: int_obs=0.2d0*6.d0/24.d0  !Time interval
  logical,parameter :: diag_err_obs=.true.
  double precision,parameter :: err_obs=1.d0              !Observation error [-]

  !---Forecast vs. Obs. error correlation (add_inf --> true)
  logical,parameter :: fo_cor=.true.                            !Forecast * Obs correlation for KF
!  logical,parameter :: fo_cor=.false.                          !Forecast * Obs correlation for KF
  integer,parameter :: yyyys_fo=2,dds_fo=1,hhs_fo=0             !Start date
  integer,parameter :: yyyye_fo=12,dde_fo=0,hhe_fo=0             !End date

  !System noise vs. Obs. (add_inf --> true)
  logical,parameter :: so_cor=.false.
  integer,parameter :: yyyys_so=2,dds_so=1,hhs_so=0             !Start date
  integer,parameter :: yyyye_so=12,dde_so=0,hhe_so=0             !End date
  
end module parameter
