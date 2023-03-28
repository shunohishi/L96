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
  integer,parameter :: nday=365*10     !Number of time [day]
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
  ![1:Kalman Filter, 2:3D-Var, 3:ETKF, 4:LETKF, 5:EnSRF, 6:PO, 7:PF, 8:LPF, 9: ETKFCC, 10:ETKFNC]
  integer,parameter :: execute_da=10
  integer,parameter :: execute_pf=1                             !1: SU, 2: MN, 3: Residual
  integer,parameter :: nens_exp=1                               !Ensemble experiment
  character(6),parameter :: dir(10)=(/"KF    ","3DVAR ","ETKF  ","LETKF ","EnSRF ","PO    ","PF    ","LPF   ",&
       & "ETKFCC","ETKFNC"/)
  double precision,parameter :: err_fct=5.0d0                   !Initial Forecast error[-]
  double precision,parameter :: dx_fct=1.d-3                    !For tangent linear model
  !Inflation
  !integer,parameter :: ninf=1                                   !Number of Inflation (***See make_inf***)
  !integer,parameter :: ninf=11                                  !Number of Inflation (***See make_inf***)
  integer,parameter :: ninf=23                                  !Number of Inflation (***See make_inf***)
  logical,parameter :: mult_inf=.true.                         !Multiplicative inflation
  double precision,parameter :: mult_amp=1.d-2                  !Mult. amplitude
  logical,parameter :: add_inf=.false.                           !Additive inflation
  double precision,parameter :: add_amp=1.d-4                   !Q amplitude
  logical,parameter :: rtpp_inf=.false.                          !RTPP  
  double precision,parameter :: rtpp_intv=0.1d0                 !Interval
  logical,parameter :: rtps_inf=.false.                          !RTPP  
  double precision,parameter :: rtps_intv=0.1d0                 !Interval
  !Localization
  logical,parameter :: Pf_loc=.false.                           !Pf Localization for EnSRF/PO
  logical,parameter :: K_loc=.false.                            !K localization for EnSRF
  integer,parameter :: nsigma_loc=1
  !EnSRF
  logical,parameter :: SRF_serial=.false.                        !Serial EnSRF or not
  
  !---Observation
  integer,parameter :: yyyys_obs=1,dds_obs=1,hhs_obs=0    !Start date
  integer,parameter :: yyyye_obs=10,dde_obs=0,hhe_obs=0    !End date
  double precision,parameter :: int_obs=0.2d0*6.d0/24.d0  !Time interval
  logical,parameter :: diag_err_obs=.true.
  double precision,parameter :: err_obs=1.d0              !Observation error [-]

  !---Forecast vs. Obs. error correlation
  integer,parameter :: npco=20
  integer,parameter :: nday_spin=365 !Number of time for spin-up[day]
  !integer,parameter :: nR=2
  integer,parameter :: nR=18
  
end module parameter
