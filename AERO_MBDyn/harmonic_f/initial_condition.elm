# harmonic impossed movement
#driven: CURR_BLADE + BLADE_FLAP_OFFSET+100, reference, IC,
/*
joint: CURR_BLADE + BLADE_FLAP_OFFSET+100, total joint,
  HUB,
    position, reference, CURR_BLADE + BLADE_FLAP_OFFSET, null,
    position orientation, reference, CURR_BLADE + BLADE_FLAP_OFFSET, eye,
    rotation orientation, reference, CURR_BLADE, eye,
  CURR_BLADE,
    position, reference, CURR_BLADE + BLADE_FLAP_OFFSET, null,
    position orientation, reference, CURR_BLADE + BLADE_FLAP_OFFSET, eye,
    rotation orientation, reference, CURR_BLADE + BLADE_FLAP_OFFSET, eye,
  position constraint, 0, 0, 0, null,
  orientation constraint, 0, 0, 1,
    0., 0., 1.,
      #const,XID_0+XI_0/T_REV;
      fourier series,
        dt,
        OMEGA_100,
        1,
          0.01,
          0.005, 0.,
          0.002, 0.,
          0.0005, 0.,
          0.0008, 0.,
          0.006, 0.,
          0.0009, 0.,
          0.0002, 0.,
          forever,
        0.;
      element, HARMONIC_EXCITATION, loadable, string, "psi", 
          string, "0.001*sin(Var+HARM*(real_eval(CURR_BLADE)-1)*DELTAPSI)";
          #string, "0.001*sin(Var+DELTAPSI)";
      #step,T_REV-dt,-XI_0/T_REV,XID_0+XI_0/T_REV;
      #string, "real_eval(XI_0D+XI_0/T_REV*(1-step(Time-T_REV+0.01)))";
      */


couple:CURR_BLADE + BLADE_FLAP_OFFSET+100, absolute,
  CURR_BLADE,
  position , null,
  0.,0.,1., 
    fourier series,
        -(BLADE-1)*DELTAPSI/OMEGA_100,
        OMEGA_100,
        HARMONIC,
          2000.,
          2400., 0.,
          3400, 0.,
          1500., 0.,
          1000., 0.,
          2200., 0.,
          800., 0.,
          250., 0.,
          forever,
        0.;
  
  /* element, HARMONIC_EXCITATION, loadable, string, "psi", 
          string, "0.001*sin( Var+HARM*(real_eval(CURR_BLADE)-1)*DELTAPSI)";
  */

