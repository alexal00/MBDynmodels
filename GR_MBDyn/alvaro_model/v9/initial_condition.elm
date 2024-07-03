# Initial lead-lag angle
driven: CURR_BLADE + BLADE_FLAP_OFFSET+100, reference, IC,
/*
joint: CURR_BLADE + BLADE_FLAP_OFFSET+100, axial rotation,
  HUB,
    position, reference, CURR_BLADE + BLADE_FLAP_OFFSET, null,
    orientation, reference, CURR_BLADE, eye,
  CURR_BLADE,
    position, reference, CURR_BLADE + BLADE_FLAP_OFFSET, null,
    orientation, reference, CURR_BLADE + BLADE_FLAP_OFFSET, eye,
      #const,XID_0+XI_0/T_REV;
      step5,T_REV-6*dt,XID_0+XI_0/(T_REV-6*dt),T_REV-dt,XID_0;
      #step,T_REV-dt,-XI_0/T_REV,XID_0+XI_0/T_REV;
      #string, "real_eval(XI_0D+XI_0/T_REV*(1-step(Time-T_REV+0.01)))";
*/

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
  orientation constraint, 0, 0, angular velocity,
    0., 0., 1.,
      const,XID_0+XI_0/T_REV;
      #step5,T_REV-6*dt,XID_0+XI_0/(T_REV-6*dt),T_REV-dt,XID_0;
      #step,T_REV-dt,-XI_0/T_REV,XID_0+XI_0/T_REV;
      #string, "real_eval(XI_0D+XI_0/T_REV*(1-step(Time-T_REV+0.01)))";
  