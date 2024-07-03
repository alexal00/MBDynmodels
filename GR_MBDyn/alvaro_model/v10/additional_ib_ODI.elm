#set: A_BAR = 0.;
#set: B_BAR = 0.;
#set: C_BAR = X_BLADE_FLAP_OFFSET/4;
  
joint: CURR_BLADE + BLADE_FLAP_OFFSET+1, rod,
  CURR_BLADE,
    position, reference, CURR_BLADE + A_NOD, null,
  NXT_BLADE,
    position, reference, NXT_BLADE + B_NOD, null,
  from nodes,
  linear viscoelastic,
    0.,
    0.;
