
joint: CURR_BLADE + BLADE_FLAP_OFFSET+1, rod,
  CURR_BLADE,
    position, reference, CURR_BLADE + A_NOD, null,
  NXT_BLADE,
    position, reference, NXT_BLADE + B_NOD, null,
  from nodes,
  nlsf viscous,
  0.,		# Linear damping
  "nlvisc_rod";	# Hiperbolic tangent non-linear damping
