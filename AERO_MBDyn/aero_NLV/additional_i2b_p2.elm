# Once all bars are defined the viscoelastic cylindrical damper is defined

# Spherical bearing (lead-lag visco-elastic connector)

joint: CURR_BLADE + BLADE_FLAP_OFFSET+1, deformable hinge,
  CURR_BLADE+D_NOD,
    position, reference, CURR_BLADE + C_NOD, null,
    orientation, reference, CURR_BLADE + D_NOD,eye,
  CURR_BLADE+F_NOD,
    position, reference, CURR_BLADE + C_NOD, null,
    orientation, reference, CURR_BLADE + F_NOD,eye,
  reference, CURR_CL_BLADE_ROOT;

#           o  Di+1
#            \Ci+1
#     o-------o====o
#     |      / 3   | Ei+1
#     | Di  o Fi+1 |  7
#  Fo | o --|------o Bi+1
#    \|/  5 |     6
#     o Ci  | 2
#     ||    |
#     || 4  |
#     o----o
#     Ei 1  Ai

# The combination of a total hinge, whose only allowed movement is the lead-lag
# and the deformable hinge with only K_z and C_z creates the lead-lag damper
# union with the hub
