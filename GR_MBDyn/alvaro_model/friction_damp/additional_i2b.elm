# Sketch of the blade damping mechanism
#           o  Di+1
#            \Ci+1
#     o-------o====o      [INT_BLADE]
#     |      / 3   | Ei+1
#     | Di  o Fi+1 |  7
#  Fo | o --|------o Bi+1
#    \|/ 5  |   6
#     o Ci  | 
#     ||    | 2
#   4 ||    |
#     o-----o
#     Ei 1  Ai
#   [CURR_BLADE]

body: CURR_BLADE+AF_ROD, CURR_BLADE+AF_ROD,
  1e-5,
  reference,CURR_BLADE+AF_ROD,null,
  diag,1e-5,1e-5,1e-5;
body: CURR_BLADE+D_NOD, CURR_BLADE+D_NOD,
  1e-5,
  reference,CURR_BLADE+D_NOD,null,
  diag,1e-5,1e-5,1e-5;
body: INT_BLADE+BD_ROD, INT_BLADE+BD_ROD,
  1e-5,
  reference, INT_BLADE+BD_ROD,null,
  diag,1e-5,1e-5,1e-5;
body: INT_BLADE+F_NOD, INT_BLADE+F_NOD,
  1e-5,
  reference,INT_BLADE+F_NOD,null,
  diag,1e-5,1e-5,1e-5;

# 1.

joint: CURR_BLADE + AF_ROD, spherical hinge,
  CURR_BLADE,
    position, reference, CURR_BLADE + A_NOD , null,
    orientation, reference, CURR_BLADE + A_NOD, eye,
  CURR_BLADE+AF_ROD,
    position, reference, CURR_BLADE + A_NOD, null,
    orientation, reference, CURR_BLADE + A_NOD, eye;

/*
joint: CURR_BLADE + AF_ROD, total joint,
  CURR_BLADE,
    position, reference, CURR_BLADE + A_NOD , null,
    position orientation, reference, CURR_BLADE + A_NOD, eye,
    rotation orientation, reference, CURR_BLADE + A_NOD, eye,
  CURR_BLADE+AF_ROD,
    position, reference, CURR_BLADE + A_NOD, null,
    position orientation, reference, CURR_BLADE + A_NOD, eye,
    rotation orientation, reference, CURR_BLADE + A_NOD, eye,
  position constraint, 1, 1, 1,null,
  orientation constraint, 0, 0, 0, null;
*/
# 2.
/*
joint: INT_BLADE + F_NOD, revolute hinge,
  HUB,  # Static node C is substituted by an offset from the hub node
    position, reference, INT_BLADE + C_NOD , null,
    orientation, reference, INT_BLADE + C_NOD, eye,
  INT_BLADE+F_NOD,
    position, reference, INT_BLADE + C_NOD, null,
    orientation, reference, INT_BLADE + C_NOD, eye;
*/

joint: INT_BLADE + F_NOD, total joint,
  HUB,  # Static node C is substituted by an offset from the hub node
    position, reference, INT_BLADE + C_NOD , null,
    position orientation, reference, INT_BLADE + C_NOD, eye,
    rotation orientation, reference, INT_BLADE + C_NOD, eye,
  INT_BLADE+F_NOD,
    position, reference, INT_BLADE + C_NOD, null,
    position orientation, reference, INT_BLADE + C_NOD, eye,
    rotation orientation, reference, INT_BLADE + C_NOD, eye,
  position constraint, 1, 1, 1, null,
  orientation constraint, 1, 1, 0,null;

# 3.

joint: INT_BLADE + AF_LINK, spherical hinge,
  CURR_BLADE+AF_ROD,  # Static node C is substituted by an offset from the hub node
    position, reference, INT_BLADE + F_NOD+1 , null,
    orientation, reference, INT_BLADE + F_NOD+1, eye,
  INT_BLADE+F_NOD,
    position, reference, INT_BLADE + F_NOD+1, null,
    orientation, reference, INT_BLADE + F_NOD+1, eye;

/*
joint: CURR_BLADE + AF_LINK, total joint,
  CURR_BLADE+AF_ROD,  # Static node C is substituted by an offset from the hub node
    position, reference, INT_BLADE + F_NOD+1 , null,
    position orientation, reference, INT_BLADE + F_NOD+1, eye,
    rotation orientation, reference, INT_BLADE + F_NOD+1, eye,
  INT_BLADE+F_NOD,
    position, reference, INT_BLADE + F_NOD+1, null,
    position orientation, reference, INT_BLADE + F_NOD+1, eye,
    rotation orientation, reference, INT_BLADE + F_NOD+1, eye,
  position constraint, 1, 1, 1, null,
  orientation constraint, 0, 0, 0,null;
*/
/*
joint: CURR_BLADE + AF_LINK, distance,
  CURR_BLADE+AF_ROD,  # Static node C is substituted by an offset from the hub node
    position, reference, CURR_BLADE+A_NOD , null,
  INT_BLADE+D_NOD,
    position, reference, INT_BLADE+D_NOD+1, null,
  from nodes;
*/
# 4.
/*
joint: CURR_BLADE + D_NOD, revolute hinge,
  HUB,
    position, reference, CURR_BLADE + C_NOD , null,
    orientation, reference, CURR_BLADE + C_NOD, eye,
  CURR_BLADE+D_NOD,
    position, reference, CURR_BLADE + C_NOD, null,
    orientation, reference, CURR_BLADE + C_NOD, eye;
*/

joint: CURR_BLADE + D_NOD, total joint,
  HUB,
    position, reference, CURR_BLADE + C_NOD , null,
    position orientation, reference, CURR_BLADE + C_NOD, eye,
    rotation orientation, reference, CURR_BLADE + C_NOD, eye,
  CURR_BLADE+D_NOD,
    position, reference, CURR_BLADE + C_NOD, null,
    position orientation, reference, CURR_BLADE + C_NOD, eye,
    rotation orientation, reference, CURR_BLADE + C_NOD, eye,
  position constraint, 1, 1, 1,null,
  orientation constraint, 1, 1, 0, null;

# 5.

joint: INT_BLADE + B_NOD, spherical hinge,
  INT_BLADE,
    position, reference, INT_BLADE + B_NOD , null,
    orientation, reference, INT_BLADE + B_NOD, eye,
  INT_BLADE+BD_ROD,
    position, reference, INT_BLADE + B_NOD, null,
    orientation, reference, INT_BLADE + B_NOD, eye;

/*
joint: INT_BLADE + B_NOD, total joint,
  INT_BLADE,
    position, reference, INT_BLADE + B_NOD , null,
    position orientation, reference, INT_BLADE + B_NOD, eye,
    rotation orientation, reference, INT_BLADE + B_NOD, eye,
  INT_BLADE+BD_ROD,
    position, reference, INT_BLADE + B_NOD, null,
    position orientation, reference, INT_BLADE + B_NOD, eye,
    rotation orientation, reference, INT_BLADE + B_NOD, eye,
  position constraint, 1, 1, 1,null,
  orientation constraint, 0, 0, 0, null;
*/

# 6.
/*
joint: INT_BLADE + BD_LINK, total joint,
  INT_BLADE+BD_ROD,  # Static node C is substituted by an offset from the hub node
    position, reference, CURR_BLADE+D_NOD+1 , null,
    position orientation, reference, CURR_BLADE+D_NOD+1, eye,
    rotation orientation, reference, CURR_BLADE+D_NOD+1, eye,
  CURR_BLADE+D_NOD,
    position, reference, CURR_BLADE+D_NOD+1, null,
    position orientation, reference, CURR_BLADE+D_NOD+1, eye,
    rotation orientation, reference, CURR_BLADE+D_NOD+1, eye,
  position constraint, 1, 1, 1, null,
  orientation constraint, 0, 0, 0,null;
*/

joint: INT_BLADE + BD_LINK, spherical hinge,
  INT_BLADE+BD_ROD,  # Static node C is substituted by an offset from the hub node
    position, reference, CURR_BLADE+D_NOD+1 , null,
    orientation, reference, CURR_BLADE+D_NOD+1, eye,
  CURR_BLADE+D_NOD,
    position, reference, CURR_BLADE+D_NOD+1, null,
    orientation, reference, CURR_BLADE+D_NOD+1, eye;

/*
joint: INT_BLADE + BD_LINK, distance,
  INT_BLADE+BD_ROD,  # Static node C is substituted by an offset from the hub node
    position, reference, CURR_BLADE+B_NOD , null,
  CURR_BLADE+D_NOD,
    position, reference, CURR_BLADE+D_NOD+1, null,
  from nodes;
*/

  

