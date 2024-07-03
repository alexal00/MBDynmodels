# thesis
This git repo contains the main scripts used inside the WSL environment where MBDyn was installed.
It spans several folders each of which is used in a different chapt of the thesis titled:
A Comparative Study of the Impact of Innovative Lead-lag Damping Configurations on Helicopter Rotor Stability and Loads
## Folders
1. GR_MBDyn
    Chapter 1: Ideal multi-body simulation of ground resonance
    alvaro_model/
        v9 : Initial conditions prescribed by means of constant value, NOTE: Non-zero velocity values
        v10 : Initial conditions prescribed by means of step function that goes to zero at a t = 2pi/Omega
        v11 : Improved MB model of the i2b TODO
    cassoni_model/
        reference GR MB model
2. GR_AERO
    Chapter 3: Influence of aerodynamics and out-of-plane motions on ground resonance stability
    Same as GR_MBDyn model with aerodynamics and swashplate
3. AERO_MBDyn
    Chapter 4: Loads in forward flight
    * KINCUP
        TODO
2. FRICTION
    Verification of the correct implementation of the non-linear CL used in
    Chapter 5: Implementation and analysis of non-linear elastomeric dampers

