# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2017
# 
# Pierangelo Masarati	<masarati@aero.polimi.it>
# Paolo Mantegazza	<mantegazza@aero.polimi.it>
# 
# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
# via La Masa, 34 - 20156 Milano, Italy
# http://www.aero.polimi.it
# 
# Changing this copyright notice is forbidden.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (version 2 of the License).
# 
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

aerodynamic body: CURR_BLADE,
	CURR_BLADE, # blade node
	induced velocity, 99,
	reference, CURR_BLADE, (R + R_CUTOUT)/2, 0., 0.,
	reference, CURR_BLADE,
		1, 0., 1., 0.,	# local axis 1 is trailing egde to leading edge
		2, 0., 0., 1.,	# axis 2 is "up" (?)...
		# 3, 1., 0., 0., # which means that local axis 3 is along the blade axis
	# span (const)
	R - R_CUTOUT,
	# chord distribution (shape 1D)
	### constant chord
	const, CHORD,
	# location of the aerodynamic center
	### aerodynamic center at the feathering axis
	const, 0.,
	# location of boundary condition evaluation point
	## angle of attack evaluated at the aerodynamic center
	const, 0.,
	# built-in twist (positive nose up)
	### no twist
	# const, 0.,
	### linear, negative twist; nominal pitch at 75 % of the span
	linear, 3./4.*BLADE_TWIST, -1./4.*BLADE_TWIST,
	# number of Gauss points
	N_GAUSS_POINTS_AERO,
	# airfoil data
	c81, NACA0012;

# vim:ft=mbd
