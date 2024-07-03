# fixed swashplate
	joint: SWASHPLATE_FIXED, total joint,
		AIRFRAME,
			position, reference, SWASHPLATE_FIXED, null,
			position orientation, reference, SWASHPLATE_FIXED, eye,
			rotation orientation, reference, SWASHPLATE_FIXED, eye,
		SWASHPLATE_FIXED,
			position, reference, SWASHPLATE_FIXED, null,
			position orientation, reference, SWASHPLATE_FIXED, eye,
			rotation orientation, reference, SWASHPLATE_FIXED, eye,
		position constraint, 1, 1, 1,
			0., 0., COLLECTIVE_GEAR_RATIO, # gear ratio
				# displacement controls collective pitch
				reference, D_PITCH_INPUT_COLLECTIVE,
		orientation constraint, 1, 1, 1,
			component,
				# rotation controls cyclic pitch
				# lateral
					drive, linear, 0., CYCLIC_GEAR_RATIO, reference, D_PITCH_INPUT_CYCLIC_LATERAL,
				# longitudinal
					drive, linear, 0., CYCLIC_GEAR_RATIO, reference, D_PITCH_INPUT_CYCLIC_LONGITUDINAL,
				inactive;

	# fixed to rotating swashplate
	joint: SWASHPLATE_ROTATING, total joint,
		SWASHPLATE_FIXED,
			position, reference, SWASHPLATE_FIXED, null,
			position orientation, reference, SWASHPLATE_FIXED, eye,
			rotation orientation, reference, SWASHPLATE_FIXED, eye,
		SWASHPLATE_ROTATING,
			position, reference, SWASHPLATE_ROTATING, null,
			position orientation, reference, SWASHPLATE_ROTATING, eye,
			rotation orientation, reference, SWASHPLATE_ROTATING, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 1, 0, null;

	# rotating swashplate to hub
	joint: SWASHPLATE_ROTATING + 1, total joint,
		HUB,
			position, reference, SWASHPLATE_ROTATING, null,
			position orientation, reference, SWASHPLATE_ROTATING, eye,
			rotation orientation, reference, SWASHPLATE_ROTATING, eye,
		SWASHPLATE_ROTATING,
			position, reference, SWASHPLATE_ROTATING, null,
			position orientation, reference, SWASHPLATE_ROTATING, eye,
			rotation orientation, reference, SWASHPLATE_ROTATING, eye,
		position constraint, 0, 0, 0, null,
		orientation constraint, 0, 0, 1, null;