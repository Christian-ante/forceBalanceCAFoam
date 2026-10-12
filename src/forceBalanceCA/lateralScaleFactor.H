static scalarField theta_dyn_OLD(patch().size(), theta0_);			// Theta scalar field of the previous timestep

scalarField factor(patch().size(), 0.0);							// Scale factor of evaluated contact angle (to consider lateral surface of the advancing droplet) (LINEAR TREND AT THE MOMENT)
scalarField& factor_ref = factor;									// Reference scalar field to access factor scalar field

scalarField uwall(patch().size());									// Velocity of the water front (equal to specific droplet velocity in this simulation)
scalarField& uwall_ref = uwall;										// Reference scalar field to access uwall scalar field

for (int dropIdx = 0; dropIdx < dropNumber; dropIdx++)				// For every droplet in the domain
{
	if (mesh.time().value() != 0)									// Skip the first timestep
	{
		double max_theta = 0, min_theta = 180;						// Max/min values in the previous theta scalar field

		forAll(faceCells, faceI)									// Cycle through all the patch faces
		{
			label cellI = faceCells[faceI];							// Cell to which the face belong

			if (vofID[cellI] == Labels[dropIdx])					// Evaluate max and min applied contact angle for each droplet
			{
				if (theta_dyn_OLD[faceI] > max_theta)
				{
					max_theta = theta_dyn_OLD[faceI];
				}
				else if (theta_dyn_OLD[faceI] < min_theta)
				{
					min_theta = theta_dyn_OLD[faceI];
				}
			}
		}
		if (Pstream::parRun()) // return TRUE if parallel simulation
		{
			reduce(max_theta, maxOp<double>());
			reduce(min_theta, minOp<double>());
		}

		double n_x_max = 1 * sin(max_theta * deg2rad);				// Max (positive) projection of the normal versor to the interface on X-dir
		double n_x_min = 1 * sin(min_theta * deg2rad);				// Min (negative) projection of the normal versor to the interface on X-dir

		forAll(faceCells, faceI)									// Cycle through all the patch faces
		{
			label cellI = faceCells[faceI];							// Cell to which the face belong

			if (vofID[cellI] == Labels[dropIdx])
			{
				if (n_Inter[cellI][0] <= 0.0)						// Advancing front (the normal to the interface is ENTERING the droplet)
				{
					factor_ref[faceI] = std::abs(n_Inter[cellI][0] / n_x_max);		// Scale factor of the advancing contact angle
					uwall_ref[faceI] = -velocityMean[dropIdx][0];					// Contact line velocity of the advancing front (front-half droplet)
				}
				else												// Receding front
				{
					factor_ref[faceI] = std::abs(n_Inter[cellI][0] / n_x_min);		// Scale factor of the receding contact angle
					uwall_ref[faceI] = velocityMean[dropIdx][0];					// Contact line velocity of the receding front (back-half droplet)
				}

				// Limit the value of factor_ref to 1
				if (factor_ref[faceI] > 1)
				{
					factor_ref[faceI] = 1;
				}
			}
		}
	}
}

theta_dyn_OLD = theta0_;		// Reinitialize theta_dyn_OLD for next timestep; values in the field will be modified when evaluating contact angle of droplets