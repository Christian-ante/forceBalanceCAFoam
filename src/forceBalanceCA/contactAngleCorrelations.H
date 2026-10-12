scalarField Ca = mu_water * uwall / sigma;		// Capillary number scalar field

for (int dropIdx = 0; dropIdx < dropNumber; dropIdx++)		// For every droplet in the domain
{
	if (ratio[dropIdx] > 1)									// Adhesion force is lower than drag force, need to calculate the proper advancing and receding angles
	{
		forAll(faceCells, faceI)								// Cycle through all the patch faces
		{
			label cellI = faceCells[faceI];						// Cell to which the face belong
			if (vofID[cellI] == Labels[dropIdx])
			{
				const scalar& u = uwall[faceI];					// contact velocity of the i-th boundary cell

				if (u >= 0)	// Dynamic Receding Contact Angle model
				{
					double theta_rad = thetaR_ * deg2rad;
					double cubedValue = std::pow(theta_rad, 3);
					double innerValue = std::cbrt(cubedValue - 72.0 * Ca[faceI]);
					double theta_dyn_temp = innerValue * rad2deg;
					double delta = theta0_ - theta_dyn_temp;
					theta_dyn_ref[faceI] = theta0_ - delta * factor[faceI];
				}
				else 		// Kistler's (advancing angle) model
				{
					double theta_rad = thetaA_ * deg2rad;
					double cubedValue = std::pow(theta_rad, 3);
					double f_hoff_inv = cubedValue / 72;
					double x_hoff = Ca[faceI] + f_hoff_inv;
					double temp = acos(1.0 - 2.0 * tanh(5.16 * pow((x_hoff / (1.0 + 1.31 * pow(x_hoff, 0.99))), 0.706)));
					double theta_dyn_temp = temp * rad2deg;
					double delta = theta_dyn_temp - theta0_;
					theta_dyn_ref[faceI] = theta0_ + delta * factor[faceI];
				}
				theta_dyn_OLD[faceI] = theta_dyn_ref[faceI];								// Pass theta scalar field to next timestep
			}
		}
		dropAlreadySlidingCheck[dropIdx] = true;											// The droplet is sliding --> update the flag!
	}
}

return theta_dyn;