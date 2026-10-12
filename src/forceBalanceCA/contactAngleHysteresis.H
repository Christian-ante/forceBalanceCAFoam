Foam::tmp<Foam::scalarField> theta_dyn(new Foam::scalarField(patch().size()));		// Output scalar field
scalarField& theta_dyn_ref = theta_dyn.ref();										// Scalar field to access to variable stored in theta_dyn in the for cycle
theta_dyn_ref = theta0_;															// Initialize field value to equilibrium contact angle

static std::vector<bool> dropAlreadySlidingCheck(dropNumber, false);				// Flag to determine if the droplet is already sliding in the simulation

for (int dropIdx = 0; dropIdx < dropNumber; dropIdx++)		// For every droplet in the domain
{
	std::string fileName_check = "Adhesion_Check_" + std::to_string(dropIdx + 1) + ".txt";
	std::string fileName_dec = "dropDeceleration_" + std::to_string(dropIdx + 1) + ".txt";

	if (Pstream::master())
	{
		// Cancel txt files (will be written again on next iteration if still in adhesion condition)
		std::remove(fileName_check.c_str());
		std::remove(fileName_dec.c_str());
	}

	if (ratio[dropIdx] <= 1)			// Adhesion force is higher than drag force, need to calculate the proper advancing and receding angles
	{
		if (Pstream::master())
		{
			std::ofstream outFile_adh;
			outFile_adh.open(fileName_check, std::ios::out);		// Write this txt only in adhesion condition -> it works as a flag for the modified solver
			outFile_adh.close();

			double dropDeceleration = std::abs((totalDragForce[dropIdx][0] - adhesionForce[dropIdx])) / dropMass[dropIdx];		// Droplet deceleration due to adhesion on surface
			std::ofstream outFile_adh_dec;
			outFile_adh_dec.open(fileName_dec, std::ios::out);																	// Pass droplet deceleration to the solver
			if (dropAlreadySlidingCheck[dropIdx] == false)			// If the droplet was not moving 
			{
				outFile_adh_dec << 0;
			}
			else													// If the droplet was moving
			{
				outFile_adh_dec << dropDeceleration;
			}
			outFile_adh_dec.close();
		}

		forAll(faceCells, faceI)								// Cycle through all the patch faces
		{
			label cellI = faceCells[faceI];						// Cell to which the face belong
			if (vofID[cellI] == Labels[dropIdx])
			{
				const scalar& uu = uwall[faceI];				// Contact velocity of the i-th boundary cell
				
				if (uu >= 0)									// Receding front
				{
					theta_dyn_ref[faceI] = theta0_ - ratio[dropIdx] * (theta0_ - thetaR_) * factor[faceI];		// Linear trend in between contact angle hysteresis
				}
				
				else											// Advancing front
				{
					theta_dyn_ref[faceI] = theta0_ + ratio[dropIdx] * (thetaA_ - theta0_) * factor[faceI];		// Linear trend in between contact angle hysteresis
				}

				theta_dyn_OLD[faceI] = theta_dyn_ref[faceI];	// Pass theta scalar field to next timestep
			}
		}
	}
}