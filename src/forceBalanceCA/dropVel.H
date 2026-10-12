volVectorField& U = const_cast<volVectorField&>(this->db().lookupObject<volVectorField>("U"));	// read velocity vector field at each timestep

std::vector<vector> velocityMean(dropNumber, vector(0, 0, 0));	// Mean velocity for each droplet


for (int dropIdx = 0; dropIdx < dropNumber; dropIdx++)			// For every droplet in the domain
{
	vector velocityMeanSum(0, 0, 0);	// Sum of velocity for water cell (and single droplet)
	int countSum = 0;					// Water droplet counter (for single droplet)

	forAll(mesh.cells(), celli)			// Cycle through all the cells in the domain
	{
		if (alpha1[celli] > treshold_alpha && vofID[celli] == Labels[dropIdx])		// If water cell and it is part of the dropIdx droplet
		{
			countSum++;
			velocityMeanSum += U[celli];
		}
	}

	if (Pstream::parRun())				// Return TRUE if parallel simulation
	{
		reduce(countSum, sumOp<int>());
		reduce(velocityMeanSum, sumOp<vector>());
	}

	if (countSum > 0)					// If there is water in the domain
	{
		velocityMeanSum /= countSum;	// Mean velocity
	}

	velocityMean[dropIdx] = velocityMeanSum;
}