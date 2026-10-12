std::vector<double> dropMass(dropNumber, 0.0);							// Droplet(s) Mass
std::vector<vector> GCenter(dropNumber, vector(0, 0, 0));				// Droplet(s) center of mass position
static std::vector<vector> GCenterOld(dropNumber, vector(0, 0, 0));		// Droplet(s) center of mass position of previous simulation (check concistency in droplet numbering)

for (int dropIdx = 0; dropIdx < dropNumber; dropIdx++)		// For every droplet in the domain
{
	int waterCell = 0;							// Water cell counter
	double dropMass_s = 0;						// Water droplet mass [kg]
	vector cellCenter(0, 0, 0);					// Temporary vector to store sum of water cell coordinates
	vector GCenter_s(0, 0, 0);					// Droplet center of mass position

	forAll(mesh.cells(), celli)					// Cycle through all the cells in the domain
	{
		if (alpha1[celli] > treshold_alpha && vofID[celli] == Labels[dropIdx])		// If water cell and it is part of the dropIdx droplet
		{
			waterCell++;
			cellCenter += mesh.C()[celli];    
			dropMass_s += rho_water * mesh.V()[celli] * alpha1[celli];
		}
	}

	if (Pstream::parRun())                    // Return TRUE if parallel simulation
	{
		reduce(waterCell, sumOp<int>());      // Sum the number of water cell over all partitions
		reduce(cellCenter, sumOp<vector>());  // Sum the cellCenter vector over all partitions
		reduce(dropMass_s, sumOp<double>());    // Sum the water mass over all partitions
	}

	// Evaluate water droplet center of mass coordinates
	if (waterCell > 0)
	{
		GCenter_s = cellCenter / waterCell;
	}

	// Save data to vectors
	dropMass[dropIdx] = dropMass_s;
	GCenter[dropIdx] = GCenter_s;
}



/*
// Dati i centri di massa, voglio che la prima goccia sia sempre la prima e la seconda sia sempre la seconda
// Il centro di massa della prima goccia deve essere quindi a massimo una cella di distanza (courant = 0.5 nella simulazione) rispetto alla prima goccia dell'iterazione precedente
// Questo deve essere valido per tutte le N (=dropMass) gocce


Info << mesh.deltaCoeffs()[0] << endl;
for (int i = 0; i < dropNumber; i++)
{

	for (int j = 0; j < dropNumber; j++)
	{
		scalar dist = mag(GCenter[i] - GCenterOld[j]);

		if (dist < mesh.deltaCoeffs()[0])
		{
			std::swap(Labels[i], Labels[j]);
			std::swap(dropMass[i], dropMass[j]);
			std::swap(GCenter[i], GCenter[j]);

			break;
		}
	}
}

GCenterOld = GCenter;
*/