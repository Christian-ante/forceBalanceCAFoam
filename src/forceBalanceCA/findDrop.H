// Flood-fill algorithm

const volScalarField& alpha1 = mesh.lookupObject<volScalarField>("alpha.water");	// Read alpha of the primary phase (water) in the mesh at each timestep
const scalar treshold_alpha = 0.02;													// Volume fraction threshold to define water cells

// Warn (once) if this boundary condition is applied to more than one patch: per-patch state is kept
// separately, but the droplet exchange files (dropNumber.txt, Adhesion_Check_<i>.txt, ...) are shared
static bool multiPatchChecked = false;
if (!multiPatchChecked && mesh.time().timeIndex() > mesh.time().startTimeIndex())
{
	label nPatchesWithBC = 0;
	forAll(alpha1.boundaryField(), patchI)
	{
		if (alpha1.boundaryField()[patchI].type() == typeName)
		{
			nPatchesWithBC++;
		}
	}
	if (nPatchesWithBC > 1)
	{
		WarningInFunction
			<< typeName << " is applied to " << nPatchesWithBC << " patches." << nl
			<< "    Forces and contact angles are evaluated per patch, but the exchange files read by the" << nl
			<< "    solver are shared: merge the walls into a single patch (e.g. with createPatch)." << endl;
	}
	multiPatchChecked = true;
}

const fvPatch& patch_ = this->patch();					// Reference to the patch where the dynamic angle boundary is applied
const labelUList& faceCells = patch_.faceCells();		// Reference to the faces of the patch

volScalarField vofID						// vofID field to store different labels for each droplet
(
	IOobject
	(
		"vofID",                           
		mesh.time().timeName(),             
		mesh,                               
		IOobject::NO_READ,
		IOobject::AUTO_WRITE    
	),
	mesh,
	dimensionedScalar("zero", dimless, 0)   // Initialize to 0
);

globalIndex globalNumbering(mesh.nCells());		// Global cell index to avoid confusion between processors


forAll(mesh.cells(), cellI)						// Cycle through all the cells in the domain
{
	if (alpha1[cellI] > treshold_alpha && vofID[cellI] == 0.0)		// Cell is considered liquid and not labelled yet
	{
		label globalCellID = globalNumbering.toGlobal(cellI);		// Convert local processor cellIDs to global cellIDs (for parallel simulations)

		vofID[cellI] = globalCellID;								// Assign to the cell its ID

		DynamicList<label> cells(1);								// Create dynamic list to store all neighbours
		cells.append(mesh.cellCells()[cellI]);						// Get neighbours of cellI

		while (!cells.empty())						// Loop until cells vector is empty
		{
			DynamicList<label> neighbours(0);		// Create dynamic list to store all neighbourood neighbours

			forAll(cells, i)						// Loop over cells vector
			{
				label cellI_ = cells[i];

				if (alpha1[cellI_] > treshold_alpha && vofID[cellI_] == 0.0)		// Check if neighbour cell is liquid and has no ID yet
				{
					vofID[cellI_] = globalCellID;									// Assign to the cell cellI ID (since continuous liquid volume)
					neighbours.append(mesh.cellCells()[cellI_]);					// Add neighbours to neighbourood list
				}
			}
			cells.clear();							// Empty cells
			cells.transfer(neighbours);				// Move neighbourood list to cells
		}
	}
}
vofID.correctBoundaryConditions();					// Update vofID on all processors


// Correct for parallel processing
if (Pstream::parRun())			// Return TRUE if parallel simulation
{
	bool corr_check = 0;		// Flag to exit (early) from the cycle
	for (int i = 0; i < Pstream::nProcs() - 1; i++)
	{
		forAll(mesh.boundary(), patchI)								// Loop over all partition patches
		{
			const fvPatch& pPatch = mesh.boundary()[patchI];		// Current patch

			if (pPatch.coupled())									// Check if patch is coupled, e.g. processor patch
			{
				forAll(pPatch, faceI)								// Loop over all faces of processor patch
				{
					label cellI = pPatch.faceCells()[faceI];												// Cell label that "own" faceI on current core (NB: the face is shared on both processors!)
					const scalar ID_own = vofID.boundaryField()[patchI].patchInternalField()()[faceI];		// vofID field value for the current cell
					const scalar ID_nei = vofID.boundaryField()[patchI].patchNeighbourField()()[faceI];		// vofID field value for the neighbor cell

					if (ID_own > ID_nei && ID_nei > 0)												// Both cell are water cell, but ID_own is bigger (higher "globalCellID") than ID_nei
					{
						vofID[cellI] = ID_nei;														// Update vofID
						corr_check = 1;																// Update correction flag (still correcting vofID in this cycle)

						DynamicList<label> cells_corr(1);											// Create dynamic list to store all neighbours
						cells_corr.append(mesh.cellCells()[cellI]);									// Get neighbours of cellI
						
						while (!cells_corr.empty())													// Loop until cells_corr vector is empty
						{
							DynamicList<label> neighbours_corr(1);									// Create dynamic list to store all neighbourood neighbours

							forAll(cells_corr, i)													// Loop over cells vector
							{
								if (vofID[cells_corr[i]] > vofID[cellI])							// Check if vofID of the current cell is bigger than "ID_nei" (value from neighboor processor)
								{
									vofID[cells_corr[i]] = vofID[cellI];							// Update vofID

									neighbours_corr.append(mesh.cellCells()[cells_corr[i]]);		// Add neighbours to neighbourood list
								}
							}
							cells_corr.clear();														// Empty cells_corr
							cells_corr.transfer(neighbours_corr);									// Move neighbourood list to cells_corr
						}
					}			
				}
			}
		}
		vofID.correctBoundaryConditions();				// Update vofID on all processors

		reduce(corr_check, maxOp<bool>());				// Synchronize corr_check value over all processors
		if (corr_check == 0)							// If vofID was not updated in this cycle
		{
			break;										// Break the for loop
		}
		corr_check = 0;									// Refresh corr_check value for all processors
	}
}

vofID.write();

// Synchronize labels for all processors
int dropNumber;					// Droplet deposited on this boundary
List<label> Labels;
if (Pstream::parRun())			// Return TRUE if parallel simulation
{
	DynamicList<label> localLabels;		// Local list of unique vofID label for each processor
	forAll(vofID, cellI)
	{
		label id = vofID[cellI];
		if (id > 0 && findIndex(localLabels, id) == -1)		// New id must be positive and must not be already present in the localLabels list
		{
			localLabels.append(id);							// Append it to the list
		}
	}

	// Pstream::gatherList seems to not like DynamicList, thus we create a List variable and copy localLabels values
	List<label> localLabelsList(localLabels.size());
	forAll(localLabels, i)
	{
		localLabelsList[i] = localLabels[i];
	}

	List<List<label>> gatheredLabels(Pstream::nProcs());			// List of label lists to gather data from all processors
	gatheredLabels[Pstream::myProcNo()] = localLabelsList;			// Each processor fills its line
	Pstream::gatherList(gatheredLabels);							// Complete lList of label lists is passed to master processor

	DynamicList<label> globalLabels;		// Dynamic list to store unique value of vofID
	dropNumber = 0;							// Droplet number counter
	if (Pstream::master())					// Execute only in master processor
	{
		// Check for unique value in gatheredLabels list and remove duplicates
		forAll(gatheredLabels, gatI)
		{
			forAll(gatheredLabels[gatI], i)
			{
				label id = gatheredLabels[gatI][i];
				if (findIndex(globalLabels, id) == -1)				// Append to globalLabels only if not already present in the list
				{
					globalLabels.append(id);
					dropNumber++;
				}
			}
		}
	}
	reduce(dropNumber, sumOp<int>());			// Synchronize droplet number on all processors
				
	Labels.setSize(dropNumber);					// List to synchronize unique vofID in all processors
	if (Pstream::master())						// Execute only in master processor
	{
		for (int i = 0; i < dropNumber; i++)
		{
			Labels[i] = globalLabels[i];
		}
	}
	else										// Execute only in slave processors
	{
		for (int i = 0; i < dropNumber; i++)
		{
			Labels[i] = 0;
		}
	}
	reduce(Labels, sumOp<List<label>>());		// The sum will be equal to the value stored in master processor, other cores will always have 0!
}
else	// NOT parallel simulation
{			
	dropNumber = 0;						// Droplet number counter
	forAll(vofID, cellI)
	{
		label id = vofID[cellI];
		if (id > 0 && findIndex(Labels, id) == -1)		// New id must be positive and must not be already present in the localLabels list
		{
			Labels.append(id);							// Append it to the list
			dropNumber++;
		}
	}
}

// Pass to the solver the number of found droplets and respective labels
if (Pstream::master())
{
	static int dropNumberOld = 0;

	if (dropNumber != dropNumberOld)						// Rewrite the file only if droplet number changes (break up/coalescence)
	{
		std::remove("dropNumber.txt");

		std::ofstream drop_info;
		drop_info.open("dropNumber.txt", std::ios::out);
		drop_info << dropNumber;
		drop_info.close();
	}
	dropNumberOld = dropNumber;

	// Labels
	forAll(Labels, i)
	{
		std::string fileName_label = "droplet_" + std::to_string(i + 1) + "_label.txt";
		std::remove(fileName_label.c_str());

		std::ofstream outFile_label;
		outFile_label.open(fileName_label, std::ios::out);
		outFile_label << Labels[i];
		outFile_label.close();
	}
}