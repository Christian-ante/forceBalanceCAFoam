int dropNumber = 0;							// Droplet number (evaluated from boundary)
double treshold_alpha = 0.01;				// Treshold to define water cells

// Droplet number (from the boundary)
std::ifstream infile("dropNumber.txt");
infile >> dropNumber;
infile.close();

// List of droplet label (from the boundary)
List<label> Labels(dropNumber, 0);
forAll(Labels, i)
{
	std::string fileName_label = "droplet_" + std::to_string(i + 1) + "_label.txt";
	std::ifstream inFile_label(fileName_label);
	inFile_label >> Labels[i];
	inFile_label.close();
}

// Retrieve vofID field (evaluated by boundary!)
volScalarField vofID
(
	IOobject
	(
		"vofID",
		mesh.time().timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	),
	mesh
);


for (int dropIdx = 0; dropIdx < dropNumber; dropIdx++)								// For every droplet in the domain
{
	std::string fileName_check = "Adhesion_Check_" + std::to_string(dropIdx + 1) + ".txt";
	std::string fileName_dec = "dropDeceleration_" + std::to_string(dropIdx + 1) + ".txt";

	std::ifstream infile(fileName_check);
	if (infile.is_open())															// Open "Adhesion_Check.txt" file
	{
		double dropDeceleration;
		std::ifstream infile_dec(fileName_dec);
		infile_dec >> dropDeceleration;												// Read and set dropDeceleration from the value calculated in user boundary
		infile_dec.close();

		if (dropDeceleration == 0)													// The droplet was not moving -> pin in position
		{
			forAll(mesh.boundaryMesh(), patchi)										// Loop through all boundary patches
			{
				const polyPatch& boundaryPatch = mesh.boundaryMesh()[patchi];

				if (boundaryPatch.name() == "defaultFaces")							// Skip "defaultFaces" (for 2D simulations)
				{
					continue;
				}

				forAll(boundaryPatch, faceI)										// Loop over all faces in the boundary
				{
					for (label faceI = 0; faceI < boundaryPatch.size(); ++faceI)
					{
						label faceID = boundaryPatch.start() + faceI;				// Global index of the face (= Start index of the boundary + specific index of the boundary)
						label cellI = mesh.faceOwner()[faceID];						// Global index of the cell

						if (alpha1[cellI] > treshold_alpha && vofID[cellI] == Labels[dropIdx])			// If the cell is a water or interface cell and it is the correct droplet
						{
							U[cellI] = vector::zero;								// Pin the droplet
						}
					}
				}
			}
		}
		else																		// The droplet was sliding -> decelerate
		{
			forAll(mesh.cells(), celli)												// Loop over all cells
			{
				if (alpha1[celli] > treshold_alpha && vofID[celli] == Labels[dropIdx])					// If the cell is a water or interface cell and it is the correct droplet
				{
					U[celli][0] = std::max(0.0, U[celli][0] - dropDeceleration * runTime.deltaTValue()); 		// Reduce the velocity of water cell (if in adhesion condition)
				}
			}
		}
	}
}

// Remove extra saved timesteps							NEED TO CHANGE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if (!mesh.time().outputTime())
{
	word lastTimeDir = mesh.time().path() / mesh.time().timeName() / "vofID";
	std::remove(lastTimeDir.c_str());
}
