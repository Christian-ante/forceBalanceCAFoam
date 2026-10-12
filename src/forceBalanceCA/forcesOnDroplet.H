constexpr double deg2rad = M_PI / 180.0;		// From angles in degrees to radians
constexpr double rad2deg = 180.0 / M_PI;		// From angles in radiants to degrees
			
const volVectorField gradAlpha1 = fvc::grad(alpha1);							// Vector field of gradient of volume fraction in the domain
dimensionedScalar smallDim("smallDim", gradAlpha1.dimensions(), SMALL);			// Avoid division by 0!!
const volVectorField n_Inter = gradAlpha1 / (mag(gradAlpha1) + smallDim);		// Normal versor to the interface
const volTensorField gradU = fvc::grad(U);										// Tensor field of gradient of velocity in the domain

using Point = std::vector<double>;						// Declaration for a point
using Edge = std::pair<Point, Point>;					// Declaration for a couple of points
using EdgeList = std::vector<Edge>;						// Declaration for a list of couples of points (edges)
using PointCloud = std::vector<Point>;					// Declaration for a list of points (cloud)
using VectorList = std::vector<std::vector<double>>;	// Declaration for a list of vectors

std::vector<double> ratio(dropNumber, 0);							// Ratio between total drag and adhesion to check if the droplet is in adhesion or slide/roll condition
std::vector<double> adhesionForce(dropNumber, 0);					// Adhesion force for each droplet
std::vector<vector> viscousForce(dropNumber, vector(0, 0, 0));		// Sum of viscous force acting on the water-air interface for each droplet
std::vector<vector> pressureForce(dropNumber, vector(0, 0, 0));		// Sum of pressure forces acting on the water-air interface for each droplet
std::vector<double> g_force_X(dropNumber, 0);						// Gravitational force [N] in X-direction
std::vector<vector> totalDragForce(dropNumber, vector(0, 0, 0));	// Total drag force (viscous+pressure) acting on each droplet

std::vector<double> dA_TOT(dropNumber, 0);				// Total interface area [m2] for each droplet

for (int dropIdx = 0; dropIdx < dropNumber; dropIdx++)			// For every droplet in the domain
{
	double dA = 0;							// Water-air interface area of the cell
	double D = 0;							// Parameter of the cutting plane that approximates the water-air interface

	forAll(mesh.cells(), celli)
	{
		if (alpha1[celli] > treshold_alpha && alpha1[celli] < 0.99 && vofID[celli] == Labels[dropIdx])		// Only if water-air interface cells and current droplet (Labels[dropIdx])
		{
			// Coefficient of the cutting plane in the form: ax + by + cz = d
			double a = n_Inter[celli][0];
			double b = n_Inter[celli][1];
			double c = n_Inter[celli][2];

			// NB: CELL MUST BE A CUBE (EVEN IN 2D CASE)
			double cellVolume = mesh.V()[celli];	// Cell volume [m3]
			double l = cbrt(cellVolume);			// Edge of the cube [m]

			// The plane is expressed in the form "y = Ax + D" or "z = Ax + By + D", we need to find the proper D to obtain the target alpha for each cell
			// Parameters for the research of the correct D
			double D_min = -10;			// Lower limit for D
			double D_max = 10;			// Higher limit for D
			int maxIterations = 100;	// Max # of iterations
			double tolerance = 0.01;	// Difference between alpha and alpha_target (in terms of volume or area fraction)

			if (std::abs(c) < 0.05)		// c ~ 0 -> The plane is expressed as "y = Ax + D"
			{
				if (std::abs(a) > 0.95 || std::abs(b) > 0.95)	// (a ~ 0 || b ~ 0) -> The plane is expressed as "x = D" ("y = D")
				{
					dA = l * l;			// The plane is parallel to one of the cube face
				}
				else
				{
					// y = Ax + D
					double A = -a / b;

					// Edges of the square (as couples of vertices)
					EdgeList edges = {
					{{0, 0}, {l, 0}},
					{{0, 0}, {0, l}},
					{{l, 0}, {l, l}},
					{{0, l}, {l, l}}
					};

					// Create points cloud dividing the cube numPoints times along each direction
					int numPoints = 10;							// Number of divisions
					PointCloud points;							// Points cloud declaration
					points.reserve(numPoints * numPoints);		// Reserve space to improve efficiency

					// Points cloud definition
					double step = l / (numPoints - 1);			// Distance between subsequent points
					for (int i = 0; i < numPoints; ++i)
					{
						for (int j = 0; j < numPoints; ++j)
						{
							double x = i * step;
							double y = j * step;
							points.push_back({ x, y });
						}
					}

					D = findInterfacePlane(l, A, D_min, D_max, maxIterations, tolerance, alpha1[celli], points);		// Function to find proper D
					VectorList intersections = findIntersection(A, D, edges);											// Function to find intersection between plane and cube
					dA = intersectionArea(intersections, l);															// Function to evaluate interface area dA

					// Additional control to check if there is a wrong number of intersection for a 2D case
					if (intersections.size() != 2)
					{
						std::stringstream fileName;
						fileName << "WARNING_2D_wrongIntersectionNumber=" << intersections.size() << "_time=" << mesh.time().value() << "[s].txt";

						std::ofstream outFileW;								// Create void txt file
						outFileW.open(fileName.str(), std::ios::out);
						outFileW.close();
					}
				}
			}
			else
			{
				// z = Ax + By + D
				double A = -a / c;
				double B = -b / c;

				// Edges of the cube (as couples of vertices)
				EdgeList edges = {
				{{0, 0, 0}, {l, 0, 0}},
				{{0, 0, 0}, {0, l, 0}},
				{{0, 0, 0}, {0, 0, l}},
				{{l, 0, 0}, {l, l, 0}},
				{{l, 0, 0}, {l, 0, l}},
				{{0, l, 0}, {l, l, 0}},
				{{0, l, 0}, {0, l, l}},
				{{0, 0, l}, {l, 0, l}},
				{{0, 0, l}, {0, l, l}},
				{{l, l, 0}, {l, l, l}},
				{{l, 0, l}, {l, l, l}},
				{{0, l, l}, {l, l, l}}
				};

				// Create points cloud dividing the cube numPoints times along each direction
				int numPoints = 10;										// Number of divisions
				PointCloud points;										// Points cloud declaration
				points.reserve(numPoints * numPoints * numPoints);		// Reserve space to improve efficiency

				// Points cloud definition
				double step = l / (numPoints - 1);			// Distance between subsequent points
				for (int i = 0; i < numPoints; ++i)
				{
					for (int j = 0; j < numPoints; ++j)
					{
						for (int k = 0; k < numPoints; ++k)
						{
							double x = i * step;
							double y = j * step;
							double z = k * step;
							points.push_back({ x, y, z });
						}
					}
				}

				D = findInterfacePlane(l, A, B, D_min, D_max, maxIterations, tolerance, alpha1[celli], points);				// Function to find proper D
				VectorList intersections = findIntersection(A, B, D, edges);												// Function to find intersection between plane and cube
				dA = intersectionArea(intersections);																		// Function to evaluate interface area dA

				// Additional control to check if there is a wrong number of intersection for a 3D case
				if (intersections.size() < 3 || intersections.size() > 6)
				{
					std::stringstream fileName;
					fileName << "WARNING_3D_wrongIntersectionNumber=" << intersections.size() << "_time=" << mesh.time().value() << "[s].txt";

					std::ofstream outFileW;														// Create void txt file
					outFileW.open(fileName.str(), std::ios::out);
					outFileW.close();
				}
			}

			// NB: PRESSURE AND VISCOUS DRAG FORCES ARE EVALUATED AT A CERTAIN DISTANCE (PARAMETER "k") FROM THE INTERFACE CELL
			int k = 3;															// Number of cells from the interface
			vector p_dyn_point = mesh.C()[celli] - k * n_Inter[celli] * l;		// Position to evaluate dynamic pressure and viscous forces contributes

			int closestCell_p = -1;				// ID of the closest cell to reference cell
			double MinDistance_p = GREAT;		// Variable to store the minimum distance between cell centroid and reference cell
			double distance = 0;				// Variable to store distance between each cell and reference one
			double p_dyn = 0;					// Dynamic pressure
			vector p_dyn_point_cell;			// Position of the closest cell to reference cell
			vector p_dyn_point_cell_vel;		// Velocity of the closest cell to reference cell
			vector cellC;						// Centroid of the cell in the cycle
			vector localPressureForce;			// Local pressure drag force
			vector localViscousForce;			// Local viscous drag force
			tensor viscousStress;				// Local viscous stress

			forAll(mesh.cells(), celli_p) 		// Cycle through all the cell in the domain, looking to the closest one to "p_dyn_point"
			{
				cellC = mesh.C()[celli_p];					// Center of celli_p cell
				distance = mag(p_dyn_point - cellC);		// Distance form selected point and celli_p center
				if (distance < MinDistance_p)				// Check if the new distance is lower than the previously saved one
				{
					MinDistance_p = distance;				// Update the minimum distance from the reference cell
					closestCell_p = celli_p;				// ID of the closest cell to reference cell
				}
			}

			p_dyn_point_cell = mesh.C()[closestCell_p];					// Position of the closest cell to reference cell
			p_dyn_point_cell_vel = U[closestCell_p];					// Velocity of the closest cell to reference cell
			viscousStress = mu_water * gradU[closestCell_p];			// Local viscous stress

			p_dyn = 0.5 * rho_air * magSqr(p_dyn_point_cell_vel);				// Dynamic pressure = 1/2 * rho_air * U^2
			localPressureForce = n_Inter[celli] * p_dyn * dA;					// Local pressure drag force
			pressureForce[dropIdx] += localPressureForce;								// Total pressure force

			localViscousForce = viscousStress & (n_Inter[celli] * dA);			// Local viscous drag force
			viscousForce[dropIdx] += localViscousForce;									// Total viscous force

			dA_TOT[dropIdx] += dA;
		}
	}
	totalDragForce[dropIdx] = pressureForce[dropIdx] + viscousForce[dropIdx];			// Sum the total drag force for each interface cell

	// Need to calculate the droplet wet area on the wall for the adhesion force calculation
	double dA_adh = 0;										// Droplet wet area on the surface
	forAll(faceCells, faceI)								// Cycle through all the patch faces
	{
		label cellI = faceCells[faceI];						// Cell to which the face belong
		if (alpha1[cellI] > treshold_alpha && vofID[cellI] == Labels[dropIdx])							// If the cell is a water or interface cell
		{
			double cellVolume_adh = mesh.V()[cellI];		// Cell volume [m3]
			double l_adh = cbrt(cellVolume_adh);			// Edge of the cube [m]
			dA_adh += l_adh * l_adh;						// Contact area on the surface (to evaluate adhesion force)
		}
	}

	if (Pstream::parRun())									// Return TRUE if parallel simulation
	{
		reduce(pressureForce[dropIdx], sumOp<vector>());	// Evaluate the sum of pressure force over all the partitions
		reduce(viscousForce[dropIdx], sumOp<vector>());		// Evaluate the sum of viscous force over all the partitions
		reduce(totalDragForce[dropIdx], sumOp<vector>());	// Evaluate the sum of drag force over all the partitions
		reduce(dA_adh, sumOp<double>());					// Evaluate the total contact area
		reduce(dA_TOT[dropIdx], sumOp<double>());			// Evaluate total interface area for each droplet				
	}

	// ElSherbini, A.I., &Jacobi, A.M. (2006).Retention forces and contact angles for critical liquid drops on non - horizontal surfaces.Journal of colloid and interface science, 299(2), 841 - 849.
	double D_adh = std::sqrt(4 * dA_adh / M_PI);																					// Equivalent diameter of the circle with the same area as the wetted area
	adhesionForce[dropIdx] = kk * 24.0 * sigma * D_adh * (cos(thetaR_ * deg2rad) - cos(thetaA_ * deg2rad)) / (std::pow(M_PI, 3));		// Maximum adhesion force given the receding and advancing angles

	// Gravitational force [N] in X-direction
	g_force_X[dropIdx] = dropMass[dropIdx] * g.value().x();

	// Ratio between total drag and adhesion to check if the droplet is in adhesion or slide/roll condition
	if (adhesionForce[dropIdx] != 0)
	{
		ratio[dropIdx] = std::abs((totalDragForce[dropIdx][0]+ g_force_X[dropIdx]) / adhesionForce[dropIdx]);
	}
	else	// Droplet is not in contact with wall!
	{
		ratio[dropIdx] = 999;
	}

	//Info << "Interface evaluated for Droplet# " << dropIdx + 1 << ":\t" << dA_TOT[dropIdx] * pow(10,6) << "\t [mm2]" << endl;

}