// FUNCTION TO FIND THE PROPER "D" COEFFICIENT FOR THE PLANE THAT APPROXIMATES THE WATER-AIR INTERFACE
// THE SAME FUNCTION WORKS FOR 2D OR 3D CASES WITH DIFFERENT NUMBER OF ARGUMENTS


using Point = std::vector<double>;						// Declaration for a point
using PointCloud = std::vector<Point>;					// Declaration for a list of points (cloud)


// 2D CASE: y = Ax + D
double findInterfacePlane(double l, double A, double D_min, double D_max, int maxIterations, double tolerance, double targetAlpha, const PointCloud& points)
{
	double D = (D_max + D_min) / 2;			// Initial guess for D
	double alpha_temp;						// Alpha of the cell to compare with target

	// Binary research to find a D close enough to the target
	for (int iteration = 0; iteration < maxIterations; ++iteration) 
	{
        int pointsBelow = 0;				// Points below the cutting plane/line
		int pointsTotal = 0;				// Total number of points

        for (const auto& point : points)	// Cycle through all the points in the cloud
		{ 
			// Coordinates of the points cloud
           	double xx = point[0];			// X component of the point
			double yy = point[1];			// Y component of the point
			double y_plane = A * xx + D;	// Value of y_plane as a function of D
			pointsTotal++;

			if (yy < y_plane)				// If the point is under the plane
			{
	     		pointsBelow++;				// Increase counter
            }
        }

        alpha_temp = static_cast<double>(pointsBelow) / pointsTotal;	// Evaluate the fraction of points under the plane (over all the points)
        double err = alpha_temp - targetAlpha;							// Difference between alpha_temp (is a function of D) and target (volume fraction of the cell, retrieved from simulation)
        
        if (err > 0)		// Overshoot	->	Cut the top half of the domain
		{
			D_max = D;
        }
		else				// Undershoot	->	Cut the lower half of the domain
		{
            D_min = D;
        }

		if (std::abs(err) <= tolerance)		// Break the for cycle and return D if the error is lower than the required tolerance
		{
            break;
		}

		D = (D_max + D_min) / 2;			// Refresh D value over the new domain

	}

    return D;
}


// 3D CASE: z = Ax + By + D
double findInterfacePlane(double l, double A, double B, double D_min, double D_max, int maxIterations, double tolerance, double targetAlpha, const PointCloud& points)
{
	double D = (D_max + D_min) / 2;			// Initial guess for D
	double alpha_temp;						// Alpha of the cell to compare with target

	// Binary research to find a D close enough to the target
	for (int iteration = 0; iteration < maxIterations; ++iteration) 
	{
		int pointsBelow = 0;				// Points below the cutting plane/line
		int pointsTotal = 0;				// Total number of points

        for (const auto& point : points)	// Cycle through all the points in the cloud
		{ 
			// Coordinates of the points cloud
            double xx = point[0];			// X component of the point
            double yy = point[1];			// Y component of the point
            double zz = point[2];			// Z component of the point
            double z_plane = A * xx + B * yy + D;		// Value of z_plane as a function of D
			pointsTotal++;

            if (zz < z_plane)		// If the point (from points cloud) is under the plane
			{
                pointsBelow++;		// Increase counter
            }
        }

		alpha_temp = static_cast<double>(pointsBelow) / pointsTotal;	// Evaluate the fraction of points under the plane (over all the points)
		double err = alpha_temp - targetAlpha;							// Difference between alpha_temp (is a function of D) and target (volume fraction of the cell, retrieved from simulation)

		if (err > 0)		// Overshoot	->	Cut the top half of the domain
		{
			D_max = D;
		}
		else				// Undershoot	->	Cut the lower half of the domain
		{
			D_min = D;
		}

		if (std::abs(err) <= tolerance)		// Break the for cycle and return D if the error is lower than the required tolerance
		{
			break;
		}

		D = (D_max + D_min) / 2;			// Refresh D value over the new domain

	}

	return D;
}