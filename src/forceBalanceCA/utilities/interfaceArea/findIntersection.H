// FUNCTION TO FIND THE INTERSECTION BETWEEN INTERFACE PLANE AND CELL


using Point = std::vector<double>;						// Declaration for a point
using Edge = std::pair<Point, Point>;					// Declaration for a couple of points
using EdgeList = std::vector<Edge>;						// Declaration for a list of couples of points (edges)
using VectorList = std::vector<std::vector<double>>;	// Declaration for a list of vectors


// Utility function to check if intersection points are duplicate (cutting plane passing through a vertex!)
bool arePointsEqual(const std::vector<double>& p1, const std::vector<double>& p2)
{
	return (std::fabs(p1[0] - p2[0]) < Foam::small) && (std::fabs(p1[1] - p2[1]) < Foam::small);
} 

// Utility function to check if intersection points are duplicate (cutting plane passing through a vertex!) in 3D
bool arePointsEqual3D(const std::vector<double>& p1, const std::vector<double>& p2)
{
	return (std::fabs(p1[0] - p2[0]) < Foam::small) && (std::fabs(p1[1] - p2[1]) < Foam::small) && (std::fabs(p1[2] - p2[2]) < Foam::small);
}

// 2D CASE: y = Ax + D
VectorList findIntersection(double A, double D, const EdgeList edges)
{
	VectorList intersections;					// Vector to return intersection points

	for (auto &edge : edges)					// Cycle through all edges
	{
		const Point& p1 = edge.first;			// First vertex of the edge
        const Point& p2 = edge.second;			// Second vertex of the edge

		// Every point on the edge can be expressed as: y1 + t * (y2 - y1) = A * (x1 + t * (x2 - x1)) + D		(0 <= t <= 1)
		// Then solve the equation and find t (= t_num / t_den)
    	double t_num = A * p1[0] + D - p1[1];
    	double t_den = (p2[1] - p1[1]) - A * (p2[0] - p1[0]);

		if (t_den == 0)		// No intersections
		{
            continue;		// Skip to next edge
        }

        double t = t_num / t_den;

		if (t >= 0 && t <= 1) 
		{
			std::vector<double> intersect(2);
        	intersect[0] = p1[0] + t * (p2[0] - p1[0]);
        	intersect[1] = p1[1] + t * (p2[1] - p1[1]);

			// Check for duplicate intersection (plane cutting at the vertex)
			bool alreadyExists = false;
			for (const auto& existingIntersect : intersections)
			{
				if (arePointsEqual(existingIntersect, intersect))
				{
					alreadyExists = true;
					break;
				}
			}

			if (!alreadyExists)							// If the point is not dublicate
			{
				intersections.push_back(intersect);		// Add it to intersections
			}
		}
	}
	return intersections;
}

// 3D CASE: z = Ax + By + D
VectorList findIntersection(double A, double B, double D, const EdgeList edges)
{
	VectorList intersections;		// Vector to return intersection points

	for (auto &edge : edges) 
	{
        const Point& p1 = edge.first;		// First vertex of the edge
        const Point& p2 = edge.second;		// Second vertex of the edge

		// Every point on the edge can be expressed as:  z1 + t * (z2 - z1) = A * (x1 + t * (x2 - x1)) + B * (y1 + t * (y2 - y1)) + D		(0 <= t <= 1)
		// Then solve the equation and find t (= t_num / t_den)
    	double t_num = A * p1[0] + B * p1[1] + D - p1[2];
    	double t_den = (p2[2] - p1[2]) - A * (p2[0] - p1[0]) - B * (p2[1] - p1[1]);

		if (t_den == 0)		// No intersections
		{
            continue;		// Skip to next edge
        }

        double t = t_num / t_den;

		if (t >= 0 && t <= 1) 
		{
			std::vector<double> intersect(3);
        	intersect[0] = p1[0] + t * (p2[0] - p1[0]);
        	intersect[1] = p1[1] + t * (p2[1] - p1[1]);
        	intersect[2] = p1[2] + t * (p2[2] - p1[2]);

			// Check for duplicate intersection (plane cutting at the vertex)
			bool alreadyExists = false;
			for (const auto& existingIntersect : intersections)
			{
				if (arePointsEqual3D(existingIntersect, intersect))
				{
					alreadyExists = true;
					break;
				}
			}

			if (!alreadyExists)							// If the point is not dublicate
			{
				intersections.push_back(intersect);		// Add it to intersections
			}
		}
	}
	return intersections;
}