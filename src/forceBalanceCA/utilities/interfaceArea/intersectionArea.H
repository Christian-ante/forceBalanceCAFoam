// FUNCTION TO CALCULATE THE WATER-AIR INTERSECTION AREA FOR EACH CELL


using VectorList = std::vector<std::vector<double>>;	// Declaration for a list of vectors


// 2D CASE
double intersectionArea (VectorList intersections, double l)
{
	double dA = 0;
	double height = l;
	double base = 	std::sqrt(std::pow(intersections[0][0] - intersections[1][0], 2) + std::pow(intersections[0][1] - intersections[1][1], 2));
	dA = base * height;

	return dA;
}


// 3D CASE
// The idea behind the 3D implementation is to divide the geoemtric shape in triangles and then calculate the interface area as a sum of triangles.
// To do so, the first step is to correctly order the vertexes in clockwise direction.
// Then the area is simply the sum of a varying number of triangles.

// Triangle area given the vertexes coordinate
double dA_(const std::vector<double>& A_, const std::vector<double>& B_, const std::vector<double>& C_)
{
	// Find AB and AC vectors (two edges of the triangle) 
	std::vector<double> AB = { B_[0] - A_[0], B_[1] - A_[1], B_[2] - A_[2] };
	std::vector<double> AC = { C_[0] - A_[0], C_[1] - A_[1], C_[2] - A_[2] };

	// Cross product between AB and AC
	std::vector<double> crossProduct = {
		AB[1] * AC[2] - AB[2] * AC[1],
		AB[2] * AC[0] - AB[0] * AC[2],
		AB[0] * AC[1] - AB[1] * AC[0]
	};

	// Magnitude of the resulting vector
	double crossMagnitude = std::sqrt(crossProduct[0] * crossProduct[0] +
										crossProduct[1] * crossProduct[1] +
										crossProduct[2] * crossProduct[2]);

	// Triangle area
	return 0.5 * crossMagnitude;
}

// Evaluate polar angle between a reference point and each vertex
double polarAngle(const std::vector<double>& p, const std::vector<double>& refPoint)
{
	double dx = p[0] - refPoint[0];
	double dy = p[1] - refPoint[1];
	return std::atan2(dy, dx);
}

// Function to compare (and sort) the vertexes
bool comparePolarAngles(const std::vector<double>& p1, const std::vector<double>& p2, const std::vector<double>& refPoint)
{
	if (p1 == refPoint)		// I want the first vertex of the list to be the refPoint! From there, the other vertexes will be in anti-clockwise order
	{
		return 1;
	}
	else if (p2 == refPoint)
	{
		return 0;
	}
	else
	{
		return polarAngle(p1, refPoint) < polarAngle(p2, refPoint);
	}
}

// Calculate interface area
double intersectionArea(std::vector<std::vector<double>>& intersections)
{
	double dA = 0;

	if (intersections.size() == 3)	// Cut area: triangle
	{
		// Find the vertexes of triangle
		std::vector<double> A_ = intersections[0];
		std::vector<double> B_ = intersections[1];
		std::vector<double> C_ = intersections[2];

		dA = dA_(A_, B_, C_);
	}
	else
	{
		// Reference point to evaluate polar angles, the one with lower X component
		std::vector<double> refPoint = intersections[0];
		for (const auto& point : intersections)
		{
			if (point[0] < refPoint[0])
				refPoint = point;
		}

		// Sort points based on polar angles relative to refPoint
		std::sort(intersections.begin(), intersections.end(), [&](const std::vector<double>& p1, const std::vector<double>& p2) {
			return comparePolarAngles(p1, p2, refPoint);
		});

		if (intersections.size() == 4)	// Cut area: quadrilateral
		{
			std::vector<double> A_ = intersections[0];
			std::vector<double> B_ = intersections[1];
			std::vector<double> C_ = intersections[2];
			std::vector<double> D_ = intersections[3];

			dA = dA_(A_, B_, D_) + dA_(B_, C_, D_);
		}
		else if (intersections.size() == 5)	// Cut area: pentagon
		{
			std::vector<double> A_ = intersections[0];
			std::vector<double> B_ = intersections[1];
			std::vector<double> C_ = intersections[2];
			std::vector<double> D_ = intersections[3];
			std::vector<double> E_ = intersections[4];

			dA = dA_(A_, B_, E_) + dA_(B_, C_, D_) + dA_(B_, D_, E_);
		}
		else if (intersections.size() == 6)	// Cut area: hexagon
		{
			std::vector<double> A_ = intersections[0];
			std::vector<double> B_ = intersections[1];
			std::vector<double> C_ = intersections[2];
			std::vector<double> D_ = intersections[3];
			std::vector<double> E_ = intersections[4];
			std::vector<double> F_ = intersections[5];

			dA = dA_(A_, B_, C_) + dA_(A_, C_, D_) + dA_(A_, D_, F_) + dA_(D_, E_, F_);
		}
	}

	return dA;
}