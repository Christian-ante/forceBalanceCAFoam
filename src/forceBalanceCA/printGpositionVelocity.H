if (Pstream::master()) // True only for the root processor
{
	static double lastWrittenTime = -1.0; // Keep track of the last written timestep
	double currentTime = mesh.time().value();

	// Write only if currentTime = 0 or output time of the simulation, and if it has not already been written
	if ((mesh.time().outputTime() || currentTime == 0) && currentTime != lastWrittenTime)
	{
		lastWrittenTime = currentTime; // Update last written timestep

		// Create or append to txt file for each droplet
		for (int dropIdx = 0; dropIdx < dropNumber; dropIdx++) // Loop for each droplet
		{
			std::string fileName = "CenterOfMass_" + std::to_string(dropIdx + 1) + ".txt";
			std::ofstream outFile(fileName, std::ios::out | std::ios::app); // Open file in append mode
			if (outFile.is_open())
			{
				if (outFile.tellp() == 0) // Check if file is empty (first-time write)
				{
					outFile << "Time [s]\tCenter Of Mass Position [m]\tCenter Of Mass Velocity [m/s]\n";
				}
				outFile << std::fixed << std::setprecision(6)
					<< currentTime << "\t"
					<< vectorToString(GCenter[dropIdx]) << "\t"
					<< vectorToString(velocityMean[dropIdx]) << "\n";
				outFile.close();
			}
		}
	}
}