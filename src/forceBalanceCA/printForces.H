if (Pstream::master())											// True only for the root processor
{
	static double lastWrittenTime_2 = -1.0;						// Keep track of the last written timestep
	double currentTime = mesh.time().value();
	
	if ((mesh.time().outputTime() || currentTime == 0) && currentTime != lastWrittenTime_2)		// Write only if currentTime = 0 or output time of the simulation, and if it has not already been written
	{
		lastWrittenTime_2 = currentTime;						// Update last written timestep

		for (int dropIdx = 0; dropIdx < dropNumber; dropIdx++)	// Create or append to txt file for each droplet
		{
			std::string fileName = "ForcesOnDroplet_" + std::to_string(dropIdx + 1) + ".txt";
			std::ofstream outFile(fileName, std::ios::out | std::ios::app);		// Open file in append mode
			if (outFile.is_open())
			{
				if (outFile.tellp() == 0)						// Check if file is empty (first-time write --> add header)
				{
					outFile << "Time [s]\tPressure Force [N]\tViscous Force [N]\tG Force X [N]\tTotal Drag X-dir [N]\tAdhesion X-dir [N]\tRatio[-]\n";
				}
				outFile << std::fixed << std::setprecision(8)
					<< currentTime <<								"\t"
					<< vectorToString(pressureForce[dropIdx]) <<	"\t"
					<< vectorToString(viscousForce[dropIdx]) <<		"\t"
					<< g_force_X[dropIdx] <<						"\t"
					<< totalDragForce[dropIdx][0] <<				"\t"
					<< adhesionForce[dropIdx] <<					"\t"
					<< ratio[dropIdx] <<							"\n";
				outFile.close();
			}
		}
	}
}