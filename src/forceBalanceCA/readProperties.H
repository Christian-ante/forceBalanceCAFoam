const fvMesh& mesh = this->patch().boundaryMesh().mesh();

IOdictionary transportProperties	// Read constant/transportProperties file
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

uniformDimensionedVectorField g     // Read gravity field (usually constant, but should also be able to handle rotating domain)
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


double def = -999;																				// Default value if the code cannot read "transportProperties"
double sigma = transportProperties.lookupOrDefault<scalar>("sigma", def);						// Read superficial tension from transportProperties
double rho_air = transportProperties.subDict("air").lookupOrDefault<scalar>("rho", def);		// Read air density from transportProperties
double rho_water = transportProperties.subDict("water").lookupOrDefault<scalar>("rho", def);	// Read water density from transportProperties
double nu_air = transportProperties.subDict("air").lookupOrDefault<scalar>("nu", def);			// Read air kinematic viscosity from transportProperties
double nu_water = transportProperties.subDict("water").lookupOrDefault<scalar>("nu", def);		// Read water kinematic viscosity from transportProperties

if (sigma == def || rho_air == def || rho_water == def || nu_air == def || nu_water == def)
{
    FatalErrorInFunction
        << "Cannot read material properties (water/air rho and nu, sigma)!"
        << nl << "Please check your transportProperties file."
        << exit(FatalError);
}

double mu_water = rho_water * nu_water;		// Water dynamic viscosity [Pa s] (used for the capillary number)
double mu_air = rho_air * nu_air;			// Air dynamic viscosity [Pa s] (used for the viscous drag, evaluated on the air side of the interface)
