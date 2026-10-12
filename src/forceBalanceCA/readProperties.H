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
double mu_water = transportProperties.subDict("water").lookupOrDefault<scalar>("mu", def);		// Read dynamic viscosity from transportProperties
double sigma = transportProperties.lookupOrDefault<scalar>("sigma", def);						// Read superficial tension from transportProperties
double rho_air = transportProperties.subDict("air").lookupOrDefault<scalar>("rho", def);		// Read air density from transportProperties
double rho_water = transportProperties.subDict("water").lookupOrDefault<scalar>("rho", def);	// Read water density from transportProperties

if (mu_water == def || sigma == def || rho_air == def || rho_water == def)
{
    FatalErrorInFunction
        << "Cannot read material properties (water mu/rho, air rho and sigma)!"
        << nl << "Please check your transportProperties file."
        << exit(FatalError);
}
