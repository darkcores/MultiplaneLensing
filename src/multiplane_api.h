#pragma once

#include <vector>

extern "C" {

#ifndef MULTIPLANE_INTERNAL_USE
typedef void* MultiPlaneCtx;

// Gewoon iets voor een X/Y positie. Hier worden floats
// gebruikt, wat wil zeggen dat die waarden uitgedrukt worden
// in de opgegeven 'angularUnit' (zie verder)
struct MPCUXY
{
	float X, Y;
};

// Vaste parameters voor de plummer basisfuncties. Aan de massa
// zullen we verderop dan iets veranderen door een extra factor
// op te geven
struct MPCUInitialPlummerInfo
{
	MPCUXY position;
	float width; // ook weer uitgedrukt in 'angularUnit' (dus bvb X*angularUnit (als double) geeft de echte waarde)
	double initialMass; // dit is gewoon in kg uitgedrukt.
};
#endif

// De functies van de API; een 'int' als returnwaarde is negatief bij een
// error, en 0 als alles ok ging

// Initialiseer een MultiPlaneCtx (wordt via laatste parameter teruggegeven)
int mpcuInitMultiPlaneCalculation(
		// Hoeken (dingen in XY struct) zijn uitgedrukt in deze eenheid,
		// dwz als veelvouden van deze waarde
		double angularUnit,

		// Kosmologie parameters	
		double h, double W_m, double W_r, double W_v, double w,

		// De roodverschuivingen van de lensvlakken
		const std::vector<float> &lensRedshifts,

		// Voor elk lensvlak, de startparameters van alle plummers.
		const std::vector<std::vector<MPCUInitialPlummerInfo>> &fixedPlummerParameters, 

		// De roodverschuivingen van alle bronnen
		const std::vector<float> &sourceRedshifts,

		// Per bron, de theta-vectoren (ook uitgedrukt in 'angularUnit')
		const std::vector<std::vector<MPCUXY>> &theta, 

		// Als initialisatie werkte
		MultiPlaneCtx *pCtx);

// Voor een geinitialiseerde context, pas de initiele massa's aan volgens 
// massFactors en bereken de beta-vectoren voor alle bronnen
int mpcuCalculateSourcePositions(MultiPlaneCtx ctx, const std::vector<std::vector<float>> &massFactors);

// Voor de bron met index 'srcIdx', geef de berekende beta-posities terug
const std::vector<MPCUXY> &mpcuGetSourcePositions(MultiPlaneCtx ctx, int srcIdx);

// Alles weer opruimen
void mpcuClearContext(MultiPlaneCtx ctx);

} // extern "C"

