#include "context.h"

typedef MultiPlaneContext* MultiPlaneCtx;
typedef Vector2D<float> MPCUXY;
typedef PlummerParams MPCUInitialPlummerInfo;

#define MULTIPLANE_INTERNAL_USE
#include "multiplane_api.h"

using namespace std;

int mpcuInitMultiPlaneCalculation(
		double angularUnit,
		double h, double W_m, double W_r, double W_v, double w,
		const vector<float> &lensRedshifts,
		const vector<vector<MPCUInitialPlummerInfo>> &fixedPlummerParameters, 
		const vector<float> &sourceRedshifts,
		const vector<vector<MPCUXY>> &theta, 
		MultiPlaneCtx *pCtx)
{
	MultiPlaneContext *pContext = new MultiPlaneContext(angularUnit, Cosmology(h, W_m, W_r, W_v, w));
	int status = pContext->init(lensRedshifts, fixedPlummerParameters, sourceRedshifts);
	if (status < 0)
	{
		delete pContext;
		return status;
	}

	pContext->setThetas(theta);

	*pCtx = pContext;
	return 0;
}

int mpcuCalculateSourcePositions(MultiPlaneCtx ctx, const vector<vector<float>> &massFactors)
{
	return ctx->calculatePositions(massFactors);
}

const vector<MPCUXY> &mpcuGetSourcePositions(MultiPlaneCtx ctx, int srcIdx)
{
	return ctx->getSourcePositions(srcIdx);
}

void mpcuClearContext(MultiPlaneCtx ctx)
{
	delete ctx;
}

