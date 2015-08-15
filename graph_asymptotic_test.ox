#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#import <modelbase>
#import <simula>

main()
{	decl vReturns;

	decl iSize = 1000	;

	vReturns = rann(iSize,20);
	decl vMin, vMax, vMean;
	vMin = minc(vReturns')';
	vMean = meanc(vReturns')';
	vMax = maxc(vReturns')';

//	SetDrawWindow("GRAPH TEST");
//	Draw(0, (vMin~vMean~vMax)',0,1);	
//	ShowDrawWindow();

	SetDrawWindow("DISTRIBUTION GRAPH TEST ");
	DrawDensity(0, (vReturns[1000-1][]), {"Density at 1000"});
	DrawDensity(1, (vReturns[100-1][]), {"Density at 100"});
	DrawDensity(2, (vReturns[10-1][]), {"Density at 10"});
	DrawDensity(3, (vReturns[1-1][]), {"Density at 1"});
	ShowDrawWindow();


}
