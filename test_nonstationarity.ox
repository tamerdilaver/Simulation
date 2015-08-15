#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <maximize.h>
#import <modelbase>
#import <simula>
#include <oxfloat.h>

static decl iSIZE;

fSimGARCH(const dAlpha, const dBeta, const dOmega, const dGamma, const avReturns){
	decl vTemp, vH;
	vTemp = vH =  zeros(iSIZE+1, 1);

	vH[0]= dGamma;		//by definition
	
	for(decl i = 0; i < iSIZE; i++){	
		vTemp[i] =  sqrt(vH[i])*rann(1,1);
		vH[i+1] = dOmega+ dBeta*vH[i] + dAlpha*sqr(vTemp[i]) ;
	}

	vTemp = dropr(vTemp,iSIZE);
	vH = dropr(vH,iSIZE);
	
	avReturns[0] = vTemp;
	return 1;
}

main()
{
	decl vR;

	iSIZE = 100000;

	fSimGARCH(0.1, 0.904, 0.01, 0.1, &vR);

	SetDrawWindow("test non-stationarity");
	Draw(0, vR');
	ShowDrawWindow();
}
