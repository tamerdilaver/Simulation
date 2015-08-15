	#include <oxstd.h>
	#include <oxdraw.h>
	#include <oxprob.h>
	#import <modelbase>
	#import <simula>
	#include <oxfloat.h>

main(){
	decl vX, vY, mZ;


	//decl iCorrect, dSigma2_eta, dMu, dPhi ;

	//dMu = 0.0506437;		//we fix dMu because   log((81/(77+81))/(1-(81/(77+81)))) = 0.0506437 

	//decl  dPhi_LB, dPhi_Steps, dPhi_UB;
	//dPhi_LB = 0;
	//dPhi_Steps = 0.1;
	//dPhi_UB = 1;

//	decl  dSigma2_Eta_LB, dSigma2_Eta_Steps, dSigma2_Eta_UB;
//	dSigma2_Eta_LB = 0.1;
//	dSigma2_Eta_Steps = 0.1;
//	dSigma2_Eta_UB = 10; 

decl iN = 100;
	vX = zeros(iN,1); //vX are sigma2_eta
	vY = zeros(iN,1);
	
	mZ = zeros(sizerc(vX),sizerc(vY));
	
	for(decl i = 0; i<iN;i++){
		vX[i]  = i;
		
		for(decl w = 0; w<iN ; w++){
			vY[w] = w;
			mZ[i][w] = log(50^2-(i-50)^2)+log(50^2-fabs(w-50)^2);
		}
		 
	}

//	decl iMaxCorrect = maxc(maxc(mZ)');
//
//	print("\nMax Correct estimate: ", iMaxCorrect);
//	print("\nWith sigma2_eta = ",vX[maxcindex(maxc( mZ' )')], " dPhi = ", vY[maxcindex(maxc(mZ)')]);

	SetDrawWindow("3d plot");
	DrawXYZ( 0,vY, vX, mZ,0 , "vY", "vx", "vX+vY");
	ShowDrawWindow();
}
