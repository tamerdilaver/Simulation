#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <maximize.h>
#import <modelbase>
#import <simula>
#include <oxfloat.h>

static decl iB;	 					//Repeats
static decl iSIZE;					//Size of time series
static decl iSTEPS;					//#Steps to divide the size
static decl iSIMS;					//# of Zt ~ N(0,1)
static decl iPLOT;
static decl iDIST;
static decl dALPHA;
static decl dBETA;
static decl dOMEGA;
static decl dGAMMA;
static decl iPARS;					//number of parameters
static decl vSTD_ALD;				// Zt ~ ALD(0, 1, p)
static decl vSTD_AST;				// Zt ~ AST(0, 1, ,nu, p)
static decl vRANT;				// Zt ~ ALD(0, b, p)
static decl s_vY; 					//Simulated returns
//static decl	bSTD_ERROR;				//0 or 1 boalean



findicator(const vX){
	decl vReturn, iN;
	iN = sizerc(vX);

	if(rows(vX)!=1){
		vReturn = zeros(iN,1);
	}else{
		vReturn = zeros(1,iN);
	}
	
//	print(vReturn);
	
	for(decl i=0;i<iN;i++){
		if(vX[i]<=0){
			vReturn[i]=1;
		}
	}
	return vReturn; 
}

fSign(const vX){
	decl vReturn, iN;
	iN = sizerc(vX);

	if(rows(vX)!=1){
		vReturn = ones(iN,1);
	}else{
		vReturn = ones(1,iN);
	}
	
//	print(vReturn);
	
	for(decl i=0;i<iN;i++){
		if(vX[i]<=0){
			vReturn[i]=-1;
		}
		if(vX[i]==0){
			vReturn[i]=0;
		}
	}
	return vReturn;
}

fK(const dLambda){
	decl dK;

	if(dLambda>320){
		return 1/sqrt(2*M_PI);	   //too large calculation this is the limit
	}

	dK = gammafact((dLambda+1)/2)/(gammafact(dLambda/2)*sqrt(dLambda*M_PI));

	return dK;
}

fSimALGAS(const dAlpha, const dBeta, const dOmega, const dGamma, const dP, const avReturns){
	decl vTemp, vH;
	vTemp = vH =  zeros(iSIZE+1, 1);

	vH[0]= dGamma;		//by definition
	
	for(decl i = 0; i < iSIZE; i++){	
		vTemp[i] =  sqrt(vH[i])*vSTD_ALD[i];
		vH[i+1] = dOmega + dBeta*vH[i] + dAlpha*((dP-findicator(vTemp[i]))*(2*sqrt(1-2*dP+2*sqr(dP))/((1-dP)*dP))*vTemp[i]*sqrt(vH[i])-2*vH[i]) ;
	}

	vTemp = dropr(vTemp,iSIZE);
	vH = dropr(vH,iSIZE);
	
	avReturns[0] = vH;
	return 1;
}

fSimASTGAS(const dAlpha, const dBeta, const dOmega,  const dGamma, const dLambda, const dP, const avReturns){
	decl vTemp, vH;
	vTemp = vH =  zeros(iSIZE+1, 1);

	vH[0]= dGamma; //by definition
	
	for(decl i = 0; i < iSIZE; i++){	
		vTemp[i] =  sqrt(vH[i])*vSTD_AST[i];
		vH[i+1] = dOmega + dAlpha*(dLambda+3)/dLambda*((dLambda+1)*sqrt(    4*((dP^3+(1-dP)^3)*(dLambda/(dLambda-2)))-16*fK(dLambda)^2*((1-2*dP)*(dLambda/(dLambda-1)))^2 	)/((dP-findicator(-vTemp[i]))^2*(4*dLambda))*(1+(    4*((dP^3+(1-dP)^3)*(dLambda/(dLambda-2)))-16*fK(dLambda)^2*((1-2*dP)*(dLambda/(dLambda-1)))^2 	)*sqr(vTemp[i])/(((dP-findicator(-vTemp[i]))^2*(4*dLambda))*vH[i]))^(-1)*sqr(vTemp[i])-vH[i]) + dBeta*vH[i];
	}					  																																					 //sqr(p-findicator(-vNews)) 

	vTemp = dropr(vTemp,iSIZE);
	vH = dropr(vH,iSIZE);
	
	avReturns[0] = vH;
	return 1;
}

main()
{

//	decl vRanExp, vRanExp2, vRanBern, vRanALD, vRanN;
//	decl dB;  //dB > 0
//	decl dP;  //0 < dP < 1
//	
//	dB = 1;
//	dP = 0.3;
//	iDIST = 10000000;
//	iSIZE = 1000;
//	iPLOT = 1000;
//	vRanExp = ranexp(iDIST, 1, 1);
//	vRanExp2 = ranexp(iDIST, 1, 1);
////	vRanBern = ranbinomial(iDIST, 1, 1, dP);
//	vRanALD = vRanExp/dP-vRanExp2/(1-dP);	 //~ALD(0,1,p)
//	vRanALD = dB*vRanALD;					 //~ALD(0,b,p)
//
//	vSTD_ALD = (1-dP)*dP/(dB*sqrt(1-2*dP+2*sqr(dP)))*vRanALD;
//	print(variance(vSTD_ALD));
//	vRanN = rann(iDIST,1);

	decl vRant1, vRant2, vRanu, vRanAST, vRann;

	iDIST = 10000000;
	iSIZE = 1000;
	iPLOT = 1000;

	decl dLambda;	// V>0	  (or V>2) to ensure variance!!!)
	decl dP;		// 0<alpha<1

	dLambda = 1000000;
	dP = 0.5;
	decl dOmega, dAlpha, dBeta, dGamma;
	dAlpha = .1;
	dBeta = 1;
	dOmega = 0.01;
	dGamma = 0.1;
	
	
	vRanu = ranu(iDIST, 1);
	vRann = rann(iDIST,1);
	vRant1 = rant(iDIST, 1, dLambda);
	vRant2 = rant(iDIST,1, dLambda);
	//print(vRant2[0:iPLOT]);

	vRanAST = dP* fabs(vRant1).*(fSign(vRanu-dP)-1) + (1-dP)*fabs(vRant2).*(fSign(vRanu-dP)+1) ;
	//print(vRanAST[0:iPLOT]);

	//print(fK(dLambda)^2);
	vSTD_AST = vRanAST/sqrt(    4*((dP^3+(1-dP)^3)*(dLambda/(dLambda-2)))-16*fK(dLambda)^2*((1-2*dP)*(dLambda/(dLambda-1)))^2 	);
	//print(vSTD_AST[0:iPLOT]); //hier kan hij niet verder dan 300
	
	//print(variance(vSTD_AST));
	
	decl vReturns;
	//fSimALGAS(0.1, 0.99, 0.01, 0.1, 0.45, &vReturns);
	//print(vReturns);

	fSimASTGAS(dAlpha, dBeta, dOmega, dGamma, dLambda, dP, &vReturns);

	SetDrawWindow("0_SIM_asymmetric_student's_t-GAS");
	DrawDensity(0, (vRann)', {"Density zt~ N(0,1)"});
	DrawDensity(0, (vSTD_AST)', {"Density zt~AST(0,p,nu,1)"});
	Draw(1, (vSTD_AST[0:iPLOT])');
	Draw(2, (vReturns)');
	ShowDrawWindow();
	//print(vRanALD);
}
