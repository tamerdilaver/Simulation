/*
**	Master Thesis Econometrics
**
**  Purpose:
**  	For some fixed parameter values simulate and estimate GARCH model parameters
**		with Maximum Likelikhood many times. s.t. Elog(alpha_0 z_t^2 + beta_0) = 0.
**
**  Date:
**    	30/05/2015
**
**  Author:
**	  	Tamer Dilaver
**
**	Supervisor:
**		Fransisco Blasques
**
*/



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
static decl dALPHA;
static decl dBETA;
static decl dOMEGA;
static decl dGAMMA;
static decl dLAMBDA;
static decl dPROB;
static decl iEXP;					//What constraint? GARCH/ t-GAS / al-GAS
static decl iPARS;					//number of parameters
static decl vSTD_NORM;				// Zt ~ N(0,1)
static decl vSTD_STUDENT_T;			// Zt ~ T(0,1,dLambda)
static decl vSTD_ALD;				//Zt ~ ALD(0,1,p)
static decl vSTD_AST;				//Zt ~ AST(0,1,b, p)
static decl s_vY; 					//Simulated returns
static decl	bSTD_ERROR;				//0 or 1 boalean


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


/*
**  Function:	Transform (start)parameters
**
**  Input: 		vTheta [parametervalues]
**
**  Output: 	vThetaStar
*/
fTransform2(const avThetaStar, const vTheta){
	avThetaStar[0]=		vTheta;
	
	avThetaStar[0][0] = log(vTheta[0]);

	return 1;
}

/*
**  Function:	Transform parameters back
**
**  Input: 		vThetaStar
**
**  Output: 	vTheta [parametervalues]
*/
fTransformBack2(const avTheta, const vThetaStar){
	avTheta[0]=		vThetaStar;

	avTheta[0][0] = exp(vThetaStar[0]);
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adBeta, vTheta
**
**  Output: 	1 
*/
fGetPars2(const adBeta, const vTheta){

	adBeta[0] = exp(vTheta[0]);
	return 1;
}

/*
**  Function:	Calculate targetvalue -[ E log(alpha_0 z_t^2 + beta_0) ]^2 for given parameters
**
**  Input:		vTheta [parametervalues], adFunc [adres functionvalue], avScore [the score],  amHessian [hessianmatrix]
**
**  Output:		1
**
*/
fExpectation_GARCH(const vTheta, const adFunc, const avScore, const amHessian){
	decl dBeta;
	fGetPars2( &dBeta, vTheta);

	//NOTICE: Since maxBFGS() is a function that maximimises I have squared the target function
	//and then I have put a negative sign in front of it. This ensures that maximising will
	//effectively search for the value for when the target function is zero.
	adFunc[0] = - (meanc(log(fabs(dALPHA *( vSTD_NORM.^2) + dBeta))))^2; 						
	return 1;
}

fExpectation_tGAS(const vTheta, const adFunc, const avScore, const amHessian){
	decl dBeta;
	fGetPars2( &dBeta, vTheta);

	//NOTICE: Since maxBFGS() is a function that maximimises I have squared the target function
	//and then I have put a negative sign in front of it. This ensures that maximising will
	//effectively search for the value for when the target function is zero.
	adFunc[0] = - (meanc(log(fabs(   dALPHA*((dLAMBDA+3)/dLAMBDA) *( (dLAMBDA+1)/sqr(dLAMBDA-2)*((1+(vSTD_STUDENT_T.^2)/(dLAMBDA-2)).^(-2)).* (vSTD_STUDENT_T.^4) ) +dBeta))))^2;
	return 1;
}

fExpectation_alGAS(const vTheta, const adFunc, const avScore, const amHessian){
	decl dBeta;
	fGetPars2( &dBeta, vTheta);

	//NOTICE: Since maxBFGS() is a function that maximimises I have squared the target function
	//and then I have put a negative sign in front of it. This ensures that maximising will
	//effectively search for the value for when the target function is zero.
	adFunc[0] = - (meanc(log(fabs(dALPHA*( 2*sqrt(1-2*dPROB+2*sqr(dPROB))/((1-dPROB)*dPROB)*(dPROB-findicator(vSTD_ALD)).*vSTD_ALD)+dBeta))))^2; 						
	return 1;
}

fExpectation_atGAS(const vTheta, const adFunc, const avScore, const amHessian){
	decl dBeta;
	fGetPars2( &dBeta, vTheta);

	//NOTICE: Since maxBFGS() is a function that maximimises I have squared the target function
	//and then I have put a negative sign in front of it. This ensures that maximising will
	//effectively search for the value for when the target function is zero.
	adFunc[0] = - (meanr(log(fabs(dALPHA*(dLAMBDA+3)/dLAMBDA*(     (dLAMBDA+1)./(4*dLAMBDA*sqr(dPROB-findicator(-vSTD_AST')) ) *            (( 4*(   (dPROB^3+(1-dPROB)^3)*(dLAMBDA/(dLAMBDA-2)))  -  16*fK(dLAMBDA)^2*((1-2*dPROB)*(dLAMBDA/(dLAMBDA-1)))^2 	) )          .* (1+1*sqr(vSTD_AST')*(( 4*(   (dPROB^3+(1-dPROB)^3)*(dLAMBDA/(dLAMBDA-2)))  -  16*fK(dLAMBDA)^2*((1-2*dPROB)*(dLAMBDA/(dLAMBDA-1)))^2 	) ) ./((4*dLAMBDA*sqr(dPROB-findicator(-vSTD_AST')) )*1)).^(-1).*(1*sqr(vSTD_AST')))+dBeta))))^2;
	return 1;					  //dALPHA*(dLAMBDA+3)/dLAMBDA*((dLAMBDA+1)*sqrt(    4*((dP^3+(1-dP)^3)*(dLAMBDA/(dLAMBDA-2)))-16*fK(dLAMBDA)^2*((1-2*dP)*(dLAMBDA/(dLAMBDA-1)))^2 	)/((dP-findicator(-vSTD_AST))^2*(4*dLAMBDA))*(1+sqrt(    4*((dP^3+(1-dP)^3)*(dLAMBDA/(dLAMBDA-2)))-16*fK(dLAMBDA)^2*((1-2*dP)*(dLAMBDA/(dLAMBDA-1)))^2 	)*sqr(vSTD_AST)/(((dP-findicator(-vSTD_AST))^2*(4*dLAMBDA))*1))^(-1)*sqr(vSTD_AST)-1) + dBeta;
}		

/*
**  Function:	Get the value for beta subject to Elog(alpha_0 z_t^2 + beta_0) = 0 for given parameters	and simulated Zt's
**
**  Input:		iSims, adBeta  
**
**  Output:		1
**
*/
fGetBeta(const iSims, const adBeta)
{
	decl vTheta, vThetaStart, vThetaStar, dFunc, iA;

	//initialise startparameter(s)
	vTheta = zeros(1,1);
	vTheta = <1.1>;			// dBeta
	vThetaStart = vTheta;

	//transform startparameter(s)
	fTransform2(&vThetaStar, vTheta);

	if(iEXP==0){
		//maximise
		iA=MaxBFGS(fExpectation_GARCH, &vThetaStar, &dFunc, 0, TRUE);
	}else if(iEXP==1){
		//maximise
		iA=MaxBFGS(fExpectation_tGAS, &vThetaStar, &dFunc, 0, TRUE);
	}else if(iEXP==2){
		//maximise
		iA=MaxBFGS(fExpectation_alGAS, &vThetaStar, &dFunc, 0, TRUE);
	}else if(iEXP==3){
		//maximise
		iA=MaxBFGS(fExpectation_atGAS, &vThetaStar, &dFunc, 0, TRUE);
	}
	
	//Transform thetasStar back
   	fTransformBack2(&vTheta, vThetaStar);

//	print("\nStart & Optimal parameter(s) with alpha fixed at ",dALPHA," and ",iSIMS," simulations such that we get I(1). \n",
//          "%r", { "dBeta"},
//          "%c", {"thetaStart","theta"}, vThetaStart~vTheta);

	if(fabs(dFunc)<1/1000000){
		adBeta[0] = vTheta[0];
  	}else{
		adBeta[0] = .NaN;
	}
	return 1;
}



/*
**				MAIN PROGRAM
**
**  Purpose:	Simulate GARCH returns for alpha, omega, gamma and beta many times.
**				Estimate GARCH parameters alpha, beta, omega and gamma.
**
**  Input: 		dALPHA, dBETA, dOMEGA, dGAMMA, iB, iSIZE, iSIMS, iSTEPS
**
**  Output: 	Figures
*/
main()
{
	//SET PARAMETERS
//	dALPHA = 0.1;
// 	dBETA = 0.89;
	dOMEGA = 0.01;
	dGAMMA = 0.1;
//	dGAMMA = dOMEGA/(1-dALPHA-dBETA); //doesn't work becomes negative
	iPARS = 4;


	//SET # OF SIMULATIONS 
	iB = 500; 			//max 5000
	iSIZE = 100;		//max 5000
	iSIMS = iB*iSIZE;
	vSTD_NORM = rann(iSIMS,1);
	



/*
** ..................................................................................	
**	 		SE - REGION
**	Get Betas for different alphas
**..................................................................................
*/


	decl dMin, dMax, dAmountOfDots, vNews, vBeta, mBeta;
	dMin = 0.01;
	dMax = 4.00;
	dAmountOfDots = 100;
	decl dBeta_0;  


	//GARCH
	vNews = vBeta = zeros(dAmountOfDots,1);
	vNews[0] = dMin;
	dALPHA = vNews[0];
	iEXP = 0;
	//fGetBeta(iSIMS,  &dBeta_0);
	vBeta[0] = 1;
	for(decl i = 1; i<dAmountOfDots; i++){
		vNews[i] = vNews[i-1]+(dMax-dMin)/dAmountOfDots;
		dALPHA = vNews[i];
		//fGetBeta(iSIMS,  &dBeta_0);

		if(1>0){
			vBeta[i] = 1;	  //
	   	}else{
			vBeta[i] = .NaN;
		}
	}
	mBeta = vBeta;
	print("\nfirst set finished \n");

	SetDrawWindow("SE II B.T. first set");
	decl iArea = 0;
	decl mYt = mBeta';    //decl mYt = (mBeta-vNews*(zeros(1,1)~ones(1,2)))';
	//print(vBeta);
	decl vX = vNews';
	DrawX(iArea, mYt, vX);
	ShowDrawWindow();

	//t-GAS
//	dLAMBDA = 2.1;
//	vSTD_STUDENT_T = sqrt((dLAMBDA-2)/dLAMBDA)*rant(iSIMS,1,dLAMBDA);	 //draw standardized student's t distributed 
//
//	vNews[0] = dMin;
//	dALPHA = vNews[0];
//	fGetBeta(iSIMS,  &dBeta_0);
//	vBeta[0] = dBeta_0;
//	iEXP = 1;
//	for(decl i = 1; i<dAmountOfDots; i++){
//		vNews[i] = vNews[i-1]+(dMax-dMin)/dAmountOfDots;
//		dALPHA = vNews[i];
//		fGetBeta(iSIMS,  &dBeta_0);
//		vBeta[i] = dBeta_0;
//	}
//	mBeta = mBeta~(vBeta+vNews*3/dLAMBDA);
//	//print(vBeta~vNews);
//	print("\nsecond set finished \n");

		//t-GAS
	dLAMBDA = 8;
	vSTD_STUDENT_T = sqrt((dLAMBDA-2)/dLAMBDA)*rant(iSIMS,1,dLAMBDA);	 //draw standardized student's t distributed 

	vNews[0] = dMin;
	dALPHA = vNews[0];
	iEXP = 1;
	fGetBeta(iSIMS,  &dBeta_0);
	vBeta[0] = dBeta_0;
	//print(vBeta[0] );
	
	for(decl i = 1; i<dAmountOfDots; i++){
		vNews[i] = vNews[i-1]+(dMax-dMin)/dAmountOfDots;
		dALPHA = vNews[i];
		fGetBeta(iSIMS,  &dBeta_0);
		vBeta[i] = dBeta_0;
	}
	mBeta = (vBeta+vNews*3/dLAMBDA);
	//print(vBeta~vNews);
	print("\nthird set finished \n");

	SetDrawWindow("SE II B.T. third set");
	 iArea = 0;
	 mYt = mBeta';    //decl mYt = (mBeta-vNews*(zeros(1,1)~ones(1,2)))';

	 print(vBeta~vNews);
	 
	//print(vBeta);
	 vX = vNews';
	DrawX(iArea, mYt, vX);
	ShowDrawWindow();

//	//al-GAS
//	dPROB = 0.45;
//	decl  vRanExp, vRanExp2, vRanALD;
//	vRanExp = ranexp(iSIMS, 1, 1);
//	vRanExp2 = ranexp(iSIMS, 1, 1);
//	vRanALD = vRanExp/dPROB-vRanExp2/(1-dPROB);	 //~ALD(0,1,p)
//	vSTD_ALD = (1-dPROB)*dPROB/(sqrt(1-2*dPROB+2*sqr(dPROB)))*vRanALD;
//	
//	vNews[0] = dMin;
//	dALPHA = vNews[0];
//	fGetBeta(iSIMS,  &dBeta_0);
//	vBeta[0] = dBeta_0;
//	iEXP = 2;
//	for(decl i = 1; i<dAmountOfDots; i++){
//		vNews[i] = vNews[i-1]+(dMax-dMin)/dAmountOfDots;
//		dALPHA = vNews[i];
//		fGetBeta(iSIMS,  &dBeta_0);
//		vBeta[i] = dBeta_0;
//	}
//	mBeta = (vBeta+vNews);
//	print("\nfourth set finished \n");
//
//	SetDrawWindow("Stationarity Regions B.T. fourth set");
//	 iArea = 0;
//	 mYt = mBeta';    //decl mYt = (mBeta-vNews*(zeros(1,1)~ones(1,2)))';
//	//print(vBeta);
//	 vX = vNews';
//	DrawX(iArea, mYt, vX);
//	ShowDrawWindow();
//
////	//al-GAS
////	dPROB = 0.65;
//////	decl  vRanExp, vRanExp2, vRanALD;
////	vRanExp = ranexp(iSIMS, 1, 1);
////	vRanExp2 = ranexp(iSIMS, 1, 1);
////	vRanALD = vRanExp/dPROB-vRanExp2/(1-dPROB);	 //~ALD(0,1,p)
////	vSTD_ALD = (1-dPROB)*dPROB/(sqrt(1-2*dPROB+2*sqr(dPROB)))*vRanALD;
////	
////	vNews[0] = dMin;
////	dALPHA = vNews[0];
////	fGetBeta(iSIMS,  &dBeta_0);
////	vBeta[0] = dBeta_0;
////	iEXP = 2;
////	for(decl i = 1; i<dAmountOfDots; i++){
////		vNews[i] = vNews[i-1]+(dMax-dMin)/dAmountOfDots;
////		dALPHA = vNews[i];
////		fGetBeta(iSIMS,  &dBeta_0);
////		vBeta[i] = dBeta_0;
////	}
////	mBeta = mBeta~(vBeta+vNews);
////	print("\nfifth set finished \n");
//
//	//at-GAS
//	dPROB = 0.45;
//	dLAMBDA = 8;
//	decl vRant1, vRant2, vRanu, vRanAST;
//	
//	vRanu = ranu(iSIMS, 1);
//	vRant1 = rant(iSIMS, 1, dLAMBDA);
//	vRant2 = rant(iSIMS,1, dLAMBDA);
//
//	vRanAST = dPROB* fabs(vRant1).*(fSign(vRanu-dPROB)-1) + (1-dPROB)*fabs(vRant2).*(fSign(vRanu-dPROB)+1) ;
//	vSTD_AST = vRanAST/sqrt(4*((dPROB^3+(1-dPROB)^3)*(dLAMBDA/(dLAMBDA-2)))-16*fK(dLAMBDA)^2*((1-2*dPROB)*(dLAMBDA/(dLAMBDA-1)))^2);
//	
//	vRanExp = ranexp(iSIMS, 1, 1);
//	vRanExp2 = ranexp(iSIMS, 1, 1);
//	vRanALD = vRanExp/dPROB-vRanExp2/(1-dPROB);	 //~ALD(0,1,p)
//	vSTD_ALD = (1-dPROB)*dPROB/(sqrt(1-2*dPROB+2*sqr(dPROB)))*vRanALD;
//	
//	vNews[0] = dMin;
//	dALPHA = vNews[0];
//	fGetBeta(iSIMS,  &dBeta_0);
//	vBeta[0] = dBeta_0;
//	iEXP = 3;
//	for(decl i = 1; i<dAmountOfDots; i++){
//		vNews[i] = vNews[i-1]+(dMax-dMin)/dAmountOfDots;
//		dALPHA = vNews[i];
//		fGetBeta(iSIMS,  &dBeta_0);
//		vBeta[i] = dBeta_0;
//	}
//	mBeta = (vBeta+vNews*3/dLAMBDA);
//	print("\nsixth set finished \n");
//	
//
//	//GET BETA SUCH THAT WE GET INEQUALITY "Elog(alpha_0 z_t^2 + beta_0) > 0".
//	
//	SetDrawWindow("Stationarity Regions B.T. sixth set");
//	 iArea = 0;
//	 mYt = mBeta';    //decl mYt = (mBeta-vNews*(zeros(1,1)~ones(1,2)))';
//	//print(vBeta);
//	 vX = vNews';
//	DrawX(iArea, mYt, vX);
//	ShowDrawWindow();
//

 

}

