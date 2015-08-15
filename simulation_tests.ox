/*
**	Master Thesis Econometrics
**
**  Purpose:
**  	For some fixed parameter values simulate and estimate GARCH model parameters
**		with Maximum Likelikhood many times.
**
**  Date:
**    	14/07/2014
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
static decl iSIMS;					//#of Zt~N(0,1)
static decl dALPHA;
static decl dOMEGA;
static decl dGAMMA;
static decl vSTD_NORM;
static decl s_vY; 					//Simulated returns
//static decl dBETA = 0.9081;		//brute force (grid) search 5min (temporary solution)

fTransform(const avThetaStar, const vTheta){
	avThetaStar[0]=		vTheta;
	
	avThetaStar[0][0] = log(vTheta[0]);

	return 1;
}

fTransformBack(const avTheta, const vThetaStar){
	avTheta[0]=		vThetaStar;

	avTheta[0][0] = exp(vThetaStar[0]);
	return 1;
}

fGetPars(const adBeta, const vTheta){

	adBeta[0] = exp(vTheta[0]);
	return 1;
}

fExpectation(const vTheta, const adFunc, const avScore, const amHessian){
	decl dBeta;
	fGetPars( &dBeta, vTheta);

	adFunc[0] = - (meanc(log(fabs(dALPHA * vSTD_NORM.^2 + dBeta))))^2; 						
	return 1;
}


fGetBeta(const iSims, const adBeta)
{
	decl vTheta, vThetaStart, vThetaStar, dFunc, iA;

	//initialise startparameter(s)
	vTheta = zeros(1,1);
	vTheta = <0.9>;			// dBeta
	vThetaStart = vTheta;

	//transform startparameter(s)
	fTransform(&vThetaStar, vTheta);

	//Calculate functionvalue at startparameter(s)
//	fExpectation(vThetaStar, &dFunc, 0, 0);
//	println("\nfunctionvalue at startparameter(s):", dFunc);

	iA=MaxBFGS(fExpectation, &vThetaStar, &dFunc, 0, TRUE);
	
	//Transform thetasStar back
   	fTransformBack(&vTheta, vThetaStar);

	//print optimal parameters, functionvalue at final parameter(s), time and return adBeta.
//	print("\n",MaxConvergenceMsg(iA));
//	println("\nfunctionvalue at final parameter(s):", dFunc);
	print("\nStart & Optimal parameter(s) with alpha fixed at ",dALPHA," and ",iSIMS," simulations such that we get I(1). \n",
          "%r", { "dBeta"},
          "%c", {"thetaStart","theta"}, vThetaStart~vTheta);
	adBeta[0] = vTheta[0];
//   	print("\n",time());
	
	//extra check 
//	print("\nmeanc(log(fabs(dALPHA * vSTD_NORM.^2 + dBeta))) = ",meanc(log(fabs(dALPHA * vSTD_NORM.^2 + vTheta[0]))));

	
	//draw some densities
//	SetDrawWindow("Density tests");
//	DrawDensity(0, vSTD_NORM', "test");
//	DrawDensity(1, (log(fabs(dALPHA * vSTD_NORM.^2 + vTheta[0])))', "test");
//	ShowDrawWindow();
}

fSimGARCH(const dAlpha, const dBeta, const dOmega, const dGamma, const avReturns, const iIteration){
	decl vTemp, vH;
	vTemp = vH =  zeros(iSIZE+1, 1);

	vH[0]= dGamma;		//by definition
	
	for(decl i = 0; i < iSIZE; i++){	
		vTemp[i] =  sqrt(vH[i])*vSTD_NORM[(i + (iIteration*iSIZE))];
		vH[i+1] = dOmega + dAlpha*sqr(vTemp[i]) + dBeta*vH[i];
	}

	vTemp = dropr(vTemp,iSIZE);
	vH = dropr(vH,iSIZE);
	
	avReturns[0] = vTemp;
	return 1;
}


fTransform2(const avThetaStar, const vTheta){
	avThetaStar[0]=		vTheta;
	
	avThetaStar[0][0] = log(vTheta[0]);
	avThetaStar[0][1] = log(vTheta[1]);
	return 1;
}

fGetPars2(const adAlpha, const adBeta, const vTheta){

	adAlpha[0] = exp(vTheta[0]);
	adBeta[0] = exp(vTheta[1]);
	return 1;
}

fLogLike_Garch(const vTheta, const adFunc, const avScore, const amHessian){

	decl dAlpha, dBeta;
	fGetPars2( &dAlpha,  &dBeta, vTheta);

	decl dS2 = dGAMMA;					 											//initial condition by definition
	decl vEta = zeros(sizerc(s_vY), 1);

	for(decl i = 0; i < sizerc(s_vY); ++i){
	
			//likelihood contribution
			vEta[i] = 1 / (sqrt(M_2PI * dS2)) * exp(-(s_vY[i])^2 / (2 * dS2));		//Gaussian

			//GARCH recursion
			dS2 = dOMEGA + dBeta * dS2 + dAlpha * s_vY[i]^2;

	}
	
	adFunc[0] = sumc(log(vEta))/sizerc(s_vY); 									 	//Average
	return 1;
}

fTransformBack2(const avTheta, const vThetaStar){
	avTheta[0]=		vThetaStar;

	avTheta[0][0] = exp(vThetaStar[0]);
	avTheta[0][1] = exp(vThetaStar[1]);
	return 1;
}

fEstimateGarch(const vReturns, const adAlpha_hat, const adBeta_hat, const dBeta_0){

	//initialise parametervalues
	decl vTheta = zeros(2,1);
	vTheta = <0.1 ; 0.9>;			// Alpha, Beta Startingvalues
	decl vThetaStart = vTheta;

	//globalalize returns and vectorize true pars
	s_vY = vReturns;
	decl vTheta_0  = zeros(2,1);
	vTheta_0[0] = dALPHA;
	vTheta_0[1] = dBeta_0;
	
	decl vThetaStar; 
	//transform parameters
	fTransform2(&vThetaStar, vTheta);

	decl dFunc;
//	//Calculate LogLikelihood-value at startparameter(s)
//	fLogLike_Garch(vThetaStar, &dFunc, 0, 0);
//	println("\nLoglikelihood-value at startparameter(s):", dFunc);

	decl iA;
	//Maximize the LL
	iA=MaxBFGS(fLogLike_Garch, &vThetaStar, &dFunc, 0, TRUE);

	//Transform thetasStar back
  	fTransformBack2(&vTheta, vThetaStar);

	//print optimal parameters, LL optimal parameters, time.
//	print("\n",MaxConvergenceMsg(iA));
//	println("\nLoglikelihood-value at final parameter(s):", dFunc);
//	print("\nTrue-, Start- & Optimal parameters \n",
//          "%r", { "dAlpha",  "dBeta"},
//          "%c", {"theta_0","thetaStart","theta"}, vTheta_0~vThetaStart~vTheta);
//   	print("\n",time());

	//return alpha and beta
	adAlpha_hat[0] = vTheta[0];
	adBeta_hat[0] = vTheta[1];

	return 1;
}

fMonteCarlo(const amMonteCarlo, dBeta_0){
	decl mTemp;
	mTemp = zeros(iB,2);

	for(decl i = 0; i<iB ; i++){
		decl vReturns;
		fSimGARCH(dALPHA, dBeta_0, dOMEGA, dGAMMA, &vReturns, i);

//		SetDrawWindow("Simulated Returns");
//		Draw(0, (vReturns)',0,1);	
//		ShowDrawWindow();

		decl dAlpha_hat, dBeta_hat;
		fEstimateGarch(vReturns, &dAlpha_hat, &dBeta_hat, dBeta_0);	 //Omega and Gamma fixed at true values

		mTemp[i][0] =  dAlpha_hat;
		mTemp[i][1]	=  dBeta_hat;
	}
	amMonteCarlo[0] = mTemp;

}

fMonteCarlo2(const amAlpha, amBeta, dBeta_0){
	decl mTemp, mTempAlpha, mTempBeta;
	mTempAlpha = mTempBeta = zeros((iSIZE/iSTEPS),3);
	mTemp = zeros(iB,2);

	decl iSize = iSIZE;

	for(decl j = 0; j<(iSize/iSTEPS) ; j++){
		iSIZE = ((iSTEPS)*(j+1));
	//	print("\n",iSIZE);
		for(decl i = 0; i<iB ; i++){
			decl vReturns;
			fSimGARCH(dALPHA, dBeta_0, dOMEGA, dGAMMA, &vReturns, i);
	
	//		SetDrawWindow("Simulated Returns");
	//		Draw(0, (vReturns)',0,1);	
	//		ShowDrawWindow();
	
			decl dAlpha_hat, dBeta_hat;
			fEstimateGarch(vReturns, &dAlpha_hat, &dBeta_hat, dBeta_0);	 //Omega and Gamma fixed at true values
	
			mTemp[i][0] =  dAlpha_hat;
			mTemp[i][1]	=  dBeta_hat;
		}
		// vMin, vMean, vMax;
	//	print(minc(mTemp[][])'[0]);
		mTempAlpha[j][0] = minc(mTemp[][])'[0];
		mTempAlpha[j][1] = meanc(mTemp[][])'[0];
		mTempAlpha[j][2] = maxc(mTemp[][])'[0];
	
		mTempBeta[j][0] = minc(mTemp[][])'[1];
		mTempBeta[j][1] = meanc(mTemp[][])'[1];
		mTempBeta[j][2] = maxc(mTemp[][])'[1];
	}

	amAlpha[0] = mTempAlpha;
	amBeta[0] = mTempBeta;
}


main()
{
	//Set Parameters
	dALPHA = 0.1;
	dOMEGA = 0.1;
	dGAMMA = 0.1;

	//Set number of simulations 
	iB = 100;
	iSIZE = 1000;
	iSIMS = iB*iSIZE;
	vSTD_NORM = rann(iSIMS,1);

	//Get beta such that we get equality 
	decl dBeta_0;
	fGetBeta(iSIMS,  &dBeta_0);
	
	//Get distributions of alpha and beta (to check for asymptotic normality)
	decl mMonteCarlo;
	fMonteCarlo(&mMonteCarlo, dBeta_0);	  

	SetDrawWindow("Distribution Graph alpha and beta (to show asymptotic normality)");
	DrawDensity(0, (mMonteCarlo[][0])', {"Density alpha"});
	DrawDensity(1, (mMonteCarlo[][1])', {"Density beta"});
	ShowDrawWindow();

	//Check consistency for alpha and beta
	decl mAlpha, mBeta;
	iSTEPS = iSIZE/10;				 //this takes a while (fast things up? decrease iB)
	fMonteCarlo2(&mAlpha, &mBeta, dBeta_0);

	SetDrawWindow("Graph Parameter Consistency");
	Draw(0, (mAlpha)',iSTEPS,iSTEPS);
	Draw(1, (mBeta)',iSTEPS,iSTEPS);
	DrawTitle(0,"(i) alpha");	
	DrawTitle(1,"(ii) beta");
	ShowDrawWindow();
}

