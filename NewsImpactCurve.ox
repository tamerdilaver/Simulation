/*
**	Thesis MSc Econometrics 	
**
**  Purpose:
**  	Plot a news-impact curve
**
**  Date:
**    	27/05/2015
**
**  Author:
**	  	Tamer Dilaver
**
**	Supervisor:
**		F. Blasques & S.J. Koopman
**
*/

#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <maximize.h>
#import <modelbase>
#import <simula>
#include <oxfloat.h>

static decl dALPHA;			 		
static decl dBETA;
static decl dOMEGA;
static decl dGAMMA;
static decl iSIZE; 
static decl vSTD_NORM;


/*
**  Function:	Simulate GARCH returns for given parameters
**
**  Input:		dAlpha, dBeta, dOmega, dGamma, avReturns, iIteration [to get different Zt's]
**
**  Output:		1
**
*/

//fSimGARCH(const dAlpha, const dBeta, const dOmega, const dGamma, const avReturns, const iIteration){
//	decl vTemp, vH;
//	vTemp = vH =  zeros(iSIZE+1, 1);
//
//	vH[0]= dGamma;		//by definition
//	
//	for(decl i = 0; i < iSIZE; i++){	
//		vTemp[i] =  sqrt(1)*vSTD_NORM[(i + (iIteration*iSIZE))];
//		vH[i+1] = dOmega+ dBeta*1 + dAlpha*sqr(vTemp[i]) ;
//	}
//
//	vTemp = dropr(vTemp,iSIZE);
//	vH = dropr(vH,iSIZE);
//	
//	avReturns[0] = vTemp;
//	return 1;
//}

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

fK(const dLambda){
	decl dK;

	if(dLambda>320){
		return 1/sqrt(2*M_PI);	   //too large calculation this is the limit
	}

	dK = gammafact((dLambda+1)/2)/(gammafact(dLambda/2)*sqrt(dLambda*M_PI));

	return dK;
}

//
//         thetaStart        theta    std.error
//alpha       0.10000     0.075437     0.015393
//beta        0.85000      0.91110     0.018862
//omega      0.010000     0.010267    0.0040376

main()
{
	dALPHA 	= 0.075437;
	dBETA 	= 0.91110;
	dOMEGA	= 0.010267;
//	dGAMMA	= dOMEGA/(1-dALPHA-dBETA);
//	iSIZE	= 100;
//	vSTD_NORM = rann(iSIZE,1);
//
//	decl vReturns;
//	fSimGARCH(dALPHA, dBETA, dOMEGA, dGAMMA, &vReturns, 0);
//	print(vReturns);

	decl dMin, dMax, dAmountOfDots, vNews;
	dMin = -5;
	dMax = 5;
	dAmountOfDots = 1000;

	vNews = zeros(dAmountOfDots,1);
	vNews[0] = dMin;
	for(decl i = 1; i<dAmountOfDots; i++){
		vNews[i] = vNews[i-1]+(dMax-dMin)/dAmountOfDots;
	}

	

//take dGAMMA = 1 (= h_{t-1}) 
//	SetDrawWindow("NIC_GARCH_AND_ROBUST_GARCH");
	decl iArea = 0;
//	decl mYt = dOMEGA + dBETA*dGAMMA + dALPHA*dGAMMA*sqr(vNews)'|dOMEGA + dBETA*dGAMMA + dALPHA*sqr(vNews)';
//	decl mYt = 1*sqr(vNews)'|1*sqr(vNews)'./(1+1*sqr(vNews)');  //garch and Robust GARCH
//	decl vX = vNews';
//	DrawX(iArea, mYt, vX);
//	ShowDrawWindow();
//
//
//	decl dL = 4; //dLambda
	SetDrawWindow("NIC_GARCH_AND_t-GAS");
//	iArea = 0;
//	decl mYt = dOMEGA + dBETA*dGAMMA + dALPHA*dGAMMA*sqr(vNews)'|dOMEGA + dBETA*dGAMMA + dALPHA*sqr(vNews)';
	//decl mYt = dOMEGA + dALPHA*1*sqr(vNews)'+dBETA;  //garch and GAS

	
//                thetaStart        theta    std.error
//alpha              0.10000     0.042204     0.011147
//beta               0.99000      0.99650    0.0050631
//omega            0.0010000    0.0037325    0.0027661
//gamma              0.10000      0.31763      0.18743
//lambda              10.000       5.4912      0.80216

	
	dALPHA 	= 0.042204;
	dBETA 	= 0.99650;
	dOMEGA	= 0.0037325;
	decl dL = 5.4912;
	
//	decl mYt = dOMEGA + dALPHA*(dL+3)/dL*((dL+1)/(dL-2)*(1+1*sqr(vNews)'/((dL-2)*1)).^(-1).*(1*sqr(vNews)')-1)+(dBETA);  //GAS
//	dL = 10;
//	mYt = mYt|dOMEGA + dALPHA*(dL+3)/dL*((dL+1)/(dL-2)*(1+1*sqr(vNews)'/((dL-2)*1)).^(-1).*(1*sqr(vNews)')-1)+(dBETA+dALPHA);  //GAS
//	dL = 10000;
//	mYt = mYt|dOMEGA + dALPHA*(dL+3)/dL*((dL+1)/(dL-2)*(1+1*sqr(vNews)'/((dL-2)*1)).^(-1).*(1*sqr(vNews)')-1)+(dBETA+dALPHA);  //GAS
//	decl vX = vNews';
//	DrawX(iArea, mYt, vX);
//	ShowDrawWindow();
//
//	
//	SetDrawWindow("NIC_GARCH_AND_l-GAS");
//	iArea = 0;
////	decl mYt = dOMEGA + dBETA*dGAMMA + dALPHA*dGAMMA*sqr(vNews)'|dOMEGA + dBETA*dGAMMA + dALPHA*sqr(vNews)';
//	mYt = dOMEGA + dALPHA*1*sqr(vNews)'+dBETA;  //garch 
//	mYt = mYt|dOMEGA + dALPHA*(2*sqrt(2)*fabs(vNews')-2)+(dBETA+dALPHA);  //Laplace-GAS
//	dL = 10;
//	mYt = mYt|dOMEGA + dALPHA*(dL+3)/dL*((dL+1)/(dL-2)*(1+1*sqr(vNews)'/((dL-2)*1)).^(-1).*(1*sqr(vNews)')-1)+(dBETA+dALPHA);  //GAS
//	vX = vNews';
//	DrawX(iArea, mYt, vX);
//	ShowDrawWindow();


	SetDrawWindow("NIC_GARCH_AND_Asymmetric_l-GAS");
	iArea = 0;
//	decl mYt = dOMEGA + dBETA*dGAMMA + dALPHA*dGAMMA*sqr(vNews)'|dOMEGA + dBETA*dGAMMA + dALPHA*sqr(vNews)';
//	mYt = dOMEGA + dALPHA*1*sqr(vNews)'+dBETA;  //garch 
//	mYt = mYt|dOMEGA + dALPHA*(2*sqrt(2)*fabs(vNews')-2)+(dBETA+dALPHA);  //Laplace-GAS

	//instable
//	dL = 10;
//	mYt = mYt|dOMEGA + dALPHA*(dL+3)/dL*((dL+1)/(dL-2)*(1+1*sqr(vNews)'/((dL-2)*1)).^(-1).*(1*sqr(vNews)')-1)+(dBETA+dALPHA);  //GAS

	//stable
//	dL = 2;
//	mYt = mYt|dOMEGA + dALPHA*(dL+3)/dL*(      (dL+1)*(vNews').^2'./(  (dL-2)*1+(vNews').^2'  )        -1)+(dBETA+dALPHA);  //GAS
//
//         thetaStart        theta    std.error
//alpha       0.10000     0.050859     0.013610
//beta        0.99000      0.99457    0.0067606
//omega      0.010000    0.0057160    0.0040196
//gamma       0.10000      0.33675      0.24357
//p           0.50000      0.46657    0.0090673

	dALPHA 	= 0.050859;
	dBETA 	= 0.99457;
	dOMEGA	= 0.0057160;
	decl p	= 0.46657;

//	mYt = mYt|dOMEGA + dALPHA*(2*sqrt(1-2*p+2*p^2)/ ((1-p)*p) *  (p-findicator( sqrt(1-2*p+2*p^2)*vNews'/ ((1-p)*p*sqrt(1)))).*vNews'*sqrt(1) -2)+(dBETA+dALPHA); //Asymetric Laplace-GAS
//	p= 0.5;
decl	mYt = dOMEGA + dALPHA*(2*sqrt(1-2*p+2*p^2)/ ((1-p)*p) *  (p-findicator( sqrt(1-2*p+2*p^2)*vNews'/ ((1-p)*p*sqrt(1)))).*vNews'*sqrt(1) -2)+(dBETA); //Asymetric Laplace-GAS
//	p= 0.7;
//	mYt = mYt|dOMEGA + dALPHA*(2*sqrt(1-2*p+2*p^2)/ ((1-p)*p) *  (p-findicator( sqrt(1-2*p+2*p^2)*vNews'/ ((1-p)*p*sqrt(1)))).*vNews'*sqrt(1) -2)+(dBETA+dALPHA); //Asymetric Laplace-GAS
//	dL = 10;
//	mYt = mYt|dOMEGA + dALPHA*(dL+3)/dL*((dL+1)/(dL-2)*(1+1*sqr(vNews)'/((dL-2)*1)).^(-1).*(1*sqr(vNews)')-1)+(dBETA+dALPHA);  //GAS
//	vX = vNews';
//	DrawX(iArea, mYt, vX);
//	ShowDrawWindow();

//	SetDrawWindow("NIC_GARCH_AND_Asymmetric_Student's_t-GAS");
//	iArea = 0;
//	decl mYt = dOMEGA + dBETA*dGAMMA + dALPHA*dGAMMA*sqr(vNews)'|dOMEGA + dBETA*dGAMMA + dALPHA*sqr(vNews)';
//	mYt = dOMEGA + dALPHA*1*sqr(vNews)'+dBETA;  //garch 
//	mYt = mYt|dOMEGA + dALPHA*(2*sqrt(2)*fabs(vNews')-2)+(dBETA+dALPHA);  //Laplace-GAS

//                thetaStart        theta    std.error
//alpha              0.10000     0.045741     0.011592
//beta               0.99000      0.99553    0.0055227
//omega            0.0010000    0.0044806    0.0030331
//gamma              0.10000      0.30720      0.19150
//lambda              10.000       5.2598      0.73202
//p                  0.50000      0.46554    0.0086044
	//decl p;
											 
	dALPHA 	= 0.045741;
	dBETA 	= 0.99553;
	dOMEGA	= 0.0044806;
	p= 0.46554;
	dL = 5.2598;
//	print((dL+1)*(4*dL*sqr(p-findicator(-vNews')) ).^(-1));
//	mYt = mYt|dOMEGA + dALPHA*(dL+3)/dL*((dL+1)/(dL-2)*(1+1*sqr(vNews)'/((dL-2)*1)).^(-1).*(1*sqr(vNews)')-1)+(dBETA+dALPHA);  //AST-GAS
	mYt = mYt|dOMEGA + dALPHA*(dL+3)/dL*(     (dL+1)./(4*dL*sqr(p-findicator(-vNews')) ) *            (( 4*(   (p^3+(1-p)^3)*(dL/(dL-2)))  -  16*fK(dL)^2*((1-2*p)*(dL/(dL-1)))^2 	) )          .* (1+1*sqr(vNews')*(( 4*(   (p^3+(1-p)^3)*(dL/(dL-2)))  -  16*fK(dL)^2*((1-2*p)*(dL/(dL-1)))^2 	) ) ./((4*dL*sqr(p-findicator(-vNews')) )*1)).^(-1).*(1*sqr(vNews'))-1)+(dBETA);  //AST-GAS
//	print();
//	mYt = mYt|dOMEGA + dALPHA*(dL+3)/dL         *(               (dL+1)  *  ( 4*(   (p^3+(1-p)^3)*(dL/(dL-2)))  -  16*fK(dL)^2*((1-2*p)*(dL/(dL-1)))^2 	)         /(sqr(p-findicator(-vNews)) *4*dL)   ).*(1+(    4*((p^3+(1-p)^3)*(dL/(dL-2)))-16*fK(dL)^2*((1-2*p)*(dL/(dL-1)))^2 	)*sqr(vNews')./(sqr(p-findicator(-vNews)) *4*dL*1)).^(-1).*sqr(vNews')-1) +(dBETA+dALPHA);
//	print(rows(mYt));

	decl vX = vNews';
	DrawX(iArea, mYt, vX);
	ShowDrawWindow();

}
