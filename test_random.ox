#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <maximize.h>
#import <modelbase>
#import <simula>
#include <oxfloat.h>

main()
{
	decl mTEST;
//
//	mTEST = rann(2.9,1.2);
//	print(mTEST);
//
//	print(1*mTEST[9][0]);
//	print(meanc(rann(100*1000000,1)));

//	print(ceil(0.1));
decl mData_2 = loadmat("ReturnsCloseToClose.csv");
decl mData_1 = loadmat("sbux_cropped.csv");
	decl vReturns_1 = mData_1[:][6];
	vReturns_1 = reversec(vReturns_1); 
	vReturns_1 = 100*(log(vReturns_1)-log(lag0(vReturns_1,1)));
	vReturns_1 = dropr(vReturns_1,0);

	savemat("test_sbux.csv", mData_2~vReturns_1);
	
}
