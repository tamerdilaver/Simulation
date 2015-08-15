#include <oxstd.h>

findicator(const vX){
	decl vReturn, iN;
	iN = sizerc(vX);

	if(rows(vX)!=1){
		vReturn = zeros(iN,1);
	}else{
		vReturn = zeros(1,iN);
	}
	
	print(vReturn);
	
	for(decl i=0;i<iN;i++){
		if(vX[i]<=0){
			vReturn[i]=1;
		}
	}
	return vReturn; 
}
main()
{
	print(<-3,2>);
	print(findicator(<-3,2>));
}
