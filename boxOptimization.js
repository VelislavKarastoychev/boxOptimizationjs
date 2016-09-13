/*
 * The complex M.J.Box algorithm,
 * translated from Fortran 77.
 * Autor: Velislav S. Karastoychev,
 * email:exel_mmm@abv.bg.
 * Date: 5-Aug-2016.
 * Description: Finding the minimum of a 
 * given multivariable function under double constraints
 * Purpose: this program uses the complex method of M.J.Box
 * and is initialy created by Joel A. Richardson and J.L. Kuester
 * in Fortran 77.To use the program, the user must provide the specifical
 * functions jfunc and jconst1,which describe the objective function 
 * and the constraints.An example how to use the program is the follow.
 * '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
 * ' the minimized function                                                                          '
 * ' function jfunc(i,n,m,k,l,x){                                                                    '
 * '      var f = (9 - Math.pow(x[i][0] - 3,2))*Math.pow(x[i][1],3)/(27*Math.sqrt(3));               '
 * '	  return f;                                                                                  '
 * ' }                                                                                               '
 * ' the constranits:                                                                                '
 * ' function jconst1(n,m,k,x,i,l){                                                                  '
 * '    var g = new Array(),h = new Array();                                                         '
 * '    x[i][2] = x[i][0] + Math.sqrt(3)*x[i][1];                                                    '
 * '    g[0] = 0.0;                                                                                  '
 * '    g[1] = 0.0;                                                                                  '
 * '    g[2] = 0;                                                                                    '
 * '    h[0] = 100;                                                                                  '
 * '    h[1] = x[i][0]/Math.sqrt(3);                                                                 '
 * '    h[2] = 6;                                                                                    '
 * '    return {g:g,h:h,chi:x};                                                                      '
 * ' }                                                                                               '
 * ' the parameters:                                                                                 '
 * ' var x = [],alpha = 1.3,beta = 0.001,gamma = 5,delta = 0.0001,                                   '
 * ' itmax = 100,n = 2,m = 3,l = 3,k = 4;                                                            '
 * ' for(i = 0;i < k;i++){                                                                           '
 * '	 x[i] = [];                                                                                  '
 * '	 for(j = 0;j < l;j++){                                                                       '
 * '		 x[i][j] = i === 0 && j === 0?1:                                                     '
 * '	     i === 0 && j === 1?0.5:                                                                 '
 * '		 i === 0 && j === 2?x[0][0] + Math.sqrt(3)*x[0][1]:0;                                '
 * '	 }                                                                                           '
 * ' }                                                                                               '     
 * ' executing:                                                                                      '     
 * ' programming.Box(n,m,l,k,x,alpha,beta,gamma,delta,itmax,jconst1,jfunc);                           '     
 * ' output:=> { x: [ 3.0003939306264016, 1.731787690628987, 5.99993819871818 ],                     '     
 * ' fmin: 0.9995443200842825,                                                                       '    
 * ' iterations: 90 }                                                                                '     
 * '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
 */
var programming = new Object();
(function(exports){
  /*
   * Create an k x n matrix with elements of 
   * random numbers.We use the John Burkardt routine
   * writen in Fortran (modified 23-Jan-2006).The name
   * "rn" of the original routine is renamed with "random".
   */
	function random(m,n,seed){
		var i,j,k,rn = new Array();
                for(i = 0;i < m;i++){
        	    rn[i] = [];
        	    for(j = 0;j < n;j++){
        		seed = parseInt(seed);
        		k = parseInt(seed/127773);
            	        seed = parseInt(16807*(seed - k*127773)- k*2836);
            	        if ( seed < 0 )seed += 2147483647;
        		rn[i][j] = seed*4.656612875e-10;
        	    }
               }
        return rn;
	}
	function jconsx(n,m,l,k,x,alpha,beta,gamma,delta,itmax,jconst1,jfunc){
		var f = new Array(),/*The current function values in points 0,1,2,3...k*/
		k1,/*for loop limit*/
		kode,/*integer --> 1:there is implicit variables,0: no implicit variables*/
		kount,/*counter for the gamma convergence criteria*/
		i,/*counter for the functions jfunc,jconst1 and other*/
		j,ii,jj,/*counters of the for loop*/
		iev1,iev2,/*minimum and maximum function value indices*/
		icm,/*counter of the min-max for loop*/
		g,/*lower constraint array*/
		h,/*upper constraint array*/
		xc,/*centroid array*/
		it,/*current iteration*/
		r = random(k,n,12345);/*random numbers array*/
		/*
		 * define the iterations counter and check if the problem
		 * has inplicit variables and 
		 * zero the (k - 1) rows of the x matrix if is not zero 
		 */
		 it = 1;
		 kode = m - n <= 0?0:1; 
		 for(ii = 1;ii < k;ii++){
		 	for(j = 0;j < n;j++)x[ii][j] = 0.0;
		 }
		//console.log('Initial x matrix:' + x);
		/* Calculate the complex points and check
		 * against constraints:
		 */
		 for(ii = 1;ii < k;ii++){
		 	for(j = 0;j < n;j++){
		 		i = ii;
		 		g = jconst1(n,m,k,x,i,l).g;
		 		h = jconst1(n,m,k,x,i,l).h;
		 		x = jconst1(n,m,k,x,i,l).chi;
		 		x[ii][j] = g[j] + r[ii][j]*(h[j] - g[j]);
		 	}
		 	k1 = ii + 1;
		 	x[ii][j] = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).chi[j];
		 	g = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).g;
		 	h = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).h;
		 	//xc = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).centroid;
		 }
		 k1 = k;
		 /*
		  * Calculate the f points and define the convergence
		  * factor kount.
		  */
		  for(i = 0;i < k;i++)f[i] = jfunc(i,n,m,k,l,x);
		  kount = 1;
		  /*
		   * While the current iteration it is smaller than 
		   * the maximum iterations number itmax do the follow:
		   */
		   do{
		   	  	/* find the points with lowest and hightest function value */
		   	  	iev1 = 0,iev2 = 0;
		 		for(icm = 1;icm < k;icm++){
		 			if(f[iev1] > f[icm])iev1 = icm;
		 			if(f[iev2] < f[icm])iev2 = icm;
		 		}
		 		/*check convergence criteria*/
		 		if(f[iev2] - f[iev1] >= beta){
		 			kount = 1;
		 		}
		 		else{
		 			kount += 1;
		 			if(kount >= gamma)break;
		 		}
		 		/* replace point with lowest function value */
		 		xc = jcent(n,m,k,iev1,x,k1);
		 		for(j = 0;j < n;j++){
		 			x[iev1][j] *= -alpha;
		 			x[iev1][j] += (1 + alpha)*xc[j];
		 		}
		 		i = iev1;
		 		x[i] = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).chi;
		 		g = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).g;
		 		h = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).h;
		 		xc = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).centroid
		 		f[i] = jfunc(i,n,m,k,l,x);
		 		/*
		 		 * While a point repeats with lowest function value, 
		 		 * replace the hightest function value if 
		 		 * repeats as lowest function value:
		 		 */
		 		 do{
		 		 	iev2 = 0;
		 		 	for(icm = 1;icm < k;icm++){
		 		 		if(f[iev2] > f[icm])iev2 = icm;
		 		 	}
		 		 	if(iev2 != iev1)break;
		 		 	for(jj = 0;jj < n;jj++){
		 		 		x[iev1][jj] += xc[jj];
		 		 		x[iev1][jj] /= 2;
		 		 	}
		 		 	i = iev1;
		 		 	x[i] = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).chi;
		 			g = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).g;
		 			h = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).h;
		 			xc = jcek1(i,n,m,k,x,kode,delta,l,k1,xc).centroid;
		 			f[i] = jfunc(i,n,m,k,l,x);
		 		}while(true);
		 		it += 1; 	
		   }while(it <= itmax);
		   return {x:x[iev1],fmin:f[iev1],iterations:it};
	}
	function jcek1(i,n,m,k,x,kode,delta,l,k1,xc){
		var g,h,xc,kt,jc,nn;
		do{
			kt = 0;
			jc = jconst1(n,m,k,x,i,l);
			g = jc.g;
		 	h = jc.h;
		 	x = jc.chi;
		 	/*
		 	 * check against the explict variablws
		 	 */
		 	for(j = 0;j < n;j++){
		 		if(x[i][j] <= g[j])x[i][j] = g[j] + delta;
		 		else{
		 			if(h[j] <= x[i][j])x[i][j] = h[j] - delta;
		 		}
		 	}
			if(kode === 0)break;
			/*
			 * check against the implicit constraints
			 */
			 nn = n;
			 for(j = nn;j < m;j++){
			 	jc = jconst1(n,m,k,x,i,l);
			 	g = jc.g;
		 	 	h = jc.h;
		 	 	x = jc.chi;
		 	 	if((x[i][j] < g[j]) || ((x[i][j] >= g[j]) && (h[j] < x[i][j]))){
		 	 		iev1 = i;
		 	 		kt = 1;
		 	 		xc = jcent(n,m,k,iev1,x,k1);
		 	 		for(jj = 0;jj < n;jj++){
		 	 			x[i][jj] +=xc[jj];
		 	 			x[i][jj] /= 2;
		 	 		}
		 	 	}
			 }
		}while(kt > 0);
		return{chi:x[i],g:g,h:h,centroid:xc};
	}
	function jcent(n,m,k,iev1,x,k1){
		var j,il,rk,xc = new Array();
		for(j = 0;j < n;j++){
			xc[j] = 0.0;
			for(il = 0;il < k1;il++)xc[j] += x[il][j];
			rk = k1 - 1;
			xc[j] -= x[iev1][j];
			xc[j] /= rk;
		}
		return xc;
	}
	exports.Box = jconsx;
	
})(programming); 
