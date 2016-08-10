# boxOptimizationjs
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
