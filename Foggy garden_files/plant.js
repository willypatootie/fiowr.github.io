var pt, ind, off = 0, pi = 0;
function stem(x,y,z, tx,ty,tz){
   var len = 6;
   var t = t0 = off/len, tm = t + 6;
   for(var i = 0; i < 5; i++ ){
       ind[pi++] = t++;  ind[pi++] = t;  ind[pi++] = tm;
       ind[pi++] = tm++; ind[pi++] = t;  ind[pi++] = tm;
   }
   ind[pi++] = t;  ind[pi++] = t0;  ind[pi++] = tm;
   ind[pi++] = tm; ind[pi++] = t0;  ind[pi++] = t0 + 6;
   var b0 = [.05,0,0, .025,0,.044, -.025,0,.044,  -.05,0,0, -.025,0,-.044, .025,0,-.044];

   var u = 0,  g = .7;
   for(var i = 0; i < 6; i++ ){
     pt[off++] = b0[u++]+tx;  pt[off++] = b0[u++]+ty;  pt[off++] = b0[u++]+tz;
     pt[off++] = 0;  pt[off++] = g;  pt[off++] = 0;}
   u = 0;
   for(var i = 0; i < 6; i++ ){
     pt[off++] = b0[u++] + x+tx;  pt[off++] = b0[u++] + y+ty;  pt[off++] = b0[u++] + z+tz;
     pt[off++] = 0;  pt[off++] = g+.2;  pt[off++] = 0;}
}
function splinePatch(m, n, sm, sn2, Bm, Bn, cp, strips,  saw){
   var len = 6,  sn = 2*sn2 + 1;
   var t = off/len, tm = t + sm;
   if( saw){
     for(var i = 0; i < sm-1; i++ ){
       ind[pi++] = tm++; ind[pi++] = ++t;  ind[pi++] = tm;
     }
     t++;  tm++;
     for(var j = 1; j < sn2; j++ ){
       for(var i = 0; i < sm-1; i++ ){
         ind[pi++] = t++;  ind[pi++] = t;  ind[pi++] = tm;
         ind[pi++] = tm++; ind[pi++] = t;  ind[pi++] = tm;
       }
       t++;  tm++;
     }
     for(var j = 1; j < sn2; j++ ){
       for(var i = 0; i < sm-1; i++ ){
         ind[pi++] = tm++; ind[pi++] = t;  ind[pi++] = tm;
         ind[pi++] = t++;  ind[pi++] = t;  ind[pi++] = tm;
       }
       t++;  tm++;
     }
//     t--;
     tm++;
     for(var i = 0; i < sm-1; i++ ){
       ind[pi++] = t++;  ind[pi++] = t;  ind[pi++] = tm++;
     }}
   else{
     for(var j = 0; j < sn2; j++ ){
       for(var i = 0; i < sm-1; i++ ){
         ind[pi++] = t++;  ind[pi++] = t;  ind[pi++] = tm;
         ind[pi++] = tm++; ind[pi++] = t;  ind[pi++] = tm;
       }
       t++;  tm++;
     }
     for(var j = 0; j < sn2; j++ ){
       for(var i = 0; i < sm-1; i++ ){
         ind[pi++] = tm++; ind[pi++] = t;  ind[pi++] = tm;
         ind[pi++] = t++;  ind[pi++] = t;  ind[pi++] = tm;
       }
       t++;  tm++;
     }
   }
   var pn = 0;
   for(var v = 0; v < sn; v++ ){
    var pm = 0;
    for(var u = 0; u < sm; u++ ){
     var pc = 0;
     for(var l = 0; l < len; l++ ){
       var su = 0;
       for(var j = 0; j < n; j++ ){
         var s = 0;
         for(var i = 0; i < m; i++ ) s += Bm[pm+i]*cp[pc++];
         su += Bn[pn+j]*s;}
       if(l != 4) pt[off++] = su;
       else switch( strips ){
         case 0:  pt[off++] = su;  break;
         case 1:
           if((l == 4) && (((v-u)&1)||(v==sn2)) ) pt[off++] = su - .3;
           else pt[off++] = su;  break;
         case 2: pt[off++] = su - .3*(u & 1); break;
         case 3: pt[off++] = su - .3*(v & 1); break;
       }
     }
     pm += m;}
    pn += n;
   }
}
function BezierFunc(n, sn, Bn){
   var B = new Float32Array(n);
   var dt = 1/(sn - 1), t = 0, t1, i,j,k, p = 0;
   for (k = 0; k < sn; k++){
    B[0] = 1;  t1 = 1 - t;
    for (j = 0; j < n-1; j++)
     for (i = j; i >= 0; i--){
      B[i+1] += t*B[i];
      B[i] = t1*B[i];}
    for (j = 0; j < n; j++){
      Bn[p++] = B[j];
      B[j] = 0;}
    t += dt;
  }
}
function rotX(fi, mn, cp){
   var c = Math.cos(fi),  s = Math.sin(fi),  mn2 = mn + mn;
   for(var i = 0; i < mn; i++ ){
      var y = cp[i + mn];
      cp[i + mn] = y*c + cp[i + mn2]*s;
      cp[i + mn2] = -y*s + cp[i + mn2]*c;}
}
function rotY(fi, mn, cp){
   var c = Math.cos(fi),  s = Math.sin(fi),  mn2 = mn + mn;
   for(var i = 0; i < mn; i++ ){
      var x = cp[i];
      cp[i] = x*c + cp[i + mn2]*s;
      cp[i + mn2] = -x*s + cp[i + mn2]*c;}
}
function rotZ(fi, mn, cp){
   var c = Math.cos(fi),  s = Math.sin(fi);
   for(var i = 0; i < mn; i++ ){
      var x = cp[i];
      cp[i] = x*c + cp[i + mn]*s;
      cp[i + mn] = -x*s + cp[i + mn]*c;}
}
function translate(x, y, z, mn, cp){
   var mn2 = mn + mn;
   for(var i = 0; i < mn; i++ ){
      cp[i] += x;  cp[i + mn] += y;  cp[i + mn2] += z;}
}
function scale(sc, mn, cp){
   var mn2 = mn + mn;
   for(var i = 0; i < mn; i++ ){
      cp[i] *= sc;  cp[i + mn] *= sc;  cp[i + mn2] *= sc;}
}
function rotX_p(fi, off0, off, pt){
   var c = Math.cos(fi),  s = Math.sin(fi);
   for(var i = off0 + 1; i < off; i += 6 ){
      var y = pt[i];
      pt[i] = y*c + pt[i + 1]*s;
      pt[i + 1] = -y*s + pt[i + 1]*c;}
}
function rotY_p(fi, off0, off, pt){
   var c = Math.cos(fi),  s = Math.sin(fi);
   for(var i = off0; i < off; i += 6 ){
      var x = pt[i];
      pt[i] = x*c + pt[i + 2]*s;
      pt[i + 2] = -x*s + pt[i + 2]*c;}
}
function rotZ_p(fi, off0, off, pt){
   var c = Math.cos(fi),  s = Math.sin(fi);
   for(var i = off0; i < off; i += 6 ){
      var x = pt[i];
      pt[i] = x*c + pt[i + 1]*s;
      pt[i + 1] = -x*s + pt[i + 1]*c;}
}
function translate_p(x, y, z, off0, off, pt){
   for(var i = off0; i < off; i += 6 ){
      pt[i] += x;  pt[i + 1] += y;  pt[i + 2] += z;}
}
function scale_p(sc, off0, off, pt){
   for(var i = off0; i < off; i += 6 ){
      pt[i] *= sc;  pt[i + 1] *= sc;  pt[i + 2] *= sc;}
}
function mirror_x(mn, cp){
   var mn2 = mn + mn;
   for(var i = 0; i < mn; i++ ) cp[i] = -cp[i]
}
