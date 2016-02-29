/***********************************************************************/
/*                                                                     */
/*              Finding MT2W                                           */
/*              Reference:  arXiv:1203.4813 [hep-ph]                   */
/*              Authors: Yang Bai, Hsin-Chia Cheng,                    */
/*                       Jason Gallicchio, Jiayin Gu                   */
/*              Based on MT2 by: Hsin-Chia Cheng, Zhenyu Han           */ 
/*              May 8, 2012, v1.00a                                    */
/*                                                                     */  
/***********************************************************************/

/*******************************************************************************
  Usage: 

  1. Define an object of type "mt2w":
     
     mt2w_bisect::mt2w mt2w_event;
 
  2. Set momenta:
 
     mt2w_event.set_momenta(pl,pb1,pb2,pmiss);
 
     where array pl[0..3], pb1[0..3], pb2[0..3] contains (E,px,py,pz), pmiss[0..2] contains (0,px,py) 
     for the visible particles and the missing momentum. pmiss[0] is not used. 
     All quantities are given in double.    
     
     (Or use non-pointer method to do the same.)

  3. Use mt2w::get_mt2w() to obtain the value of mt2w:

     double mt2w_value = mt2w_event.get_mt2w();       
          
*******************************************************************************/ 
              
#include "babymaker/bmaker/interface/mt2w_bisect.hh"
#include <cmath>
#include <iostream>

using namespace std;

namespace mt2w_bisect
{

  mt2w::mt2w(double upper_bound, double error_value, double scan_step)
  {
    solved_ = false;
    momenta_set_ = false;
    mt2w_b_ = 0.;  // The result field.  Start it off at zero.
    upper_bound_ = upper_bound;  // the upper bound of search for MT2W, default value is 500 GeV 
    error_value_ = error_value;  // if we couldn't find any compatible region below the upper_bound, output mt2w = error_value;
    scan_step_ = scan_step;    // if we need to scan to find the compatible region, this is the step of the scan
  }

  double mt2w::get_mt2w()
  {
    if (!momenta_set_)
      {
        cout <<" Please set momenta first!" << endl;
        return error_value_;
      }
        
    if (!solved_) mt2w_bisect();
    return mt2w_b_;
  }

  void mt2w::set_momenta(double *pl, double *pb1, double *pb2, double* pmiss)
  {
    // Pass in pointers to 4-vectors {E, px, py, px} of doubles.  
    // and pmiss must have [1] and [2] components for x and y.  The [0] component is ignored.
    set_momenta(pl[0],  pl[1],  pl[2],  pl[3],
                pb1[0], pb1[1], pb1[2], pb1[3],
                pb2[0], pb2[1], pb2[2], pb2[3],
                pmiss[1], pmiss[2]);
  }

  void mt2w::set_momenta(double El,  double plx,  double ply,  double plz,
                         double Eb1, double pb1x, double pb1y, double pb1z,
                         double Eb2, double pb2x, double pb2y, double pb2z,
                         double pmissx, double pmissy)
  {
    solved_ = false;     //reset solved tag when momenta are changed.
    momenta_set_ = true;
        
    double msqtemp;   //used for saving the mass squared temporarily

    //l is the visible lepton
        
    El_  = El;
    plx_ = plx;
    ply_ = ply;
    plz_ = plz;
        
    Elsq_ = El*El;
        
    msqtemp = El*El-plx*plx-ply*ply-plz*plz;
    if (msqtemp > 0.0) {mlsq_ = msqtemp;}
    else {mlsq_ = 0.0;}                           //mass squared can not be negative
    ml_ = sqrt(mlsq_);                             // all the input masses are calculated from sqrt(p^2)
        
    //b1 is the bottom on the same side as the visible lepton
        
    Eb1_  = Eb1;
    pb1x_ = pb1x;
    pb1y_ = pb1y;
    pb1z_ = pb1z;
        
    Eb1sq_ = Eb1*Eb1;
        
    msqtemp = Eb1*Eb1-pb1x*pb1x-pb1y*pb1y-pb1z*pb1z;
    if (msqtemp > 0.0) {mb1sq_ = msqtemp;}
    else {mb1sq_ = 0.0;}                          //mass squared can not be negative
    mb1_ = sqrt(mb1sq_);                           // all the input masses are calculated from sqrt(p^2)
        
    //b2 is the other bottom (paired with the invisible W)
        
    Eb2_  = Eb2;
    pb2x_ = pb2x;
    pb2y_ = pb2y;
    pb2z_ = pb2z;
        
    Eb2sq_ = Eb2*Eb2;
        
    msqtemp = Eb2*Eb2-pb2x*pb2x-pb2y*pb2y-pb2z*pb2z;
    if (msqtemp > 0.0) {mb2sq_ = msqtemp;}
    else {mb2sq_ = 0.0;}                          //mass squared can not be negative
    mb2_ = sqrt(mb2sq_);                           // all the input masses are calculated from sqrt(p^2)

        
    //missing pt        
        
        
    pmissx_ = pmissx; 
    pmissy_ = pmissy;
        
    //set the values of masses
        
    mv_ = 0.0;   //mass of neutrino
    mw_ = 80.385;  //mass of W-boson

    //precision?        

    if (ABSOLUTE_PRECISION > 100.*RELATIVE_PRECISION) precision_ = ABSOLUTE_PRECISION;
    else precision_ = 100.*RELATIVE_PRECISION;
  }

  void mt2w::mt2w_bisect()
  {
  
   
    solved_ = true;
    cout.precision(11);

    // In normal running, mtop_high WILL be compatible, and mtop_low will NOT.
    double mtop_high = upper_bound_; //set the upper bound of the search region
    double mtop_low;                //the lower bound of the search region is best chosen as m_W + m_b

    if (mb1_ >= mb2_) {mtop_low = mw_ + mb1_;}
    else {mtop_low = mw_ + mb2_;}
        
    // The following if and while deal with the case where there might be a compatable region
    // between mtop_low and 500 GeV, but it doesn't extend all the way up to 500.
    // 
        
    // If our starting high guess is not compatible, start the high guess from the low guess...
    if (teco(mtop_high)==0) {mtop_high = mtop_low;}
        
    // .. and scan up until a compatible high bound is found.
    //We can also raise the lower bound since we scaned over a region that is not compatible
    while (teco(mtop_high)==0 && mtop_high < upper_bound_ + 2.*scan_step_) {

      mtop_low=mtop_high;
      mtop_high = mtop_high + scan_step_;
    }
        
    // if we can not find a compatible region under the upper bound, output the error value
    if (mtop_high > upper_bound_) {
      mt2w_b_ = error_value_;
      return;
    }
        
    // Once we have an compatible mtop_high, we can find mt2w using bisection method
    while(mtop_high - mtop_low > precision_)
      {
        double mtop_mid,teco_mid;
        //bisect
        mtop_mid = (mtop_high+mtop_low)/2.;
        teco_mid = teco(mtop_mid);
      
        if(teco_mid == 0) {mtop_low  = mtop_mid;}
        else {mtop_high  = mtop_mid;}
           
      }
    mt2w_b_ = mtop_high;   //output the value of mt2w
    return;
  }

  // for a given event, teco ( mtop ) gives 1 if trial top mass mtop is compatible, 0 if mtop is not.
        
  int mt2w::teco(  double mtop)
  {
        
    //first test if mtop is larger than mb+mw   
        
    if (mtop < mb1_+mw_ || mtop < mb2_+mw_) {return 0;}

    //define delta for convenience, note the definition is different from the one in mathematica code by 2*E^2_{b2}
                
    double ETb2sq = Eb2sq_ - pb2z_*pb2z_;  //transverse energy of b2
    double delta = (mtop*mtop-mw_*mw_-mb2sq_)/(2.*ETb2sq);
        
        
    //del1 and del2 are \Delta'_1 and \Delta'_2 in the notes eq. 10,11
        
    double del1 = mw_*mw_ - mv_*mv_ - mlsq_;
    double del2 = mtop*mtop - mw_*mw_ - mb1sq_ - 2*(El_*Eb1_-plx_*pb1x_-ply_*pb1y_-plz_*pb1z_);
        
    // aa bb cc are A B C in the notes eq.15
        
    double aa = (El_*pb1x_-Eb1_*plx_)/(Eb1_*plz_-El_*pb1z_);
    double bb = (El_*pb1y_-Eb1_*ply_)/(Eb1_*plz_-El_*pb1z_);
    double cc = (El_*del2-Eb1_*del1)/(2.*Eb1_*plz_-2.*El_*pb1z_);
        
  
    //calculate coefficients for the two quadratic equations (ellipses), which are
    //
    //  a1 x^2 + 2 b1 x y + c1 y^2 + 2 d1 x + 2 e1 y + f1 = 0 ,  from the 2 steps decay chain (with visible lepton)
    //
    //  a2 x^2 + 2 b2 x y + c2 y^2 + 2 d2 x + 2 e2 y + f2 <= 0 , from the 1 stop decay chain (with W missing)
    //
    //  where x and y are px and py of the neutrino on the visible lepton chain

    a1_ = Eb1sq_*(1.+aa*aa)-(pb1x_+pb1z_*aa)*(pb1x_+pb1z_*aa);
    b1_ = Eb1sq_*aa*bb - (pb1x_+pb1z_*aa)*(pb1y_+pb1z_*bb);
    c1_ = Eb1sq_*(1.+bb*bb)-(pb1y_+pb1z_*bb)*(pb1y_+pb1z_*bb);
    d1_ = Eb1sq_*aa*cc - (pb1x_+pb1z_*aa)*(pb1z_*cc+del2/2.0);
    e1_ = Eb1sq_*bb*cc - (pb1y_+pb1z_*bb)*(pb1z_*cc+del2/2.0);
    f1_ = Eb1sq_*(mv_*mv_+cc*cc) - (pb1z_*cc+del2/2.0)*(pb1z_*cc+del2/2.0);
        
    //  First check if ellipse 1 is real (don't need to do this for ellipse 2, ellipse 2 is always real for mtop > mw+mb)
        
    double det1 = (a1_*(c1_*f1_ - e1_*e1_) - b1_*(b1_*f1_ - d1_*e1_) + d1_*(b1_*e1_-c1_*d1_))/(a1_+c1_);
        
    if (det1 > 0.0) {return 0;}
        
    //coefficients of the ellptical region
        
    a2_ = 1-pb2x_*pb2x_/(ETb2sq);
    b2_ = -pb2x_*pb2y_/(ETb2sq);
    c2_ = 1-pb2y_*pb2y_/(ETb2sq);
        
    // d2o e2o f2o are coefficients in the p2x p2y plane (p2 is the momentum of the missing W-boson)
    // it is convenient to calculate them first and transfer the ellipse to the p1x p1y plane
    d2o_ = -delta*pb2x_;
    e2o_ = -delta*pb2y_;
    f2o_ = mw_*mw_ - delta*delta*ETb2sq;
        
    d2_ = -d2o_ -a2_*pmissx_ -b2_*pmissy_;
    e2_ = -e2o_ -c2_*pmissy_ -b2_*pmissx_;
    f2_ = a2_*pmissx_*pmissx_ + 2*b2_*pmissx_*pmissy_ + c2_*pmissy_*pmissy_ + 2*d2o_*pmissx_ + 2*e2o_*pmissy_ + f2o_;
        
    //find a point in ellipse 1 and see if it's within the ellipse 2, define h0 for convenience
    double x0, h0, y0, r0;
    x0 = (c1_*d1_-b1_*e1_)/(b1_*b1_-a1_*c1_);
    h0 = (b1_*x0 + e1_)*(b1_*x0 + e1_) - c1_*(a1_*x0*x0 + 2*d1_*x0 + f1_);
    if (h0 < 0.0) {return 0;}  // if h0 < 0, y0 is not real and ellipse 1 is not real, this is a redundant check.
    y0 = (-b1_*x0 -e1_ + sqrt(h0))/c1_;
    r0 = a2_*x0*x0 + 2*b2_*x0*y0 + c2_*y0*y0 + 2*d2_*x0 + 2*e2_*y0 + f2_;
    if (r0 < 0.0) {return 1;}  // if the point is within the 2nd ellipse, mtop is compatible
        
        
    //obtain the coefficients for the 4th order equation 
    //devided by Eb1^n to make the variable dimensionless
    long double A4, A3, A2, A1, A0;

    A4 = 
      -4*a2_*b1_*b2_*c1_ + 4*a1_*b2_*b2_*c1_ +a2_*a2_*c1_*c1_ + 
      4*a2_*b1_*b1_*c2_ - 4*a1_*b1_*b2_*c2_ - 2*a1_*a2_*c1_*c2_ + 
      a1_*a1_*c2_*c2_;  
        
    A3 =
      (-4*a2_*b2_*c1_*d1_ + 8*a2_*b1_*c2_*d1_ - 4*a1_*b2_*c2_*d1_ - 4*a2_*b1_*c1_*d2_ + 
       8*a1_*b2_*c1_*d2_ - 4*a1_*b1_*c2_*d2_ - 8*a2_*b1_*b2_*e1_ + 8*a1_*b2_*b2_*e1_ + 
       4*a2_*a2_*c1_*e1_ - 4*a1_*a2_*c2_*e1_ + 8*a2_*b1_*b1_*e2_ - 8*a1_*b1_*b2_*e2_ - 
       4*a1_*a2_*c1_*e2_ + 4*a1_*a1_*c2_*e2_)/Eb1_;
        
        
    A2 =
      (4*a2_*c2_*d1_*d1_ - 4*a2_*c1_*d1_*d2_ - 4*a1_*c2_*d1_*d2_ + 4*a1_*c1_*d2_*d2_ - 
       8*a2_*b2_*d1_*e1_ - 8*a2_*b1_*d2_*e1_ + 16*a1_*b2_*d2_*e1_ + 
       4*a2_*a2_*e1_*e1_ + 16*a2_*b1_*d1_*e2_ - 8*a1_*b2_*d1_*e2_ - 
       8*a1_*b1_*d2_*e2_ - 8*a1_*a2_*e1_*e2_ + 4*a1_*a1_*e2_*e2_ - 4*a2_*b1_*b2_*f1_ + 
       4*a1_*b2_*b2_*f1_ + 2*a2_*a2_*c1_*f1_ - 2*a1_*a2_*c2_*f1_ + 
       4*a2_*b1_*b1_*f2_ - 4*a1_*b1_*b2_*f2_ - 2*a1_*a2_*c1_*f2_ + 2*a1_*a1_*c2_*f2_)/Eb1sq_;
        
    A1 =
      (-8*a2_*d1_*d2_*e1_ + 8*a1_*d2_*d2_*e1_ + 8*a2_*d1_*d1_*e2_ - 8*a1_*d1_*d2_*e2_ - 
       4*a2_*b2_*d1_*f1_ - 4*a2_*b1_*d2_*f1_ + 8*a1_*b2_*d2_*f1_ + 4*a2_*a2_*e1_*f1_ - 
       4*a1_*a2_*e2_*f1_ + 8*a2_*b1_*d1_*f2_ - 4*a1_*b2_*d1_*f2_ - 4*a1_*b1_*d2_*f2_ - 
       4*a1_*a2_*e1_*f2_ + 4*a1_*a1_*e2_*f2_)/(Eb1sq_*Eb1_);
        
    A0 =
      (-4*a2_*d1_*d2_*f1_ + 4*a1_*d2_*d2_*f1_ + a2_*a2_*f1_*f1_ + 
       4*a2_*d1_*d1_*f2_ - 4*a1_*d1_*d2_*f2_ - 2*a1_*a2_*f1_*f2_ + 
       a1_*a1_*f2_*f2_)/(Eb1sq_*Eb1sq_);
        
    /*long  double A0sq, A1sq, A2sq, A3sq, A4sq;
      A0sq = A0*A0;
      A1sq = A1*A1;
      A2sq = A2*A2;
      A3sq = A3*A3;
      A4sq = A4*A4;*/
    long double A3sq = A3*A3;
   
    long double B3, B2, B1, B0;
    B3 = 4*A4;
    B2 = 3*A3;
    B1 = 2*A2;
    B0 = A1;
   
    long double C2, C1, C0;
    C2 = -(A2/2 - 3*A3sq/(16*A4));
    C1 = -(3*A1/4. -A2*A3/(8*A4));
    C0 = -A0 + A1*A3/(16*A4);
   
    long double D1, D0;
    D1 = -B1 - (B3*C1*C1/C2 - B3*C0 -B2*C1)/C2;
    D0 = -B0 - B3 *C0 *C1/(C2*C2)+ B2*C0/C2;
   
    long double E0;
    E0 = -C0 - C2*D0*D0/(D1*D1) + C1*D0/D1;
   
    long  double t1,t2,t3,t4,t5;
    //find the coefficients for the leading term in the Sturm sequence  
    t1 = A4;
    t2 = A4;
    t3 = C2;
    t4 = D1;
    t5 = E0;
 

    //The number of solutions depends on diffence of number of sign changes for x->Inf and x->-Inf
    int nsol;
    nsol = signchange_n(t1,t2,t3,t4,t5) - signchange_p(t1,t2,t3,t4,t5);

    //Cannot have negative number of solutions, must be roundoff effect
    if (nsol < 0) nsol = 0;
        
    int out;
    if (nsol == 0) {out = 0;}  //output 0 if there is no solution, 1 if there is solution
    else {out = 1;}

    return out;
  
  }  

  inline int mt2w::signchange_n( long double t1, long double t2, long double t3, long double t4, long double t5)
  {
    int nsc;
    nsc=0;
    if(t1*t2>0) nsc++;
    if(t2*t3>0) nsc++;
    if(t3*t4>0) nsc++;
    if(t4*t5>0) nsc++;
    return nsc;
  }
  inline int mt2w::signchange_p( long double t1, long double t2, long double t3, long double t4, long double t5)
  {
    int nsc;
    nsc=0;
    if(t1*t2<0) nsc++;
    if(t2*t3<0) nsc++;
    if(t3*t4<0) nsc++;
    if(t4*t5<0) nsc++;
    return nsc;
  }

}//end namespace mt2w_bisect
