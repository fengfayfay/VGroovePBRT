#include "reflection.h"

namespace pbrt{


//optimations:  sinTheta, cosTheta, tanPhi can all be available from wo and wh
// sinTheta = wh.z, cosTheta = wh.x; tanPhi = wop.x/wop.z?
//cosPhi = wop.z

struct Hit {
    Hit(Float thetaR, int bounce, char side, Float GFactor): 
        thetaR(thetaR), GFactor(GFactor), side(side), bounce(bounce) {};

    bool isHit(Float thetaI, int bcount, char s) {
        return fabs(thetaR-thetaI) < 1e-6 && bounce == bcount && side == s;
    }
    
    Float thetaR, GFactor;
    int bounce;
    char side;
    //xmin, xmax, ratio are only used for debugging 
};


inline Float stableCeiling(Float psi, Float gAngle) {
    
    int k = int(psi/gAngle);
    CHECK(k >= 0);
    if (psi - k * gAngle < 1e-6) return k;
    return k + 1;
}

class VGroove {
  public:
    VGroove(int maxBounce, int minBounce):maxBounce(maxBounce), minBounce(minBounce), sumG(0){}
    void clear() { sumG = 0; theHits.clear(); }
    inline std::vector<Hit>& completeEval(Float thetaG, Float thetaO); 
    inline std::vector<Hit>& leftEvalOnly(Float thetaG, Float thetaO); 
    inline std::vector<Hit>& rightEvalOnly(Float thetaG, Float thetaO); 

    inline Float inverseEval(Float thetaO, Float thetaI, int bounceCount, char side, Float& thetaM);
    std::vector<Hit> theHits;
    Float sumG;
    const int maxBounce, minBounce;

  private:
    inline void leftEval(Float theta, Float phi, Float maxX, bool hasNear, Float hitrange);
    inline Float rightEval(Float theta, Float phi, Float hitrange);
    inline void addHit(Float xi, int bcount, char side, Float GFactor);
    inline static bool computeThetaM(Float thetaO, Float thetaI, int bounceCount, char side, Float& thetaM);
    inline static bool computeRightThetaM(Float thetaO, Float thetaI, int bounceCount, Float& theta);
    inline static bool computeLeftThetaM(Float thetaO, Float thetaI, int bounceCount, Float& theta);
};

bool 
VGroove::computeRightThetaM(Float thetaO, Float thetaI, int n, Float &theta) {
    Float xchi = -thetaI;
    if ((n+1)%2 == 1) xchi *= -1;
    theta = (Pi + thetaO - xchi) *.5/n;
    if (theta > thetaO && theta < .5 * Pi + 1e-6) return true;
    return false; 
}

bool
VGroove::computeLeftThetaM(Float thetaO, Float thetaI, int n, Float& theta) {
    Float xchi = -thetaI;
    if (n%2 == 1) xchi *= -1;
    theta = (Pi - thetaO - xchi) *.5/n;
    if (theta > 1e-6 && theta < .5 * Pi + 1e-6) return true;
    return false;
} 

bool
VGroove::computeThetaM(Float thetaO, Float thetaI, int bounceCount, char side, Float& thetaM) {
    if (side == 'l')  {
        return computeLeftThetaM(thetaO, thetaI, bounceCount, thetaM);
    } else {
        return computeRightThetaM(thetaO, thetaI, bounceCount, thetaM);
    }
}

Float
VGroove::inverseEval(Float thetaO, Float thetaI, int bounceCount, char side, Float &thetaM) { 

    Float GFactor = 0;

    CHECK(thetaO >= 0);
    if (thetaO + .0001 > .5 * Pi) return GFactor; 
    
    clear();
    if (side == 'r') {
        bool validGroove = computeRightThetaM(thetaO, thetaI, (int) bounceCount, thetaM);
        if (validGroove) rightEvalOnly(thetaM, thetaO);
    } else {
        bool validGroove = computeLeftThetaM(thetaO, thetaI, (int) bounceCount, thetaM);
        if (validGroove) leftEvalOnly(thetaM, thetaO);
    }
    for (int i = 0; i < theHits.size(); i++) {
        if (theHits[i].isHit(thetaI, bounceCount, side)) {
            GFactor = theHits[i].GFactor;
            return GFactor;
        }
    }
    return 0; 
}

void 
VGroove::addHit(Float xi, int bcount, char side, Float GFactor) {
   
    if (bcount >= minBounce && bcount <=maxBounce) { 
        if (GFactor > 1e-5 ) {
            GFactor = std::min((Float)1, GFactor);
            theHits.push_back( Hit(-xi, bcount, side, GFactor));
            sumG += GFactor;
            if (sumG > (1 + 1e-5)) {
                std::cout<<"unexected sumG: "<< sumG <<"\n";
                fflush(stdout);
            }
        }
    }

    //scale minInterval and maxInterval to [-1, 1] to match pat's ray tracer
    //Float minInterval = interval[0]/range * 2.0;
    //Float maxInterval = interval[1]/range * 2.0;
    //zipin ratio is still within [0, 1] range
    //Float r = (interval[1] - interval[0])/range;
} 

void 
VGroove::leftEval(Float theta, Float phi, Float maxX, bool hasNear, Float hitrange) {

    Float gAngle = 2.0 * theta;
    Float psi_min = Pi - 2.0 *(theta + phi);

    Float sinThetaPhi = sin(theta+phi);

    Float psi_max = Pi - (theta+phi);

    Float zMax = 1;

    if (hasNear == false && sinThetaPhi > 0) {
        //this is exactly Gi because when all the rays go to the left side
        //the max left side visible is Gi
        //equation 9 in Zipin Paper, hitrange is 2sinTheta

        zMax = hitrange * cos(phi)/ sinThetaPhi;
        psi_max -= asin((1.0-zMax) * sinThetaPhi);
    } 

    //print(mt.degrees(psi_min), mt.degrees(psi_max))

    int n_min = stableCeiling(psi_min, gAngle);
    int n_max = stableCeiling(psi_max, gAngle);
    
    //if n_max - n_min > 1:
    //    print(theta, phi, n_max, n_min) 

    Float xi_min = Pi - phi - n_min * gAngle;
    Float xi_max = Pi - phi - n_max * gAngle;
    
    Float x_min_intersect = -hitrange * .5;
    Float x_max_intersect = maxX;
    Float xrange = maxX - x_min_intersect;

    if (n_min%2) xi_min *= -1;
    if (n_max%2) xi_max *= -1;

    if (n_max > n_min) {
        //compute z critical intersect, length from circle center
        int k = n_min;
        Float criticalAngle = Pi - (theta+phi) - k* gAngle;
        Float z = 0;
        if (criticalAngle > 0) {
            z = sin(criticalAngle)/sinThetaPhi;
        }
        addHit(xi_max, n_max, 'l', zMax - (1-z));
        addHit(xi_min, n_min, 'l', 1-z);
    
        //Float x_critical_intersect = x_min_intersect + (1-z) / zMax  * xrange
        //addHit(hits, xi_max, n_max, 'l', zMax - (1-z), (x_critical_intersect, x_max_intersect), hitrange);
        //addHit(hits, xi_min, n_min, 'l', 1-z, (x_min_intersect, x_critical_intersect), hitrange);
    } else {
        addHit(xi_min, n_min, 'l', zMax);

        //addHit(hits, xi_min, n_min, 'l', zMax, (x_min_intersect, x_max_intersect), hitrange);
    }
}

Float
VGroove::rightEval(Float theta, Float phi, Float hitrange) {

    Float gAngle = 2.0 * theta;
   
    //Feng's correction to near psi computation from Paper 
    Float psi_max = Pi -(theta-phi);
    Float psi_min = Pi - 2.0 * (theta-phi);

    //#Zipin Paper near psi computation 
    //psi_max = mt.pi - 2.0 * phi
    //psi_min = mt.pi - theta - phi

    Float x_near_min = cos(theta) * tan(phi);

    int n_min = stableCeiling(psi_min, gAngle);
    int n_max = stableCeiling(psi_max, gAngle);

    //if (n_max - n_min) > 1:
    //    print(theta, phi, n_max, n_min) 

    Float xi_min = Pi + phi - n_min * gAngle;
    Float xi_max = Pi + phi - n_max * gAngle;

    if (n_min%2 == 0) xi_min *= -1;
    if (n_max%2 == 0) xi_max *= -1;
    
    Float x_min_intersect = x_near_min;
    Float x_max_intersect = sin(theta);

    if (n_min == n_max || theta - phi < 1e-5) {
        //addHit(hits, xi_min, n_min, 'right', 1.0, (x_min_intersect, x_max_intersect), hitrange)
        addHit(xi_min, n_min, 'r', 1.0);
    } else {
        int k = n_min;
        Float criticalAngle = Pi - (theta-phi) -  gAngle * k;
        Float z = 0;
        
        if (criticalAngle > 0) {
            z = sin(criticalAngle) /sin(theta-phi);
        }
        addHit(xi_min, n_min, 'r', 1.0-z);
        addHit(xi_max, n_max, 'r', z);
        //x_critical_intersect = x_min_intersect + (x_max_intersect - x_min_intersect) * z;
        //addHit(hits, xi_min, n_min, 'right', 1.0-z, (x_critical_intersect, x_max_intersect), hitrange)
        //addHit(hits, xi_max, n_max, 'right', z, (x_min_intersect, x_critical_intersect), hitrange)
    }
    return x_near_min;
}

std::vector<Hit> & 
VGroove::completeEval(Float thetaR, Float phiR) {
    clear();
    Float xmax = sin(thetaR);
    if (phiR > thetaR) {
        leftEval(thetaR, phiR, xmax, false, xmax*2);
    } else {
        Float  x_near_min = rightEval(thetaR, phiR, xmax*2); 
        leftEval(thetaR, phiR, x_near_min, true,  xmax*2);
    }
    return theHits;
}

std::vector<Hit> & 
VGroove::leftEvalOnly(Float thetaR, Float phiR) {
    clear();
    bool hasNear = (thetaR > phiR);
    Float xmax = sin(thetaR);
    Float farMax = hasNear? cos(thetaR) * tan(phiR) : xmax;
    leftEval(thetaR, phiR, xmax, hasNear, xmax*2);
    return theHits;
}

std::vector<Hit> & 
VGroove::rightEvalOnly(Float thetaR, Float phiR) {
    clear();
    Float xmax = sin(thetaR);
    Float  x_near_min = rightEval(thetaR, phiR, xmax*2); 
    return theHits;
}

struct EvalFrame;
class VGrooveReflection : public MicrofacetReflection {
  public:
    // MicrofacetReflection Public Methods
    VGrooveReflection(const Spectrum &R,
                      MicrofacetDistribution *distribution, Fresnel *fresnel, 
                      int maxBounce = 3, int minBounce = 1, bool uniSample = true);
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum f(const Vector3f &wo, const Vector3f &wi, Float& pdf) const;
    
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    //std::string ToString() const;
    bool testPDF() const;

  private:

    //uniform sampling for testing
    Spectrum UniSample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;

    Float microfacetReflectionWithoutG(const Vector3f& wo, const Vector3f& wi,
                   const Vector3f& wh) const;
    Float microfacetPdf(const Vector3f& wo, const Vector3f& wh) const;
    Float computeBounceBrdf(const EvalFrame& evalFrame, VGroove& vgroove, int bounce, char side,
                    Float& pdf, Spectrum& F) const;
    Spectrum eval(const EvalFrame& evalFrame, const Vector3f &wo, const Vector3f &wi, Float &pdf) const;

    Float computeGFactor(const EvalFrame& evalFrame, VGroove& vgroove, int bounce, 
                         char side, Vector3f& wm) const;

    Float computePdfIntegral(Float thetaI) const;

    

    const int maxBounce, minBounce;
    bool uniSample;
};


}//end of namespace pbrt
