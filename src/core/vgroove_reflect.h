#include "reflection.h"

namespace pbrt{


//optimations:  sinTheta, cosTheta, tanPhi can all be available from wo and wh
// sinTheta = wh.z, cosTheta = wh.x; tanPhi = wop.x/wop.z?
//cosPhi = wop.z

struct Hit {
    Hit(float thetaR, int bounce, char side, float GFactor): 
        thetaR(thetaR), GFactor(GFactor), side(side), bounce(bounce) {};

    bool isHit(float thetaI, int bcount, char s) {
        return fabs(thetaR-thetaI) < 1e-6 && bounce == bcount && side == s;
    }
    
    float thetaR, GFactor;
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
    inline std::vector<Hit>& eval(float thetaG, float thetaO); 
    inline std::vector<Hit>& leftEvalOnly(float thetaG, float thetaO); 
    inline std::vector<Hit>& rightEvalOnly(float thetaG, float thetaO); 

    inline void leftEval(float theta, float phi, float maxX, bool hasNear, float hitrange);
    inline float rightEval(float theta, float phi, float hitrange);
    inline void addHit(float xi, int bcount, char side, float GFactor);
    inline float inverseEval(float thetaO, float thetaI, int bounceCount, char side, float &thetaM);
    std::vector<Hit> theHits;

  private:
    inline static bool computeThetaM(float thetaO, float thetaI, int bounceCount, char side, float& thetaM);
    inline static bool computeRightThetaM(float thetaO, float thetaI, int bounceCount, float& theta);
    inline static bool computeLeftThetaM(float thetaO, float thetaI, int bounceCount, float& theta);
};

bool 
VGroove::computeRightThetaM(float thetaO, float thetaI, int n, float &theta) {
    float xchi = -thetaI;
    if ((n+1)%2 == 1) xchi *= -1;
    theta = (Pi + thetaO - xchi) *.5/n;
    if (theta > thetaO && theta < .5 * Pi) return true;
    return false; 
}

bool
VGroove::computeLeftThetaM(float thetaO, float thetaI, int n, float& theta) {
    float xchi = -thetaI;
    if (n%2 == 1) xchi *= -1;
    theta = (Pi - thetaO - xchi) *.5/n;
    if (theta > 1e-6 && theta < .5 * Pi) return true;
    return false;
} 

bool
VGroove::computeThetaM(float thetaO, float thetaI, int bounceCount, char side, float& thetaM) {
    if (side == 'l')  {
        return computeLeftThetaM(thetaO, thetaI, bounceCount, thetaM);
    } else {
        return computeRightThetaM(thetaO, thetaI, bounceCount, thetaM);
    }
}

float
VGroove::inverseEval(float thetaO, float thetaI, int bounceCount, char side, float &thetaM) {

    float GFactor = 0;

    assert(thetaO > 0);
    if (thetaO + .0001 > .5 * Pi) return GFactor; 
    
    theHits.clear();
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
    return GFactor; 
}

void 
VGroove::addHit(float xi, int bcount, char side, float GFactor) {
    
    if (GFactor > 1e-5) theHits.push_back( Hit(-xi, bcount, side, GFactor));

    //scale minInterval and maxInterval to [-1, 1] to match pat's ray tracer
    //float minInterval = interval[0]/range * 2.0;
    //float maxInterval = interval[1]/range * 2.0;
    //zipin ratio is still within [0, 1] range
    //float r = (interval[1] - interval[0])/range;
} 

void 
VGroove::leftEval(float theta, float phi, float maxX, bool hasNear, float hitrange) {

    float gAngle = 2.0 * theta;
    float psi_min = Pi - 2.0 *(theta + phi);

    float sinThetaPhi = sin(theta+phi);

    float psi_max = Pi - (theta+phi);

    float zMax = 1;

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

    float xi_min = Pi - phi - n_min * gAngle;
    float xi_max = Pi - phi - n_max * gAngle;
    
    float x_min_intersect = -hitrange * .5;
    float x_max_intersect = maxX;
    float xrange = maxX - x_min_intersect;

    if (n_min%2) xi_min *= -1;
    if (n_max%2) xi_max *= -1;

    if (n_max > n_min) {
        //compute z critical intersect, length from circle center
        int k = n_min;
        float criticalAngle = Pi - (theta+phi) - k* gAngle;
        float z = 0;
        if (criticalAngle > 0) {
            z = sin(criticalAngle)/sinThetaPhi;
        }
        addHit(xi_max, n_max, 'l', zMax - (1-z));
        addHit(xi_min, n_min, 'l', 1-z);
    
        //float x_critical_intersect = x_min_intersect + (1-z) / zMax  * xrange
        //addHit(hits, xi_max, n_max, 'l', zMax - (1-z), (x_critical_intersect, x_max_intersect), hitrange);
        //addHit(hits, xi_min, n_min, 'l', 1-z, (x_min_intersect, x_critical_intersect), hitrange);
    } else {
        addHit(xi_min, n_min, 'l', zMax);

        //addHit(hits, xi_min, n_min, 'l', zMax, (x_min_intersect, x_max_intersect), hitrange);
    }
}

float
VGroove::rightEval(float theta, float phi, float hitrange) {

    float gAngle = 2.0 * theta;
   
    //Feng's correction to near psi computation from Paper 
    float psi_max = Pi -(theta-phi);
    float psi_min = Pi - 2.0 * (theta-phi);

    //#Zipin Paper near psi computation 
    //psi_max = mt.pi - 2.0 * phi
    //psi_min = mt.pi - theta - phi

    float x_near_min = cos(theta) * tan(phi);

    int n_min = stableCeiling(psi_min, gAngle);
    int n_max = stableCeiling(psi_max, gAngle);

    //if (n_max - n_min) > 1:
    //    print(theta, phi, n_max, n_min) 

    float xi_min = Pi + phi - n_min * gAngle;
    float xi_max = Pi + phi - n_max * gAngle;

    if (n_min%2 == 0) xi_min *= -1;
    if (n_max%2 == 0) xi_max *= -1;
    
    float x_min_intersect = x_near_min;
    float x_max_intersect = sin(theta);

    if (n_min == n_max || theta - phi < 1e-5) {
        //addHit(hits, xi_min, n_min, 'right', 1.0, (x_min_intersect, x_max_intersect), hitrange)
        addHit(xi_min, n_min, 'r', 1.0);
    } else {
        int k = n_min;
        float criticalAngle = Pi - (theta-phi) -  gAngle * k;
        float z = 0;
        
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
VGroove::eval(float thetaR, float phiR) {
    theHits.clear();
    float xmax = sin(thetaR);
    if (phiR > thetaR) {
        leftEval(thetaR, phiR, xmax, false, xmax*2);
    } else {
        float  x_near_min = rightEval(thetaR, phiR, xmax*2); 
        leftEval(thetaR, phiR, x_near_min, true,  xmax*2);
    }
    return theHits;
}

std::vector<Hit> & 
VGroove::leftEvalOnly(float thetaR, float phiR) {
    theHits.clear();
    bool hasNear = (thetaR > phiR);
    float xmax = sin(thetaR);
    float farMax = hasNear? cos(thetaR) * tan(phiR) : xmax;
    leftEval(thetaR, phiR, xmax, hasNear, xmax*2);
    return theHits;
}

std::vector<Hit> & 
VGroove::rightEvalOnly(float thetaR, float phiR) {
    theHits.clear();
    float xmax = sin(thetaR);
    float  x_near_min = rightEval(thetaR, phiR, xmax*2); 
    return theHits;
}

struct EvalFrame;
class VGrooveReflection : public MicrofacetReflection {
  public:
    // MicrofacetReflection Public Methods
    VGrooveReflection(const Spectrum &R,
                      MicrofacetDistribution *distribution, Fresnel *fresnel, 
                      int maxBounce = 3, int minBounce = 1, bool uniSample = true) 
        : MicrofacetReflection(R, distribution, fresnel), maxBounce(maxBounce), 
          minBounce(minBounce), uniSample(uniSample) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    //std::string ToString() const;

  private:

    //uniform sampling for testing
    Spectrum UniSample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;

    float microfacetReflectionWithoutG(const Vector3f& wo, const Vector3f& wi,
                   const Vector3f& wh) const;
    float microfacetPdf(const Vector3f& wo, const Vector3f& wh) const;
    float computeBounceBrdf(const EvalFrame& evalFrame, VGroove& vgroove, int bounce, char side,
                    float& pdf) const;
    Spectrum eval(const Vector3f &wo, const Vector3f &wi, float &pdf) const;

    float computeGFactor(const EvalFrame& evalFrame, VGroove& vgroove, int bounce, char side, Vector3f& wm) const;

    int maxBounce, minBounce;
    bool uniSample;
    //VGroove theGroove;
};

}//end of namespace pbrt
