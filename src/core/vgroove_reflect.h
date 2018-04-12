#include "reflect.h"

class VGrooveReflection : public MicrofacetReflection {
  public:
    // MicrofacetReflection Public Methods
    VGrooveReflection(const Spectrum &R,
                         MicrofacetDistribution *distribution, Fresnel *fresnel)
        : MicrofacetReflection(r, distribution, fresnel) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:

    Zipin zipin;
};

struct ZipinHit {
    ZipinHit(loat thetaR, char side, char bounce, float G): 
        thetaR(thetaR), G(G), side(side), bounce(bounce) {};
    
    Float thetaR, G;
    char side;
    char bounce;
    Float G;
    //xmin, xmax, ratio are only used for debugging 
};


Float stableCeiling(Float psi, Float gAngle) {
    int k = int(psi/gAngle);
    if psi - k * gAngle < 1e-6 return k
    return k + 1;
}

class VGroove {
  public:
    std::vector<Hit> Eval(Float thetaG, Float thetaO); 
    std::vector<Hit> nearEval(Float thetaG, Float thetaO);
    std::vector<Hit> farEval(Float thetaG, Float thetaO);
};

void VGroove::addHit(std::vector<Hit>& hits, float xi, char bcount, 
    char side, float G) {

    //scale minInterval and maxInterval to [-1, 1] to match pat's ray tracer
    //float minInterval = interval[0]/range * 2.0;
    //float maxInterval = interval[1]/range * 2.0;
    //zipin ratio is still within [0, 1] range

    float r = (interval[1] - interval[0])/range;
    if (G > 1e-5) hits.push_back( ZipinHit(-xi, s, bcount, G));
} 

void VGroove::far(std::vector<Hit>& hits, float theta, float phi, float maxX, bool hasNear, 
        float hitrange):

    float gAngle = 2.0 * theta;
    float psi_min = Pi - 2.0 *(theta + phi);

    float sinThetaPhi = sin(theta+phi);

    float psi_max = Pi - (theta+phi)

    float zMax = 1;

    if (hasNear == False && sinThetaPhi > 0) {
        //this is exactly Gi because when all the rays go to the left side
        //the max left side visible is Gi
        //equation 9 in Zipin Paper, hitrange is 2sinTheta

        zMax = hitrange * cos(phi)/ sinThetaPhi 
        psi_max -= sin((1.0-zMax) * sinThetaPhi)
    } 

    //print(mt.degrees(psi_min), mt.degrees(psi_max))

    int n_min= stableCeiling(psi_min, gAngle)
    int n_max = stableCeiling(psi_max, gAngle)
    
    //if n_max - n_min > 1:
    //    print(theta, phi, n_max, n_min) 

    float xi_min = Pi - phi - n_min * gAngle
    float xi_max = Pi - phi - n_max * gAngle
    
    float x_min_intersect = -hitrange * .5
    float x_max_intersect = maxX;
    float xrange = maxX - x_min_intersect;

    if (n_min%2) xi_min *= -1;
    if (n_max%2) xi_max *= -1;

    if (n_max > n_min) {
        //compute z critical intersect, length from circle center
        int k = n_min;
        criticalAngle = Pi - (theta+phi) - k* gAngle;
        if (criticalAngle <= 0) {
            z = 0;
        } else {
            z = sin(criticalAngle)/sinThetaPhi;
        }
    
        float x_critical_intersect = x_min_intersect + (1-z) / zMax  * xrange

        //addHit(hits, xi_max, n_max, 'l', zMax - (1-z), (x_critical_intersect, x_max_intersect), hitrange);
        //addHit(hits, xi_min, n_min, 'l', 1-z, (x_min_intersect, x_critical_intersect), hitrange);
        addHit(hits, xi_max, n_max, 'l', zMax - (1-z));
        addHit(hits, xi_min, n_min, 'l', 1-z);
    } else {
        //addHit(hits, xi_min, n_min, 'l', zMax, (x_min_intersect, x_max_intersect), hitrange);
        addHit(hits, xi_min, n_min, 'l', zMax);
    }
   
    return hits;
}

void VGroove:: near(std::vector<Hits> hits, float theta, float phi, float hitrange) {

    float gAngle = 2.0 * theta;
   
    //Feng's correction to near psi computation from Paper 
    float psi_max = Pi -(theta-phi);
    float psi_min = Pi - 2.0 * (theta-phi);

    //#Zipin Paper near psi computation 
    //psi_max = mt.pi - 2.0 * phi
    //psi_min = mt.pi - theta - phi

    float x_near_min = os(theta) * tan(phi);

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

    if (n_min == n_max || theta - phi < .000001) {
        addHit(hits, xi_min, n_min, 'right', 1.0, (x_min_intersect, x_max_intersect), hitrange)
    } else {
        k = n_min
        criticalAngle = mt.pi - (theta-phi) - k * gAngle
        if criticalAngle <= 0:
            z = 0
        else:
            z = mt.sin(criticalAngle) /mt.sin(theta-phi)
        x_critical_intersect = x_min_intersect + (x_max_intersect - x_min_intersect) * z
        addHit(hits, xi_min, n_min, 'right', 1.0-z, (x_critical_intersect, x_max_intersect), hitrange)
        addHit(hits, xi_max, n_max, 'right', z, (x_min_intersect, x_critical_intersect), hitrange)
    }

def zipin(theta, phi):

    //thetaR = mt.radians(theta)
    //phiR = mt.radians(phi)
    thetaR = theta
    phiR = phi

    xmax = mt.sin(thetaR)

    if phiR > thetaR:
        return far(thetaR, phiR, mt.sin(thetaR), False, xmax*2)
    else:
        (nearhits, x_near_min) = near(thetaR, phiR, xmax*2) 
        farhits = far(thetaR, phiR, x_near_min, True,  xmax*2)
    return farhits + nearhits 
