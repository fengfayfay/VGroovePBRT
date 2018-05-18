#include<stdio.h>
#include <sampling.h>
#include "vgroove_reflect.h"

namespace pbrt {

//EvalFrame
//SampleFrame
//Jacobian


inline bool rel_eq(Float x, Float y, Float thresh = 1e-5) { return fabs(x-y) < thresh; }

inline bool rel_eq(const Vector3f&v1, const Vector3f& v2, Float thresh = 1e-3) {
    return rel_eq(v1.x, v2.x, thresh) && rel_eq(v1.y, v2.y, thresh) && rel_eq(v1.z, v2.z, thresh); 
}

//phi is always 0
float computeThetaForPhi0(const Vector3f& w) {
    float cosTheta = w.z;
    float cosPhi = 1;
    float sinPhi = 0; 
    float sinTheta = w.x/cosPhi;
    float theta = atan2(sinTheta, cosTheta);
    return theta;
}
/*
class Frame:
    def __init__(self, X, Z):
        self.Z = Z
        self.Y = self.Z.cross(X).norm()
        self.X = self.Y.cross(self.Z).norm()

    def globalToLocal(self,w):
        x = w.dot(self.X)
        y = w.dot(self.Y)
        z = w.dot(self.Z)
        return vec3.Vec3(x, y, z)

    def localToGlobal(self, w):
        g = self.X * w.x + self.Y * w.y + self.Z * w.z
        return g
*/

struct Frame {

    void createRotationZ(const Vector3f& wo, const Vector3f& wi){
        Vector3f wh = Normalize((wo + wi) *.5);
        createRotationZ(wh);
    }
    void createRotationZ(const Vector3f& wh) {
        float t = atan2(wh.y, wh.x);
        /*
        if (t < 0) {
            t += Pi * 2.0;
        }
        */
        frameTheta = -t;
        w2l = RotateZ(frameTheta * 180.0/Pi);
    }

    void createRotationY(const Vector3f& wh) {
        float t = atan2(wh.x, wh.z);
        /*
        if (t < 0) {
            t += Pi * 2.0;
        }
        */
        frameTheta = -t;
        w2l = RotateY(frameTheta * 180.0/Pi);
    }

    Vector3f worldToLocal(const Vector3f& wW) const{
        return w2l(wW);
    }
    Vector3f localToWorld(const Vector3f& wl) const {
        //this is a bit of hack to take advantage of rotation being symmetric
        Normal3f N(wl.x, wl.y, wl.z);
        return Vector3f(w2l(N));
    }
    Transform w2l;
    float frameTheta;
};
    
struct EvalFrame : public Frame {
    EvalFrame(const Vector3f& owo, const Vector3f& owi, bool autoFlip = true) : owo(owo), owi(owi) {
        wo = Normalize(owo);
        wi = Normalize(owi);
        
        createRotationZ(wo, wi);
        wo = worldToLocal(wo);
        wi = worldToLocal(wi);

        CHECK(rel_eq(wo.y, -wi.y));

        wop = Vector3f(wo.x, 0, wo.z);
        wop = Normalize(wop);
        wip = Vector3f(wi.x, 0, wi.z);
        wip = Normalize(wip);
        theta_o = computeThetaForPhi0(wop);
        theta_i = computeThetaForPhi0(wip);
        flipped = false;

        /*
        phi_o = atan2(wo.y, wo.x);
        phi_i = atan2(wi.y, wi.x);
        theta_or = computeTheta(wo, phi_o)
        theta_ir = computeTheta(wi, phi_i)
        */

        //y value and x value needs to be flipped 
        if (theta_o < 0 && autoFlip) {
            theta_o *= -1;
            theta_i *= -1;
            wo.x *= -1;
            wi.x *= -1;
            wop.x *= -1;
            wip.x *= -1;
            wo.y *= -1;
            wi.y *= -1;
            flipped = true;
        }
    }

    Vector3f wo, wi, wop, wip, owo, owi;
    float theta_o, theta_i;
    bool flipped;
    //float theta_or, phi_o, theta_ir, phi_i;

};
///////////////////////////////////SampleFrame//////////////////////////
struct SampleFrame:public Frame {

    SampleFrame(const Vector3f&  owo, const Vector3f& owh, bool autoFlip = true) {
        wo = Normalize(owo);
        createRotationZ(owh);

        wo = worldToLocal(wo);
        wh = worldToLocal(owh);

        CHECK(rel_eq(wh.y, 0));

        wop = Normalize(Vector3f(wo.x, 0, wo.z));
        theta_o = computeThetaForPhi0(wop);

        //phi_o = atan2(wo.y, wo.x)
        //theta_or = computeTheta(self.wo, self.phi_o)
        flipped = false;

        //flip frame around z axis if theta_o is negative for zipin
        //no need to flip for ray traced vgrooves
        if (theta_o < 0 && autoFlip) {
            theta_o *= -1;
            wo.x *= -1;
            wo.y *= -1;
            wh.x *= -1;
            wop.x *= -1;
            flipped = true;
        }
    }

    Vector3f constructWi(float theta_i, Vector3f& localWi) {
        Vector3f wip = Vector3f(sin(theta_i), 0, cos(theta_i)); 
        Vector3f wi(0, 0, 0);
        wi.y = -wo.y;
        float scale = sqrt(1.0 - wi.y * wi.y);

        wi.x = wip.x * scale;
        wi.z = wip.z * scale;
        float lenW = Dot(wi, wi);
        CHECK(rel_eq(lenW, 1.0));

        localWi = wi;
        //important:
        // the localWi should be computed according to local wo, wh
        //only flipped for final output
        if (flipped) {
            wi.x *= -1;
            wi.y *= -1;
        }
        return localToWorld(wi);
    }
    Vector3f wo, wh, wop;
    float theta_o;
    bool flipped;
};

struct Jacobian: public Frame{
    Jacobian(const Vector3f& owo, const Vector3f& owi, const Vector3f& owh, const Fresnel* fresnel): owo(owo), owi(owi), owh(owh), fresnel(fresnel) {

        //Vector3f owhp = Vector3f(-wh.x, 0, wh.z);
        //wih = vec3.dot(wi, wh)

        createRotationY(owh);

        N = Vector3f(0, 0, 1);
        Ng = worldToLocal(N);
        H = worldToLocal(owh);
        HP = Reflect(H, Ng);

        CHECK(rel_eq(H, Vector3f(0, 0, 1)));
        //assert(vec3.Equal(H, vec3.Vec3(0, 0, 1)))
        //assert(vec3.Equal(self.HP, self.rotateH.rotate(vec3.Vec3(-wh.x, 0, wh.z))))

        wo = worldToLocal(owo);
        wi = worldToLocal(owi);
        //assert(vec3.EqualValue(wi.z, wih))

        dxpdxa = - 1.0 + 2.0 * Ng.x * Ng.x;
        dypdya = - 1.0 + 2.0 * Ng.y * Ng.y;
        dzpdxa = 2.0 * Ng.x * Ng.z;
    }
    
    const Vector3f& getH(int bounce) const {
        return (bounce % 2) ? H: HP;
    }

    Vector2f computeDxaDya(int bounce, Vector3f&w, Spectrum& F) {
        if (bounce == 0) {
            w = -wo;
            F = Spectrum(1.);
            return Vector2f(0, 0);    
        } else {
            Vector3f wp;
            Vector2f pDxDy = computeDxaDya(bounce - 1, wp, F);
            float pDxdxa = pDxDy.x;
            float pDydya = pDxDy.y;

            Vector3f h = getH(bounce);
            float kp = Dot(wp, h);
            w = Reflect(-wp, h);
            if (fresnel) F *= fresnel->Evaluate(Dot(w, h));
            float dxdxa = 0, dydya = 0;
            if (bounce % 2 == 0) {
                float pDzdxa = -wp.x/wp.z * pDxdxa;
                dxdxa = pDxdxa - 2.0 * kp * dxpdxa;
                dxdxa -= 2.0 * HP.x * (pDxdxa * HP.x + wp.x * dxpdxa + pDzdxa * HP.z + wp.z * dzpdxa);
                dydya = pDydya + 2.0 * kp;
            } else {
                dxdxa = pDxdxa - 2.0 * kp;
                dydya = pDydya - 2.0 * kp;
            }
            return Vector2f(dxdxa, dydya);
        }
    }

    float computeJacobian(int bounce, Spectrum& F) {
        Vector3f wr;
        Vector2f dxy = computeDxaDya(bounce, wr, F);

        if (!rel_eq(wr, wi)) {
            Vector2f dxy = computeDxaDya(bounce, wr, F);
            std::cout << "computed outgoing direction: " << wr <<"\n";
            std::cout << "expected outgoing direction: " << wi <<"\n";
            fflush(stdout);
        }
        float denom = fabs(dxy.x * dxy.y);
        if (denom < 1e-7) return 1;
        float nom = fabs(wi.z);
        float jacobian = nom/denom;
        return jacobian;
    }

    const Fresnel* fresnel;
    Vector3f N, Ng, H, HP, wo, wi, owo, owi, owh;
    float dxpdxa, dypdya, dzpdxa;

};

//thetaM is in (0, .5pi) so cosThetaM > 0
//left ZipinNormal should be the same xsign as wop, right ZipinNormal has -xsign as wop
Vector3f 
computeZipinNormal(float thetaM, char side, const Vector3f& wop) {
    Vector3f n(cos(thetaM), 0, sin(thetaM));
    n.x *= wop.x > 0 ? 1: -1;
    if (side == 'r' ){
        n.x *= -1;
    }
    return n;
}

VGrooveReflection::VGrooveReflection(const Spectrum &R,
                      MicrofacetDistribution *distribution, Fresnel *fresnel,
                      int maxBounce, int minBounce, bool uniSample)
        : MicrofacetReflection(R, distribution, fresnel), maxBounce(maxBounce),
          minBounce(minBounce), uniSample(uniSample) { 
    srand(time(NULL)); 
}


float 
VGrooveReflection::computeGFactor(const EvalFrame& evalFrame, VGroove& vgroove, int bounce, char side, Vector3f& wm) const {

    //debug single bounce

    float thetaM = 0;
    float GFactor = vgroove.inverseEval(evalFrame.theta_o, evalFrame.theta_i, bounce, side, thetaM);
    
    if (GFactor > 0) {
        GFactor = std::min(1.0f, GFactor);
        wm = computeZipinNormal(thetaM, side, evalFrame.wop);
        if (bounce == 1) {
            //no relation between frameTheta and thetaM (frameTheta is related to phi_h not theta_h
            Vector3f wh = Normalize((evalFrame.wo + evalFrame.wi)*.5);
            float mh = Dot(wm, wh);
            CHECK(rel_eq(mh, 1));
        }
    }
    return GFactor;
}


float 
VGrooveReflection::computeBounceBrdf(const EvalFrame& evalFrame, VGroove& vgroove, int bounce, char side, 
                    float& pdf, Spectrum &F) const {
    Vector3f wm;
    float GFactor = computeGFactor(evalFrame, vgroove, bounce, side, wm);
    float brdf(0);
    if (GFactor > 0) {
        Vector3f owm = evalFrame.localToWorld(wm);
        float value = microfacetReflectionWithoutG(evalFrame.owo, evalFrame.owi, owm);
        float mpdf = microfacetPdf(evalFrame.owo, owm);
        Jacobian jacobian(evalFrame.wo, evalFrame.wi, wm, fresnel);
        //float J = jacobian.computeJacobian(bounce, F) * Dot(wo, wm) * 4;
        float J = jacobian.computeJacobian(bounce, F) * Dot(evalFrame.wo, wm) * 4;
        if (bounce == 1) {
            if (!rel_eq(J, 1, 1e-3)) {
                std::cout << "J is different from 1: "<< J << "\n";
                fflush(stdout);
            }
        }
        brdf = value * J * GFactor;
        pdf = mpdf * J * GFactor;     
    }
    //brdf = MicrofacetReflection::f(evalFrame.wo, evalFrame.wi);
    return brdf;
}

Spectrum 
VGrooveReflection::eval(const Vector3f &wo, const Vector3f &wi, float& pdf) const {

    //pdf = .5/Pi;
    //return MicrofacetReflection::f(wo, wi);
    if (!SameHemisphere(wo, wi)) return Spectrum(0.f);
    EvalFrame evalFrame(wo, wi);

    VGroove vgroove;
    Spectrum brdf(0);
    pdf = 0;
    for (int n = minBounce; n<=maxBounce; n++) {
        float tpdf; 
        Spectrum F;
        float tbrdf = computeBounceBrdf(evalFrame, vgroove, n, 'r', tpdf, F);
        brdf += tbrdf * F;
        pdf += tpdf;
        if (n == 1 && tbrdf > 0) continue;
        tbrdf = computeBounceBrdf(evalFrame, vgroove, n, 'l', tpdf, F);
        brdf += tbrdf * F;
        pdf += tpdf;
    }
    return Spectrum(brdf);
    //return R*brdf;
}

int weightedRandomChoice(std::vector<Hit> hits, int maxBounce, int minBounce, float& prob)
{
    //float u = .5;
    prob = 0;

    int hitCount = hits.size();
    if (hitCount == 1) {
        prob = 1;
        return 0;
    }

    if (hitCount == 2) {
        float u1 = ((float) rand()/(RAND_MAX));
        if (u1 < hits[0].GFactor) {
            prob = hits[0].GFactor;
            return 0;
        } else {
            prob = hits[1].GFactor;
            return 1;
        }
    } 

    return -1;
}

Spectrum 
VGrooveReflection::Sample_f(const Vector3f &owo, Vector3f *wi, const Point2f &u,
                       Float *pdf, BxDFType *sampledType) const {
    Vector3f wo = Normalize(owo); 

    if (wo.z <=  1e-6) {
        if (pdf) *pdf = 0;
        return Spectrum(0);
    }

    if (uniSample) {    
        return UniSample_f(wo, wi, u, pdf, sampledType);
    } else { 
   
    Vector3f wh = distribution->Sample_wh(wo, u);
    //*wi = Reflect(wo, wh);
    //return eval(wo, *wi, *pdf);
    
    SampleFrame sampleFrame(wo, wh);
    float grooveTheta = asin(wh.z);
    VGroove vgroove;
    if (sampleFrame.wh.x >0) {
        vgroove.leftEvalOnly(grooveTheta, sampleFrame.theta_o);
    } else {
        vgroove.rightEvalOnly(grooveTheta, sampleFrame.theta_o);
    }
    float prob;

    //need to make sure the function below is doing the right thing
    //should i only choose from one side of the groove? since wh is only
    //on one side I think so.
    if (vgroove.theHits.size() > 0) {
        int choice = weightedRandomChoice(vgroove.theHits, 3, 1, prob);
        if (choice >= 0) {
            Hit hit = vgroove.theHits[choice];
            if (hit.bounce >= minBounce && hit.bounce <= maxBounce) {
            Vector3f sampleFrameWi;
            *wi = sampleFrame.constructWi(hit.thetaR, sampleFrameWi); 
            float brdf = microfacetReflectionWithoutG(wo, *wi, wh);
            float tpdf = microfacetPdf(wo, wh);
            //Jacobian jacobian(wo, *wi, wh, fresnel);
            Jacobian jacobian(sampleFrame.wo, sampleFrameWi, sampleFrame.wh, fresnel);
            Spectrum F;
            //float J = jacobian.computeJacobian(hit.bounce, F) * Dot(wo, wh) * 4;
            float J = jacobian.computeJacobian(hit.bounce, F) * Dot(sampleFrame.wo, sampleFrame.wh) * 4;
            brdf *= hit.GFactor * J;
            tpdf *= prob * J;

            if (pdf) *pdf = tpdf;
            return R*F*brdf;
            }
        }
    }
    if (pdf) *pdf = 0;
    return Spectrum(0);
    }
}

Spectrum 
VGrooveReflection::UniSample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                       Float *pdf, BxDFType *sampledType) const {

    if (pdf) *pdf = .5/Pi;
    *wi = UniformSampleHemisphere(u);
    float pdfstub;
    return eval(wo, *wi, pdfstub); 
}


Spectrum 
VGrooveReflection:: f(const Vector3f &wo, const Vector3f &wi) const {
    float pdf = 0;
    return eval(wo, wi, pdf);
}

Float 
VGrooveReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {

    return .5/Pi;

    //always return uniform pdf for MIS for now
    if (uniSample) {
        return .5/Pi;
    } else {
        //real pdf computation
        float pdf = 0;
        Spectrum brdf = eval(wo, wi, pdf);
        return pdf;
    }
}

float
VGrooveReflection::microfacetReflectionWithoutG(const Vector3f& wo, const Vector3f& wi, 
                   const Vector3f& wh) const {

    Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
    // Handle degenerate cases for microfacet reflection
    if (cosThetaI == 0 || cosThetaO == 0) return 0;
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) return 0;
    //Spectrum F = fresnel->Evaluate(Dot(wi, wh));
    return distribution->D(wh) / (4 * cosThetaI * cosThetaO);
}

float
VGrooveReflection::microfacetPdf(const Vector3f& wo, const Vector3f& wh) const {
    float factor = Dot(wo, wh);
    if (factor < 1e-5) return 0;
    return distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));        
}

/*
void
VGroove::saveOutHits(float theta, float phi) {
}

void 
VGroove::VGrooveTest() {

    float phi_offset = 0.6;
    float theta_offset = 0.3;
    for (int i = 0; i < 90; i++) {
        float theta = radians(theta_offset + i);
        for (int j = 0; j < 90; j++) {
            float phi = radians(phi_offset + j);
            leftHitOnly(theta, phi);
            saveOut(theta, phi, theHits);  
            rightHitOnly(theta, phi);
            saveOut(theta, phi, theHits);  
        }
    }
}
*/       

} //end namespace
