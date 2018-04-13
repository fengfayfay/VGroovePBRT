#include "vgroove_reflect.h"

namespace pbrt {

//EvalFrame
//SampleFrame
//Jacobian

//phi is always 0
float computeThetaForPhi0(const Vector3f& w) {
    float cosTheta = w.z;
    float cosPhi = 1;
    float sinPhi = 0; 
    float sinTheta = w.x/cosPhi;
    float theta = math.atan2(sinTheta, cosTheta);
    return theta;
}

struct Frame {

    void createRotationZ(const Vector3f& wo, const Vector3f& wi){
        Vector3f wh = Normalize((wo + wi) *.5);
        createRotationZ(wh);
    }
    void createRotationZ(const Vector3f& wh) {
        float t = math.atan2(wh.y, wh.x);
        /*
        if (t < 0) {
            t += Pi * 2.0;
        }
        */
        float theta = -t;
        w2l = RotateZ(theta * 180.0/Pi);
    }

    void createRotationY(const Vector& wh) {
        float t = math.atan2(wh.x, wh.z);
        /*
        if (t < 0) {
            t += Pi * 2.0;
        }
        */
        float theta = -t;
        w2l = RotateY(theta * 180.0/Pi);
    }

    Vector3f worldToLocal(const Vector3f& wW) {
        return w2l(wW);
    }
    Vector3f localToWorld(const Vector3f& wl) {
        //this is a bit of hack to take advantage of rotation being symmetric
        Normal3f N(wl.x, wl.y, wl.z);
        return l2w(N);
    }
    Transform w2l;
};
    
struct EvalFrame : public Frame {
    EvalFrame(const Vector3f& owo, const Vector3f& owi, bool autoFlip = true) {
        wo = owo.norm()
        wi = owi.norm()
        
        w2l = createRotationZ(wo, wi);
        wo = worldToLocal(wo);
        wi = worldToLocal(wi);

        wop = Vector3f(wo.x, 0, wo.z);
        wop = Normalize(wop);
        wip = Vector3f(wi.x, 0, wi.z);
        wip = Normalize(wip);
        theta_o = computeThetaForPhi0(wop);
        theta_i = computeThetaForPhi0(wip);
        flipped = False

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

    Vector3f wo, wi, wop, wip;
    float theta_o, theta_i;
    bool flipped;
    //float theta_or, phi_o, theta_ir, phi_i;

};
///////////////////////////////////SampleFrame//////////////////////////
struct SampleFrame::public Frame {

    SampleFrame(const Vector3f&  owo, const Vector3f& owh, bool autoFlip = true) {
        wo = Normalize(owo);
        w2l = createRotationZ(owh);

        wo = worldToLocal(wo);
        wh = worldToLoal(wh);
        assert(math.fabs(wh.y) < .000001);

        wop = Normalize(Vector3f(wo.x, 0, wo.z);
        theta_o = computeThetaForPhi0(self.wop, 0);

        //phi_o = atan2(self.wo.y, self.wo.x)
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

    Vector3f constructWi(theta_i) {
        Vector3f wip = Vector3f(sin(theta_i), 0, cos(theta_i)); 
        Vector3f wi(0, 0, 0);
        wi.y = -wo.y;
        float scale = sqrt(1.0 - wi.y * wi.y)

        wi.x = wip.x * scale;
        wi.z = wip.z * scale;
        //lenW = wi.dot(wi);
        //assert(math.fabs(lenW - 1.0) < .0001)

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
    Jacobian(const Vector3f& wo, const Vector3f& wi, const Vector3f& wh): owo(wo), owi(wi), owh(wh) {

        Vector3f owhp = Vector3f(-wh.x, 0, wh.z);
        N = Vector3f(0, 0, 1);
        //wih = vec3.dot(wi, wh)

        w2l = createRotationY(wh);
        Ng = worldToLocal(N);
        H = worldToLocal(wh);
        HP = Reflect(H, Ng)

        //assert(vec3.Equal(self.H, vec3.Vec3(0, 0, 1)))
        //assert(vec3.Equal(self.HP, self.rotateH.rotate(vec3.Vec3(-wh.x, 0, wh.z))))

        wo = worldToLocal(wo);
        wi = worldToLocal(wi);
        //assert(vec3.EqualValue(self.wi.z, self.wih))

        dxpdxa = - 1.0 + 2.0 * Ng.x * Ng.x
        dypdya = - 1.0 + 2.0 * Ng.y * Ng.y
        dzpdxa = 2.0 * Ng.x * Ng.z
    }
    
    const Vector3f& getH(int bounce) const {
        return (bounce % 2) ? H: HP;
    }

    Vector2f computeDxaDya(int bounce, Vector3f&w) {
        if (bounce == 0) {
            w = -wo;
            return Vector2f(0, 0);    
        } else {
            Vector3f wp;
            pDxDy = computeDxaDya(bounce - 1, wp);
            float pDxdxa = pDxDy.x;
            float pDydya = pDxDy.y;

            Vector3f h = getH(bounce);
            float kp = Dot(wp, h);
            w = Reflect(-wp, h);
            float dxdxa = 0, dydya = 0;
            if (bounce % 2 == 0 {
                float pDzdxa = -wp.x/wp.z * pDxdxa;
                dxdxa = pDxdxa - 2.0 * kp * dxpdxa;
                dxdxa -= 2.0 * HP.x * (pDxdxa * HP.x + wp.x * dxpdxa + pDzdxa * HP.z + wp.z * dzpdxa);
                dydya = pDydya + 2.0 * kp
            } else {
                dxdxa = pDxdxa - 2.0 * kp
                dydya = pDydya - 2.0 * kp
            }
            return Vector2f(dxdxa, dydya);
        }
    }

    float computeJacobian(int bounce) {
        Vector3f wr;
        Vec2f dxy = computeDxaDya(bounce, wr);

        //assert(vec3.Equal(wr, self.wi))
        float denom = fabs(dxy.x * dxy.y);
        if (denom < 1e-6) return 0;
        float nom = fabs(wi.z);
        float jacobian = nom/denom;
        return jacobian;
    }

    Vector3f N, Ng, H, HP, wo, wi;
    float dxpdxa, dypdya, dzpdxa;

};

Vector3f computeZipinNormal(float thetaM, char side, const Vector3f& wop) {
    Vector3f n(math.cos(thetaM), 0, math.sin(thetaM));
    n.x *= wop.x > 0 ? 1: -1;
    if (side == 'r' ){
        n.x *= -1;
    }
    return n;
}

float comuteGFactor(EvalFrame& evalFrame, VGroove& vgroove, int bounce, char side, Vector3f& wm)
{
    float G = 0, thetaM = 0;
    bool valid = vgroove.inverseEval(evalFrame.theta_o, evalFrame.theta_i, bounce, side, thetaM, G);    
    if (valid) {
        wm = computZipinNormal(thetaM, side, evalFrame.wop);
    }
    return G;
}


float 
VGrooveReflection::computeBounceBrdf(const EvalFrame& evalFrame, VGroove& vgroove, int bounce, char side, 
                    float& pdf) const {
    Vector3f wm;
    float G = computeGFactor(evalFrame, vgroove, bounce, side, wm);
    float brdf = 0;
    if (G > 0) {
        float value = microfaceReflectionWithoutG(evalFrame.wo, evalFrame.wi, wm);
        float mpdf = microfacetReflectionPdf(evalFrame.wo, wh);
        Jacobian jacobian(evalFrame.wo, evalFrame.wi, wh);
        float J = jacobian.computeJacobian(bounce);

        //4 is due to the single bounce jacobian accountment in the microfacetReflection* functions
        float brdf = value * J * G * 4;
        pdf = mpdf * J * G * 4;     
    }
    return brdf;
}

Spectrum 
VGrooveReflection::Eval(const Vector3f &wo, const Vector3f &wi, int maxBounce, int minBounce, float& pdf) const {

    VGroove vgroove;
    EvalFrame evalFrame(wo, wi);
    float brdf = 0;
    pdf = 0;
    for (int n = minBounce; n<=maxBounce; n++) {
        float tpdf; 
        float tbrdf = computeBounceBrdf(evalFrame, vgroove, n, 'l', tpdf);
        brdf += tbrdf;
        pdf += tpdf;
        if (n == 1 && brdf > 0) continue;
        tbrdf = computeBounceBrdf(evalFrame, vgroove, n, 'r', tpdf);
        brdf += tbrdf;
        pdf += tpdf;
    }
    return R*brdf;
}

int weightedRandomChoice(std::vector<Hit> hits, int maxBounce, int minBounce, float& prob)
{
    //float u = .5;
    if (hits.size() >0) {
        return 0;
    }
    return -1;
}

Spectrum 
VGrooveReflection::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                       Float *pdf, BxDFType *sampledType) const {
    wo = wo.norm();

    if (wo.z <=  1e-6) {
        if (pdf) *pdf = 0;
        return Spectrum(0);
    }
    Vector3f wh = distribtion->Sample_wh(wo, u);
    SampleFrame sampleFrame(wo, wh);
    float grooveTheta = asin(wh.z);
    VGroove vgroove;
    vgroove.eval(grooveTheta, sampleFrame.theta_o);
    float prob;

    //need to make sure the function below is doing the right thing
    //should i only choose from one side of the groove? since wh is only
    //on one side I think so.
    int choice = weightedRandomChoice(vgroove.theHits, 3, 1, prob);
    if (choice >= 0) {
        Hit hit = vgroove.theHits[choice];
        float brdf = microfacetValueWithoutG(wo, wh);
        float tpdf = microfacetPdf(wo, wh);
        Vector3f wi = sampleFrame.constructWi(hit.thetaR); 
        Jacobian jacobian(wo, wi, wh);
        Float J = jacobian.computeJacobian(hit.n);
        brdf *= hit.G * J * 4;
        tpdf *= prob * J * 4;

        if (pdf) *pdf = tpdf;
        return R*brdf;
    } else {
        if (pdf) *pdf = 0;
        return Spectrum(0);
    }
}

Spectrum 
VGrooveReflection:: f(const Vector3f &wo, const Vector3f &wi) const {
    float pdf = 0;
    return Eval(wo, wi, 3, 1, pdf);
}

Float 
VGrooveReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {

    float pdf = 0;
    Spectrum brdf = eval(wo, wi, 3, 1, pdf);
    return pdf;
}

float
VGrooveReflection::microfacetReflectionWithoutG(const Vector3f& wo, const Vector3f& wi, 
                   const Vector3f& wh) const {

    Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
    // Handle degenerate cases for microfacet reflection
    if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0.);
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
    wh = Normalize(wh);
    //Spectrum F = fresnel->Evaluate(Dot(wi, wh));
    return distribution->D(wh) / (4 * cosThetaI * cosThetaO);
}

float
VGrooveReflection::microfacetPdf(const Vector3f& wo, const Vector3f& wh) const {
    float factor = Dot(wo, wh);
    if (factor < 1e-5) return 0;
    return distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));        
}

} //end namespace
