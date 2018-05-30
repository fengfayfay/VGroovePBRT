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
Float computeThetaForPhi0(const Vector3f& w) {
    Float cosTheta = w.z;
    Float cosPhi = 1;
    Float sinPhi = 0; 
    Float sinTheta = w.x/cosPhi;
    Float theta = atan2(sinTheta, cosTheta);
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

Transform MyRotateY(Float theta) {
    Float sinTheta = std::sin(theta);
    Float cosTheta = std::cos(theta);
    Matrix4x4 m(cosTheta, 0, sinTheta, 0, 0, 1, 0, 0, -sinTheta, 0, cosTheta, 0,
                0, 0, 0, 1);
    return Transform(m, Transpose(m));
}

Transform MyRotateZ(Float theta) {
    Float sinTheta = std::sin(theta);
    Float cosTheta = std::cos(theta);
    Matrix4x4 m(cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0, 0, 0, 0, 1, 0,
                0, 0, 0, 1);
    return Transform(m, Transpose(m));
}

struct Frame {

    void createRotationZ(const Vector3f& wh) {
        Float t = atan2(wh.y, wh.x);
        /*
        if (t < 0) {
            t += Pi * 2.0;
        }
        */
        frameTheta = -t;
        w2l = MyRotateZ(frameTheta);
        l2w = MyRotateZ(-frameTheta);
    }

    void createRotationY(const Vector3f& wh) {
        Float t = atan2(wh.x, wh.z);
        /*
        if (t < 0) {
            t += Pi * 2.0;
        }
        */
        frameTheta = -t;
        w2l = MyRotateY(frameTheta);
        l2w = MyRotateY(-frameTheta);
    }

    Vector3f worldToLocal(const Vector3f& wW) const{
        return w2l(wW);
    }
    Vector3f localToWorld(const Vector3f& wl) const {
        //this is a bit of hack to take advantage of rotation being symmetric

        Vector3f t(wl);
        if (flipped) {
            t.x = -t.x;
            //N.y = -N.y;
        }
        return Vector3f(l2w(t));
    }
    Transform w2l, l2w;
    Float frameTheta;
    bool flipped;
};
    
struct EvalFrame : public Frame {
    EvalFrame() {flipped = false; };
    
    EvalFrame(const Vector3f& owo, const Vector3f& owi, bool autoFlip = true) : owo(owo), owi(owi) {
        wo = Normalize(owo);
        wi = Normalize(owi);

        owh = Normalize((wo+wi));
        
        createRotationZ(owh);
        wo = worldToLocal(wo);
        wi = worldToLocal(wi);
        wh = worldToLocal(owh); 

        CHECK(rel_eq(wo.y, -wi.y));

        wop = Vector3f(wo.x, 0, wo.z);
        wop = Normalize(wop);
        wip = Vector3f(wi.x, 0, wi.z);
        wip = Normalize(wip);
        theta_o = computeThetaForPhi0(wop);
        theta_i = computeThetaForPhi0(wip);
        flipped = false;

        //y value and x value needs to be flipped 
        if (theta_o < 0 && autoFlip) {
            theta_o *= -1;
            theta_i *= -1;
            wo.x *= -1;
            wi.x *= -1;
            wop.x *= -1;
            wip.x *= -1;
            //wo.y *= -1;
            //wi.y *= -1;
            flipped = true;
        }
    }
    Vector3f wo, wi, wh, wop, wip, owo, owi, owh;
    Float theta_o, theta_i;
};

struct SampleFrame:public EvalFrame {
    SampleFrame(const Vector3f&  oWo, const Vector3f& oWh, bool autoFlip = true) {
        owo = oWo;
        owh = oWh;
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
            //wo.y *= -1;
            wh.x *= -1;
            wop.x *= -1;
            flipped = true;
        }
    }

    Vector3f constructWi(Float t_i) {
        theta_i = t_i;
        wip = Vector3f(sin(theta_i), 0, cos(theta_i)); 
        wi = Vector3f(0, 0, 0);
        wi.y = -wo.y;
        Float scale = sqrt(1.0 - wi.y * wi.y);

        wi.x = wip.x * scale;
        wi.z = wip.z * scale;
        Float lenW = Dot(wi, wi);
        CHECK(rel_eq(lenW, 1.0));

        //important:
        // the localWi should be computed according to local wo, wh
        //only flipped for final output
        return owi =  localToWorld(wi);
    }
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
            Float pDxdxa = pDxDy.x;
            Float pDydya = pDxDy.y;

            Vector3f h = getH(bounce);
            Float kp = Dot(wp, h);
            w = Reflect(-wp, h);
            if (fresnel) F *= fresnel->Evaluate(Dot(w, h));
            Float dxdxa = 0, dydya = 0;
            if (bounce % 2 == 0) {
                Float pDzdxa = -wp.x/wp.z * pDxdxa;
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

    Float computeJacobian(int bounce, Spectrum& F) {
        Vector3f wr;
        Vector2f dxy = computeDxaDya(bounce, wr, F);

        if (!rel_eq(wr, wi)) {
            Vector2f dxy = computeDxaDya(bounce, wr, F);
            std::cout << "computed outgoing direction: " << wr <<"\n";
            std::cout << "expected outgoing direction: " << wi <<"\n";
            fflush(stdout);
        }
        Float denom = fabs(dxy.x * dxy.y);
        if (denom < 1e-6) return 0;
        Float nom = fabs(wi.z);
        Float jacobian = nom/denom;
        return jacobian;
    }

    const Fresnel* fresnel;
    Vector3f N, Ng, H, HP, wo, wi, owo, owi, owh;
    Float dxpdxa, dypdya, dzpdxa;

};

//thetaM is in (0, .5pi) so cosThetaM > 0
//left ZipinNormal should be the same xsign as wop, right ZipinNormal has -xsign as wop
Vector3f 
computeZipinNormal(Float thetaM, char side, const Vector3f& wop) {
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


Float 
VGrooveReflection::computeGFactor(const EvalFrame& evalFrame, VGroove& vgroove, int bounce, 
                   char side, Vector3f& wm) const {

    //debug single bounce

    Float thetaM = 0;
    Float GFactor = vgroove.inverseEval(evalFrame.theta_o, evalFrame.theta_i, bounce, side, thetaM);
    
    if (GFactor > 0) {
        GFactor = std::min((Float)1.0, GFactor);
        wm = computeZipinNormal(thetaM, side, evalFrame.wop);
        if (bounce == 1) {
            //no relation between frameTheta and thetaM (frameTheta is related to phi_h not theta_h
            Vector3f wh = Normalize((evalFrame.wo + evalFrame.wi)*.5);
            Float mh = Dot(wm, wh);
            if (!rel_eq(mh, 1.f)){
                std::cout<<"wm: "<< wm << " wh: "<<wh<<"\n";
                fflush(stdout);
            }
        }
    }
    return GFactor;
}


Float 
VGrooveReflection::computeBounceBrdf(const EvalFrame& evalFrame, VGroove& vgroove, int bounce, char side, 
                   Float& pdf, Spectrum &F, Float& weight) const {
    
    Vector3f wm;
    Float GFactor = computeGFactor(evalFrame, vgroove, bounce, side, wm);
    Float NGFactor = vgroove.sumG > 1e-6? GFactor/vgroove.sumG : 0;

    pdf = 0;
    weight = 0;
    F = Spectrum(1.f);
    Float brdf(0);
    if (GFactor > 0) {
        Vector3f owm = evalFrame.localToWorld(wm);
        if (Dot(evalFrame.wo, wm) < 1e-6) return brdf;
        Float value = microfacetReflectionWithoutG(evalFrame.owo, evalFrame.owi, owm);
        Float mpdf = microfacetPdf(evalFrame.owo, owm);
        //Float value = microfacetReflectionWithoutG(evalFrame.wo, evalFrame.wi, wm);
        //Float mpdf = microfacetPdf(evalFrame.wo, wm);
        Jacobian jacobian(evalFrame.wo, evalFrame.wi, wm, fresnel);
        Float Jac = jacobian.computeJacobian(bounce, F);
        if (Jac > 0) {
            Float J = Jac * Dot(evalFrame.wo, wm) * 4; 
            if (bounce == 1) {
                if (!rel_eq(J, 1, 1e-3)) {
                    std::cout << "J is different from 1: "<< J << "\n";
                    fflush(stdout);
                    //Float J = jacobian.computeJacobian(bounce, F) * Dot(evalFrame.wo, wm) * 4;
                }
            }
            brdf = value * J * GFactor;
            if (NGFactor + 1e-4 < GFactor) {
                std::cout << "NG: "<< NGFactor << " G:"<< GFactor<< "\n";
                fflush(stdout);
            }
            //pdf = distribution->Pdf(evalFrame.owo, owm) * Jac * NGFactor;
            pdf = mpdf * J * NGFactor;
            weight = 1;
            return brdf;
        }
    }
    return brdf;
}

Spectrum 
VGrooveReflection::eval(const EvalFrame& evalFrame, const Vector3f &wo, const Vector3f &wi, Float& pdf) const {

    //pdf = .5/Pi;
    //return MicrofacetReflection::f(wo, wi);
    Spectrum brdf(0);
    pdf = 0;
    if (!SameHemisphere(wo, wi)) return brdf;
    if (evalFrame.theta_o < 1e-6) return brdf;

    VGroove vgroove(maxBounce, minBounce);

    Float weightSum = 0; 
    
    for (int n = minBounce; n<=maxBounce; n++) {
        Float tpdf = 0,  weight = 0; 
        Spectrum F;
        Float tbrdf = computeBounceBrdf(evalFrame, vgroove, n, 'r', tpdf, F, weight);
        brdf += tbrdf * F * weight;
        pdf += tpdf;
        weightSum += weight;
        if (n == 1 && tbrdf > 0) continue;
        tbrdf = computeBounceBrdf(evalFrame, vgroove, n, 'l', tpdf, F, weight);
        brdf += tbrdf * F * weight;
        pdf += tpdf;
        weightSum += weight;
    }
    return Spectrum(brdf);
    //return R*brdf;
}

int weightedRandomChoice(std::vector<Hit> hits, Float sumG, Float& prob)
{
    //Float u = .5;
    prob = 0;
     
    int validHits = hits.size();
    if (validHits == 1) {
        prob = 1;
        return 0;
    }

    if (validHits == 2) { 
        Float weight = hits[0].GFactor/sumG;
        Float u1 = (((Float) rand())/(RAND_MAX));
        if (u1 < weight) {
            prob = weight;
            return 0;
        } else {
            prob = 1.0 - weight;
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
    if (sampleFrame.theta_o < 1e-6) {
        if (pdf) *pdf = 0;
        return Spectrum(0);
    }

    Float grooveTheta = asin(wh.z);
    VGroove vgroove(maxBounce, minBounce);
    if (sampleFrame.wh.x * sampleFrame.wo.x >0) {
        vgroove.leftEvalOnly(grooveTheta, sampleFrame.theta_o);
    } else {
        vgroove.rightEvalOnly(grooveTheta, sampleFrame.theta_o);
    }
    if (vgroove.theHits.size() > 0) {
        Float prob = 0;
        int choice = weightedRandomChoice(vgroove.theHits, vgroove.sumG, prob);
        if (choice >= 0) {
            Hit hit = vgroove.theHits[choice];
            *wi = sampleFrame.constructWi(hit.thetaR);
            Float tmppdf = 0;
            Spectrum brdf = eval(sampleFrame, wo, *wi, tmppdf);
            if (pdf) *pdf = tmppdf;
            return brdf;
            /*
            Float brdf = microfacetReflectionWithoutG(wo, *wi, wh);
            Float tpdf = microfacetPdf(wo, wh);
            Jacobian jacobian(sampleFrame.wo, sampleFrame.wi, sampleFrame.wh, fresnel);
            Spectrum F(1.f);
            Float J = jacobian.computeJacobian(hit.bounce, F) * Dot(sampleFrame.wo, sampleFrame.wh) * 4;
            brdf *= hit.GFactor * J;
            tpdf *= prob * J;

            if (pdf) *pdf = tpdf;
            return R*F*brdf;
            */
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
    Float pdfstub;
    EvalFrame evalFrame(wo, *wi);
    return eval(evalFrame, wo, *wi, pdfstub); 
}


Spectrum 
VGrooveReflection:: f(const Vector3f &wo, const Vector3f &wi) const {
    Float pdf = 0;
    EvalFrame evalFrame(wo, wi);
    return eval(evalFrame, wo, wi, pdf);
}

Float 
VGrooveReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {

    //return .5/Pi;

    //always return uniform pdf for MIS for now
    if (uniSample) {
        return .5/Pi;
    } else {
        //real pdf computation
        Float pdf = 0;
        EvalFrame evalFrame(wo, wi);
        Spectrum brdf = eval(evalFrame, wo, wi, pdf);
        return pdf;
    }
}

Float
VGrooveReflection::microfacetReflectionWithoutG(const Vector3f& wo, const Vector3f& wi, 
                   const Vector3f& wh) const {

    Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
    // Handle degenerate cases for microfacet reflection
    if (cosThetaI == 0 || cosThetaO == 0) return 0;
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) return 0;
    //Spectrum F = fresnel->Evaluate(Dot(wi, wh));
    return distribution->D(wh) / (4 * cosThetaI * cosThetaO);
}

Float
VGrooveReflection::microfacetPdf(const Vector3f& wo, const Vector3f& wh) const {
    return distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));        
}

Float
VGrooveReflection::computePdfIntegral(Float thetaI) const {
    int phiCount = 100;
    int muCount = 200;
    Float sum = 0;

    Vector3f wo(sin(thetaI), 0, cos(thetaI));
    for (int i = 0; i < phiCount; i++) {
        Float phi = i;
        phi = Pi * 2.0 * phi /phiCount;
        Float cosPhi = cos(phi);
        Float sinPhi = sin(phi);

        for (int j = 0; j < muCount; j++) {
            Float mu = j + 0.5;
            mu /= muCount;
            Float sinTheta = sqrt(1 - mu * mu);
            Vector3f wi(sinTheta * cosPhi, sinTheta * sinPhi, mu);
            Float pdf = Pdf(wo, wi);
            sum += pdf;
        }
    }
    return sum*Pi*2/(muCount * phiCount);
}

bool
VGrooveReflection::testPDF() const {
    for (int i = 0; i < 90; i++){
        Float thetaI = Radians(i - 0.5);
        Float pdfIntegral = computePdfIntegral(thetaI);
        if (pdfIntegral -1e-2 > 1.0 || 1.0 - pdfIntegral > 1e-3) {
            std::cout << "pdfIntegral problem: " << pdfIntegral << " thetaI" << i << "\n";
            return false;
        }
    }
    return true;
}
/*
void
VGroove::saveOutHits(Float theta, Float phi) {
}

void 
VGroove::VGrooveTest() {

    Float phi_offset = 0.6;
    Float theta_offset = 0.3;
    for (int i = 0; i < 90; i++) {
        Float theta = radians(theta_offset + i);
        for (int j = 0; j < 90; j++) {
            Float phi = radians(phi_offset + j);
            leftHitOnly(theta, phi);
            saveOut(theta, phi, theHits);  
            rightHitOnly(theta, phi);
            saveOut(theta, phi, theHits);  
        }
    }
}
*/       

} //end namespace
