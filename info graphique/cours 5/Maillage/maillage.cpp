#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.14159265359

#include <algorithm>
#include <complex>
#include <random>
#include <string>
#include <iostream>
#include <stdio.h>

static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform(0, 1);

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};


class Vecteur{
    public:
        explicit Vecteur(double x=0, double y=0, double z=0) {
            coords[0] = x;
            coords[1] = y;
            coords[2] = z;
        };
        double operator[](int i) const {return coords[i]; };
        double &operator[](int i) {return coords[i]; };
        double sqrNorm() {
            return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
        }
        Vecteur get_normalized(){
            double n = sqrt(sqrNorm());
            return Vecteur(coords[0]/n, coords[1]/n, coords[2]/n);
        }
        Vecteur& operator+=(const Vecteur& a){
            coords[0] += a[0];
            coords[1] += a[1];
            coords[2] += a[2];
            return *this;
        }

    private:
        double coords[3];
    };

Vecteur operator+(const Vecteur& a, const Vecteur& b){
    return Vecteur(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vecteur operator-(const Vecteur& a, const Vecteur& b){
    return Vecteur(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vecteur operator-(const Vecteur& a){
    return Vecteur(-a[0], -a[1], -a[2]);
}
Vecteur operator*(double a, const Vecteur& b){
    return Vecteur(a*b[0], a*b[1], a*b[2]);
}
Vecteur operator*(const Vecteur& a, double b){
    return Vecteur(b*a[0], b*a[1], b*a[2]);
}
Vecteur operator*(const Vecteur& a, const Vecteur& b){
    return Vecteur(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vecteur operator/(const Vecteur& a, double b){
    return Vecteur(a[0]/b, a[1]/b, a[2]/b);
}
Vecteur cross(const Vecteur& a, const Vecteur& b){
    return Vecteur(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}
double prodScal(const Vecteur& a, const Vecteur& b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
double sqr(const float a){
    return a*a;
}



Vecteur random_cos(const Vecteur &N){
    double u1 = uniform(engine);
    double u2 = uniform(engine);
    double x = cos(2*M_PI*u1)*sqrt(1-u2);
    double y = sin(2*M_PI*u1)*sqrt(1-u2);
    double z = sqrt(u2);
    Vecteur T1;
    if (N[0] < N[1] && N[0] < N[2]){
        T1 = Vecteur(0, N[2], -N[1]);
    }
    else{
        if (N[1] < N[2] && N[1] < N[0]){
            T1 = Vecteur(N[2], 0, -N[0]);
        }
        else{
            T1 = Vecteur(N[1], -N[0], 0);
        }

    }
    T1 = T1.get_normalized();
    Vecteur T2 = cross(N, T1);
    return z*N + x*T1 + y*T2;
}


class Ray
{
    public:
        explicit Ray(const Vecteur& C, const Vecteur& u): C(C), u(u) {
        };
        Vecteur C, u;
};

class Object{
public:
    Object(){};
    virtual bool intersect(const Ray& r, Vecteur& P, Vecteur& N, double& t) = 0;
    bool isMirror, isTransp, isCreux;
    double n, R;
    Vecteur albedo, O;
};

class TriangleMesh : public Object{
public:
  ~TriangleMesh() {}
    TriangleMesh(const Vecteur& albedo, bool mirror = false, bool transp=false, double n=1.4, bool creux=false) {
    this->albedo = albedo;
    this->isMirror = mirror;
    this->isCreux = creux;
    this->isTransp = transp;
    this->n = n;
    };
    void readOBJ(const char* obj) {

        char matfile[255];
        char grp[255];

        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ') {
                Vecteur vec;

                Vecteur col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);

                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vecteur vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vecteur vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char* consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }

            }

        }
        fclose(f);

    }

    bool intersect(const Ray& r, Vecteur& P, Vecteur& normal, double& t){
            t = 1E9;
            bool has_inter = false;
            for (int i = 0; i < indices.size(); i++){
                const Vecteur &A = vertices[indices[i].vtxi];
                const Vecteur &B = vertices[indices[i].vtxj];
                const Vecteur &C = vertices[indices[i].vtxk];
                Vecteur e1 = B - A;
                Vecteur e2 = C - A;
                Vecteur N = cross(e1, e2);
                Vecteur AO = r.C - A;
                Vecteur AO_cross_u = cross(AO, r.u);
                double invU_cross_N = 1./prodScal(r.u, N);
                double beta = - prodScal(e2, AO_cross_u)*invU_cross_N;
                double gamma = prodScal(e1, AO_cross_u)*invU_cross_N;
                double alpha = 1 - beta - gamma;
                double localt = -prodScal(AO, N)*invU_cross_N;
                if (beta >= 0 && gamma >= 0 && beta <= 1 && gamma <= 1 && alpha >= 0 && localt > 0){
                    has_inter = true;
                    if (localt < t){
                        t = localt;
                        normal = -N.get_normalized();
                        P = r.C + t * r.u;
                    }
                }
            }
            return has_inter;

        };

    std::vector<TriangleIndices> indices;
    std::vector<Vecteur> vertices;
    std::vector<Vecteur> normals;
    std::vector<Vecteur> uvs;
    std::vector<Vecteur> vertexcolors;

};

class Sphere : public Object {
    public:
        explicit Sphere(const Vecteur& O, double R, const Vecteur& albedo, bool isMirror = false, bool isTransp = false, double n = 1.4, bool isCreux = false) {
            this->albedo = albedo;
            this->isMirror = isMirror;
            this->isTransp = isTransp;
            this->isCreux = isCreux;
            this->n = n;
            this->R = R;
            this->O = O;
        };
        bool intersect(const Ray& r, Vecteur& P, Vecteur& N, double& t){
            double a = 1;
            double b = 2*prodScal(r.u, r.C-O);
            double c = (r.C-O).sqrNorm() - R*R;
            double delta = b*b - 4*a*c;
            if (delta < 0) return false;

            double sqDelta = sqrt(delta);
            double t2 = (-b + sqDelta)/ (2*a);
            if (t2 < 0) return false;

            double t1 = (-b - sqDelta) / (2*a);
            if (t1 > 0)
                t = t1;
            else
                t = t2;

            P = r.C + t*r.u;
            N = (P - O).get_normalized();
            return true;
        };
};


class Scene{
public:
    Scene() {};
    bool intersect(const Ray& r, Vecteur& P, Vecteur& N, Vecteur &albedo, bool &mirror, bool &transp, double &n, double &t, int &objectId){
        t = 1E10;
        bool has_inter = false;
        for(int i = 0; i < objects.size(); i++){
            Vecteur localP, localN;
            double localt;
            if (objects[i]->intersect(r, localP, localN, localt) && localt < t){
                has_inter = true;
                t = localt;
                P = localP;
                N = localN;
                if (objects[i]->isCreux) {
                    N = -N;
                }
                albedo = objects[i]->albedo;
                mirror = objects[i]->isMirror;
                transp = objects[i]->isTransp;
                n = objects[i]->n;
                objectId = i;
            }
        }
        return has_inter;
     }

     Vecteur getColor(const Ray& r, int rebond, bool last_diffu) {
            Vecteur P, N, albedo;
            double t, n2;
            bool mirror, transp;
            int objectId;
            bool inter = intersect(r, P, N, albedo, mirror, transp, n2, t, objectId);
            Vecteur color(0, 0, 0);
            if (rebond > 10) return color;

            if (inter) {
                if (objectId == 0){
                    if (rebond == 0 || !last_diffu)
                        return Vecteur(I, I, I)/(4*M_PI*M_PI*objects[0]->R*objects[0]->R);
                    else
                        return Vecteur(0., 0., 0.);
                        }
                if (mirror) {
                    Vecteur reflectedDir = r.u - 2*prodScal(r.u, N)*N;
                    Ray reflectedRay(P + 1E-5*N, reflectedDir);
                    return getColor(reflectedRay, rebond + 1, false);
                } else {
                    if (transp) {
                        double n1 = 1;
                        Vecteur N2 = N;
                        if (prodScal(r.u, N) > 0) {
                            std::swap(n1, n2);
                            N2 = -N;
                        }
                        Vecteur Tt = n1 / n2 * (r.u - prodScal(r.u, N2)*N2);
                        double rad = 1 - sqr(n1 / n2) * (1 - sqr(prodScal(r.u, N2)));
                        if (rad < 0) {
                            Vecteur reflectedDir = r.u - 2*prodScal(r.u, N)*N;
                            Ray reflectedRay(P + 1E-5*N, reflectedDir);
                            return getColor(reflectedRay, rebond + 1, false);

                        }
                        Vecteur Tn = -sqrt(rad)*N2;
                        Vecteur refractedDir = Tt + Tn;
                        Ray refractedRay(P - 1E-5*N2, refractedDir);
                        return getColor(refractedRay, rebond + 1, false);
                    }
                    else {
                        /*Vecteur PL = L - P;
                        double d = sqrt(PL.sqrNorm());
                        Vecteur shadowP, shadowN, shadowAlbedo;
                        double shadowt, shadown2;
                        Ray shadowRay(P+1E-5*N, PL/d);
                        bool shadowMirror, shadowTransp;
                        int shadowId;
                        bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadown2, shadowt, shadowId);
                        if (shadowInter && shadowt < d) {
                        } else {
                            color = I/(4*M_PI*d*d)*albedo/M_PI*std::max(0., prodScal(N, PL/d));
                        }*/
                        // Eclairage direct
                        Vecteur PL = L - P;
                        PL = PL.get_normalized();
                        Vecteur random_w = random_cos(-PL);
                        Vecteur xprime = random_w*objects[0]->R + objects[0]->O;
                        Vecteur Pxprime = xprime - P;
                        double distance = sqrt(Pxprime.sqrNorm());
                        Pxprime = Pxprime/distance;

                        Vecteur shadowP, shadowN, shadowAlbedo;
                        double shadowt, shadown2;
                        Ray shadowRay(P+1E-4*N, Pxprime);
                        bool shadowMirror, shadowTransp;
                        int shadowId;
                        bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadown2, shadowt, shadowId);
                        if (shadowInter && shadowt < distance - 2E-4) {
                            color = Vecteur(0., 0., 0.);
                        } else {
                            double Jacobien = std::max(1E-8, prodScal(random_w, -Pxprime))/(distance*distance);
                            double probabilite = std::max(1E-8, prodScal(-PL, random_w))/(M_PI*objects[0]->R*objects[0]->R);
                            color = I/(4*M_PI*M_PI*objects[0]->R*objects[0]->R)*albedo/M_PI*std::max(0., prodScal(N, Pxprime))*Jacobien/probabilite;
                        }

                        // Eclairage indirect
                        Vecteur wi = random_cos(N);
                        Ray Rayw_inc(P + 1E-5*N, wi);
                        color += albedo*getColor(Rayw_inc, rebond + 1, true);

                    }
                }
            };

            return color;
        }

    std::vector<Object*> objects;
    Vecteur L;
    double I;
};


int main() {
    int W = 128;
    int H = 128;
    int nbrays = 1;

    Scene scene;
    Vecteur C(0, 0, 55);
    double dist_mise_pt = 55.;
    scene.L = Vecteur(-10, 20, 40);
    Sphere Lumiere(scene.L, 15, Vecteur(1., 1., 1.));
    Sphere S1(Vecteur(-12, 0, 0), 10, Vecteur(1., 1., 1.), true);
    Sphere S2(Vecteur(8, 4, 0), 10, Vecteur(1., 0., 0.), false, true, 1.4);
    Sphere S3(Vecteur(8, 4, 0), 9.8, Vecteur(1., 0., 0.), false, true, 1.4, true);
    Sphere S4(Vecteur(0, -5, 15), 5, Vecteur(1., 1., 0.));
    Sphere Ssol(Vecteur(0, -1000, 0), 990, Vecteur(1., 1., 1.));
    Sphere Smur1(Vecteur(-1000, 0, 0), 940, Vecteur(1., 0., 0.));
    Sphere Smur2(Vecteur(1000, 0, 0), 940, Vecteur(1., 0., 0.));
    Sphere Smur3(Vecteur(0, 0, 1000), 940, Vecteur(0., 0., 0.));
    Sphere Smur4(Vecteur(0, 0, -1000), 940, Vecteur(0., 1., 0.));
    Sphere Splafond(Vecteur(0, 1000, 0), 940, Vecteur(0., 0., 1.));
    TriangleMesh m(Vecteur(1., 1., 1.));
    m.readOBJ("cat_object.obj");
    double angle_rot = 130.;
    double theta = angle_rot/360.*2*M_PI;
    for(int i=0; i<m.vertices.size(); i++){
        m.vertices[i][2] -= 15;
        m.vertices[i][0] /= 2;
        m.vertices[i][1] /= 2;
        m.vertices[i][2] /= 2;
        double old_vertices = m.vertices[i][0];
        m.vertices[i][0] = cos(theta)*m.vertices[i][0] - sin(theta)*m.vertices[i][1];
        m.vertices[i][1] = sin(theta)*old_vertices + cos(theta)*m.vertices[i][1];
        m.vertices[i][1] += 20;
        m.vertices[i][2] -= 3;
        std::swap(m.vertices[i][1], m.vertices[i][2]);
    }
    scene.objects.push_back(&Lumiere);
    /*scene.objects.push_back(&S1);
    scene.objects.push_back(&S2);
    scene.objects.push_back(&S3);
    scene.objects.push_back(&S4);*/
    scene.objects.push_back(&m);

    scene.objects.push_back(&Ssol);
    scene.objects.push_back(&Splafond);
    scene.objects.push_back(&Smur1);
    scene.objects.push_back(&Smur2);
    scene.objects.push_back(&Smur3);
    scene.objects.push_back(&Smur4);
    double fov = 60 * M_PI / 180;
    scene.I = 5E9;

    std::vector<unsigned char> image(W*H * 3, 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        std::cout << i << '\n';
        for (int j = 0; j < W; j++) {

            Vecteur color(0,0,0);
            for (int k=0; k<nbrays; k++){
                double u1 = uniform(engine);
                double u2 = uniform(engine);
                double x1 = 0.25 * cos(2*M_PI*u1)*sqrt(-2*log(u2));
                double x2 = 0.25 * sin(2*M_PI*u1)*sqrt(-2*log(u2));

                double y1 = 1 * cos(2*M_PI*u1)*sqrt(-2*log(u2));
                double y2 = 1 * sin(2*M_PI*u1)*sqrt(-2*log(u2));

                Vecteur u(j - W/2 + x2 + 0.5, i - H/2 + x1 + 0.5, - W / (2.*tan(fov/2)));
                u = u.get_normalized();

                Vecteur Cible = C + dist_mise_pt*u;
                Vecteur Cprime = C + Vecteur(y1, y2, 0);
                Vecteur uprime = (Cible - Cprime).get_normalized();

/*                Ray r(Cprime, uprime);
*/                Ray r(C, u);
                color += scene.getColor(r, 0, false);
            }
            color = color / nbrays;
            image[((H - i - 1)* W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));
            image[((H - i - 1)* W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));
            image[((H - i - 1)* W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));
            }
        }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}
