#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"



class Vecteur{ //Definition d'une classe pour les vecteurs
public:
    explicit Vecteur(double x=0, double y=0, double z=0){
        coordonnees[0] = x;
        coordonnees[1] = y;
        coordonnees[2] = z;
    };
    double operator[](int i) const {return coordonnees[i]; };
    double &operator[](int i) {return coordonnees[i]; };
    double NormeAuCarre() { //norme au carre
        return coordonnees[0]*coordonnees[0]+ coordonnees[1]*coordonnees[1]+ coordonnees[2]*coordonnees[2];
    }
    Vecteur normaliser() { //pour normaliser un vecteur
        double norme = sqrt(NormeAuCarre());
        return Vecteur(coordonnees[0] / norme,coordonnees[1] / norme,coordonnees[2] / norme);
    }
private:
    double coordonnees[3];
};

Vecteur operator+(const Vecteur& a, const Vecteur& b) { //somme de deux vecteurs
    return Vecteur(a[0] + b[0],a[1] + b[1],a[2] + b[2]);
}

Vecteur operator-(const Vecteur& a, const Vecteur& b) { //difference de deux vecteurs
    return Vecteur(a[0] - b[0],a[1] - b[1],a[2] - b[2]);
}

Vecteur operator*(double a, const Vecteur& b) { //multiplication d'un vecteur pour un double
    return Vecteur(a*b[0],a*b[1],a*b[2]);
}

Vecteur operator*(const Vecteur& a, double b) { //idem
    return Vecteur(a[0]*b,a[1]*b,a[2]*b);
}

Vecteur operator/(const Vecteur& a, double b) { //division d'un vecteur pour un double
    return Vecteur(a[0]/b,a[1]/b,a[2]/b);
}


double prodScal(const Vecteur& a, const Vecteur& b) { //produit scalaire entre 2 vecteurs
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


class Rayon{ //Classe pour les rayons (caracterises par leur centre et leur direction)
public:
    Rayon(const Vecteur& C, const Vecteur& u) : C(C), u(u){
    }
    Vecteur C,u;
};

class Sphere{ //classe pour les spheres
public:
    Sphere(const Vecteur& O, double R, const Vecteur &albedo): O(O), R(R), albedo(albedo) {
    }
    bool intersection(const Rayon& r, Vecteur& P, Vecteur& N, double &t){
        // resout a*te + b*t + c = 0
        double a=1;
        double b=2*prodScal(r.u, r.C - O);
        double c= (r.C-O).NormeAuCarre() - R*R;
        double delta = b*b - 4 * a*c;
        if (delta < 0) return false; //pas de solution, donc pas d'intersection

        double racinedelta = sqrt(delta);
        double t2 = (-b + racinedelta)/(2*a);
        if (t2 < 0) return false; //pas de solution positive, donc pas d'intersection

        double t1 = (-b - racinedelta)/(2*a);
        if (t1 > 0) // t est la plus petite valeur positive entre t1 et t2
            t=t1;

        else
            t=t2;

        P = r.C + t*r.u; // Point d'intersection
        N = (P-O).normaliser(); //Normale au point d'intersection

        return true;

    };
    Vecteur O;
    double R;
    Vecteur albedo;
};

class Scene{ // classe de la scene, qui peut comporter plusieurs spheres
public:
    Scene() {};
    bool intersection(const Rayon& r, Vecteur& P, Vecteur& N, Vecteur &albedo){
        double t =1E9;
        bool intersecte = false;
        for (int i = 0; i < objets.size(); i++) { // pour chacune des spheres de la scene
            Vecteur Pobjet, Nobjet;
            double tobjet;
            if (objets[i].intersection(r, Pobjet, Nobjet, tobjet) && tobjet<t){ // s'il y a une intersection plus proche que celles existantes, alors on prend en compte celle-ci et pas les autres
                intersecte = true;
                t=tobjet;
                P=Pobjet;
                N=Nobjet;
                albedo = objets[i].albedo;
            }

        }
        return intersecte;
    }
    std::vector<Sphere> objets;

};

int main() {
	int W = 512;
	int H = 512;

	 Vecteur C(0, 0, 55);
    Scene scene;
    Sphere S1(Vecteur(0, 0, 0), 10, Vecteur(1., 0.3, 0.2));
    Sphere S2(Vecteur(10, 10, 0), 4, Vecteur(0., 1., 0.));
    Sphere Sgauche(Vecteur(0, 0, 1000), 940, Vecteur(1., 1., 1.));
    Sphere Sdroite(Vecteur(0, 0, -1000), 940, Vecteur(1., 1., 1.));
    Sphere Shaut(Vecteur(0, 1000, 0), 940, Vecteur(1., 0., 0.));
    Sphere Sbas(Vecteur(0, -1000, 0), 990, Vecteur(0., 0., 1.));
    Sphere Splafond(Vecteur(1000, 0, 0), 940, Vecteur(1., 1., 1.));
    Sphere Ssol(Vecteur(-1000, 0, 0), 940, Vecteur(1., 1., 1.));

    scene.objets.push_back(S1);
    scene.objets.push_back(S2);
    scene.objets.push_back(Sgauche);
    scene.objets.push_back(Sdroite);
    scene.objets.push_back(Shaut);
    scene.objets.push_back(Sbas);
    scene.objets.push_back(Splafond);
    scene.objets.push_back(Ssol);

    double fov = 60 * M_PI / 180;
    double I =5E9; //intensite lumineuse
    Vecteur L(-10,20,40); //source lumineuse

	std::vector<unsigned char> image(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

            Vecteur u(j-W/2,i-H/2, -W / (2.*tan(fov / 2)));
            u=u.normaliser();

            Rayon r(C, u);
            Vecteur P,N,albedo ;
            Vecteur couleur(0,0,0);
            if (scene.intersection(r,P,N,albedo)) { // en cas d'intersection avec un objet de la scene
                double normePL = sqrt((L-P).NormeAuCarre()); //norme de PL
				double prov = std::max(0., prodScal(N,(L-P)/normePL));
                couleur = I/(4*M_PI*normePL*normePL) *albedo/M_PI *prov;
            }

            couleur[0] = std::pow(couleur[0],0.45); //application du gamma
            couleur[1] = std::pow(couleur[1],0.45);
            couleur[2] = std::pow(couleur[2],0.45);

			image[((H- i - 1)*W + j) * 3 + 0] = std::min(255., couleur[0]);
			image[((H- i - 1)*W + j) * 3 + 1] = std::min(255., couleur[1]);
			image[((H- i - 1)*W + j) * 3 + 2] = std::min(255., couleur[2]);
		}
	}
	stbi_write_png("image_gamma.png", W, H, 3, &image[0], 0);

	return 0;
}
