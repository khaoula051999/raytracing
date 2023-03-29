#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

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
    Sphere(const Vecteur& O, double R): O(O), R(R) {
    }
    bool intersection(const Rayon& r){
        // resout a*te + b*t + c = 0, et renvoie s'il y a une solution
        double a=1;
        double b=2*prodScal(r.u, r.C - O);
        double c= (r.C-O).NormeAuCarre() - R*R;
        double delta = b*b - 4 * a*c;
        if (delta >= 0) return true;

        return false;
    };
    Vecteur O;
    double R;
};


int main() {
	int W = 512;
	int H = 512;

	Vecteur C(0, 0, 55);
    Vecteur O(0,0,0);
    double R = 10; //rayon de la sphere
    Sphere S(O,R);
    double fov = 60 * M_PI / 180;

	std::vector<unsigned char> image(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

            Vecteur u(j-W/2,i-H/2, -W / (2.*tan(fov / 2)));
            u=u.normaliser();
            Rayon r(C, u);
            Vecteur couleur(0,0,0);
            if (S.intersection(r)) { //en cas d'intersection avec la sphere, on met le pixel blanc
                couleur = Vecteur(255,255,255);
            }

			image[(i*W + j) * 3 + 0] = couleur[0];
			image[(i*W + j) * 3 + 1] = couleur[1];
			image[(i*W + j) * 3 + 2] = couleur[2];
		}
	}
	stbi_write_png("image_blanche.png", W, H, 3, &image[0], 0);

	return 0;
}
