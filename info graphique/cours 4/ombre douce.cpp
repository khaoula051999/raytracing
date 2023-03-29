#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <random>


static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0, 1);

class Vecteur { //Definition d'une classe pour les vecteurs
public:
    explicit Vecteur(double x = 0, double y = 0, double z = 0) {
        coordonnees[0] = x;
        coordonnees[1] = y;
        coordonnees[2] = z;
    };
    double operator[](int i) const { return coordonnees[i]; };
    double& operator[](int i) { return coordonnees[i]; };
    double NormeAuCarre() { //norme au carre
        return coordonnees[0] * coordonnees[0] + coordonnees[1] * coordonnees[1] + coordonnees[2] * coordonnees[2];
    }
    Vecteur normaliser() { //pour normaliser un vecteur
        double norme = sqrt(NormeAuCarre());
        return Vecteur(coordonnees[0] / norme, coordonnees[1] / norme, coordonnees[2] / norme);
    }
private:
    double coordonnees[3];
};

Vecteur operator+(const Vecteur& a, const Vecteur& b) { //somme de deux vecteurs
    return Vecteur(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vecteur operator-(const Vecteur& a, const Vecteur& b) { //difference de deux vecteurs
    return Vecteur(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vecteur operator-(const Vecteur& a) { //inverse d'un vecteur
    return Vecteur(-a[0], -a[1], -a[2]);
}

Vecteur operator*(double a, const Vecteur& b) { //multiplication d'un vecteur pour un double
    return Vecteur(a * b[0], a * b[1], a * b[2]);
}


Vecteur operator*(const Vecteur& a, double b) { //idem
    return Vecteur(a[0] * b, a[1] * b, a[2] * b);
}

Vecteur operator/(const Vecteur& a, double b) { //division d'un vecteur pour un double
    return Vecteur(a[0] / b, a[1] / b, a[2] / b);
}

Vecteur prodVect(const Vecteur& a, const Vecteur& b) { //produit vectoriel entre 2 vecteurs
    return Vecteur(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);

}

double prodScal(const Vecteur& a, const Vecteur& b) { //produit scalaire entre 2 vecteurs
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vecteur prodTermeaterme(const Vecteur& a, const Vecteur& b) { //produit terme e terme entre 2 vecteur
    return Vecteur(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}



Vecteur aleatoire(const Vecteur& N) {//Renvoie un vecteur z*N+x*T1+y*T2, avec (N,T1,T2) un repere, et x,y,z aleatoires suivant une loi de probabilite cosinus
    double u1 = uniform(engine); // nombre aleatoire entre 0 et 1
    double u2 = uniform(engine);
    double x = cos(2 * M_PI * u1) * sqrt(1 - u2);
    double y = sin(2 * M_PI * u1) * sqrt(1 - u2);
    double z = sqrt(u2);
    Vecteur T1;
    if (N[0] < N[1] && N[0] < N[2]) {
        T1 = Vecteur(0, N[2], -N[1]);
    }
    else {
        if (N[1] < N[2] && N[1] < N[0]) {
            T1 = Vecteur(N[2], 0, -N[0]);
        }
        else {
            T1 = Vecteur(N[1], -N[0], 0);
        }
    }

    T1 = T1.normaliser();
    Vecteur T2 = prodVect(N, T1);
    return x * T1 + y * T2 + z * N;

}



class Rayon { //Classe pour les rayons (caracterises par leur centre et leur direction)
public:
    Rayon(const Vecteur& C, const Vecteur& u) : C(C), u(u) {
    }
    Vecteur C, u;
};

class Sphere { //classe pour les spheres (rajout d'un attribut pour savoir si c'est un miroir ou pas, et un autre pour la transparence)
public:
    Sphere(const Vecteur& O, double R, const Vecteur& albedo, bool Miroir = false, bool Transparent = false) : O(O), R(R), albedo(albedo), Miroir(Miroir), Transparent(Transparent) {
    }
    bool intersection(const Rayon& r, Vecteur& P, Vecteur& N, double& t) {
        // resout a*te + b*t + c = 0
        double a = 1;
        double b = 2 * prodScal(r.u, r.C - O);
        double c = (r.C - O).NormeAuCarre() - R * R;
        double delta = b * b - 4 * a * c;
        if (delta < 0) return false; //pas de solution, donc pas d'intersection

        double racinedelta = sqrt(delta);
        double t2 = (-b + racinedelta) / (2 * a);
        if (t2 < 0) return false; //pas de solution positive, donc pas d'intersection

        double t1 = (-b - racinedelta) / (2 * a);
        if (t1 > 0) // t est la plus petite valeur positive entre t1 et t2
            t = t1;

        else
            t = t2;

        P = r.C + t * r.u; // Point d'intersection
        N = (P - O).normaliser(); //Normale au point d'intersection

        return true;

    };
    Vecteur O;
    double R;
    Vecteur albedo;
    bool Miroir, Transparent;
};

class Scene { // classe de la scene, qui peut comporter plusieurs spheres
public:
    Scene() {};
    bool intersection(const Rayon& r, Vecteur& P, Vecteur& N, Vecteur& albedo, double& t, bool& miroir, bool& transp, int& id) {
        t = 1E9;
        bool intersecte = false;
        for (int i = 0; i < objets.size(); i++) {// pour chacune des spheres de la scene
            Vecteur Pobjet, Nobjet;
            double tobjet;
            if (objets[i].intersection(r, Pobjet, Nobjet, tobjet) && tobjet < t) { //s'il y a une intersection plus proche que celles existantes, alors on prend en compte celle-ci et pas les autres
                intersecte = true;
                t = tobjet;
                P = Pobjet;
                N = Nobjet;
                albedo = objets[i].albedo;
                miroir = objets[i].Miroir;
                transp = objets[i].Transparent;
				id=i;
            }

        }
        return intersecte;
    }

    Vecteur obtenirCouleur(const Rayon& r, int rebond) { //pour obtenir la couleur
        if (rebond > 5) return Vecteur(0., 0., 0.); // si on depasse 5 rebonds, on renvoit la couleur noire
        else {
            double t;
            bool miroir, transp;
            Vecteur P,N,albedo;
            Vecteur couleur(0,0,0);
			int id;
            if (intersection(r, P, N, albedo, t, miroir, transp,id)) {// en cas d'intersection avec un objet de la scene
				if (id == 0) { //la 1ere sphere correspond e la lumiere
                    return Vecteur(I, I, I) / (4 * M_PI * M_PI * objets[0].R * objets[0].R);
                }

                if (miroir) { // si c'est un miroir
                    Vecteur Directionreflechie = r.u - 2 * prodScal(r.u, N) * N; // direction de la reflexion
                    Rayon Rayonreflechi(P + 0.00001 * N, Directionreflechie); //rayon reflechi, partant du point d'intersection et allant dans la direction reflechie
                    return obtenirCouleur(Rayonreflechi, rebond + 1);
                }
                else {
                    if (transp) { // si  il est transparent

                        double n1 = 1., n2 = 1.4; //indices optiques
                        Vecteur N2 = N;
                        if (prodScal(r.u, N) > 0) { //si on sort de la sphere, on inverse les n et la normale
                            std::swap(n1, n2);
                            N2 = -N;
                        }
                        double angle = 1 - n1 * n1 / (n2 * n2) * (1 - prodScal(r.u, N2) * prodScal(r.u, N2));
                        if (angle < 0) { // si l'angle est plus petit que 0, il y a reflexion
                            Vecteur Directionreflechie = r.u - 2 * prodScal(r.u, N) * N;// direction de la reflexion
                            Rayon Rayonreflechi(P + 0.00001 * N, Directionreflechie);//rayon reflechi, partant du point d'intersection et allant dans la direction reflechie
                            return obtenirCouleur(Rayonreflechi, rebond + 1);
                        }
                        Vecteur Tt = n1/n2 * (r.u -prodScal(r.u,N2)*N2); //Composante tangentielle de la direction de refraction
                        Vecteur Tn=-sqrt(angle)*N2; //Composante normale de la direction de refraction
                        Vecteur Directionrefractee = Tt + Tn; //direction de refraction
                        Rayon Rayonrefracte(P-0.0001*N2,Directionrefractee); //rayon refracte, partant du point d'intersection et allant dans la direction refractee
                        return obtenirCouleur(Rayonrefracte, rebond + 1);
                    }
                    else {
                        //eclairage indirect
                        Vecteur wi = aleatoire(N);
                        Rayon Rayonwi(P + 0.00001 * N, wi);//rayon aleatoire
                        couleur = couleur + prodTermeaterme(albedo, obtenirCouleur(Rayonwi, rebond + 1));


                    }
                    }
                }
                return couleur;
            }
        }
        std::vector<Sphere> objets;
        Vecteur L;
        double I;

    };


int main() {
	int W = 512;
	int H = 512;

	Vecteur C(0, 0, 55);
	Scene scene;
	scene.L = Vecteur(-10, 20, 40);//source lumineuse	
	Sphere Slumiere(scene.L, 15, Vecteur(1., 0.3, 0.2));
    Sphere S1(Vecteur(0, 0, 0), 10, Vecteur(1., 0.3, 0.2));
    Sphere S2(Vecteur(10, 10, 0), 4, Vecteur(0., 1., 0.));
    Sphere Sgauche(Vecteur(0, 0, 1000), 940, Vecteur(1., 1., 1.));
    Sphere Sdroite(Vecteur(0, 0, -1000), 940, Vecteur(1., 1., 1.));
    Sphere Shaut(Vecteur(0, 1000, 0), 940, Vecteur(1., 0., 0.));
    Sphere Sbas(Vecteur(0, -1000, 0), 990, Vecteur(0., 0., 1.));
    Sphere Splafond(Vecteur(1000, 0, 0), 940, Vecteur(1., 1., 1.));
    Sphere Ssol(Vecteur(-1000, 0, 0), 940, Vecteur(1., 1., 1.));
    scene.objets.push_back(Slumiere);
    scene.objets.push_back(S1);
    scene.objets.push_back(S2);
    scene.objets.push_back(Sgauche);
    scene.objets.push_back(Sdroite);
    scene.objets.push_back(Shaut);
    scene.objets.push_back(Sbas);
    scene.objets.push_back(Splafond);
    scene.objets.push_back(Ssol);

	double fov = 60 * M_PI / 180;
	scene.I = 5E9; //intensite lumineuse
	int nmbrayon = 100;

	std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {


			Vecteur u(j - W / 2, i - H / 2, -W / (2. * tan(fov / 2)));
			u = u.normaliser();

			Rayon r(C, u);
			Vecteur couleur(0, 0, 0);
            for (int k = 0; k < nmbrayon; k++) {
                double u1 = uniform(engine);
                double u2 = uniform(engine);
                double x1 = 0.25 * cos(2 * M_PI * u1) * sqrt(-2 * log(u2));
                double x2 = 0.25 * sin(2 * M_PI * u1) * sqrt(-2 * log(u2));
                Vecteur u(j - W / 2 + x2 + 0.5, i - H / 2 + x1 + 0.5, -W / (2. * tan(fov / 2))); // antialiasing : on ne vise plus le centre des pixels, mais un endroit aleatoire proche de ce centre
                u = u.normaliser();

                Rayon r(C, u);

                couleur = couleur + scene.obtenirCouleur(r, 0);
            }
			couleur = couleur / nmbrayon;


			couleur[0] = std::pow(couleur[0], 0.45);//application du gamma
			couleur[1] = std::pow(couleur[1], 0.45);
			couleur[2] = std::pow(couleur[2], 0.45);


			image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., couleur[0]);
			image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., couleur[1]);
			image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., couleur[2]);
		}
	}
	stbi_write_png("ombredouce.png", W, H, 3, &image[0], 0);

	return 0;
}

