#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include "ConfigFile.hpp"
#include "ConfigFile.tpp"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std;
typedef vector<complex<double> > vec_cmplx;

// Fonction résolvant le système d'équations A * solution = rhs
// où A est une matrice tridiagonale
template <class T> void triangular_solve(vector<T> const& diag,
                                         vector<T> const& lower,
                                         vector<T> const& upper,
                                         vector<T> const& rhs,
                                         vector<T>& solution)
{
  vector<T> new_diag = diag;
  vector<T> new_rhs = rhs;
  
  // forward elimination
  for(int i(1); i<diag.size(); ++i) 
  {
    T pivot = lower[i-1] / new_diag[i-1];
    new_diag[i] -= pivot * upper[i-1];
    new_rhs[i] -= pivot * new_rhs[i-1];
  }
  
  solution.resize(diag.size());
  
  // solve last equation
  solution[diag.size()-1] = new_rhs[diag.size()-1] / new_diag[diag.size()-1];
  
  // backward substitution
  for(int i = diag.size() - 2; i >= 0; --i) 
  {
    solution[i] = (new_rhs[i] - upper[i] * solution[i+1]) / new_diag[i];
  }
}


// Potentiel V(x) :
double V(double const& x, double const& omega, double const& delta,double const& x0)
{
    //return .5*omega*omega*min((x-delta)*(x-delta),(x+delta)*(x+delta));
    
    if ((x < x0-delta*0.5)or(x > x0+delta*0.5)){
        return 0;
    }else{
        return 1.;
    }
    
    //return 0;
}


// Déclaration des diagnostiques de la particule d'après sa fonction d'onde psi :
//  - prob calcule la probabilité de trouver la particule entre les points nL.dx et nR.dx,
//  - E calcule son énergie,
//  - xmoy calcule sa position moyenne,
//  - x2moy calcule sa position au carré moyenne,
//  - pmoy calcule sa quantité de mouvement moyenne,
//  - p2moy calcule sa quantité de mouvement au carré moyenne.
// Les définitions de ces fonctions sont en dessous du main.
double prob(vec_cmplx const& psi, int nL, int nR, double dx);
double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx);
double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double pmoy(vec_cmplx const& psi, double const& dx);
double p2moy(vec_cmplx const& psi, double const& dx);


int main(int argc,char **argv)
{
  complex<double> complex_i = complex<double> (0,1); // Nombre imaginaire i
  
  string inputPath = "configuration.in"; // Fichier d'input par défaut
  if(argc>1) // Fichier d'input spécifié par l'utilisateur ("./Onde config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath);

  for(int i(2); i<argc; ++i) // Input complémentaires ("./Onde config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);
  
  // Paramètres physiques :
  double hbar    = 1.;
  double m       = 1.;
  double tfin    = configFile.get<double>("tfin");
  double xL      = configFile.get<double>("xL");
  double xR      = configFile.get<double>("xR");
  double omega   = configFile.get<double>("omega");
  double delta   = configFile.get<double>("delta");
  double x0      = configFile.get<double>("x0");
  double k0      = 2. * M_PI * configFile.get<int>("n") / (xR-xL);
  double sigma0  = configFile.get<double>("sigma_norm") * (xR-xL);
  
  // Paramètres numériques :
  double dt      = configFile.get<double>("dt");
  int Ninters    = configFile.get<int>("Ninters");
  int Npoints    = Ninters + 1;
  double dx      = (xR-xL) / Ninters;
  
  // Maillage :
  vector<double> x(Npoints);
  for(int i(0); i<Npoints; ++i)
    x[i] = xL + i*dx;

  // Initialisation et normalisation de la fonction d'onde :
  vec_cmplx psi(Npoints);
  // Initialiser le paquet d'onde, équation (4.107) du cours
  // Normaliser la fonction d'onde à l'aide de la fonction prob()
  for(int i(0);i<Npoints;i++){
	  psi[i]=exp(complex_i*k0*x[i]-0.5*((x[i]-x0)/sigma0)*((x[i]-x0)/sigma0));
  }
  double norme(sqrt(prob(psi,0,Npoints-1,dx)));
  
  for(int i(0);i<Npoints;i++){
	  psi[i]/=norme;
  }
  // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
  vec_cmplx dH(Npoints), aH(Ninters), cH(Ninters); // matrice Hamiltonienne
  vec_cmplx dA(Npoints), aA(Ninters), cA(Ninters); // matrice du membre de gauche de l'équation (4.90)
  vec_cmplx dB(Npoints), aB(Ninters), cB(Ninters); // matrice du membre de droite de l'équation (4.90)
  
  complex<double> a, b; // Coefficients complexes de l'équation (4.91)
  double a_1(hbar/(4.0*m*dx*dx));
  a=complex_i*dt*a_1;
  
  for(int i(0);i<Npoints;i++){
	    double b_1(V(x[i],omega,delta,x0)/(hbar*2.0));
		b=complex_i*dt*b_1;
		// Calculer les éléments des matrices A, B et H. 
		// Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales supérieures et inférieures
		dH[i]=2.0*hbar*(2.0*a_1+b_1);
		dA[i]=1.0+2.0*a+b;
		dB[i]=1.0-2.0*a-b;
	}

  for(int i(0);i<Ninters;i++){
	  aH[i]=-2.0*hbar*a_1;
	  cH[i]=aH[i];
	  aA[i]=-a;
	  cA[i]=aA[i];
	  aB[i]=a;
	  cB[i]=aB[i];
  }
  // Conditions aux limites: psi nulle aux deux bords
  //Modifier les matrices A et B pour satisfaire les conditions aux limites
	dA[0]=1.;
	cA[0]=0.;
	aA[Ninters-1]=0.;
	dA[Ninters]=1.;
	
	dB[0]=0.;
	cB[0]=0.;
	aB[Ninters-1]=0.;
	dB[Ninters]=0.;
	

  // Fichiers de sortie :
  string output = configFile.get<string>("output");
  
  ofstream fichier_potentiel(output + "_pot.dat");
  fichier_potentiel.precision(15);
  for(int i(0); i<Npoints; ++i)
    fichier_potentiel << x[i] << " " << V(x[i], omega, delta,x0) << endl;
  fichier_potentiel.close();
  
  ofstream fichier_psi(output + "_psi2.dat");
  fichier_psi.precision(15);
  
  ofstream fichier_observables(output + "_obs.dat");
  fichier_observables.precision(15);
  
  // Boucle temporelle :
  double t;
  for(t=0.; t<tfin; t+=dt)
  {
    // Ecriture de |psi|^2 :
    for(int i(0); i<Npoints; ++i)
      fichier_psi << abs(psi[i]) * abs(psi[i]) << " ";
    fichier_psi << endl;
    
    // Ecriture des observables :
    fichier_observables << t << " "
                        << prob(psi,0,Ninters/2,dx) << " "       // probabilité que la particule soit à gauche
                        << prob(psi,Ninters/2,Ninters,dx) << " " // probabilité que la particule soit à droite
                        << E(psi,dH,aH,cH,dx) << " "  // Energie
                        << xmoy(psi,x,dx) << " "      // Position moyenne
                        << x2moy(psi,x,dx) << " "     // Position^2 moyenne
                        << pmoy(psi,dx) << " "      // Quantité de mouvement moyenne
                        << p2moy(psi,dx) << " "     // (Quantité de mouvement)^2 moyenne
                        << sqrt(x2moy(psi,x,dx)-pow(xmoy(psi,x,dx),2)) << " "
					    << sqrt(p2moy(psi,dx)-pow(pmoy(psi,dx),2))
                        << endl;   
    
    // Calcul du membre de droite psi_tmp = B * psi :
    vec_cmplx psi_tmp(Npoints,0.);
    for(int i(0); i<Npoints; ++i)
      psi_tmp[i] = dB[i] * psi[i];
    for(int i(0); i<Ninters; ++i) 
    {
      psi_tmp[i] += cB[i] * psi[i+1];
      psi_tmp[i+1] += aB[i] * psi[i];
    }
    
    // Résolution de A * psi = psi_tmp :
    triangular_solve(dA, aA, cA, psi_tmp, psi);
  
  } // Fin de la boucle temporelle
  
  for(int i(0); i<Npoints; ++i)
    fichier_psi << abs(psi[i]) * abs(psi[i]) << " ";
  
  fichier_observables << t << " "
                      << prob(psi,0,Ninters/2,dx) << " "
                      << prob(psi,Ninters/2,Ninters,dx) << " "
                      << E(psi,dH,aH,cH,dx) << " "
                      << xmoy(psi,x,dx) << " "
                      << x2moy(psi,x,dx) << " "
                      << pmoy(psi,dx) << " "
                      << p2moy(psi,dx) << " "
					  << sqrt(x2moy(psi,x,dx)-pow(xmoy(psi,x,dx),2)) << " "
					  << sqrt(p2moy(psi,dx)-pow(pmoy(psi,dx),2))
					  << endl;
  fichier_observables.close();
  fichier_psi.close();
}


double prob(vec_cmplx const& psi, int nL, int nR, double dx)
{
  // Calculer la probabilite de trouver la particule entre les points nL.dx et nR.dx
  double temp(0.0);
  for(int i(nL);i<nR;i++){
	  temp+=real(psi[i]*conj(psi[i]) + psi[i+1]*conj(psi[i+1]));
  }
  return 0.5*dx*temp;
}


double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx)
{
  //Calculer la moyenne de l'Hamiltonien

  // H(psi)
  // On utilise la matrice H calculée plus haut
  //...
  // Integrale de psi* H(psi) dx
  //...

  int psi_sz(psi.size());
  vec_cmplx psi_temp(psi_sz);
  complex<double> energie(0.0);
  psi_temp[0]=diagH[0]*psi[0]+upperH[0]*psi[1];
  for(int i(1);i<psi_sz-1;i++){
	  psi_temp[i]=lowerH[i-1]*psi[i-1]+diagH[i]*psi[i]+upperH[i]*psi[i+1];
  }
  psi_temp[psi_sz-1]=lowerH[psi_sz-2]*psi[psi_sz-2]+diagH[psi_sz-1]*psi[psi_sz-1];
  for(int i(0);i<psi_sz-1;i++){
	  energie+=conj(psi[i+1])*psi_temp[i+1]+conj(psi[i])*psi_temp[i];
  }
  return real(0.5*dx*energie);
}


double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
  //Calculer la moyenne de la position
  int x_size(x.size());
  double x_temp(0.0);
  for(int i(0);i<x_size;i++){
	  x_temp+=real(conj(psi[i+1])*x[i+1]*psi[i+1]+conj(psi[i])*x[i]*psi[i]);
  }
  return 0.5*dx*x_temp;
}


double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
  //Calculer la moyenne du x^2
  int x2_size(x.size());
  double x2_temp(0.0);
  for(int i(0);i<x2_size;i++){
	  x2_temp+=real(conj(psi[i+1])*x[i+1]*x[i+1]*psi[i+1]+conj(psi[i])*x[i]*x[i]*psi[i]);
  }
  return 0.5*dx*x2_temp;
}


double pmoy(vec_cmplx const& psi, double const& dx)
{
  //Calculer la moyenne de p 
  int psi_size(psi.size());
  double hbar(1.);
  complex<double> complex_i(0,1);
  vector<complex<double>> ppsi(psi_size-2);
  vector<complex<double>> psistar(psi_size-2);
  for(int i(0);i<psi_size-2;i++){
	  ppsi[i]=-complex_i*hbar*(psi[i+2]-psi[i])/(2.0*dx);
	  psistar[i]=conj(psi[i+1]);
  }
  complex<double> mean_p (0.0);
  for(int i(0);i<psi_size-3;i++){
	  mean_p+=psistar[i]*ppsi[i]+psistar[i+1]*ppsi[i+1];
  }
  return real(0.5*dx*mean_p);
}


double p2moy(vec_cmplx const& psi, double const& dx)
{
  // Calculer la moyenne du p^2
  int psi_size(psi.size());
  double hbar(1.);
  complex<double> complex_i(0,1);
  vector<complex<double>> p2psi(psi_size-2);
  vector<complex<double>> psistar(psi_size-2);
  for(int i(0);i<psi_size-2;i++){
	  p2psi[i]=-hbar*hbar*(psi[i+2]-2.0*psi[i+1]+psi[i])/(dx*dx);
	  psistar[i]=conj(psi[i+1]);
  }
  complex<double> mean_p2 (0.0);
  for(int i(0);i<psi_size-3;i++){
	  mean_p2+=psistar[i]*p2psi[i]+psistar[i+1]*p2psi[i+1];
  }
  return real(0.5*dx*mean_p2);
}



