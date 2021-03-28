#include "PhysEngine.hpp"
#include "ConfigFile.hpp"
#include "ConfigFile.tpp"

#define Xt 0
#define Yt 1
#define Xa 2
#define Ya 3
#define Xl 4
#define Yl 5
#define Vxt 6
#define Vyt 7
#define Vxa 8
#define Vya 9
#define Vxl 10
#define Vyl 11

#define terre 1
#define asteroide 2
#define lune 3
#define asteroide_lune 4

using namespace std;

typedef valarray<double> Vecteur;


/** Constructor
 */
PhysEngine::PhysEngine(const ConfigFile configFile) :  ra({configFile.get<double>("ra0x"),configFile.get<double>("ra0y")}), rt({configFile.get<double>("rt0x"),configFile.get<double>("rt0y")}),rl({configFile.get<double>("rl0x"),configFile.get<double>("rl0y")}),va({configFile.get<double>("va0x"),configFile.get<double>("va0y")}),vt({configFile.get<double>("vt0x"),configFile.get<double>("vt0y")}),vl({configFile.get<double>("vl0x"),configFile.get<double>("vl0y")}){
  outputFile=NULL;
  try{    
    dt=configFile.get<double>("dt");
    tfin=configFile.get<double>("tfin");
    G=configFile.get<double>("G");
    Mt=configFile.get<double>("Mt");
    Ma=configFile.get<double>("Ma");
    Ml= configFile.get<double>("Ml");
    Rt=configFile.get<double>("Rt");
    Rl=configFile.get<double>("Rl");
    t=0;    
    alpha=configFile.get<double>("alpha");
    epsilon=configFile.get<double>("epsilon");
    
    R=sqrt((ra[0]-rt[0])*(ra[0]-rt[0])+(ra[1]-rt[1])*(ra[1]-rt[1]));
	  V=sqrt(va[0]*va[0]-vt[0]+va[1]*va[1]);
	  Emec=0.5*Ma*V*V-G*Ma*Mt/R;
  

  
      /*ra[0] = configFile.get<double>("ra0x");
      ra[1] = configFile.get<double>("ra0y");
      Vecteur rt(2);
      rt[0] = configFile.get<double>("rt0x");
      rt[1] = configFile.get<double>("rt0y");
      
      Vecteur rl(2);
      rl[0] = configFile.get<double>("rl0x");
      rl[1] = configFile.get<double>("rl0y");
            Vecteur va(2);

      va[0] = configFile.get<double>("va0x");
      va[1] = configFile.get<double>("va0y");
            Vecteur vt(2);

      vt[0] = configFile.get<double>("vt0x");
      vt[1] = configFile.get<double>("vt0y");
            Vecteur vl(2);

      vl[0] = configFile.get<double>("vl0x");
      vl[1] = configFile.get<double>("vl0y");*/

      Jeanbob = configFile.get<bool>("adapt");




      
    
    //manque l'inerge potentielle
    
    std::string path=configFile.get<std::string>("outputPath");    
    outputFile = new std::ofstream();
    outputFile->open( path.c_str());
    if (!outputFile->is_open())
    {
      delete outputFile;
      outputFile=NULL;
    }
    
    //cout << "dt="<< dt << "G=" << G << ", mt=" << Mt << ", ma=" << Ma ;
    //cout << ", ml=" << Ml << ", rt1=" << rt[0] << ", rl1=" << rl[0] << endl;
    //cout << "r0=" << r0 << ", v1=" << v1 ;
    //cout << ", alpha=" << alpha << "Pas de temps adaptatif=" << Jeanbob endl;
    
  }catch(std::string e){
    cerr << "Exception: " << e.c_str() << endl;
    if(outputFile!=NULL){outputFile->close();delete outputFile;}
  }
}

/** Destructor
 */
PhysEngine::~PhysEngine(){
  outputFile->close();
   cout << "delete" <<endl;

 delete outputFile;
 cout << "fin destructeur" <<endl;
}
  
int PhysengineRungeKutta4::run(){
  
  //Si il y a eu un probleme a la construction, on sort.
  if(outputFile==NULL){
    return -1;
  }
    
    
    cout << "Start the physics engine"<< " Ma=" << Ma <<"rt1=" << rt[1] << ", rl1=" << rl[1] << endl;
    
    

    //cout << "Mt=" << Mt << endl;

    //cout << "step : t=" <<  t << " dt=" << dt << "rt[O] " << rt[0] << "rt[1]" << rt[1] << endl;
    
    //cout <<"ra[0]" << ra[0] << "ra[1]" << ra[1] ;
    
    double n(0);

    double h(get_distance(asteroide));


    //adaptateur();

  while(t <= tfin-0.5*dt){
      
		

      
      ++n;
     // cout << "dans adaptateur avant : dt=" << dt << endl;

      step();

      adaptateur();
      
      R=sqrt((ra[0]-rt[0])*(ra[0]-rt[0])+(ra[1]-rt[1])*(ra[1]-rt[1]));
	  V=sqrt(va[0]*va[0]-vt[0]+va[1]*va[1]);
	  Emec=0.5*Ma*V*V-G*Ma*Mt/R;
	  
      h= get_distance(asteroide);

      
      printout(h);

      

      //cout << "sortie" << endl;
      //cout << "dans adaptateur apres : dt=" << dt << endl;
      
      //cout << "t=" << t << endl;
      

  }
  cout << "Stop the physics engine" << endl;
    //cout << "n=" << n << endl;
  
return 0;
} 

/** Print the data in the output file
 * */
void PhysEngine::printout(double h){
  *outputFile << setprecision(14) << t << " " << rt[0] << " " << rt[1] ;
  *outputFile << " " << ra[0] << " " << ra[1] << " "<< rl[0] << " " <<rl[1] << " " <<vt[0] << " " <<vt[1] << " " ;
  *outputFile <<va[0] << " " <<va[1]<< " " <<vl[0]<< " " <<vl[1] << " " << dt << " " << h << " "<< Emec << " "<< get_distance(asteroide_lune) <<endl;
}

/** Return the angle acceleration
 */
void PhysEngine::acceleration(Vecteur& A, Vecteur r0a,Vecteur r0t,Vecteur r0l,Vecteur v0a, Vecteur v0t, Vecteur v0l){
    double rta (sqrt((r0t[0] - r0a[0])*(r0t[0] - r0a[0])+(r0t[1] - r0a[1])*(r0t[1] - r0a[1])));
    //cout << "rta = " << rta << endl;
    //cout << "r0a[0] - r0t[0]=" << r0a[0] - r0t[0] << endl;
    double rla (sqrt((r0l[0] - r0a[0])*(r0l[0] - r0a[0])+(r0l[1] - r0a[1])*(r0l[1] - r0a[1])));

    double rtl (sqrt(pow((r0t[0] - r0l[0]),2)+pow((r0t[1] - r0l[1]),2)));
    //cout << "rtl = " << rtl << endl;
    
    
    A[Xt] = v0t[0];
    A[Yt] = v0t[1];
    A[Xa] = v0a[0];
    A[Ya] = v0a[1];
    A[Xl] = v0l[0];
    A[Yl] = v0l[1];

    
    A[Vxt] = -((Ma*G)/(rta*rta*rta))*(r0t[0] - r0a[0])  -((Ml*G)/(pow(rtl,3)))*(r0t[0] - r0l[0]);
    A[Vyt] = -((Ma*G)/(rta*rta*rta))*(r0t[1] - r0a[1])  -((Ml*G)/(pow(rtl,3)))*(r0t[1] - r0l[1]);
    
    
    A[Vxa] = -(r0a[0] - r0t[0])*((Mt*G)/(rta*rta*rta)) -((Ml*G)/(pow(rla,3)))*(r0a[0] - r0l[0]);
    A[Vya] = -(r0a[1] - r0t[1])*((Mt*G)/(rta*rta*rta)) -((Ml*G)/(pow(rla,3)))*(r0a[1] - r0l[1]);
    
    
    A[Vxl] = -((Mt*G)/(pow(rtl,3)))*(r0l[0] - r0t[0]) -((Ma*G)/(pow(rla,3)))*(r0l[0] - r0a[0]);
    A[Vyl] = -((Mt*G)/(pow(rtl,3)))*(r0l[1] - r0t[1])  -((Ma*G)/(pow(rla,3)))*(r0l[1] - r0a[1]);

} 


///////////////////////////////////////////////
//
// Implementation of the virtual step() method
//
//////////////////////////////////////////////


void PhysengineRungeKutta4::adaptateur(){

    if(Jeanbob){
        
       
        Vecteur ra_p (ra);
        Vecteur rt_p (rt);
        Vecteur rl_p (rl);
        Vecteur va_p (va);
        Vecteur vt_p (vt);
        Vecteur vl_p (vl);

        double t_p (t);

        long double d1(0);
        long double d2(0);
        long double d(0);
        
        do{
        

            double dt_temp(dt);
        
        step();
        
            
        
        Vecteur ra_1 (ra);
        Vecteur rt_1 (rt);
        Vecteur rl_1 (rl);
        Vecteur va_1 (va);
        Vecteur vt_1 (vt);
        Vecteur vl_1 (vl);
            
            
        ra = ra_p;
        rt = rt_p;
        rl = rl_p;
        va = va_p;
        vt = vt_p;
        vl = vl_p;
        t = t_p;
        
        
        dt = dt*0.5;
        step();
        step();
            
            cout << "rl_1[0]=" << rl_1[0] << endl;
            cout << "rl[0]=" << rl[0] <<  endl;
        d1 = sqrt((ra_1[0] - ra[0])*(ra_1[0] - ra[0]) + (ra_1[1] - ra[1])*(ra_1[1] - ra[1]));
        d2 = sqrt((rl_1[0] - rl[0])*(rl_1[0] - rl[0]) + (rl_1[1] - rl[1])*(rl_1[1] - rl[1]));
            
            
        ra = ra_p;
        rt = rt_p;
        rl = rl_p;
        va = va_p;
        vt = vt_p;
        vl = vl_p;
        t = t_p;
        dt = dt_temp;

            //cout << "d1=" << d1 << endl;
            //cout << "d2=" << d2 << endl;


            cout  << "epsilon=" << epsilon << endl;
            
            d = d1+d2;
            //d = d1;
        if (abs(d) > epsilon){
            
            dt = 0.95*dt*pow((epsilon/(abs(d))),0.2);
            cout << "dt=" << dt << endl;
            
        }else if (d <= epsilon){
            dt = dt*pow((epsilon/(abs(d))),0.2);
            
        }}while(d > epsilon);
        }
    t+= dt;
}



void PhysengineRungeKutta4::step(){

    Vecteur A(12);
   
    //cout << "Terre :"<<  rt[0] << " " << rt[1] << "Asteroide : " <<  ra[0] << " " << ra[1] << endl;

    
   acceleration (A,ra,rt,rl,va,vt,vl);
    //cout << "A=" <<  A[Xt] << " " << A[Yt] << " " << A[Xa] << " " << A[Ya] << " " << A[Xl] << " " << A[Yl] << " " << A[Vxt] << " " << A[Vyt] << " " << A[Vxa] << " " << A[Vya] << endl;

    
    
    Vecteur A1 (12);
    A1 = A*dt;
    
    Vecteur rt1({A1[Xt],A1[Yt]});
    Vecteur ra1({A1[Xa],A1[Ya]});
    Vecteur rl1({A1[Xl],A1[Yl]});
    Vecteur vt1({A1[Vxt],A1[Vyt]});
    Vecteur va1({A1[Vxa],A1[Vya]});
    Vecteur vl1({A1[Vxl],A1[Vyl]});

    //cout << "A1="<< A1[Xt] << " " << A1[Yt] << " " << A1[Xa] << " " << A1[Ya] << " " << A1[Xl] << " " << A1[Yl] << " " << A1[Vxt] << " " << A1[Vyt] << " " << A1[Vxa] << " " << A1[Vya] << endl;
    
    
    Vecteur A2 (12);

    acceleration(A2,ra + 0.5*ra1,rt + 0.5*rt1,rl + 0.5*rl1, va + 0.5*va1, vt + 0.5*vt1, vl + 0.5*vl1);
    
    A2 = dt*A2;

       // cout <<"A2=" <<  A2[Xt] << " " << A2[Yt] << " " << A2[Xa] << " " << A2[Ya] << " " << A2[Xl] << " " << A2[Yl] << " " << A2[Vxt] << " " << A2[Vyt] << " " << A2[Vxa] << " " << A2[Vya] << endl;


    Vecteur rt2({A2[Xt],A2[Yt]});
     Vecteur ra2({A2[Xa],A2[Ya]});
     Vecteur rl2({A2[Xl],A2[Yl]});
     Vecteur vt2({A2[Vxt],A2[Vyt]});
     Vecteur va2({A2[Vxa],A2[Vya]});
     Vecteur vl2({A2[Vxl],A2[Vyl]});

    
    
    Vecteur A3(12);
    acceleration(A3,ra + 0.5*ra2,rt + 0.5*rt2,rl + 0.5*rl2, va + 0.5*va2, vt + 0.5*vt2, vl + 0.5*vl2);
    A3 = dt*A3;
    
    
    //cout <<"A3=" <<  A3[Xt] << " " << A3[Yt] << " " << A3[Xa] << " " << A3[Ya] << " " << A3[Xl] << " " << A3[Yl] << " " << A3[Vxt] << " " << A3[Vyt] << " " << A3[Vxa] << " " << A3[Vya] << endl;
    
    
    Vecteur A4(12);

    Vecteur rt3({A3[Xt],A3[Yt]});
    Vecteur ra3({A3[Xa],A3[Ya]});
    Vecteur rl3({A3[Xl],A3[Yl]});
    Vecteur vt3({A3[Vxt],A3[Vyt]});
    Vecteur va3({A3[Vxa],A3[Vya]});
    Vecteur vl3({A3[Vxl],A3[Vyl]});
    
    acceleration(A4,ra + ra3,rt + rt3,rl + rl3, va + va3, vt + vt3, vl + vl3);
    A4 = dt*A4;
    
    
   // cout <<"A4=" <<  A4[Xt] << " " << A4[Yt] << " " << A4[Xa] << " " << A4[Ya] << " " << A4[Xl] << " " << A4[Yl] << " " << A4[Vxt] << " " << A4[Vyt] << " " << A4[Vxa] << " " << A4[Vya] << endl;
    
    Vecteur rt4({A4[Xt],A4[Yt]});
    Vecteur ra4({A4[Xa],A4[Ya]});
    Vecteur rl4({A4[Xl],A4[Yl]});
    Vecteur vt4({A4[Vxt],A4[Vyt]});
    Vecteur va4({A4[Vxa],A4[Vya]});
    Vecteur vl4({A4[Vxl],A4[Vyl]});

    
    rt += (1.0/6)*(rt1 + 2.0*rt2 + 2.0*rt3 + rt4);
    ra += (1.0/6)*(ra1 + 2.0*ra2 + 2.0*ra3 + ra4);
    rl += (1.0/6)*(rl1 + 2.0*rl2 + 2.0*rl3 + rl4);
    vt += (1.0/6)*(vt1 + 2.0*vt2 + 2.0*vt3 + vt4);
    va += (1.0/6)*(va1 + 2.0*va2 + 2.0*va3 + va4);
    vl += (1.0/6)*(vl1 + 2.0*vl2 + 2.0*vl3 + vl4);

  
    
}


// Getters, Setters

void PhysEngine::set_dt(const double& d){
    dt = d;
}

double PhysEngine::get_distance(const unsigned int& astre) const{

    if (astre == terre){
        return sqrt(rt[0]*rt[0] + rt[1]*rt[1]);}
    else if (astre == asteroide){
    
        return sqrt(pow(ra[0] - rt[0],2) + pow(ra[1] - rt[1],2));}
    else if (astre == lune){
    
        return sqrt(pow(rl[0] - rt[0],2) + pow(rl[1] - rt[1],2));
        cout << "lune" << endl;
    }else if(astre == asteroide_lune){
    
        return sqrt(pow(ra[0] - rl[0],2) + pow(ra[1] - rl[1],2));
    }
    cout << "erreur" << endl;
    return 1;
}

double PhysEngine::get_dt() const{return dt;}





//constructeur copie










