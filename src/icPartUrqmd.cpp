#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cfloat>

#include "eos.h"
#include "eoChiral.h"
#include "rmn.h"
#include "fld.h"
#include "icPartUrqmd.h"

using namespace std;

IcPartUrqmd::IcPartUrqmd(Fluid* f, char* filename, double _Rgt, double _Rgz,
                         double _tau0) {
 nx = f->getNX();
 ny = f->getNY();
 nz = f->getNZ();
 dx = f->getDx();
 dy = f->getDy();
 dz = f->getDz();
 xmin = f->getX(0);
 xmax = f->getX(nx - 1);
 ymin = f->getY(0);
 ymax = f->getY(ny - 1);
 zmin = f->getZ(0);
 zmax = f->getZ(nz - 1);

 tau0 = _tau0;
 Rgx = _Rgt;
 Rgy = _Rgt;
 Rgz = _Rgz;
 nsmoothx = (int)(6.0 * Rgx / dx);  // smoothly distribute to +- this many cells
 nsmoothy = (int)(6.0 * Rgy / dy);
 nsmoothz = (int)(1.5 * Rgz / dz);

 T00 = new double**[nx];
 T0x = new double**[nx];
 T0y = new double**[nx];
 T0z = new double**[nx];
 QB = new double**[nx];
 QE = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  T00[ix] = new double*[ny];
  T0x[ix] = new double*[ny];
  T0y[ix] = new double*[ny];
  T0z[ix] = new double*[ny];
  QB[ix] = new double*[ny];
  QE[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   T00[ix][iy] = new double[nz];
   T0x[ix][iy] = new double[nz];
   T0y[ix][iy] = new double[nz];
   T0z[ix][iy] = new double[nz];
   QB[ix][iy] = new double[nz];
   QE[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    T00[ix][iy][iz] = 0.0;
    T0x[ix][iy][iz] = 0.0;
    T0y[ix][iy][iz] = 0.0;
    T0z[ix][iy][iz] = 0.0;
    QB[ix][iy][iz] = 0.0;
    QE[ix][iy][iz] = 0.0;
   }
  }
 }
#ifdef TSHIFT
 tau0 += tshift;
#endif
 // ---- read the events
 nevents = 0;
 ifstream fin(filename);
 if (!fin.good()) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }
 int np = 0;  // particle counter
 string line;
 istringstream instream;
 while (!fin.eof()) {
  getline(fin, line);
  instream.str(line);
  instream.seekg(0);
  instream.clear();
  instream >> Tau[np] >> X[np] >> Y[np] >> Eta[np] >> Mt[np] >> Px[np] >>
      Py[np] >> Rap[np] >> Id[np] >> Charge[np];
#ifdef TSHIFT
  Eta[np] = TMath::ATanH(Tau[np] * sinh(Eta[np]) /
                         (Tau[np] * cosh(Eta[np]) + tshift));
  Tau[np] += tshift;
#endif
  if (!instream.fail())
   np++;
  else if (np > 0) {
   // cout<<"readF14:instream: failure reading data\n" ;
   // cout<<"stream = "<<instream.str()<<endl ;
   if (nevents % 100 == 0) {
    cout << "event = " << nevents << "  np = " << np << "\r";
    cout << flush;
   }
   makeSmoothTable(np);
   np = 0;
   nevents++;
   // if(nevents>10000) return ;
  }
  if (np > NP - 1) cout << "ERROR: increase NP constant\n";
 }
 if (nevents > 1)
  cout << "++ Warning: loaded " << nevents << "  initial UrQMD events\n";
}

IcPartUrqmd::~IcPartUrqmd() {
 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   delete[] T00[ix][iy];
   delete[] T0x[ix][iy];
   delete[] T0y[ix][iy];
   delete[] T0z[ix][iy];
   delete[] QB[ix][iy];
   delete[] QE[ix][iy];
  }
  delete[] T00[ix];
  delete[] T0x[ix];
  delete[] T0y[ix];
  delete[] T0z[ix];
  delete[] QB[ix];
  delete[] QE[ix];
 }
 delete[] T00;
 delete[] T0x;
 delete[] T0y;
 delete[] T0z;
 delete[] QB;
 delete[] QE;
}

void IcPartUrqmd::makeSmoothTable(int npart) {
 // To add Baryon contributions to net-baryon number
 // PDG IDs corresponding to Baryons included in SMASH-1.7
 int PDG_Codes_Baryons[142] = {2112, 2212, 12112, 12212, 1214, 2124, 22112,
                           22212, 32112, 32212, 2116, 2216, 12116, 12216, 21214,
                           22124, 42112, 42212, 31214, 32124, 9902114, 9902214,
                           9952112, 9952212, 9962112, 9962212, 9912114, 9912214,
                           9902118, 9902218, 9922116, 9922216, 9922114, 9922214,
                           9972112, 9972212, 9932114, 9932214, 1218, 2128,
                           19922119, 19922219, 19932119, 19932219,
                           1114, 2114, 2214, 2224, 31114, 32114, 32214,
                           32224,1112, 1212, 2122, 2222, 11114, 12114,
                           12214, 12224, 11112, 11212, 12122, 12222, 1116,
                           1216, 2126, 2226, 21112, 21212, 22122, 22222,
                           21114, 22114, 22214, 22224, 11116, 11216,
                           12126, 12226, 1118, 2118, 2218, 2228,
                           3122, 13122, 3124, 23122, 23122, 13124,
                           43122, 53122, 3126, 13126, 23124, 3128,
                           23126, 19903129,
                           3112, 3212, 3222, 3114, 3214, 3224, 13112,
                           13212, 13222, 13114, 13214, 13224, 23112,
                           23212, 23222, 3116, 3216, 3226, 13116, 13216,
                           13226, 23114, 23214, 23224, 3118, 3218, 3228,
                           9903118, 9903218, 9903228,
                           3312, 3322, 3314, 3324, 203312, 203322, 13314,
                           13324, 103316, 103326, 203316, 203326,
                           3334, 203338};


 for (int ip = 0; ip < npart; ip++) {  // particle loop
  int ixc = (int)round((X[ip] - xmin) / dx);
  int iyc = (int)round((Y[ip] - ymin) / dy);
  int izc = (int)round((Eta[ip] - zmin) / dz);
  // finding the norm
  double norm_gauss = 0.0;
  const double gammaz = cosh(Rap[ip] - Eta[ip]);
  for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++)
   for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++)
    for (int iz = izc - nsmoothz; iz < izc + nsmoothz + 1; iz++)
     if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {
      const double xdiff = X[ip] - (xmin + ix * dx);
      const double ydiff = Y[ip] - (ymin + iy * dy);
      const double zdiff = Eta[ip] - (zmin + iz * dz);
      if (fabs(Rgz) < 1e-5)
       norm_gauss +=
           exp(-xdiff * xdiff / Rgx / Rgx - ydiff * ydiff / Rgy / Rgy);
      else
       norm_gauss +=
           exp(-xdiff * xdiff / Rgx / Rgx - ydiff * ydiff / Rgy / Rgy -
               zdiff * zdiff / Rgz / Rgz * gammaz * gammaz * tau0 *
                   tau0);  // this term may become big, >800 !
     }

  for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++)
   for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++)
    for (int iz = izc - nsmoothz; iz < izc + nsmoothz + 1; iz++)
     if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {
      const double xdiff = X[ip] - (xmin + ix * dx);
      const double ydiff = Y[ip] - (ymin + iy * dy);
      const double zdiff = Eta[ip] - (zmin + iz * dz);
      double weight;
      if (fabs(Rgz) < 1e-5)  // jun19
       weight = 1.0 / norm_gauss *
                exp(-xdiff * xdiff / Rgx / Rgx - ydiff * ydiff / Rgy / Rgy);
      else
       weight = 1.0 / norm_gauss *
                exp(-xdiff * xdiff / Rgx / Rgx - ydiff * ydiff / Rgy / Rgy -
                    zdiff * zdiff / Rgz / Rgz * gammaz * gammaz * tau0 *
                        tau0);  // -zdiff*zdiff*...*tau*tau (Jun 16)
      if (weight != weight || fabs(weight) > DBL_MAX) {
       weight = 0.0;
      }
      T00[ix][iy][iz] += weight * Mt[ip] * cosh(Rap[ip] - Eta[ip] + zdiff);
      T0x[ix][iy][iz] += weight * Px[ip];
      T0y[ix][iy][iz] += weight * Py[ip];
      T0z[ix][iy][iz] += weight * Mt[ip] * sinh(Rap[ip] - Eta[ip] + zdiff);
      // if(ix==50 && iy==50 && iz==17){ // debug output
      // foutd<<setw(14)<<Mt[ip]<<setw(14)<<Rap[ip]<<setw(14)<<Eta[ip]-zdiff<<setw(14)<<weight<<endl
      // ;
      //}
      // Add baryons and antibaryons
      for (auto &code : PDG_Codes_Baryons) {
        if (Id[ip] == code) {
          QB[ix][iy][iz] += weight;  // correct expression (only weight)!
        } else if (Id[ip] == (-1 * code)) {
          QB[ix][iy][iz] -= weight;
        }
      }
      // Add electrice charge quantum number
      QE[ix][iy][iz] += Charge[ip] * weight;
      //}
     }
 }  // end particle loop
}

void IcPartUrqmd::setIC(Fluid* f, EoS* eos) {
 double E = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0, Nb = 0.0, S = 0.0;
 double Q[7], e, p, nb, nq, ns, vx, vy, vz;
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    Q[T_] = T00[ix][iy][iz] / nevents / dx / dy / dz / tau0;  // /tau for Milne
    Q[X_] = T0x[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q[Y_] = T0y[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q[Z_] = T0z[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q[NB_] = QB[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q[NQ_] = QE[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q[NS_] = 0.0;
    if (ix == nx / 2 && iy == ny / 2 && iz == nz / 2)
     cout << "IcUrqmd, center: " << xmin + ix * dx << "  " << zmin + iz * dz
          << "  " << Q[T_] << "  " << Q[Z_] << endl;
    transformPV(eos, Q, e, p, nb, nq, ns, vx, vy, vz);
    if (e < 1e-7 || fabs(f->getX(ix)) > 10. || fabs(f->getY(iy)) > 10. ||
        fabs(f->getZ(iz)) > 5.) {
     e = nb = nq = 0.0;
     vx = vy = vz = 0.0;
    }
    Cell* c = f->getCell(ix, iy, iz);
    c->setPrimVar(eos, tau0, e, nb, nq, ns, vx, vy, vz);
    if (e > 0.) c->setAllM(1.);
    const double gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
    double u[4] = {gamma, gamma * vx, gamma * vy, gamma * vz};
    double eta = zmin + iz * dz;
    E += tau0 * ((e + p) * u[0] * (u[0] * cosh(eta) + u[3] * sinh(eta)) -
                 p * cosh(eta)) *
         dx * dy * dz;
    // if(zmin+iz*dz>0.)
    Pz += tau0 * ((e + p) * u[0] * (u[0] * sinh(eta) + u[3] * cosh(eta)) -
                  p * sinh(eta)) *
          dx * dy * dz;
    Px += tau0 * (e + p) * u[1] * u[0] * dx * dy * dz;
    Py += tau0 * (e + p) * u[2] * u[0] * dx * dy * dz;
    Nb += tau0 * nb * u[0] * dx * dy * dz;
    S += tau0 * eos->s(e, nb, nq, ns) * u[0] * dx * dy * dz;
   }
 cout << "hydrodynamic E = " << E << "  Pz = " << Pz << "  Nbar = " << Nb
      << endl
      << "  Px = " << Px << "  Py = " << Py << endl;
 cout << "initial_entropy S_ini = " << S << endl;
}
