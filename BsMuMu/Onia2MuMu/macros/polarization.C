////////// constants

// beam energy in GeV
const double pbeam = 5000.;

// masses
const double Mprot = 0.9382720;
const double Mmu = 0.10566;

// etc
double gPI = TMath::Pi();


double Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot );

TLorentzVector beam1_LAB( 0., 0., pbeam, Ebeam );
TLorentzVector beam2_LAB( 0., 0., -pbeam, Ebeam );


      // in the event loop

      // input variables (to be defined):
      // laboratory 4-momenta of positive and negative muon

      TLorentzVector muplus_LAB;

      TLorentzVector muminus_LAB;

      TLorentzVector qqbar_LAB = muplus_LAB + muminus_LAB;


      // calculation of decay angular parameters

      // boost beams and positive muon into the q-qbar rest frame:

      TVector3 LAB_to_QQBAR = -qqbar_LAB.BoostVector();

      TLorentzVector beam1_QQBAR = beam1_LAB;
      beam1_QQBAR.Boost( LAB_to_QQBAR );

      TLorentzVector beam2_QQBAR = beam2_LAB;
      beam2_QQBAR.Boost( LAB_to_QQBAR );

      TLorentzVector muplus_QQBAR = muplus_LAB;
      muplus_QQBAR.Boost( LAB_to_QQBAR );


      // reference directions in the Jpsi rest frame:

      TVector3 beam1_direction     = beam1_QQBAR.Vect().Unit();
      TVector3 beam2_direction     = beam2_QQBAR.Vect().Unit();
      TVector3 qqbar_direction     = qqbar_LAB.Vect().Unit();
      TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();

      // all polarization frames have the same Y axis = the normal to the plane formed by
      // the directions of the colliding hadrons
      TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();


      /////////////////////////////////////////////////////////////////////
      // CS frame

      TVector3 newZaxis = beam1_beam2_bisect;
      TVector3 newYaxis = Yaxis;
      TVector3 newXaxis = newYaxis.Cross( newZaxis );

      TRotation rotation;
      rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
      rotation.Invert();   // transforms coordinates from the "xyz" system
                           // to the "new" (rotated) system having the polarization axis
                           // as z axis

      TVector3 muplus_QQBAR_rotated(muplus_QQBAR.Vect());

      muplus_QQBAR_rotated.Transform( rotation );

      double costh_CS = muplus_QQBAR_rotated.CosTheta();

      double phi_CS = muplus_QQBAR_rotated.Phi() * 180. / gPI;
      if ( phi_CS < 0. ) phi_CS = 360. + phi_CS;      // phi defined in degrees from 0 to 360


      /////////////////////////////////////////////////////////////////////
      // HELICITY frame

      newZaxis = qqbar_direction;
      newYaxis = Yaxis;
      newXaxis = newYaxis.Cross( newZaxis );

      rotation.SetToIdentity();
      rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
      rotation.Invert();

      muplus_QQBAR_rotated = muplus_QQBAR.Vect();

      muplus_QQBAR_rotated.Transform( rotation );

      double costh_HX = muplus_QQBAR_rotated.CosTheta();

      double phi_HX = muplus_QQBAR_rotated.Phi() * 180. / gPI;
      if ( phi_HX < 0. ) phi_HX = 360. + phi_HX;

