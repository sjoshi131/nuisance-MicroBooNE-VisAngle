// Copyright 2024 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

/*******************************************************************************
 *    This file is part of NUISANCE.
 *
 *    NUISANCE is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    NUISANCE is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    
 *    You should have received a copy of the GNU General Public License
 *    along with NUISANCE.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#include "MicroBooNE_CC1Mu1p_VisAngle_XSec_1D_nu.h"
#include "MicroBooNE_SignalDef.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TVectorD.h"

//********************************************************************
MicroBooNE_CC1Mu1p_VisAngle_XSec_1D_nu::MicroBooNE_CC1Mu1p_VisAngle_XSec_1D_nu(
    nuiskey samplekey) {
  //********************************************************************
  fSettings = LoadSampleSettings(samplekey);

  std::string name = fSettings.GetS("name");
  std::string objSuffix;

  // work out which sample you need, and set axii
  if (!name.compare("MicroBooNE_CC1Mu1p_VisAngle_XSec_1DThetaVis_nu")) {
    fDist = kThetaVis;
    objSuffix = "ThetaVis";
    fSettings.SetXTitle("#theta_{vis} (deg)");
    fSettings.SetYTitle("d#sigma/d#theta_{vis} (cm^{2}/deg/^{40}Ar)");
  } else if (!name.compare("MicroBooNE_CC1Mu1p_VisAngle_XSec_1DCosThetaVis_nu")) {
    fDist = kCosThetaVis;
    objSuffix = "CosThetaVis";
    fSettings.SetXTitle("cos#theta_{vis}");
    fSettings.SetYTitle("d#sigma/dcos#theta_{vis} (cm^{2}/^{40}Ar)");
  } else if (!name.compare("MicroBooNE_CC1Mu1p_VisAngle_XSec_1DPMiss_nu")) {
    fDist = kPMiss;
    objSuffix = "PMiss";
    fSettings.SetXTitle("p_{miss} (GeV/c)");
    fSettings.SetYTitle("d#sigma/dp_{miss} (cm^{2}/(GeV/c)/^{40}Ar)");
  } else if (!name.compare("MicroBooNE_CC1Mu1p_VisAngle_XSec_1DMuonCosTheta_nu")) {
    fDist = kMuonCosTheta;
    objSuffix = "MuonCosTheta";
    fSettings.SetXTitle("cos#theta_{#mu}");
    fSettings.SetYTitle("d#sigma/dcos#theta_{#mu} (cm^{2}/^{40}Ar)");
  } else {
    NUIS_ABORT(
        "MicroBooNE_CC1Mu1p_VisAngle_XSec_1D_nu: Didn't get a valid name: " << name);
  }

  // Sample overview ---------------------------------------------------
  std::string descrip = name + " sample.\n"
                               "Target: Ar\n"
                               "Flux: BNB FHC numu\n"
                               "Signal: CC1Mu1p (VisAngle Analysis)\n";
  fSettings.SetDescription(descrip);
  fSettings.SetTitle(name);
  fSettings.SetAllowedTypes("FULL,DIAG/FREE,SHAPE,FIX/SYSTCOV/STATCOV",
                            "FIX/FULL");
  fSettings.SetEnuRange(0.0, 6.8);
  fSettings.DefineAllowedTargets("Ar");
  fSettings.DefineAllowedSpecies("numu");
  FinaliseSampleSettings();

  // Load data ---------------------------------------------------------
  // Path for the new data file
  std::string inputFile = FitPar::GetDataBase() +
              "/MicroBooNE/CC1Mu1p/New_All_XSecs_Combined.root";
  
  // Naming convention from newdata_dump.txt: 
  // XSecReco{objSuffix}Plot
  // UnfCov{objSuffix}Plot
  // Ac{objSuffix}Plot
  
  SetDataFromRootFile(inputFile, "RecoFullUnc" + objSuffix + "Plot");
  ScaleData(1E-38);

  fScaleFactor = GetEventHistogram()->Integral("width") / fNEvents * 1.0E-38 /
                 TotalIntegratedFlux() * 40;

  SetCovarFromRootFile(inputFile, "UnfCov" + objSuffix + "Plot");

  // The covariance matrix in the ROOT file is stored in integrated form (not divided by bin widths).
  // NUISANCE expects the covariance matrix to be in differential form.
  // We divide each element (i, j) by bin_width(i) * bin_width(j).
  if (fFullCovar) {
    for (int i = 0; i < fFullCovar->GetNrows(); i++) {
      double w_i = fDataHist->GetBinWidth(i + 1);
      for (int j = 0; j < fFullCovar->GetNcols(); j++) {
        double w_j = fDataHist->GetBinWidth(j + 1);
        (*fFullCovar)(i, j) /= (w_i * w_j);
      }
    }
    // Reset the inverse and decomposition matrices so they are re-calculated in FinaliseMeasurement()
    if (covar) {
      delete covar;
      covar = NULL;
    }
    if (fDecomp) {
      delete fDecomp;
      fDecomp = NULL;
    }
  }

  TFile *inputRootFile = TFile::Open(inputFile.c_str());
  assert(inputRootFile && inputRootFile->IsOpen());
  TH2D *hsmear =
      (TH2D *)inputRootFile->Get(("Ac" + objSuffix + "Plot").c_str());
  assert(hsmear);

  int nrows = hsmear->GetNbinsX();
  int ncols = hsmear->GetNbinsY();
  fSmearingMatrix = new TMatrixD(nrows, ncols);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      fSmearingMatrix->operator()(i, j) = hsmear->GetBinContent(i + 1, j + 1);
    }
  }

  // Final setup ------------------------------------------------------
  FinaliseMeasurement();
}

bool MicroBooNE_CC1Mu1p_VisAngle_XSec_1D_nu::isSignal(FitEvent *event) {
  // Reuse existing CC1Mu1p signal definition
  return SignalDef::MicroBooNE::isCC1Mu1p(event, EnuMin, EnuMax);
}

void MicroBooNE_CC1Mu1p_VisAngle_XSec_1D_nu::FillEventVariables(FitEvent *event) {

  if (!isSignal(event)) {
    fXVar = -999;
    return;
  }

  // Identify leading proton
  int signal_proton_index = SignalDef::MicroBooNE::getCC1Mu1pProtonIndex(event);
  FitParticle *signal_proton = event->GetParticle(signal_proton_index);

  TVector3 vpmu = event->GetHMFSParticle(13)->P3();
  TVector3 vp = signal_proton->P3();

  if (fDist == kThetaVis) {
    TVector3 vSum = vpmu + vp;
    fXVar = vSum.Theta() * 180. / TMath::Pi();
  }
  else if (fDist == kCosThetaVis) {
    TVector3 vSum = vpmu + vp;
    fXVar = vSum.CosTheta();
  }
  else if (fDist == kPMiss) {
    // Under the VisAngle analysis (arXiv:2504.17758), the missing momentum p_miss is defined as:
    // p_miss = E_reco - b
    // where b = |p_mu + p_p| (magnitude of visible momentum)
    // and E_reco = E_mu + K_p + removal_energy (30.9 MeV for Argon)
    TVector3 vSum = vpmu + vp;
    double b = vSum.Mag();
    double E_mu = event->GetHMFSParticle(13)->E();
    double K_p = signal_proton->KE();
    double B = 30.9; // Nucleon removal energy for Argon in MeV
    double E_reco = E_mu + K_p + B;
    fXVar = (E_reco - b) / 1000.0; // Convert to GeV/c
  }
  else if (fDist == kMuonCosTheta) {
    fXVar = vpmu.CosTheta();
  }
}

void MicroBooNE_CC1Mu1p_VisAngle_XSec_1D_nu::ConvertEventRates() {
  int n = fMCHist->GetNbinsX();

  TVectorD v(n), e(n);
  for (int i = 0; i < n; i++) {
    v(i) = fMCHist->GetBinContent(i + 1);
    e(i) = std::pow(fMCHist->GetBinError(i + 1), 2);
  }

  TVectorD vs = (*fSmearingMatrix) * v;
  TVectorD es = (*fSmearingMatrix) * e;

  for (int i = 0; i < n; i++) {
    fMCHist->SetBinContent(i + 1, vs(i));
    fMCHist->SetBinError(i + 1, std::sqrt(es(i)));
  }

  Measurement1D::ConvertEventRates();
}
