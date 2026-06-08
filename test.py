TFitResultPtr ptr2 = 0;
int statusFit = 1;

start1oPFit = 0;
for (int b = 1; b <= p->GetNbinsX(); b++) {
    if (p->GetBinContent(b) > 0) {
        start1oPFit = p->GetBinCenter(b);
        break;
    }
}

if (useFitP) {
    const float endFracs[5] = {0.9f, 0.8f, 0.7f, 0.6f, 0.5f};
    float peak = p->GetBinCenter(p->GetMaximumBin());

    for (incrFit = 0; incrFit < 5; incrFit++) {

        if (useOld1oPFit)  end1oPFit = endFracs[incrFit] * peak
        else end1oPFit = 0.6 * peak;

        if (end1oPFit <= start1oPFit) { statusFit = 1; continue; }

        ptr2 = p->Fit(&f_p, "QRS", "", start1oPFit, end1oPFit);
        statusFit = ptr2.Get() ? ptr2->Status() : 1;

        if (statusFit == 0) break;  // good fit, we keep it
    }

    if (statusFit == 0) {
        TF1* f_p2 = &f_p;
        ROOT::Math::IntegratorOneDim intOneDim_p(*f_p2, ROOT::Math::IntegrationOneDim::kGAUSS);
        float intFp = intOneDim_p.Integral(start1oPFit, end1oPFit);
        if (intFp <= 0) std::cout << "ERROR > INTEGRAL FIT P IS <= 0.   ITG = " << intFp << std::endl;

        float intP = p->Integral(p->FindBin(start1oPFit), p->FindBin(end1oPFit), "width");
        SFp = (intFp > 0) ? intP / intFp : -1;
        if (SFp < 0) useFitP = false;
        f_p3 = f_p2;
    }
}

if (useFitP && statusFit != 0) {      // Bad fit
    if (saveFits) { OutputHisto->cd(); p->Write(); }
    std::cout << "  Bad fit 1/p, eta = " << eta->GetBinCenter(i) << " " << p->GetName() << "  " << std::to_string(workerID) << std::endl;
    useFitP = false;
} else {                              // Good fit
    if (saveFits) { OutputHisto->cd(); p->Write(); }
}