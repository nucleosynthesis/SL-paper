double integrateExp(double s, double a, double b){

    // return integral_a^b  e^{sx} dx 

    double B = (1./s)*TMath::Exp(s*b);
    double A = (1./s)*TMath::Exp(s*a);
    return B-A;
}


TH1F *vectorToHisto(std::string name, double *v, int n, int nmin=-1, int nmax=-1){

   TH1F *h; 
   if (nmin>=0 && nmax >=0){
     h = new TH1F(name.c_str(),"",nmax-nmin,0,nmax-nmin);
     int count=1;
     for (int i=nmin;i<nmax;i++){
       h->SetBinContent(count,v[i-1]);
       h->GetXaxis()->SetBinLabel(count,Form("%d",i));
       h->GetXaxis()->SetBinLabel(count,Form("%d",i));
       count++;
     }


   } else{ 
     h = new TH1F(name.c_str(),"",n,0,n);

     for (int i=1;i<=n;i++){
       h->SetBinContent(i,v[i-1]);
       h->GetXaxis()->SetBinLabel(i,Form("%d",i));
       h->GetXaxis()->SetBinLabel(i,Form("%d",i));
     }
   }
   h->GetXaxis()->SetTitle("Signal Region");
   h->GetYaxis()->SetTitle("Number of events");

   h->SetMarkerSize(0.8);
   h->SetMarkerStyle(20);
   h->SetLineWidth(1);
   h->SetLineColor(1);

   return h;

}
TH1F *vectorToHisto(std::string name, double *v, double *e, int n){
   TH1F * h = vectorToHisto(name,v,n);
   for (int i=1;i<=n;i++){
       std::cout << i << " " << v[i-1] << " " << e[i-1] << std::endl;
     h->SetBinError(i,e[i-1]);
   }
   return h;
}

void makeSimpleLikelihoodToy(){

    gStyle->SetOptStat(0);

    
    RooRealVar JES("JES","",0);

    double N1 = 2000;
    double N2 = 500;
    double N3 = 100;
    
    double S1 = 200;
    double S2 = 70;
    double S3 = 50;

    // say JES has 10%,-20%,20%  uncertainty on the normalisations lnN ~ 1.1, 1./1.2, 1.2!
    double lnN1 = 1.05;
    double lnN2 = 1./1.05;
    double lnN3 = 1.1;

    double lnNIsr1 = 1/1.05;
    double lnNIsr2 = 1.05;
    double lnNIsr3 = 1.1;

    // background slopes
    double s1 = -0.05;
    double s2 = -0.03;
    double s3 = -0.01;

    // signal slopes
    double si1 = -0.02;
    double si2 = -0.01;
    double si3 = -0.004;

    double isr1 = 0.04;
    double isr2 = 0.06;
    double isr3 = 0.05;

    double es1 = 0.001;
    double es2 = -0.002;
    double es3 = 0.0008;
    
    // bounds will be between 100 and 1000, 30 bins per cat (30 GeV)
    double XMIN = 10;
    double XMAX = 100;
    double DX   = 3;

    // Now lets make the nominal, up and down for some imaginary JES scale!
   
    double *nominal_bkg = new double[90];
    double *JES_up_bkg  = new double[90];
    double *JES_dn_bkg  = new double[90];
    double *ISR_up_bkg  = new double[90];
    double *ISR_dn_bkg  = new double[90];

    double *data          = new double[90];
    double *mcUnc          = new double[90];
    double *signal        = new double[90];



    TRandom3 *r = new TRandom3();
    TRandom3 *rMC = new TRandom3();

    double mcUncMult = 20.;

    int j=0;

    for (int i=0; i<90; i++){
      
      if (i==30 || i ==60) j = 0;

      double a = XMIN+j*DX;
      double b = a+DX;

      double nominal, jesu,jesd,isru,isrd, s, d;
      
      if (i<30) { 
      	nominal = N1*integrateExp(s1,a,b)/integrateExp(s1,XMIN,XMAX);
      	jesu    = lnN1*N1*integrateExp(s1+es1,a,b)/integrateExp(s1+es1,XMIN,XMAX);
      	jesd    = (1./lnN1)*N1*integrateExp(s1-es1,a,b)/integrateExp(s1-es1,XMIN,XMAX);
      	isru    = lnNIsr1*N1*integrateExp(s1+isr1*isr1,a,b)/integrateExp(s1+isr1*isr1,XMIN,XMAX);
      	isrd    = (1./lnNIsr1)*N1*integrateExp(s1-isr1*isr1,a,b)/integrateExp(s1-isr1*isr1,XMIN,XMAX);
	s       = S1*integrateExp(si1,a,b)/integrateExp(si1,XMIN,XMAX);

      }
      else if (i>=30 && i < 60){
      	nominal = N2*integrateExp(s2,a,b)/integrateExp(s2,XMIN,XMAX);
      	jesu    = lnN2*N2*integrateExp(s2+es2,a,b)/integrateExp(s2+es2,XMIN,XMAX);
      	jesd    = (1./lnN2)*N2*integrateExp(s2-es2,a,b)/integrateExp(s2-es2,XMIN,XMAX);
      	isru    = lnNIsr2*N2*integrateExp(s2+isr2*isr2,a,b)/integrateExp(s2+isr2*isr2,XMIN,XMAX);
      	isrd    = (1./lnNIsr2)*N2*integrateExp(s2-isr2*isr2,a,b)/integrateExp(s2-isr2*isr2,XMIN,XMAX);
	s       = S2*integrateExp(si2,a,b)/integrateExp(si2,XMIN,XMAX);

      }
      else{
      	nominal = N3*integrateExp(s3,a,b)/integrateExp(s3,XMIN,XMAX);
      	jesu    = lnN3*N3*integrateExp(s3+es3,a,b)/integrateExp(s3+es3,XMIN,XMAX);
      	jesd    = (1./lnN3)*N3*integrateExp(s3-es3,a,b)/integrateExp(s3-es3,XMIN,XMAX);
      	isru    = lnNIsr3*N3*integrateExp(s3+isr3*isr3,a,b)/integrateExp(s3+isr3*isr3,XMIN,XMAX);
      	isrd    = (1./lnNIsr3)*N3*integrateExp(s3-isr3*isr3,a,b)/integrateExp(s3-isr3*isr3,XMIN,XMAX);
	s       = S3*integrateExp(si3,a,b)/integrateExp(si3,XMIN,XMAX);
      }
    //MC Unc
      double nominal_fluc = (double) rMC->Poisson(nominal*mcUncMult)/mcUncMult;
      mcUnc[i] = nominal*TMath::Sqrt(nominal*mcUncMult)/(nominal*mcUncMult);
      double msc = nominal_fluc/nominal;

      if (j<15) s*=double(j)/15;
      d = (double)r->Poisson(nominal);
      nominal_bkg[i]=nominal*msc;
      JES_up_bkg[i]=jesu*msc;
      JES_dn_bkg[i]=jesd*msc;
      ISR_up_bkg[i]=isru*msc;
      ISR_dn_bkg[i]=isrd*msc;
      data[i]=d;
      signal[i]=s;

      j++;
    }

    TCanvas *can = new TCanvas("canvas","",900,400);
    can->SetLogy();

    TH1F *h_signal = (TH1F*)vectorToHisto("signal",signal,90); h_signal->SetLineColor(2);
    TH1F *h_data   = (TH1F*)vectorToHisto("data",data,90);
    TH1F *h_bkg    = (TH1F*)vectorToHisto("nominal_bkg",nominal_bkg,mcUnc,90); h_bkg->SetLineColor(4); h_bkg->SetFillColor(kAzure+2); h_bkg->SetMarkerSize(0);
    
    TH1F *h_bkg_JES_up_bkg    = (TH1F*)vectorToHisto("nominal_bkg_JESUp",JES_up_bkg,90); h_bkg_JES_up_bkg->SetLineColor(1);
    TH1F *h_bkg_JES_dn_bkg    = (TH1F*)vectorToHisto("nominal_bkg_JESDown",JES_dn_bkg,90); h_bkg_JES_dn_bkg->SetLineColor(1);
    TH1F *h_bkg_ISR_up_bkg    = (TH1F*)vectorToHisto("nominal_bkg_ISRUp",ISR_up_bkg,90); h_bkg_ISR_up_bkg->SetLineColor(1),h_bkg_ISR_up_bkg->SetLineStyle(2);
    TH1F *h_bkg_ISR_dn_bkg    = (TH1F*)vectorToHisto("nominal_bkg_ISRDown",ISR_dn_bkg,90); h_bkg_ISR_dn_bkg->SetLineColor(1),h_bkg_ISR_dn_bkg->SetLineStyle(2) ;
    
    TH1F *h_bkg_MCUNC         = (TH1F*)vectorToHisto("mc_unc",mcUnc,90); h_bkg_MCUNC->SetLineColor(1),h_bkg_MCUNC->SetLineStyle(2) ;

    for (int i=0;i<h_bkg->GetNbinsX();i++){
	h_bkg->SetBinError(i+1,h_bkg_MCUNC->GetBinContent(i+1));
    }
    
    h_data->GetXaxis()->LabelsOption("v");
    //h_data->SetMinimum(0.001);
    h_data->Draw("P");
    h_bkg->Draw("E2L");
    h_bkg->Draw("Lsame");
    h_bkg_JES_dn_bkg->Draw("histsame");
    h_bkg_JES_up_bkg->Draw("histsame");
    h_bkg_ISR_dn_bkg->Draw("histsame");
    h_bkg_ISR_up_bkg->Draw("histsame");
    h_signal->Draw("histsame");
    h_data->Draw("Psame");
   
    TLegend *leg = new TLegend(0.65,0.5,0.89,0.89);
    leg->SetBorderSize(0);
    leg->AddEntry(h_data,"Observed data","P");
    leg->AddEntry(h_bkg,"Nominal background (#pm stat unc.)","LFE");
    leg->AddEntry(h_bkg_JES_up_bkg,"Energy scale up/down","L");
    leg->AddEntry(h_bkg_ISR_up_bkg,"Theory uncertainty up/down","L");
    leg->AddEntry(h_signal,"BSM signal","L");
    leg->Draw();

    TLine l1(30,h_data->GetMinimum(),30,0.75*h_data->GetMaximum());l1.SetLineColor(1);l1.SetLineStyle(3);
    TLine l2(60,h_data->GetMinimum(),60,0.75*h_data->GetMaximum());l2.SetLineColor(1);l2.SetLineStyle(3);
    l1.Draw();
    l2.Draw();
    can->SaveAs("t.pdf");

    // now lets make a ROOT file with the histograms 
    TFile *fout = new TFile("histos.root","RECREATE");
    fout->cd();
    h_data->Write();
    h_bkg->Write();
    h_signal->Write();
    h_bkg_JES_dn_bkg->Write();
    h_bkg_JES_up_bkg->Write();
    h_bkg_ISR_dn_bkg->Write();
    h_bkg_ISR_up_bkg->Write();
    h_bkg_MCUNC->Write();

// DO the same in 3 categories 
    for (int c=1;c<4;c++){
      int min = (c-1)*30 + 1;
      int max = (c-1)*30+31;
      h_signal = (TH1F*)vectorToHisto(Form("C%d_signal",c),signal,90,min,max); h_signal->SetLineColor(2);
      h_data   = (TH1F*)vectorToHisto(Form("C%d_data",c),data,90,min,max);
      h_bkg    = (TH1F*)vectorToHisto(Form("C%d_nominal_bkg",c),nominal_bkg,90,min,max); h_bkg->SetLineColor(4); h_bkg->SetFillColor(kAzure+2); h_bkg->SetMarkerSize(0);
      
      h_bkg_JES_up_bkg    = (TH1F*)vectorToHisto(Form("C%d_nominal_bkg_JESUp",c),JES_up_bkg,90,min,max); h_bkg_JES_up_bkg->SetLineColor(1);
      h_bkg_JES_dn_bkg    = (TH1F*)vectorToHisto(Form("C%d_nominal_bkg_JESDown",c),JES_dn_bkg,90,min,max); h_bkg_JES_dn_bkg->SetLineColor(1);
      h_bkg_ISR_up_bkg    = (TH1F*)vectorToHisto(Form("C%d_nominal_bkg_ISRUp",c),ISR_up_bkg,90,min,max); h_bkg_ISR_up_bkg->SetLineColor(1),h_bkg_ISR_up_bkg->SetLineStyle(2);
      h_bkg_ISR_dn_bkg    = (TH1F*)vectorToHisto(Form("C%d_nominal_bkg_ISRDown",c),ISR_dn_bkg,90,min,max); h_bkg_ISR_dn_bkg->SetLineColor(1),h_bkg_ISR_dn_bkg->SetLineStyle(2);

      TH1F *h_bkg_MCUNC         = (TH1F*)vectorToHisto(Form("mc_unc_c%d",c),mcUnc,90,min,max); h_bkg_MCUNC->SetLineColor(1),h_bkg_MCUNC->SetLineStyle(2) ;
      for (int i=0;i<h_bkg->GetNbinsX();i++){
         h_bkg->SetBinError(i+1,h_bkg_MCUNC->GetBinContent(i+1));
      }
      h_data->Write();
      h_bkg->Write();
      h_signal->Write();
      h_bkg_JES_dn_bkg->Write();
      h_bkg_JES_up_bkg->Write();
      h_bkg_ISR_dn_bkg->Write();
      h_bkg_ISR_up_bkg->Write();

    }




    fout->Close();


}
