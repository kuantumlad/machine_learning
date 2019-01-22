#include <iostream>
#include <TMarker.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <vector>
#include <map>
#include <TF1.h>
#include <TCanvas.h>
#include <cstdlib>
#include <ctime>

void PrintResults( std::vector< TF1* > results, int best_grid_point, double best_chi2 ){

  std::cout << " >> Optimal grid point " << best_grid_point << " chi2 " << best_chi2 << std::endl;
  std::cout << " >> Best Parameters: ";
  
  int f_num_parameters = results[best_grid_point]->GetNpar();
  for( int i=0; i < f_num_parameters; i++ ){
    std::cout << " " << results[best_grid_point]->GetParameter(i);
  }
  std::cout << "\n";
  
  double best_mean1 = results[best_grid_point]->GetParameter(1);
  double best_sig1 = results[best_grid_point]->GetParameter(2);

  double best_mean2 = results[best_grid_point]->GetParameter(4);
  double best_sig2 = results[best_grid_point]->GetParameter(5);
  
  double best_mean3 = results[best_grid_point]->GetParameter(7);
  double best_sig3 = results[best_grid_point]->GetParameter(8);
  
  double best_mean4 = results[best_grid_point]->GetParameter(10);
  double best_sig4 = results[best_grid_point]->GetParameter(11);

  std::cout << " Best parameters " << best_mean1 <<  " " << best_sig1 << " " << best_mean2 << " " << best_sig2 << " " << best_mean3 << " " << best_sig3 << " " << best_mean4 << " " << best_sig4 << std::endl;
}

void DrawBestResults(TF1 *f_model, std::vector<TF1*> results, int best_grid_point ){
  
  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->Divide(1,1);
  c1->cd(1);
  c1->SetTitle("Gaussian Parameter Grid Search");
  f_model->SetLineColor(kRed);
  f_model->Draw();

  for( int i = 0; i < results.size(); i++ ){
    if( i == best_grid_point ) {
      results[i]->SetLineColor(kBlue);
      results[i]->SetLineStyle(10);
      results[i]->Draw("same");
    }
  }
}

std::map<int, double> GetBestGridPoint( std::map<int, double> grid_chi ){
  
  std::map<int, double> temp_out;
  int temp_grid_point = -1;
  double temp_chi2 = 1000000;
  for( std::map<int, double>::iterator it=grid_chi.begin(); it != grid_chi.end(); ++it ){
    int counter=it->first;
    double chi2 = it->second;
    if( chi2 < temp_chi2 ){
      temp_chi2 = chi2;
      temp_grid_point = counter;
    }
  }
  temp_out[temp_grid_point]=temp_chi2;
  return temp_out;
}



int grid_search_functions(){

  std::cout << " starting example grid search program " << std::endl;
  std::string grid_search_type = "random";

  //define model range and parameters
  double min_x = 0.0;
  double max_x = 10.0;
  TF1 *f_model = new TF1("f_model","gaus(0)*gaus(3)*gaus(6)*gaus(9)",min_x, max_x);
  
  f_model->SetParameter(0,1.0);
  f_model->SetParameter(1,5.0);
  f_model->SetParameter(2,10.0);

  f_model->SetParameter(3,1.0);
  f_model->SetParameter(4,3.0);
  f_model->SetParameter(5,1.0);

  f_model->SetParameter(6,1.0);
  f_model->SetParameter(7,8.0);
  f_model->SetParameter(8,3.0);

  f_model->SetParameter(9,1.0);
  f_model->SetParameter(10,2.0);
  f_model->SetParameter(11,6.0);

  
  //set phase space limits on independent variables
  double model_x1_max = 6.5;
  double model_x1_min = 4.0;

  double model_x2_max = 4.0;
  double model_x2_min = 1.0;

  double model_x3_max = 9.5;
  double model_x3_min = 7.5;

  double model_x4_max = 2.5;
  double model_x4_min = 1.0;

  double resX1 = 0.2;
  double resX2 = 0.15;
  double resX3 = 0.2;
  double resX4 = 0.1;

  int x1_dim = (model_x1_max - model_x1_min)/resX1;
  int x2_dim = (model_x2_max - model_x2_min)/resX2;
  int x3_dim = (model_x3_max - model_x3_min)/resX3;
  int x4_dim = (model_x4_max - model_x4_min)/resX4;

  std::cout << " Binning " << std::endl;
  std::cout << " >> " <<x1_dim << " " << x2_dim << " " << x3_dim << " " << x4_dim << std::endl;
  
  std::vector<TF1*> models;
  TF1* f_dim1 = new TF1("f_{dim1}","gaus(0)",0,10);
  TF1* f_dim2 = new TF1("f_{dim2}","gaus(0)",0,10);
  TF1* f_dim3 = new TF1("f_{dim3}","gaus(0)",0,10);
  TF1* f_dim4 = new TF1("f_dim4","gaus(0)",0,10);
  
  f_dim1->SetParameter(0,1.0);
  f_dim1->SetParameter(1,5.0);
  f_dim1->SetParameter(2,10.0);
  models.push_back(f_dim1);

  f_dim2->SetParameter(0,1.0);
  f_dim2->SetParameter(1,3.0);
  f_dim2->SetParameter(2,1.0);
  models.push_back(f_dim2);

  f_dim3->SetParameter(0,1.0);
  f_dim3->SetParameter(1,8.0);
  f_dim3->SetParameter(2,3.0);
  models.push_back(f_dim3);

  f_dim4->SetParameter(0,1.0);
  f_dim4->SetParameter(1,2.0);
  f_dim4->SetParameter(2,6.0);
  models.push_back(f_dim4);

  TCanvas *c0 = new TCanvas("c0","c0",300,900);
  c0->Divide(1,3);
  c0->SetTitle("functions");  
  c0->cd(1);
  f_dim1->GetHistogram()->SetTitle("f_{1}");  
  f_dim1->GetXaxis()->SetTitle("x1");
  f_dim1->Draw();
  c0->cd(2);
  f_dim2->GetHistogram()->SetTitle("f_{2}");  
  f_dim2->GetXaxis()->SetTitle("x2");
  f_dim2->Draw();
  c0->cd(3);
  f_dim3->GetHistogram()->SetTitle("f_{3}");  
  f_dim3->GetXaxis()->SetTitle("x3");
  f_dim3->Draw();
  c0->SaveAs("models.pdf");
  
  std::cout << " >> Gauss1 Paramaters " << " " << f_dim1->GetParameter(1) << " " << f_dim1->GetParameter(2) << std::endl;
  std::cout << " >> Gauss2 Paramaters " << " " << f_dim2->GetParameter(1) << " " << f_dim2->GetParameter(2) << std::endl;
  std::cout << " >> Gauss3 Paramaters " << " " << f_dim3->GetParameter(1) << " " << f_dim3->GetParameter(2) << std::endl;

  std::vector< TF1* > grid_fits;
  std::map<int,double> grid_chi2;
  std::map<int, std::vector<double> > grid_coord;
  int grid_counter=0;
  int grid_counter2=0;
  int best_grid_point = -1;
  double best_chi2 = 10000000.0;

  TGraph2D *g_grids = new TGraph2D();
  
  if( grid_search_type == "grid" ){

    for( int dim1=0; dim1 < x1_dim; dim1++ ){
      double x1 = model_x1_min + dim1*resX1;
      double fdim1_eval = f_dim1->Eval(x1);      

      for( int dim2=0; dim2 < x2_dim; dim2++ ){
	double x2 = model_x2_min + dim2*resX2;
	double fdim2_eval = f_dim2->Eval(x2);
	
	for( int dim3=0; dim3 < x3_dim; dim3++ ){
	  double x3 = model_x3_min + dim3*resX3;
	  double fdim3_eval = f_dim3->Eval(x3);

	  g_grids->SetPoint(grid_counter2,x1,x2,x3);
	  grid_counter2++;

	  double f_guess = fdim1_eval*fdim2_eval*fdim3_eval;
	  double chi2 = pow((f_guess-1.0),2);	    
	  grid_chi2[grid_counter]=chi2;

	  std::vector<double> temp_coord;
	  temp_coord.push_back(x1);
	  temp_coord.push_back(x2);
	  temp_coord.push_back(x3);

	  grid_coord[grid_counter]=temp_coord;
	  temp_coord.clear();
	  grid_counter++;

	}
      }
    }

    std::cout << " Number of iterations " << grid_counter << std::endl;
    std::map<int, double> grid_results = GetBestGridPoint( grid_chi2 );    
    for( std::map<int,double>::iterator it = grid_results.begin(); it != grid_results.end(); ++ it){
      best_chi2 = it->second;
      best_grid_point = it->first;
    }
    
    TCanvas *c2 = new TCanvas("c2","c2",900,900);
    g_grids->SetTitle("Grid");
    g_grids->SetMarkerStyle(20);    
    g_grids->SetMarkerSize(1);
    g_grids->SetMarkerColor(kBlue);
    g_grids->Draw("PO");
    g_grids->SaveAs("graph_grid_dist.pdf");
    
    std::vector<double> best_coord = grid_coord[best_grid_point];
    std::cout << " >> BEST COORDINATES: ";
    for(int i = 0; i < best_coord.size(); i++ ){
      std::cout<< " " << best_coord[i];
    }
    std::cout << "\n";
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Execute the random grid search here for the toy model

  std::vector<double> iteration;
  std::vector<double> iteration_chi2;
  std::vector<TF1*> random_grid_fits;
  std::map<int,double> random_grid_chi;
  srand(time(0));
  int max_iter = 1000;
  if( grid_search_type == "random" ){
    for( int i = 0; i < max_iter; i++ ){
      random_grid_fits.push_back(new TF1(Form("f_rand_grid_%d_",i),"gaus(0)*gaus(3)*gaus(6)*gaus(9)", min_x, max_x) );
      
      float x1 = model_x1_min + (rand()) /((RAND_MAX/(model_x1_max - model_x1_min)));
      float x2 = model_x2_min + (rand()) /((RAND_MAX/(model_x2_max - model_x2_min)));
      float x3 = model_x3_min + (rand()) /((RAND_MAX/(model_x3_max - model_x3_min)));

      g_grids->SetPoint(grid_counter2,x1,x2,x3);
      grid_counter2++;
      
      double fdim1_eval = f_dim1->Eval(x1);
      double fdim2_eval = f_dim2->Eval(x2);
      double fdim3_eval = f_dim3->Eval(x3);

      double f_guess = fdim1_eval * fdim2_eval * fdim3_eval;
      double chi2 = pow( (f_guess-1), 2 );
      random_grid_chi[i]=chi2;
      iteration.push_back(i);
      iteration_chi2.push_back(chi2);

      std::vector<double> temp_coord;
      temp_coord.push_back(x1);
      temp_coord.push_back(x2);
      temp_coord.push_back(x3);

      grid_coord[i]=temp_coord;
      temp_coord.clear();
      //close random for loop
    }

    //print best coordinates that maximize function
    std::cout << " Number of iterations " << max_iter << std::endl;    
    std::map<int, double> grid_results = GetBestGridPoint( random_grid_chi );    
    for( std::map<int,double>::iterator it = grid_results.begin(); it != grid_results.end(); ++ it){
      best_chi2 = it->second;
      best_grid_point = it->first;
    }
    std::vector<double> best_coord = grid_coord[best_grid_point];
    std::cout << " >> BEST COORDINATES: ";
    for(int i = 0; i < best_coord.size(); i++ ){
      std::cout<< " " << best_coord[i];
    }
    std::cout << "\n";

    TCanvas *c2a = new TCanvas("c2a","c2a",900,900);
    g_grids->SetTitle("Grid");   
    g_grids->SetMarkerStyle(20);    
    g_grids->SetMarkerSize(1);
    g_grids->SetMarkerColor(kBlue);
    g_grids->Draw("PO");
    g_grids->SaveAs("graph_grid_dist.pdf");
      
       
    TCanvas *c2 = new TCanvas("c2","c2",900,900);
    c2->Divide(1,1);
    c2->cd(1);
    TGraph *g_chi2 = new TGraph(iteration.size(), &iteration[0],&iteration_chi2[0]);
    g_chi2->SetTitle("#chi^{2} per iteration");
    g_chi2->SetMarkerStyle(20);
    g_chi2->SetMarkerSize(1);
    g_chi2->SetMarkerColor(kBlue);
    g_chi2->GetXaxis()->SetTitle("Iteration");
    g_chi2->GetXaxis()->CenterTitle(); 
    g_chi2->GetYaxis()->SetTitle("#chi^{2}");
    g_chi2->GetYaxis()->CenterTitle();
    g_chi2->Draw("AP");  
    
    
    
  }
  
  
 
  std::cout << " done " << endl;  
  return 0;
}
