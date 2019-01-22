#include <iostream>
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
  
  double best_mean1 = results[best_grid_point]->GetParameter(1);
  double best_sig1 = results[best_grid_point]->GetParameter(2);

  double best_mean2 = results[best_grid_point]->GetParameter(4);
  double best_sig2 = results[best_grid_point]->GetParameter(5);
  
  double best_mean3 = results[best_grid_point]->GetParameter(7);
  double best_sig3 = results[best_grid_point]->GetParameter(8);
  
  double best_mean4 = results[best_grid_point]->GetParameter(10);
  double best_sig4 = results[best_grid_point]->GetParameter(11);

  std::cout << " Best parameters " << best_mean1 <<  " " << best_sig1 << " " << best_mean2 << " " << best_sig2 << " " << best_mean3 << " " << best_sig3 << " "  << best_mean4 << " " << best_sig4 << std::endl;
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


int grid_search(){

  std::cout << " starting example grid search program " << std::endl;
  std::string grid_search_type = "grid";

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
    
  //set generated grid parameter limits
  double model_mean1_max = 6.5;
  double model_mean1_min = 4;
  double model_sig1_max = 14.0;
  double model_sig1_min = 9.0;

  double model_mean2_max = 4.0;
  double model_mean2_min = 1.0;
  double model_sig2_max = 2.5;
  double model_sig2_min = 0.5;

  double model_mean3_max = 9.5;
  double model_mean3_min = 7.5;
  double model_sig3_max = 4.0;
  double model_sig3_min = 1.0;

  double model_mean4_max = 2.5;
  double model_mean4_min = 1.0;
  double model_sig4_max = 7.0;
  double model_sig4_min = 5.0;

  //define desired resolution to use (bin size)
  double resmean1 = 0.5;
  double resmean2 = 0.5;
  double resmean3 = 0.5;
  double resmean4 = 0.5;

  double ressig1 = 0.1;
  double ressig2 = 0.25;
  double ressig3 = 0.10;
  double ressig4 = 1.0;

  //calculate the number of grid points for each 
  double grid_mean_num1 = (model_mean1_max - model_mean1_min)/resmean1;
  double grid_mean_num2 = (model_mean2_max - model_mean2_min)/resmean2;
  double grid_mean_num3 = (model_mean3_max - model_mean3_min)/resmean3;
  double grid_mean_num4 = (model_mean4_max - model_mean4_min)/resmean4;

  double grid_sig_num1 = (model_sig1_max - model_sig1_min)/ressig1;
  double grid_sig_num2 = (model_sig2_max - model_sig2_min)/ressig2;
  double grid_sig_num3 = (model_sig3_max - model_sig3_min)/ressig3;
  double grid_sig_num4 = (model_sig4_max - model_sig4_min)/ressig4;

  std::cout << " >> grid points for mean dimension 1  " << grid_mean_num1 << std::endl;
  std::cout << " >> grid points for mean dimension 2  " << grid_mean_num2 << std::endl;
  std::cout << " >> grid points for mean dimension 3  " << grid_mean_num3 << std::endl;
  std::cout << " >> grid points for mean dimension 4  " << grid_mean_num4 << std::endl;

  std::cout << " >> grid points for sig dimension 1  " << grid_sig_num1 << std::endl;
  std::cout << " >> grid points for sig dimension 2  " << grid_sig_num2 << std::endl;
  std::cout << " >> grid points for sig dimension 3  " << grid_sig_num3 << std::endl;
  std::cout << " >> grid points for sig dimension 4  " << grid_sig_num4 << std::endl;

  double res1 = resmean1;
  double res2 = resmean2;
  double res3 = resmean3;
  double res4 = resmean4;

  double grid_num1 = grid_mean_num1;
  double grid_num2 = grid_mean_num2;
  double grid_num3 = grid_mean_num3;
  double grid_num4 = grid_mean_num4;
   
  std::vector< TF1* > grid_fits;
  std::map<int,double> grid_chi;
  int grid_counter=0;
  int grid_counter2=0;
  int best_grid_point = -1;
  double best_chi2 = 10000000.0;

  std::cout << " True Parameters " << f_model->GetParameter(1) << " " << f_model->GetParameter(2) << " " << f_model->GetParameter(4) << " " << f_model->GetParameter(5) << " " << f_model->GetParameter(7) << " " << f_model->GetParameter(8) << " " << f_model->GetParameter(10) << " " << f_model->GetParameter(11) << std::endl;


  //the dim meanN and dim sigN are for finding a random number via grid search or random search
  //the xN parameters are for scanning over the function's independent variable space to find the maximum of the function
  double temp_best_chi2 = 10000;
  if( grid_search_type == "grid" ){

    
    for( int dim1 = 0; dim1 < grid_num1; dim1++){
      double dim_mean1 = model_mean1_min + dim1*res1;
      double dim_sig1 = model_sig1_min + dim1*res1;
      
      for( int dim2 = 0; dim2 < grid_num2; dim2++){
	double dim_mean2 = model_mean2_min + dim2*res2;
	double dim_sig2 = model_sig2_min + dim2*res2;
	      
	for( int dim3 = 0; dim3 < grid_num3; dim3++){
	  double dim_mean3 = model_mean3_min + dim3*res3;
	  double dim_sig3 = model_sig3_min + dim3*res3;
	  	  
	  g_grids->SetPoint(grid_counter2,dim_mean1, dim_mean2, dim_mean3);
	  grid_counter2++;
	  
	  for( int dim4 = 0; dim4 < grid_num4; dim4++){
	    double dim_mean4 = model_mean4_min + dim4*res4;
	    double dim_sig4 = model_sig4_min + dim4*res4;
	    
	    grid_fits.push_back(new TF1(Form("f_test%d_%d_%d_%d",dim1,dim2,dim3,dim4),"gaus(0)*gaus(3)*gaus(6)*gaus(9)", min_x, max_x) );
	    
	    grid_fits[grid_counter]->SetParameter(0,1.0);
	    grid_fits[grid_counter]->SetParameter(1,dim_mean1);
	    grid_fits[grid_counter]->SetParameter(2,dim_sig1);

	    grid_fits[grid_counter]->SetParameter(3,1.0);
	    grid_fits[grid_counter]->SetParameter(4,dim_mean2);
	    grid_fits[grid_counter]->SetParameter(5,dim_sig2);

	    grid_fits[grid_counter]->SetParameter(6,1.0);
	    grid_fits[grid_counter]->SetParameter(7,dim_mean3);
	    grid_fits[grid_counter]->SetParameter(8,dim_sig3);

	    grid_fits[grid_counter]->SetParameter(9,1.0);
	    grid_fits[grid_counter]->SetParameter(10,dim_mean4);
	    grid_fits[grid_counter]->SetParameter(11,dim_sig4);
	 	    	    
	    //grid_fits.push_back(f_test);
	    double chi2 = pow((dim_mean1 - f_model->GetParameter(1) ),2) + pow((dim_sig1 - f_model->GetParameter(2) ),2) +  pow((dim_mean2 - f_model->GetParameter(4) ),2) + pow((dim_sig2 - f_model->GetParameter(5) ),2) + pow((dim_mean3 - f_model->GetParameter(7) ),2) + pow((dim_sig3 - f_model->GetParameter(8) ),2) + pow((dim_mean4 - f_model->GetParameter(10) ),2) + pow((dim_sig4 - f_model->GetParameter(11) ),2);
	    //std::cout << ">> grid counter " << grid_counter << " >> chi2 " << chi2 <<  std::endl;
	    //std::cout << " >> " << dim_mean1 << " " << dim_sig1 << " " << dim_mean2 << " " << dim_sig2 << " " << dim_mean3 << " " << dim_sig3 << " " << dim_mean4 << " " << dim_sig4 << std::endl;
	    grid_chi[grid_counter]=chi2;
	    grid_counter++;
	    	   	  	    
	  }
	}
      }
    }
    
    std::map<int, double> grid_results = GetBestGridPoint( grid_chi );    
    for( std::map<int,double>::iterator it = grid_results.begin(); it != grid_results.end(); ++ it){
      best_chi2 = it->second;
      best_grid_point = it->first;      
    }
    
    PrintResults( grid_fits, best_grid_point, best_chi2 );    
    DrawBestResults(f_model, grid_fits, best_grid_point);

    TCanvas *c2 = new TCanvas("c2","c2",900,900);
    g_grids->SetMarkerStyle(20);    
    g_grids->SetMarkerSize(1);
    g_grids->SetMarkerColor(kBlue);
    g_grids->Draw("PO");
    
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
      
      float dim_mean1 = model_mean1_min + (rand()) /((RAND_MAX/(model_mean1_max-model_mean1_min)));
      float dim_sig1 = model_sig1_min + (rand()) /((RAND_MAX/(model_sig1_max-model_sig1_min)));
    
      float dim_mean2 = model_mean2_min + (rand()) /((RAND_MAX/(model_mean2_max-model_mean2_min)));
      float dim_sig2 = model_sig2_min + (rand()) /((RAND_MAX/(model_sig2_max-model_sig2_min)));
    
      float dim_mean3 = model_mean3_min + (rand()) /((RAND_MAX/(model_mean3_max-model_mean3_min)));
      float dim_sig3 = model_sig3_min + (rand()) /((RAND_MAX/(model_sig3_max-model_sig3_min)));
    
      //      float dim_mean4 = model_mean4_min + (rand()) /((RAND_MAX/(model_mean4_max-model_mean4_min)));
      //float dim_sig4 = model_sig4_min + (rand()) /((RAND_MAX/(model_sig4_max-model_sig4_min)));
        
      random_grid_fits[i]->SetParameter(0,1.0);
      random_grid_fits[i]->SetParameter(1,dim_mean1);
      random_grid_fits[i]->SetParameter(2,dim_sig1);
    
      random_grid_fits[i]->SetParameter(3,1.0);
      random_grid_fits[i]->SetParameter(4,dim_mean2);
      random_grid_fits[i]->SetParameter(5,dim_sig2);
    
      random_grid_fits[i]->SetParameter(6,1.0);
      random_grid_fits[i]->SetParameter(7,dim_mean3);
      random_grid_fits[i]->SetParameter(8,dim_sig3);
    
      //random_grid_fits[i]->SetParameter(9,1.0);
      //random_grid_fits[i]->SetParameter(10,dim_mean4);
      //random_grid_fits[i]->SetParameter(11,dim_sig4);
    
      double chi2 = pow((dim_mean1 - f_model->GetParameter(1) ),2) + pow((dim_sig1 - f_model->GetParameter(2) ),2) +  pow((dim_mean2 - f_model->GetParameter(4) ),2) + pow((dim_sig2 - f_model->GetParameter(5) ),2) + pow((dim_mean3 - f_model->GetParameter(7) ),2) + pow((dim_sig3 - f_model->GetParameter(8) ),2);// + pow((dim_mean4 - f_model->GetParameter(10) ),2) + pow((dim_sig4 - f_model->GetParameter(11) ),2);
      random_grid_chi[i]=chi2;
      iteration.push_back(i);
      iteration_chi2.push_back(chi2);


    }

    std::map<int, double> random_grid_results = GetBestGridPoint( random_grid_chi );    
    for( std::map<int,double>::iterator it = random_grid_results.begin(); it != random_grid_results.end(); ++ it){
      best_chi2 = it->second;
      best_grid_point = it->first;      
    }
    
    PrintResults(random_grid_fits, best_grid_point, best_chi2);
    DrawBestResults(f_model, random_grid_fits, best_grid_point);
      
    TCanvas *c2 = new TCanvas("c2","c2",900,900);
    c2->Divide(1,1);
    c2->cd(1);
    TGraph *g_chi2 = new TGraph(iteration.size(), &iteration[0],&iteration_chi2[0]);
    g_chi2->SetTitle("#chi^{2} per iteration");
    g_chi2->SetMarkerSize(5);
    g_chi2->SetMarkerColor(kBlue);
    g_chi2->GetXaxis()->SetTitle("Iteration");
    g_chi2->GetXaxis()->CenterTitle(); 
    g_chi2->GetYaxis()->SetTitle("#chi^{2}");
    g_chi2->GetYaxis()->CenterTitle();
    g_chi2->Draw("AP");  
  }  
  
  return 0;
}

    

  

