const int n = 4, m = 10;
double a[n][m], b[m][n], x[m], y[n];

void xlog(int& np, double* deriv, double& f, double par[], int flag)
{
  // load problem parameters
  for (int i = 0; i<n; ++i) { 
    for (int j = 0; j<m; ++j) {
      a[i][j] = par[i*m+j];
    }
  }

  // AB = I
  double sum1 = 0;
  for (int i = 0; i<n; ++i) {
    for (int k = 0; k<n; ++k) {
      double sum = 0;
      for (int j = 0; j<m; ++j) {
	sum += a[i][j]*b[j][k];
      }
      if (i==k) sum1 += (sum-1)*(sum-1); else sum1 += sum*sum;
    }
  }

  // Ax = y
  double sum2 = 0;
  for (int i = 0; i<n; ++i) {
    double sum = 0;
    for (int j = 0; j<m; ++j) {
      sum += a[i][j]*x[j];
    }
    sum2 += (sum-y[i])*(sum-y[i]);
  }

  // min ||A||
  double sum3 = 0;
  for (int i = 0; i<n; ++i) {
    for (int j = 0; j<m; ++j) {
      sum3 += a[i][j]*a[i][j];
    }
  }

  f = (sum1+sum2)*1e6 + sum3;
}

void print(double par[])
{
  cout << "B =" << endl;
  // load problem parameters
  for (int i = 0; i<n; ++i) { 
    for (int j = 0; j<m; ++j) {
      cout << TString::Format("%8.4f",b[i][j]);
    }
    cout << endl;
  }

  cout << "x =" << endl;
  for (int j = 0; j<m; ++j) {
    cout << TString::Format("%8.4f",x[j]);
  }
  cout << endl;

  cout << "y =" << endl;
  for (int i = 0; i<n; ++i) {
    cout << TString::Format("%8.4f",y[i]);
  }
  cout << endl;

  cout << "A (solution) =" << endl;
  // load problem parameters
  for (int i = 0; i<n; ++i) { 
    for (int j = 0; j<m; ++j) {
      a[i][j] = par[i*m+j];
      cout << TString::Format("%8.4f",a[i][j]);
    }
    cout << endl;
  }

  cout << "verify A*B =" << endl;
  // AB = I
  for (int i = 0; i<n; ++i) {
    for (int k = 0; k<n; ++k) {
      double sum = 0;
      for (int j = 0; j<m; ++j) {
	sum += a[i][j]*b[j][k];
      }
      cout << TString::Format("%8.4f",sum);
    }
    cout << endl;
  }

  cout << "verify A*x-y =" << endl;
  // Ax = y
  for (int i = 0; i<n; ++i) {
    double sum = 0;
    for (int j = 0; j<m; ++j) {
      sum += a[i][j]*x[j];
    }
    cout << TString::Format("%8.4f",sum-y[i]);
  }
  cout << endl;
}

void p()
{
  TRandom3 rnd;
  // prepare the problem (B, x, y)
  for (int i = 0; i<n; ++i) {
    y[i] = rnd.Uniform(0,1);
  }
  for (int j = 0; j<m; ++j) {
    x[j] = rnd.Uniform(0,1);
  }
  for (int i = 0; i<n; ++i) {
    for (int j = 0; j<m; ++j) {
      b[j][i] = rnd.Uniform(0,1);
    }
  }

  // set up the fit procedure
  const int npar = n*m;
  TMinuit* minuit = new TMinuit(npar);
  minuit->SetFCN(xlog);
  double par[npar];
  for (int ipar = 0; ipar<npar; ++ipar) {
    minuit->DefineParameter(ipar,TString("par")+ipar,rnd.Uniform(0,1),0.1,0,0);
  }

  // minimize
  minuit->Migrad();

  // retrieve fit results
  double parerr[npar];
  for (int ipar = 0; ipar<npar; ++ipar) {
    minuit->GetParameter(ipar, par[ipar], parerr[ipar]);
  }

  print(par);
}
