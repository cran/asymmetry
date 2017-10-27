/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package asymmetry;

/**
 *
 * @author berrie
 */
public class unfolding extends slidevector {
    public unfolding(){
        System.out.println("hello from unfolding " );
    }
    public double[][] unfoldingmodel(double d[][], int ndim,boolean verb, int maxiterations,double epsi)
        {
        // double x[][] ;
         return slidevectormodel(d, ndim,verb, maxiterations,epsi);
         //return x;
         }

    public void smacof() {
        setDim(dim,nobs);
        makeV();
        System.out.println("hello from unfolding smacof " + nobj);
        Z = new double[nobj][nobs];
        B = new double[nobs][dim];
        for(int i=0;i<nobs;i++){
            Z[i][i]=1.0;
            Z[nobs+i][i]=1.0;
        }
        
        ZVZ = new double[nobs][nobs];
               for(int i=0; i<nobs; i++)
                for(int j=0; j<nobs; j++) 
                  for (int k=0; k<nobj; k++)
                      for(int l=0; l<nobj; l++)                       
                        ZVZ[i][j]+=Z[k][i]*V[k][l]*Z[l][j];
//create inverse                        
        Cholesky cd = new Cholesky(nobs);
        cd.setdouble(ZVZ);
            for(int i=0; i<nobs; i++)
                for(int j=i; j<nobs; j++) 
                    cd.addValue(i,j,1.0/nobs);
                    try {
                    cd.toUppertriangle();
                    } catch (CholeskyException CE) {;}
                    cd.invert();
             for(int i=0; i<nobs; i++)
                for(int j=i; j<nobs; j++)                    
                    cd.addValue(i,j,(-1.0/nobs));
        
        VMAT = cd.getdouble();
        BMAT=new double[nobj][nobj];
        X = new double[nobj][undim+nobj];
//random starting values        
        for (int i=0;i<nobj;i++) {
            for (int j=0;j<undim;j++) {
                X[i][j]=Math.random();
            }
        } 
        for (int i=0;i<nobj;i++)
            X[i][i+undim]=Math.random();
        BMAT = BMAT();
 /*       double XTIL[][] = new double[nobj][dim+nobj];
                     for (int i=0; i<nobj; i++)
                      for(int j=0; j<dim; j++)
                          for(int k=0;k<nobj;k++)
                          XTIL[i][j]+=BMAT[i][k]*X[k][j]; */
        
    }
public double[][] weight(){
        double XT[][] = new double[nobs][dim];
                for (int i=0; i<nobs; i++)
                    for(int j=0; j<dim; j++)      
                         B[i][j] = 0.0D;      
                    for (int i=0; i<nobs; i++)
                      for(int j=0; j<dim; j++)
                          for(int k=0;k<nobj;k++)
                          XT[i][j]+=Z[k][i]*X[k][j]; 
                   for (int i=0; i<nobs; i++)
                      for(int j=0; j<undim; j++)
                          for(int k=0;k<nobs;k++)
                          B[i][j]+=VMAT[i][k]*XT[k][j];                   
                   for (int i=0; i<nobj; i++)
                      for(int j=0; j<undim; j++)// < undim + 
                         X[i][j]=0.0;                  
                   for (int i=0; i<nobj; i++)
                      for(int j=0; j<undim; j++)
                          for(int k=0;k<nobs;k++)
                          X[i][j]+=Z[i][k]*B[k][j];        // symmetric constraints  
                     for (int i=0; i<nobj; i++)
                      for(int j=undim; j<dim; j++)
                      {
                         if(i==j-undim)
                             X[i][j]=X[i][j]/nobs;  
                         else X[i][j]=0.0; 
                      }
        return X;
    }
    public double[][] getConfiguration()
    {
        return X;
    }


    
    
    
    public void setDim(int a,int nobj)
    {
        dim=2*nobj+a;
        undim=a;
    }
int undim;
}
