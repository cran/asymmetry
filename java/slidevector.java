/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package asymmetry;
import java.math.*;

/**
 *
 * @author berrie
 */
public class slidevector { 

    /**
     * @param args the command line arguments
     */
    double X[][];
    double BMAT[][],V[][],W[][],Z[][],ZVZ[][],VMAT[][],B[][];
    double diss[][];
    int dim,nobj,nobs,maxiter,iter;
    double stress,eps;
    public void makeV()
    {
        V=new double[nobj][nobj];
        for (int i=0;i<nobs;i++)
            for(int j=nobs;j<nobj;j++){
                V[i][j]=-1.0;
                V[j][i]=-1.0;
            }
        for (int i=0;i<nobs;i++)
            for(int j=nobs;j<nobj;j++){
                V[i][i]-=V[i][j];
                V[j][j]-=V[i][j];                
            }       
    }
     public void makeV(double[][] W)
    {
        V=new double[nobj][nobj];
        for (int i=0;i<nobs;i++)
            for(int j=nobs;j<nobj;j++){
                V[i][j]=-1.0*W[i][j-nobs];
                V[j][i]=-1.0*W[j-nobs][i];
            }
        for (int i=0;i<nobs;i++)
            for(int j=nobs;j<nobj;j++){
                V[i][i]-=V[i][j];
                V[j][j]-=V[i][j];                
            }       
    }
         public double[][] slidevectormodel(double d[][], int ndim, int maxiterations,double epsi)
    {
        nobs = d[0].length;
        nobj = 2*nobs;
        maxiter=maxiterations;
        dim=ndim;
        eps=epsi;
                diss=new double[nobj][nobj];
        for (int i=0;i<nobs;i++) {
            for (int j=0;j<nobs;j++) {
                diss[i][j+nobs]=d[i][j];
                diss[j+nobs][i]=diss[i][j+nobs];
            }
        }

        smacof();
        info();
   
        return getConfiguration();
    }

    public void smacof() {
        makeV();
        Z = new double[nobj][nobs+1];
        B = new double[nobs+1][dim];
        for(int i=0;i<nobs;i++){
            Z[i][i]=1.0;
            Z[nobs+i][i]=1.0;
            Z[i][nobs]=1.0;
        }
        ZVZ = new double[nobs+1][nobs+1];
               for(int i=0; i<nobs+1; i++)
                for(int j=0; j<nobs+1; j++) 
                  for (int k=0; k<nobj; k++)
                      for(int l=0; l<nobj; l++)                       
                        ZVZ[i][j]+=Z[k][i]*V[k][l]*Z[l][j];
                        
        Cholesky cd = new Cholesky(nobs+1);
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
        X = new double[nobj][dim];
        for (int i=0;i<nobj;i++) {
            for (int j=0;j<dim;j++) {
                X[i][j]=Math.random();
            }
        }        
        BMAT = BMAT();
        double XTIL[][] = new double[nobj][dim];
                     for (int i=0; i<nobj; i++)
                      for(int j=0; j<dim; j++)
                          for(int k=0;k<nobj;k++)
                          XTIL[i][j]+=BMAT[i][k]*X[k][j];
        
    }
    public double[][] BMAT(){
        for (int i=0;i<nobj;i++) {
            for (int j=0;j<nobj;j++) {
                if(distance(X,dim,i,j)>1E-16)
                BMAT[i][j]=V[i][j]*diss[i][j]/distance(X,dim,i,j);
                else BMAT[i][j]=0.0D;
            }
        }
        for(int i=0;i<nobj;i++)
            BMAT[i][i]=0.0D;
 
        for (int i=0;i<nobj;i++) {
            for (int j=0;j<nobj;j++) {
                if(i!=j)
                BMAT[i][i]-=BMAT[i][j];           
            }
        }
        return BMAT;
        
    }
    public double[][] update(){
        double XT[][] = new double[nobj][dim];
                   for (int i=0; i<nobj; i++)
                      for(int j=0; j<dim; j++)  
                          for(int k=0; k<nobj; k++)
                              XT[i][j]+=BMAT[i][k]*X[k][j];

        return XT;
    }
    public double[][] weight(){
        double XT[][] = new double[nobs+1][dim];
                for (int i=0; i<nobs+1; i++)
                    for(int j=0; j<dim; j++)      
                         B[i][j] = 0.0D;      
                    for (int i=0; i<nobs+1; i++)
                      for(int j=0; j<dim; j++)
                          for(int k=0;k<nobj;k++)
                          XT[i][j]+=Z[k][i]*X[k][j];       
                   for (int i=0; i<nobs+1; i++)
                      for(int j=0; j<dim; j++)
                          for(int k=0;k<nobs+1;k++)
                          B[i][j]+=VMAT[i][k]*XT[k][j];
                   
                   for (int i=0; i<nobj; i++)
                      for(int j=0; j<dim; j++)
                         X[i][j]=0.0;                  
                   for (int i=0; i<nobj; i++)
                      for(int j=0; j<dim; j++)
                          for(int k=0;k<nobs+1;k++)
                          X[i][j]+=Z[i][k]*B[k][j];                  
        return X;
    }
    
    public void info()
    {
        stress=stress(diss,X,dim,nobj);
        double stressold=stress+10.0D;
        while(iter<maxiter&&(stressold-stress)>eps)
        { 
            stressold=stress;
            BMAT = BMAT();
             X = update();             
             X = weight(); 
             stress=stress(diss,X,dim,nobj);            
             System.out.println(" stress " +stress);
             iter++;
         }   
        
    }
    //aanpassen voor gewichten
    public double stress(double diss[][],double X[][],int dim, int nobj){
        double stress = 0.0D;
                   for (int i=0; i<nobs; i++)
                      for(int j=nobs; j<nobj; j++)  
                        stress+=Math.pow(diss[i][j]-distance(X,dim,i,j),2);

        return stress;
    }

    public double[][] getConfiguration()
    {
        return B;
    }
    public static double distance(double X[][],int dim,int i,int j){
        double result=0.0;  
        for(int k=0;k<dim;k++) {
            result+=Math.pow(X[i][k]-X[j][k],2);                   
        }
        return Math.sqrt(result);
    
    }
    public int getDimension()
    {
        return dim;
    }
    public int getNiter()
    {
        return iter;
    }
    public int getNobs()
    {
        return nobs;
    }
    public double getStress()
    {
        return stress;
    }
    
}
