 
package asymmetry;

public class SYMAT
{
/** Constructor: constructs a lower triangle  corresponding to the integer n of a square matrix
Sets all elements of the array to zero **/

    public SYMAT(int n)
    {
        rows = n;
        B = new double[(n * (n + 1)) / 2];
    }
/** Constructor: constructs a copy of the lower triangle Array. The integer n is the number of rows corresponding to the square matrix
All elements of the array are equal to the array passed on to the constructor **/

    public SYMAT(double A[], int n)
    {
        rows = n;
        B = new double[(n * (n + 1)) / 2];
        for(ik = 0; ik < A.length; ik++)
            B[ik] = A[ik];

    }
/** Constructor:  Creates a Sums of squares and crossproduct matrix (SSCP) of the indexed array, no check is made of varying column dimensions. The number of rows of the SSCP matrix is taken equal to the number of columns in the first row**/

    public SYMAT(double Mat[][])
    {
        nobs = Mat.length;
        rows = Mat[0].length;
        B = new double[(rows * (rows + 1)) / 2];
        for(ik = 0; ik < rows; ik++)
            for(jk = ik; jk < rows; jk++)
            {
                number = 0.0D;
                for(kk = 0; kk < Mat.length; kk++)
                    number += Mat[kk][ik] * Mat[kk][jk];

                setValue(ik, jk, number);
            }
    }
/** Constructor:  Creates a Variance Covariance Matrix (VC) of the indexed array, no check is made of varying column dimensions. The number of rows of the VC matrix is taken equal to the number of columns in the first row. The boolean b is a flag, distinguising the current constructor from the previous one**/

    public SYMAT(double Mat[][],boolean b)
    {
        nobs = Mat.length;
        rows=Mat[0].length+1;
        B = new double[(rows * (rows + 1)) / 2];

        double Mat2[][]=new double[Mat.length][Mat[0].length];

         for(int i = 0; i < Mat.length; i++)
              for(int j = 0; j < Mat[0].length; j++)
                 Mat2[i][j]=Mat[i][j];

        for(kk = 0; kk < Mat.length; kk++)
        {
            for(jk = 1; jk <= Mat[0].length; jk++)
            {
                double d = Mat2[kk][jk-1] - B[((jk + 1) * jk) / 2];
                for(ik = jk; ik <= Mat[0].length; ik++)
                B[((ik+1)*ik)/2+jk]+=((double)kk*d*(Mat2[kk][ik-1]-B[((ik + 1) * ik) / 2])) / (double)(kk + 1);

            }


            for(jk = 0; jk < Mat[0].length; jk++)
                B[((jk + 2) * (jk + 1)) / 2] = ((double)kk * B[((jk + 2) * (jk + 1)) / 2] + Mat2[kk][jk]) / (double)(kk + 1);

        }
        B[0]*=-1.0/Mat.length;
    }
/** Returns the number of rows**/

    public final int getRowDimension()
    {
        return rows;
    }
/** Returns the size of the array**/

    public final int getLength()
    {
        return B.length;
    }
/** Sets the value in the appropriate place in the array**/

    public final void setValue(int i, int j, double val)
    {
        if(i > j)
        {
            B[((i + 1) * i) / 2 + j] = val;
            return;
        }
        else
        {
            B[((j + 1) * j) / 2 + i] = val;
            return;
        }
    }
    /** Adds the value to the value in the appropriate place in the array**/

    public final void addValue(int i, int j, double val)
    {
        if(i > j)
        {
            B[((i + 1) * i) / 2 + j] += val;
            return;
        }
        else
        {
            B[((j + 1) * j) / 2 + i] += val;
            return;
        }
    }

/** Multiplies the the value in the array indexed by i and j by the specified value**/

    public final void multiply(int i, int j, double val)
    {
        if(i > j)
        {
            B[((i + 1) * i) / 2 + j] *= val;
            return;
        }
        else
        {
            B[((j + 1) * j) / 2 + i] *= val;
            return;
        }
    }
/** Sets all elemenst of the array to zero **/
    public final void setzero()
    {
        for(ik = 0; ik < B.length; ik++)
            B[ik] = 0.0D;

    }
/** Returns the array of the lower triangular part of this object**/

    public double[] getArray()
    {
        double C[] = new double[B.length];
        for(ik = 0; ik < B.length; ik++)
            C[ik] = B[ik];
        return C;

    }
/** Gets the value at the appropriate place in the array**/

    public final double getValue(int i, int j)
    {
        if(i > j)
            number = B[((i + 1) * i) / 2 + j];
        else
            number = B[((j + 1) * j) / 2 + i];
        return number;
    }
/** Returns a String representation the matrix**/

    public String toString()
    {
        String s = "";
        for(ik = 0; ik < rows; ik++)
        {
            for(jk = ik; jk < rows; jk++)
                s += getValue(ik, jk) + "  ";

            s += "\n";
        }

        return s;
    }
/** Returns the number of observations on which this lower triangle in the case of a VC or SSCP matrix is based**/

    public final int getnobs()
    {
        return nobs;
    }

    public int nobs;

    double B[];
    double number;
    int rows;
    int ik;
    int jk;
    int kk;
}
