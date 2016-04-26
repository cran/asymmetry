package asymmetry;

/*
* Cholesky Decompostion  version 1.0
* @Author: The Code Elves
* @version: 1.0
*/

/** The theory of the Cholesky decomposition is described in Thisted (1988), Elements of Statistical Computing. Other references include Schafer (1997)  Analysis of incomplete data, among many others.**/

public class Cholesky extends SYMAT
{
/**  Creates a Cholesky object with an empty lower triangle. The lower triangle is replaced with data on which a Cholesky decomposition is desired**/


    public Cholesky(int n)
    {
        super(n);
    }
/**  Creates a Cholesky decomposition of the symmetric array Mat. No check of symmetry is made, so be sure the array is symmetric.**/

    public Cholesky(double Mat[][])
        throws CholeskyException
    {
        super(Mat);
        for(i = 0; i < getRowDimension(); i++)
        {
            sum = getValue(i, i);
            for(k = 0; k <= i - 1; k++)
                sum -= getValue(k, i) * getValue(k, i);

            if(sum > 9.9999999999999995E-007D)
                setValue(i, i, Math.sqrt(sum));
            else
                throw new CholeskyException();
            for(j = i + 1; j < getRowDimension(); j++)
            {
                sum = getValue(i, j);
                for(k = 0; k <= i - 1; k++)
                    sum -= getValue(k, i) * getValue(k, j);

                setValue(i, j, sum / getValue(i, i));
            }

        }

    }
/**  Performs a Cholesky decomposition of the symmetric array. Used when a blank Cholesky object is created, and the elements are put into this object later on by the setmatrix method**/

    public void toUppertriangle()
        throws CholeskyException
    {
        for(i = 0; i < getRowDimension(); i++)
        {
            sum = getValue(i, i);
            for(k = 0; k <= i - 1; k++)
                sum -= getValue(k, i) * getValue(k, i);

            if(sum > 9.9999999999999995E-007D)
                setValue(i, i, Math.sqrt(sum));
            else
                throw new CholeskyException();
            for(j = i + 1; j < getRowDimension(); j++)
            {
                sum = getValue(i, j);
                for(k = 0; k <= i - 1; k++)
                    sum -= getValue(k, i) * getValue(k, j);

                setValue(i, j, sum / getValue(i, i));
            }

        }

    }
/** 
   * Solves the equation Ab = y, where A is a symmetric matrix 
   * @return returns a solution b
**/
    public double[] solve(double y[])
    {
        double temp2[] = new double[rows];
        double sum = 0.0D;
        temp2[0] = y[0] / getValue(0, 0);
        for(i = 1; i < rows; i++)
        {
            sum = 0.0D;
            for(j = 0; j < i; j++)
                sum += getValue(i, j) * temp2[j];

            temp2[i] = (y[i] - sum) / getValue(i, i);
        }

        double temp[] = new double[rows];
        temp[rows - 1] = temp2[rows - 1] / getValue(rows - 1, rows - 1);
        for(i = rows - 1; i >= 0; i--)
        {
            sum = 0.0D;
            for(j = rows - 1; j >= i + 1; j--)
                sum += getValue(i, j) * temp[j];

            temp[i] = (temp2[i] - sum) / getValue(i, i);
        }

        return temp;
    }

/** sets the symmetric array into this Cholesky decomposition **/
    public void setdouble(double A[][])
    {
        for(i = 0; i < rows; i++)
        {
            setValue(i, i, A[i][i]);
            for(j = i + 1; j < rows; j++)
                setValue(i, j, A[i][j]);

        }

    }
    public double [][] getdouble()
    {
        double tmp[][]=new double[rows][rows];
        for(i = 0; i < rows; i++)
        {
            tmp[i][i]=getValue(i, i);
            for(j = i + 1; j < rows; j++)
            {
               tmp[i][j]=getValue(i, j);
               tmp[j][i]=getValue(i, j);

            }

        }
    return tmp;

    }
/**  Creates  the inverse of this matrix, this object is constructed by a direct call to this object, or a sequence of calls consisting of first a call to the blank constructor and second to the toUppertriangle() method **/

    public void invert()
    {
        SYMAT G = new SYMAT(B, rows);
        double sum = 0.0D;
        for(j = rows - 1; j >= 0; j--)
        {
            setValue(j, j, 1.0D / G.getValue(j, j));
            for(k = j - 1; k >= 0; k--)
            {
                sum = 0.0D;
                for(i = k + 1; i <= j; i++)
                    sum -= G.getValue(k, i) * getValue(i, j);

                setValue(k, j, sum / G.getValue(k, k));
            }

        }

        G = new SYMAT(B, rows);
        for(i = 0; i < rows; i++)
            for(j = i; j < rows; j++)
            {
                sum = 0.0D;
                for(k = j; k < rows; k++)
                    sum += G.getValue(i, k) * G.getValue(k, j);

                setValue(i, j, sum);
            }


    }

    private int i;
    private int j;
    private int k;
    private double sum;
}
