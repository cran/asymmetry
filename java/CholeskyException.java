
package asymmetry;
public class CholeskyException extends Exception
{
/**  displays the messag
     * @return e**/


    public String toString()
    {
        return new String("Matrix not positive semidefinite");
    }
/** Constructor: constructs  a Cholesky exception which is thrown when the the matrix is not positive definite **/


    CholeskyException()
    {
    }
}
