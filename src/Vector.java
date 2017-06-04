/******************************************************************************
 *  Compilation:  javac Vector.java
 *  Execution:    java Vector
 *
 *  Implementation of a vector of real numbers.
 *
 *  This class is implemented to be immutable: once the client program
 *  initialize a Vector, it cannot change any of its fields
 *  (N or data[i]) either directly or indirectly. Immutability is a
 *  very desirable feature of a data type.
 *
 *
 *  % java Vector
 *  x        =  (1.0, 2.0, 3.0, 4.0)
 *  y        =  (5.0, 2.0, 4.0, 1.0)
 *  x + y    =  (6.0, 4.0, 7.0, 5.0)
 *  10x      =  (10.0, 20.0, 30.0, 40.0)
 *  |x|      =  5.477225575051661
 *  <x, y>   =  25.0
 *  |x - y|  =  5.0990195135927845
 *
 *  Note that java.util.Vector is an unrelated Java library class.
 *
 ******************************************************************************/

public class Vector {

    private final int n;         // length of the vector

    public double[] getData() {
        return data;
    }

    private double[] data;       // array of vector's components

    // create the zero vector of length n
    public Vector(int n) {
        this.n = n;
        this.data = new double[n];
    }

    // create a vector from an array
//    public Vector(double[] data) {
//        n = data.length;
//
//        // defensive copy so that client can't alter our copy of data[]
//        this.data = new double[n];
//        for (int i = 0; i < n; i++)
//            this.data[i] = data[i];
//    }

    // create a vector from either an array or a vararg list
    // this constructor uses Java's vararg syntax to support
    // a constructor that takes a variable number of arguments, such as
    // Vector x = new Vector(1.0, 2.0, 3.0, 4.0);
    // Vector y = new Vector(5.0, 2.0, 4.0, 1.0);

    public Vector(double... data) {
        n = data.length;

        // defensive copy so that client can't alter our copy of data[]
        this.data = new double[n];
        for (int i = 0; i < n; i++)
            this.data[i] = data[i];
    }

    // copy constructor
    public Vector(Vector vec) {
        this.n = vec.n;

        // defensive copy so that client can't alter our copy of data[]
        this.data = new double[n];
        for (int i = 0; i < n; i++)
            this.data[i] = vec.data[i];
    }

    // return the length of the vector
    public int length() {
        return n;
    }

    // return the inner product of this Vector a and b
    public double dot(Vector vec) {
        if (this.length() != vec.length())
            throw new IllegalArgumentException("dimensions disagree");
        double sum = 0.0;
        for (int i = 0; i < n; i++)
            sum = sum + (this.data[i] * vec.data[i]);
        return sum;
    }

    // return the cross product of this Vector a and b
    public Vector cross(Vector vec) {
        if (this.length() != vec.length())
            throw new IllegalArgumentException("dimensions disagree");
        double[] vecData = new double[n];

        vecData[0] = this.data[1] * vec.data[2] - vec.data[1] * this.data[2];
        vecData[1] = this.data[2] * vec.data[0] - vec.data[2] * this.data[0];
        vecData[2] = this.data[0] * vec.data[1] - vec.data[0] * this.data[1];

        return new Vector(vecData);
    }

    // return the Euclidean norm of this Vector
    public double norm() {
        return Math.sqrt(this.dot(this));
    }

    // return the Euclidean distance between this and that
    public double distanceTo(Vector vec) {
        if (this.length() != vec.length())
            throw new IllegalArgumentException("dimensions disagree");
        return this.minus(vec).norm();
    }

    // return this + that
    public Vector plus(Vector that) {
        if (this.length() != that.length())
            throw new IllegalArgumentException("dimensions disagree");
        Vector c = new Vector(n);
        for (int i = 0; i < n; i++)
            c.data[i] = this.data[i] + that.data[i];
        return c;
    }

    // return this - that
    public Vector minus(Vector vec) {
        if (this.length() != vec.length())
            throw new IllegalArgumentException("dimensions disagree");
        Vector c = new Vector(n);
        for (int i = 0; i < n; i++)
            c.data[i] = this.data[i] - vec.data[i];
        return c;
    }

    // return the corresponding coordinate
    public double cartesian(int i) {
        return data[i];
    }

    // create and return a new object whose value is (this * factor)
    public Vector scale(double factor) {
        Vector c = new Vector(n);
        for (int i = 0; i < n; i++)
            c.data[i] = factor * data[i];
        return c;
    }

    // return the corresponding unit vector
    public Vector direction() {
        if (this.norm() == 0.0)
            throw new ArithmeticException("zero-vector has no direction");
        return this.scale(1.0 / this.norm());
    }

    public Vector reflect(Vector reflectFrom) {
        // reflection direction
        double dp = this.dot(reflectFrom);
        if(dp < 0) {
            dp = 0;
        }

        // calculating reflection direction vector
        double[] result = {0.0, 0.0, 0.0};
        for (int i = 0; i < 3; i++) {
            double normalAxis = reflectFrom.cartesian(i);
            double directionAxis = this.cartesian(i);
            result[i] = directionAxis - 2*normalAxis*dp;
        }
        Vector refDirection = (new Vector(result)).direction();
        return refDirection;
    }

    public boolean isZeroVector(){
        return (this.data[0] == 0 && this.data[1] == 0 && this.data[2] == 0);
    }

    // return a string representation of the vector
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append('(');
        for (int i = 0; i < n; i++) {
            s.append(data[i]);
            if (i < n-1) s.append(", ");
        }
        s.append(')');
        return s.toString();
    }

    public static Vector multipleByCoordinates(Vector _X_, Vector _Y_){
        if(_X_.length() != _Y_.length()){
            return null;
        }
        int len = _X_.length();
        Vector _Z_ = new Vector();
        for (int i = 0; i < len; i++) {
            _Z_.getData()[i] = _X_.getData()[i] * _Y_.getData()[i];
        }
        return _Z_;
    }

    public boolean isEqualTo(Vector other){
        if(this.length() != other.length()){
            return false;
        }
        for (int i = 0; i < this.length(); i++) {
            if(this.data[i] != other.data[i]){
                return false;
            }
        }
        return true;
    }

    // test client
//    public static void main(String[] args) {
//        double[] xdata = { 1.0, 2.0, 3.0, 4.0 };
//        double[] ydata = { 5.0, 2.0, 4.0, 1.0 };
//
//        Vector x = new Vector(xdata);
//        Vector y = new Vector(ydata);
//
//        System.out.println("x        =  " + x);
//        System.out.println("y        =  " + y);
//        System.out.println("x + y    =  " + x.plus(y));
//        System.out.println("10x      =  " + x.scale(10.0));
//        System.out.println("|x|      =  " + x.norm());
//        System.out.println("<x, y>   =  " + x.dot(y));
//        System.out.println("|x - y|  =  " + x.minus(y).norm());
//        System.out.println("x dir  =  " + x.direction());
//        System.out.println("x = "+x.toString());
//    }
}