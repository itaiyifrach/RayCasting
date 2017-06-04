import java.util.Comparator;

public class Intersection implements Comparable {

    private Surface hitSurface;
    private Vector hitPoint;
    private Ray hitRay;

    public Intersection(Surface hitSurface, Vector hitPoint, Ray hitRay) {
        this.hitPoint = new Vector(hitPoint);
        this.hitSurface = hitSurface;
        this.hitRay = hitRay;
    }

    public Ray getReflectionRay() {
        Vector dir = hitRay.getDirection();

        // getting normal vector
        Vector N = getNormal();

        // reflection direction
        double dp = dir.dot(N);
        if(dp < 0) {
            dp = 0;
        }

        // calculating reflection direction vector
        double[] result = {0.0, 0.0, 0.0};
        for (int i = 0; i < 3; i++) {
            double normalAxis = N.cartesian(i);
            double directionAxis = dir.cartesian(i);
            result[i] = directionAxis - 2*normalAxis*dp;
        }
        Vector refDirection = (new Vector(result)).direction();

        return new Ray(hitRay.getStart(), refDirection);
    }

    public Vector getNormal() {
        Vector N = hitSurface.getNormalVector(hitPoint);
        //N = N.scale(-1);

        return N;
    }

    public Ray getHitRay() {
        return hitRay;
    }

    public Surface getSurface() {
        return hitSurface;
    }

    public Vector getPoint() {
        return hitPoint;
    }


    @Override
    public int compareTo(Object o) {
        return 0;
    }
}
