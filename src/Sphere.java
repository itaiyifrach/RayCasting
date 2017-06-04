

public class Sphere extends Geometry implements Surface {
    private Vector center;      // position of the sphere center (x, y, z)
    private float radius;       // radius
    private int mat_idx;        // material index

    public Sphere(Vector center, float radius, int mat_idx) {
        this.center = center;
        this.radius = radius;
        this.mat_idx = mat_idx;
    }

    public Vector findIntersection(Ray ray) {
        // using geometric method (not algebraic)
        // L = O - p0
        Vector L = center.minus(ray.getStart());

        // t_ca = L * V
        double t_ca = L.dot(ray.getDirection());
        if (t_ca < 0) {     // no intersection
            //System.out.println("no inter 1");
            return null;
        }

        double d_2 = L.dot(L) - (t_ca * t_ca);
        double r_2 = radius * radius;
        if (d_2 > r_2) {    // no intersection
            //System.out.println("no inter 2");
            return null;
        }

        // two intersection points:
        // t_hc = sqrt(r^2 - d^2)
        double t_hc = Math.sqrt(r_2 - d_2);

        // two intersection points...
        double t0 = t_ca - t_hc;
        double t1 = t_ca + t_hc;

        // checking which point is closer to p0 (camera position)
        Vector p1 = ray.getPoint(t0);
        Vector p2 = ray.getPoint(t1);
        if (p1.distanceTo(ray.getStart()) <=  p2.distanceTo(ray.getStart())) {
            return p1;
        }
        else {
            return p2;
        }
    }

    public Vector getNormalVector(Vector vec)  {
        // N = P - O
        Vector N = vec.minus(center);
        // N = (P - O) / ||P - O||
        N = N.direction();

        return N;
    }

    public int getMaterialIndex() {
        return mat_idx;
    }

    @Override
    public Material getMaterial() {
        return this.material;
    }
}
