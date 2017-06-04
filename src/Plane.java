

public class Plane extends Geometry implements Surface {
    private Vector position;
    private double offset;
    private Vector pointOnPlane;
    private int mat_idx;

    public Plane(Vector position, double offset, int mat_idx) {
        this.position = new Vector(position);
        this.offset = offset;
        this.mat_idx = mat_idx;

        double x = this.position.cartesian(0);
        double y = this.position.cartesian(1);
        double z = this.position.cartesian(2);
        if (x != 0) {
            this.pointOnPlane = new Vector(this.offset / x, 0, 0);
        }
        else if (y != 0) {
            this.pointOnPlane = new Vector(0, this.offset / y, 0);
        }
        else if (z != 0) {
            this.pointOnPlane = new Vector(0, 0, this.offset / z);
        }
    }

    public Vector findIntersection(Ray ray) {
        // Ray = p0 + tV
        Vector p0 = new Vector(ray.getStart());
        Vector V = new Vector(ray.getDirection());

        if (V.dot(position) == 0) {
            return null;
        }

        Vector temp = pointOnPlane.minus(p0);
        double temp1 = temp.dot(position);
        // V*N
        double temp2 = V.dot(position);

        double t = temp1 / temp2;
        if (t < 0) {
            return null;
        }
        else {
            return ray.getPoint(t);
        }
    }

    public Vector getNormalVector(Vector vec)  {
        return position;
    }

    public int getMaterialIndex() {
        return mat_idx;
    }

    public Vector getPointOnPlane() {
        return pointOnPlane;
    }

    @Override
    public Material getMaterial() {
        return this.material;
    }
}
