

public class Triangle implements Surface {
    private Vector[] vertices;
    private Vector N;
    private int mat_idx;

    public Triangle(Vector pos1, Vector pos2, Vector pos3, int mat_idx) {
        this.vertices = new Vector[3];
        this.vertices[0] = new Vector(pos1);
        this.vertices[1] = new Vector(pos2);
        this.vertices[2] = new Vector(pos3);
        this.mat_idx = mat_idx;

        Vector vec1 = vertices[1].minus(vertices[0]);
        Vector vec2 = vertices[2].minus(vertices[0]);
        this.N = (vec1.cross(vec2)).direction();
    }

    public Vector findIntersection(Ray ray) {
        // calculating the triangle plane
        double offset =  vertices[0].dot(this.N);
        Plane trgPlane = new Plane(this.N, offset, mat_idx);

        // finding intersection point with triangle plane
        Vector p = trgPlane.findIntersection(ray);
        if (p == null) {
            return null;
        }

        // checking intersection point with triangle
        // first edge
        Vector edge1 = vertices[1].minus(vertices[0]);
        Vector r1 = p.minus(vertices[0]);
        Vector normal = edge1.cross(r1);
        if (this.N.dot(normal) < 0) {
            return null;
        }
        // second edge
        Vector edge2 = vertices[2].minus(vertices[1]);
        Vector r2 = p.minus(vertices[1]);
        normal = edge2.cross(r2);
        if (this.N.dot(normal) < 0) {
            return null;
        }
        // third edge
        Vector edge3 = vertices[0].minus(vertices[2]);
        Vector r3 = p.minus(vertices[2]);
        normal = edge3.cross(r3);
        if (this.N.dot(normal) < 0) {
            return null;
        }

        // point intersects the triangle
        return p;
    }

    public Vector getNormalVector(Vector vec)  {
        return this.N;
    }

    public int getMaterialIndex() {
        return mat_idx;
    }
}
