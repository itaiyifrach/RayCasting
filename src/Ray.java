

public class Ray {
    Vector start;
    Vector direction;

    public Ray(Vector start, Vector direction) {
        this.start = new Vector(start);
        this.direction = new Vector(direction);
    }

    public static Ray constructRayThroughPixel(Camera camera, int i, int j, int imageWidth, int imageHeight, double pixelWidth, double pixelHeight) {
        Vector P = camera.getScreenCenterPoint();
        double moveX = (imageWidth / 2 - i) * pixelWidth;
        double moveY = (imageHeight / 2 - j) * pixelHeight;

        Vector vecX = camera.getRightDirection();
        vecX = vecX.scale(moveX);

        Vector vecY = camera.getUpDirection();
        vecY = vecY.scale(moveY);

        // calculating P new position from the center
        P = P.plus(vecX).plus(vecY);

        // calculating direction vector (P - p0)
        Vector p0 = camera.getPosition();
        Vector V = (P.minus(p0)).direction();

        return new Ray(p0, V);
    }

    /**
     *  calculating (P = p0 + t*V)
     */
    public Vector getPoint(double t) {
        return start.plus(direction.scale(t));
    }

    public Vector getStart() {
        return start;
    }

    public Vector getDirection() {
        return direction;
    }

    public void printRay() {
        System.out.println(start.toString() + "+ t" + direction.toString());
    }
}
