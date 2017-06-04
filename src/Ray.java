import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Ray {
    Vector start;
    Vector direction;

    public Ray(Vector start, Vector direction) {
        this.start = new Vector(start);
        this.direction = new Vector(direction);
    }

    public static List<Ray> constructRayThroughPixel(Camera camera, int i, int j, int imageWidth, int imageHeight, double pixelWidth, double pixelHeight, int superSamplingLvl) {
        int ss_level = superSamplingLvl / 2;
        Vector p0 = camera.getPosition();
        List<Ray> raysTrohughPixelList = new ArrayList<>();

        Vector P = camera.getScreenCenterPoint();
        double moveX = (imageWidth / 2 - i) * pixelWidth;
        double moveY = (imageHeight / 2 - j) * pixelHeight;

        Vector vecX = camera.getRightDirection();
        vecX = vecX.scale(moveX);

        Vector vecY = camera.getUpDirection();
        vecY = vecY.scale(moveY);

        // calculating P new position from the center
        P = P.plus(vecX).plus(vecY);
        // calculating direction vector (pos - p0)
        Vector V = (P.minus(p0)).direction();

        if (ss_level == 1 || ss_level == 0) {
            raysTrohughPixelList.add(new Ray(p0,V));
            return raysTrohughPixelList;
        }

        double innerPixelWidth = pixelWidth / (double)ss_level;
        double innerPixelHeight = pixelHeight / (double)ss_level;
        double randomX, randomY;
        double moveRight, moveUp;
        Vector pos;
        vecX = camera.getRightDirection();
        vecY = camera.getUpDirection();

        //build point lights
        for (int k = 0; k < ss_level; k++) {
            for (int l = 0; l < ss_level; l++) {
                randomX = innerPixelWidth * new Random().nextDouble();
                randomY = innerPixelHeight * new Random().nextDouble();
                moveRight = ((ss_level / 2 - k - 1) * innerPixelWidth) + (randomX * innerPixelWidth);
                moveUp = ((ss_level / 2 - l - 1) * innerPixelHeight) + (randomY * innerPixelHeight);
                vecX = vecX.scale(moveRight);
                vecY = vecY.scale(moveUp);

                pos = P.plus(vecX).plus(vecY);
                // calculating direction vector (pos - p0)
                V = (pos.minus(p0)).direction();

                raysTrohughPixelList.add(new Ray(p0,V));
            }
        }
        return raysTrohughPixelList;
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
}
