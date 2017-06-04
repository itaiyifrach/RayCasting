import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Ray {
    Vector start;

    public void setDirection(Vector direction) {
        this.direction = direction;
    }

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
        //System.out.print("hey from constructRayThroughPixel before do direction method on" + p0.toString());
        Vector V = (P.minus(p0)).direction();

        return new Ray(p0, V);
    }

    public static List<Ray> constructRayThroughPixel(Camera camera, int i, int j, int imageWidth, int imageHeight, double pixelWidth, double pixelHeight, int superSamplingLvl) {

        List<Ray> raysTrohughPixelList = new ArrayList<>();
        Random rnd = new Random();

        Vector P = camera.getScreenCenterPoint();
        double moveX = (imageWidth / 2 - i) * pixelWidth;
        double moveY = (imageHeight / 2 - j) * pixelHeight;

        Vector vecX = camera.getRightDirection();
        vecX = vecX.scale(moveX);

        Vector vecY = camera.getUpDirection();
        vecY = vecY.scale(moveY);

        // calculating P new position from the center
        P = P.plus(vecX).plus(vecY);

        double innerPixelWidth = pixelWidth / (double)superSamplingLvl;
        double innerPixelHeigh = pixelHeight / (double)superSamplingLvl;


        for (int k = 0; k < superSamplingLvl; k++) {
            for (int l = 0; l < superSamplingLvl; l++) {

                double innerMoveX = (pixelWidth / 2 - i) * innerPixelWidth;
                double innerMoveY = (pixelHeight / 2 - j) * innerPixelHeigh;

                double infMoveRangeX = innerMoveX/4;
                double infMoveRangeY = innerMoveX/4;

                Vector innerVecX = camera.getRightDirection();
                Vector innerVecY = camera.getUpDirection();


                innerMoveX =  infMoveRangeX * (rnd.nextDouble());
                innerVecX = innerVecX.scale(innerMoveX);

                innerMoveY =  infMoveRangeY * (rnd.nextDouble());
                innerVecY = innerVecY.scale(innerMoveY);

                // calculating P new position from the center
                P = P.plus(innerVecX).plus(innerVecY);

                // calculating direction vector (P - p0)
                Vector p0 = camera.getPosition();
                //System.out.print("hey from constructRayThroughPixel before do direction method on" + p0.toString());
                Vector V = (P.minus(p0)).direction();

                raysTrohughPixelList.add(new Ray(p0,V));

            }
        }

        return raysTrohughPixelList;
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
