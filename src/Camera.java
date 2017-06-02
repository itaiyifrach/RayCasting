

public class Camera {
    private Vector position;
    private Vector look;
    private Vector up;
    private Vector towardsDirection;
    private Vector rightDirection;
    private Vector upDirection;
    private float screenDist;
    private float screenWidth;

    public Camera(Vector position, Vector look, Vector up, float dist, float width) {
        this.position = position;
        this.look = look;
        this.up = up;
        this.screenDist = dist;
        this.screenWidth = width;

        // calculating towards, right, and up directions of camera
        this.towardsDirection = (look.minus(position)).direction();
        this.rightDirection = (towardsDirection.cross(up)).direction();
        this.upDirection = rightDirection.cross(towardsDirection);
    }

    public Vector getScreenCenterPoint() {
        // center point is (p0 + d*towards)
        Vector center = position.plus(towardsDirection.scale(screenDist));

        return center;
    }

    public Vector getPosition() {
        return position;
    }

    public Vector getLook() {
        return look;
    }

    public Vector getUp() {
        return up;
    }

    public Vector getTowardsDirection() {
        return towardsDirection;
    }

    public Vector getRightDirection() {
        return rightDirection;
    }

    public Vector getUpDirection() {
        return upDirection;
    }

    public float getScreenDist() {
        return screenDist;
    }

    public float getScreenWidth() {
        return screenWidth;
    }
}
