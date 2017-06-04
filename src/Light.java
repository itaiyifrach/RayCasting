import java.util.ArrayList;
import java.util.Random;

public class Light {
    public Vector getPosition() {
        return position;
    }

    private Vector position;
    private Color color;
    private float spec;
    private float shadow;
    private float radius;

    public Light(Vector position, Vector color, float spec, float shadow, float radius) {
        this.position = position;
        this.color = new Color(color);
        this.spec = spec;
        this.shadow = shadow;
        this.radius = radius;
    }

    public Light(Vector position, Color color, float spec, float shadow, float radius) {
        this.position = position;
        this.color = color;
        this.spec = spec;
        this.shadow = shadow;
        this.radius = radius;
    }

    public Light[] getAreaLight(Vector hitPoint, int sh_rays) {
        // get direction of light (the Normal of the area grid)
        Vector N = (hitPoint.minus(position)).direction();
        // get the offset on the Plane
        double offset = position.dot(N);
        // constructing the plane
        Plane lightPlane = new Plane(N, offset, -1);


        // TODO: compute the light grid
        Vector rightDirection = (lightPlane.getPointOnPlane().minus(position)).direction();
        Vector upDirection = rightDirection.cross(N);

        Light lights[] = new Light[sh_rays*sh_rays];
        double lightWidth = radius / sh_rays;
        double lightHeight = lightWidth;
        double randomX, randomY;

        int count = 0;
        double moveRight, moveUp;
        Vector vecX, vecY, pos;
        //build point lights
        for (int i = 0; i < sh_rays; i++) {
            for (int j = 0; j < sh_rays; j++) {
                randomX = lightWidth * new Random().nextDouble();
                randomY = lightHeight * new Random().nextDouble();
                moveRight = ((radius / 2 - i - 1) * lightWidth) + (randomX * lightWidth);
                moveUp = ((radius / 2 - j - 1) * lightHeight) + (randomY * lightHeight);
                vecX = rightDirection.scale(moveRight);
                vecY = upDirection.scale(moveUp);

                pos = position.plus(vecX).plus(vecY);

                lights[count] = new Light(pos, color, spec, shadow, radius);
                count++;
            }
        }

        return lights;
    }

    public double computeSoftShadow(Light[] lights, Intersection hit, Scene scene) {
        double count = 0;

        for (Light light: lights) {
            Ray lightRay = new Ray(light.position, (hit.getPoint().minus(light.position).direction()));  // the ray from light to hit point
            count += light.getShadowHit(lightRay, hit, scene);
        }
        return (count / lights.length);
    }

    private double getShadowHit(Ray ray, Intersection hit, Scene scene) {
        Vector tempHit;
        Vector p0 = ray.getStart();
        double tempDist, tempTransValue;
        int mat_idx;
        ArrayList<Double> hitTranspValues = new ArrayList<>();

        // distance to checked surface
        double maxDist = hit.getPoint().distanceTo(p0);

        // searching if there is an intersection before checked surface. is so, then add it to the Surfaces List
        for (Surface surface : scene.getSurfaces()) {
            if (surface == hit.getSurface()) {          // if the checked surface is the same then continue to next surface
                continue;
            }
            tempHit = surface.findIntersection(ray);
            if (tempHit != null) {
                tempDist = tempHit.distanceTo(p0);
                if (tempDist < maxDist) {               // if there is hit before the checked surface
                    mat_idx = surface.getMaterialIndex();
                    Material hitMat = scene.getMaterials().get(mat_idx - 1);
                    tempTransValue = hitMat.getTrans();
                    if (tempTransValue == 0) {       // if any of the surfaces isn't transparent, then exit
                        return (1 - this.shadow);
                    } else {                              // the surface is transparent
                        hitTranspValues.add(tempTransValue);
                    }
                }
            }
        }
        // if no surfaces found, then it is direct hit from the light
        if (hitTranspValues.size() == 0) {
            return 1;
        }
        // otherwise, mult all transparency values
        else {
            double mults = 0;
            for (Double transValue: hitTranspValues) {
                if (mults == 0) {
                    mults = transValue;
                }
                else {
                    mults *= transValue;
                }
            }
            double sum = mults + (1 - this.shadow);
            // the sum might be greater than 1. in that case, we return 1
            return Math.min(sum, 1);
        }
    }

    public Color getColor() {
        return color;
    }

    public float getSpec() {
        return spec;
    }

    public float getShadow() {
        return shadow;
    }

    public float getRadius() {
        return radius;
    }
}

