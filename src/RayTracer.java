//package RayTracing;

import java.awt.*;
import java.awt.color.*;
import java.awt.image.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import javax.imageio.ImageIO;

/**
 *  Main class for ray tracing exercise.
 */
public class RayTracer {

    private Camera camera;
    private Scene scene;
    public int imageWidth;
    public int imageHeight;
    private Vector bg_color;    // background color (r, g, b)
    private int sh_rays;     // root number of shadow rays (N^2 rays will be shot)
    private int rec_max;        // maximum number of recursions
    private int ss_level;       // super sampling level

    /**
     * Runs the ray tracer. Takes scene file, output image file and image size as input.
     */
    public static void main(String[] args) {

        try {

            RayTracer tracer = new RayTracer();

            // Default values:
            tracer.imageWidth = 500;
            tracer.imageHeight = 500;

            if (args.length < 2)
                throw new RayTracerException("Not enough arguments provided. Please specify an input scene file and an output image file for rendering.");

            String sceneFileName = args[0];
            String outputFileName = args[1];

            if (args.length > 3)
            {
                tracer.imageWidth = Integer.parseInt(args[2]);
                tracer.imageHeight = Integer.parseInt(args[3]);
            }


            // Parse scene file:
            tracer.parseScene(sceneFileName);

            // Render scene:
            tracer.renderScene(outputFileName);

//		} catch (IOException e) {
//			System.out.println(e.getMessage());
        } catch (RayTracerException e) {
            System.out.println(e.getMessage());
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }


    }

    /**
     * Parses the scene file and creates the scene. Change this function so it generates the required objects.
     */
    public void parseScene(String sceneFileName) throws IOException, RayTracerException
    {
        FileReader fr = new FileReader(sceneFileName);

        BufferedReader r = new BufferedReader(fr);
        String line = null;
        int lineNum = 0;
        System.out.println("Started parsing scene file " + sceneFileName);

        this.scene = new Scene();

        while ((line = r.readLine()) != null)
        {
            line = line.trim();
            ++lineNum;

            if (line.isEmpty() || (line.charAt(0) == '#'))
            {  // This line in the scene file is a comment
                continue;
            }
            else
            {
                String code = line.substring(0, 3).toLowerCase();
                // Split according to white space characters:
                String[] params = line.substring(3).trim().toLowerCase().split("\\s+");

                if (code.equals("cam"))
                {
                    this.camera = new Camera(new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2])),
                                    new Vector(Double.parseDouble(params[3]), Double.parseDouble(params[4]), Double.parseDouble(params[5])),
                                    new Vector(Double.parseDouble(params[6]), Double.parseDouble(params[7]), Double.parseDouble(params[8])),
                                    Float.parseFloat(params[9]), Float.parseFloat(params[10]));

                    System.out.println(String.format("Parsed camera parameters (line %d)", lineNum));
                }
                else if (code.equals("set"))
                {
                    this.bg_color = new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2]));
                    this.sh_rays = Integer.parseInt(params[3]);
                    this.rec_max = Integer.parseInt(params[4]);
                    this.ss_level = 2; // default value is 2
                    if (params.length == 6) {
                        this.ss_level = Integer.parseInt(params[5]);
                    }

                    this.scene.setmBackGroundColor(this.bg_color);

                    System.out.println(String.format("Parsed general settings (line %d)", lineNum));
                }

                else if (code.equals("mtl"))
                {
                    Material mtr = new Material(new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2])),
                            new Vector(Double.parseDouble(params[3]), Double.parseDouble(params[4]), Double.parseDouble(params[5])),
                            new Vector(Double.parseDouble(params[6]), Double.parseDouble(params[7]), Double.parseDouble(params[8])),
                            Float.parseFloat(params[9]), Float.parseFloat(params[10]));
                    this.scene.getMaterials().add(mtr);

                    System.out.println(String.format("Parsed material (line %d)", lineNum));
                }
                else if (code.equals("sph"))
                {
                    Sphere sph = new Sphere(new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2])),
                            Float.parseFloat(params[3]), Integer.parseInt(params[4]));
                    this.scene.getSurfaces().add(sph);

                    System.out.println(String.format("Parsed sphere (line %d)", lineNum));
                }
                else if (code.equals("pln"))
                {
                    Plane pln = new Plane(new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2])),
                            Integer.parseInt(params[3]), Integer.parseInt(params[4]));
                    this.scene.getSurfaces().add(pln);

                    System.out.println(String.format("Parsed plane (line %d)", lineNum));
                }
                else if (code.equals("trg"))
                {
                    Triangle trg = new Triangle(new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2])),
                            new Vector(Double.parseDouble(params[3]), Double.parseDouble(params[4]), Double.parseDouble(params[5])),
                            new Vector(Double.parseDouble(params[6]), Double.parseDouble(params[7]), Double.parseDouble(params[8])),
                            Integer.parseInt(params[9]));
                    this.scene.getSurfaces().add(trg);

                    System.out.println(String.format("Parsed cylinder (line %d)", lineNum));
                }
                else if (code.equals("lgt"))
                {
                    Light lgt = new Light(new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2])),
                            new Vector(Double.parseDouble(params[3]), Double.parseDouble(params[4]), Double.parseDouble(params[5])),
                            Float.parseFloat(params[6]), Float.parseFloat(params[7]), Float.parseFloat(params[8]));
                    this.scene.getLights().add(lgt);

                    System.out.println(String.format("Parsed light (line %d)", lineNum));
                }
                else
                {
                    System.out.println(String.format("ERROR: Did not recognize object: %s (line %d)", code, lineNum));
                }
            }
        }

        // It is recommended that you check here that the scene is valid,
        // for example camera settings and all necessary materials were defined.

        // TODO: what should we check here??
        boolean isValid = true;
        //checking if #materials == #surfaces
//        if (this.materials.size() != this.surfaces.size()) {
//            isValid = false;
//        }

        if (isValid) {
            System.out.println("Finished parsing scene file " + sceneFileName);
        }
        else {
            System.out.println("Error in parsing scene file " + sceneFileName);
        }

    }

    /**
     * Renders the loaded scene and saves it to the specified file location.
     */
    public void renderScene(String outputFileName)
    {
        long startTime = System.currentTimeMillis();

        // Create a byte array to hold the pixel data:
        byte[] rgbData = new byte[this.imageWidth * this.imageHeight * 3];

        // Put your ray tracing code here!
        //
        // Write pixel color values in RGB format to rgbData:
        // Pixel [x, y] red component is in rgbData[(y * this.imageWidth + x) * 3]
        //            green component is in rgbData[(y * this.imageWidth + x) * 3 + 1]
        //             blue component is in rgbData[(y * this.imageWidth + x) * 3 + 2]
        //
        // Each of the red, green and blue components should be a byte, i.e. 0-255

        // TODO: complete these loops...
        List<Intersection> intersections;
        Color totalColor;
        double imageRatio = imageHeight / imageWidth;
        double pixelWidth = camera.getScreenWidth() / imageWidth;
        double pixelHeight = imageRatio * pixelWidth;
        for (int i = 0; i < this.imageWidth; i++) {
            for (int j = 0; j < this.imageHeight; j++) {
//                List<Ray> mListRays = Ray.constructRayThroughPixel(camera, i, j, imageWidth, imageHeight, pixelWidth, pixelHeight, this.getSs_level());
//                // find intersection and find the closest intersection
//                this.superSamplingMehod(mListRays, i,j);

                int r = 0;
                if (i == 171 && j == 162) {
                    r = 1;
                }

                Ray ray = Ray.constructRayThroughPixel(camera, i, j, imageWidth, imageHeight, pixelWidth, pixelHeight);
                // find intersection and find the closest intersection
                intersections = findAllIntersectionOnRay(ray);
                // calculating the pixel color
                totalColor = getColor(ray, intersections, 0, 0);

                // setting the pixel color to the byte array
                fillpixelColor(rgbData, totalColor.getRgbValues(), i, j);
                System.out.println("(" + i + ", " + j +")");
            }
        }

        long endTime = System.currentTimeMillis();
        Long renderTime = endTime - startTime;

        // The time is measured for your own conveniece, rendering speed will not affect your score
        // unless it is exceptionally slow (more than a couple of minutes)
        System.out.println("Finished rendering scene in " + renderTime.toString() + " milliseconds.");

        // This is already implemented, and should work without adding any code.
        saveImage(this.imageWidth, rgbData, outputFileName);

        System.out.println("Saved file " + outputFileName);

    }


    //////////////////////// FUNCTIONS TO SAVE IMAGES IN PNG FORMAT //////////////////////////////////////////

    /*
     * Saves RGB data as an image in png format to the specified location.
     */
    public static void saveImage(int width, byte[] rgbData, String fileName)
    {
        try {

            BufferedImage image = bytes2RGB(width, rgbData);
            ImageIO.write(image, "png", new File(fileName));

        } catch (IOException e) {
            System.out.println("ERROR SAVING FILE: " + e.getMessage());
        }
    }

    /*
     * Producing a BufferedImage that can be saved as png from a byte array of RGB values.
     */
    public static BufferedImage bytes2RGB(int width, byte[] buffer) {
        int height = buffer.length / width / 3;
        ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_sRGB);
        ColorModel cm = new ComponentColorModel(cs, false, false,
                Transparency.OPAQUE, DataBuffer.TYPE_BYTE);
        SampleModel sm = cm.createCompatibleSampleModel(width, height);
        DataBufferByte db = new DataBufferByte(buffer, width * height);
        WritableRaster raster = Raster.createWritableRaster(sm, db, null);
        BufferedImage result = new BufferedImage(cm, raster, false, null);

        return result;
    }

    public static class RayTracerException extends Exception {
        public RayTracerException(String msg) {  super(msg); }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     *
     * finding the intersection surface of the Ray in the scene
     * Return Intersection object which includes the hit point and the hit surface
     * if no there is no intersection at all, then return null
     *
     */
    public Intersection getIntersection(Ray ray) {
        Vector tempHit, hitPoint = null;
        Vector p0 = ray.getStart();
        Surface hitSurface = null;
        double tempDist, minDist = Double.POSITIVE_INFINITY;

        // searching for minimum intersection point
        for (Surface surface : scene.getSurfaces()) {
            tempHit = surface.findIntersection(ray);
            if (tempHit != null) {
                tempDist = tempHit.distanceTo(p0);
                if (tempDist < minDist || hitSurface == null) {
                    minDist = tempDist;
                    hitSurface = surface;
                    hitPoint = new Vector(tempHit);
                }
            }
        }
        if (hitSurface == null) {       // no intersection...
            return null;
        }
        else {
            return new Intersection(hitSurface, hitPoint, ray);
        }
    }

    public Intersection getIntersection(Ray ray, Surface originSurface) {
        Vector tempHit, hitPoint = null;
        Vector p0 = ray.getStart();
        Surface hitSurface = null;
        double tempDist, minDist = Double.POSITIVE_INFINITY;

        // searching for minimum intersection point
        for (Surface surface : scene.getSurfaces()) {
            tempHit = surface.findIntersection(ray);
            if (tempHit != null && originSurface != originSurface) {
                tempDist = tempHit.distanceTo(p0);
                if (tempDist < minDist || hitSurface == null) {
                    minDist = tempDist;
                    hitSurface = surface;
                    hitPoint = new Vector(tempHit);
                }
            }
        }
        if (hitSurface == null) {       // no intersection...
            return null;
        }
        else {
            return new Intersection(hitSurface, hitPoint, ray);
        }
    }

    private Ray buildLightSourceRay(Light sourceLight, Intersection intersection) {
        Vector dirVector;

        // TODO: should switch direction?
        dirVector = (intersection.getPoint().minus(sourceLight.getPosition())).direction();
        Ray lightRay = new Ray(sourceLight.getPosition(), dirVector);

        return lightRay;
    }

    /**
     *
     * Getting the color with recursion method
     *
     */
    private Color getColor(Ray ray, List<Intersection> intersections, int listIndex, int recIndex) {
        if (intersections == null) {
            return new Color(this.scene.getmBackGroundColor());
        }
        // if all the intersections calculated, or no more recursions, then return scene background color
        if (listIndex >= intersections.size() || recIndex >= this.rec_max) {
            return new Color(this.scene.getmBackGroundColor());
        }

        // TODO: check transparency color (background color)
        // calculating the background color
        Color bgColor = getColor(ray, intersections, listIndex + 1, recIndex + 1);

        Intersection currIntersection = intersections.get(listIndex);
        Vector N = currIntersection.getNormal();
        Vector rayDir = ray.getDirection();
        //set normal towards source position
        if (rayDir.dot(N) > 0) {
            N = N.scale(-1);
        }
        Vector mirror = (rayDir.proj(N).scale(-2).plus(rayDir)).direction();
        // calculating all lights impact + shadow = (diffuse+specular) on CURRENT intersection
        Color pixelColor = getColorBySurface(currIntersection, N, mirror);

        // TODO: (reflection color)
        Color reflectionColor;
        Material material = this.scene.getMaterials().get(currIntersection.getSurface().getMaterialIndex() - 1);
        // if the material has reflection
        if (!material.getRefl().isZeroVector()) {
            // constructing the reflected ray
            Ray reflected = new Ray(currIntersection.getPoint(), mirror);
            // getting all the hit surfaces from this ray
            List<Intersection> reflectedIntersecteds = findAllIntersectionOnRay(reflected);
            // calculating the reflected color
            reflectionColor = getColor(reflected, reflectedIntersecteds, 0, recIndex+1);
            reflectionColor.multColor(new Color(material.getRefl()));
        }
        // if the material doesn't have reflection, then set reflected color to zero
        else {
            reflectionColor = new Color();
        }

        // (transparency) value
        double trValue = this.scene.getMaterials().get(currIntersection.getSurface().getMaterialIndex() - 1).getTrans();


        if (pixelColor == null) {   // TODO: check this option (DEBUG)
            return (new Color(this.scene.getmBackGroundColor()));
        }
        else {
            /**
             * output color = (background color)*transparency
             *                + (diffuse+specular)*(1-transparency)
             *                + (reflection color)
             */
            bgColor.multipleByScalar(trValue);              // (background color)*transparency
            pixelColor.multipleByScalar(1-trValue);         // (diffuse+specular)*(1-transparency)
            pixelColor.addColor(pixelColor);
            pixelColor.addColor(reflectionColor);           // + (reflection color)

            return pixelColor;
        }
    }

    // intersection parameter is from pixel ray!!
    private Color getColorBySurface(Intersection intersection, Vector N, Vector mirror) {
        Surface hitSurface = intersection.getSurface();
        Material hitMaterial = this.scene.getMaterials().get(hitSurface.getMaterialIndex() - 1);
        double shadowValue;
        Color totalColor;

        //Vector specularColor = new Vector(0, 0, 0);
        Vector specularColor, diffColor;
        double diffuseVal, specularVal;

        // iterating over all scene lights
        for (Light light: this.scene.getLights()) {
            specularColor = new Vector(0, 0, 0);

            // constructing the light ray to hit point
            Ray lightRay = buildLightSourceRay(light, intersection);
            diffuseVal = this.calculateDiffuseVal(lightRay, N);                    // (N*L)
            specularVal = this.calculateSpecularVal(intersection, mirror);         // (V*R)^n (n == phong coef)

            if (diffuseVal <= 0) {    // if the light is behind the object
                // TODO: check what to return in that case
                diffColor = new Vector(0,0,0);
                continue;
                //return (new Color(diffColor));
            }
            else {                  // if light in front of the object
                if (!hitMaterial.getDiff().isZeroVector() || !hitMaterial.getSpec().isZeroVector()) {
                    diffColor = (hitMaterial.getDiff()).scale(diffuseVal);          // K_d * (N*L)
                    if (specularVal > 0) {
                        specularColor = hitMaterial.getSpec().scale(light.getSpec());
                        specularColor = specularColor.scale(specularVal);           // K_s * (V*R)^n
                    }

                    // TODO: this is the shadow calculation (Nataly's):
                    /**
                     lightScreen = new Screen(L, light.position, light.lightWidth, light.lightWidth, this.shadowRays, this.shadowRays, 1);
                     p = lightScreen.howMuchIlluminated(hitPoint);
                     //						if(p==0) p=(1-light.shadowIntensity);
                     IL = light.color.cloneVector().multiplyBy((1+light.shadowIntensity*(p-1)));
                     totalColor = IL.multiplyColor(difColor.add(specularColor)).add(totalColor);
                     */

                    // TODO: this is the shadow calculation (Itai's):
                    // getting the area light grid according number of shadow rays
                    Light[] lightsGrid = light.getAreaLight(camera.getUp(), intersection.getPoint(), this.sh_rays);
                    // computing the shadow value at hit point
                    shadowValue = light.computeSoftShadow(lightsGrid, intersection, this.scene);
                    Color IL = new Color(light.getColor());
                    IL.multipleByScalar(shadowValue);
                    /////

                    totalColor = new Color(IL);
                    totalColor.multColor(new Color(diffColor.plus(specularColor)));
                    return totalColor;
                }
            }
        }
        return null;
    }

    private double calculateDiffuseVal(Ray lightRay, Vector N) {
        // switching hit direction vector...
        Vector L = (lightRay.getDirection()).scale(-1);
        double cosTeta = N.dot(L);

        return cosTeta;
    }

    private double calculateSpecularVal(Intersection intersection, Vector mirror) {
        // switching hit direction vector...
        Vector V = (intersection.getHitRay().getDirection()).scale(-1);

        double cosTeta = V.dot(mirror);
        double phong = this.scene.getMaterials().get(intersection.getSurface().getMaterialIndex() - 1).getPhong();

        double specularVal = Math.pow(cosTeta, phong);

        return specularVal;
    }

    /**
     private void superSamplingMehod (List<Ray> allRayTroughPixel, int i, int j){

     int superSamplingLvl = allRayTroughPixel.size();

     for (Ray ray : allRayTroughPixel) {
     Intersection hit = getIntersection(ray);

     if (hit == null) {       // no intersection, setting background color
     this.fillpixelColor(bg_color, i, j);
     }
     else {
     // get color of pixel (i,j) using rbgData
     Color pixelColor = this.getColor(hit);
     this.fillpixelColor(pixelColor.getRgbValues(), i, j);
     }
     }

     this.AVGPixelVal(superSamplingLvl, i, j);
     }
     */

    private void AVGPixelVal(int superSamplingLvl, int i, int j) {
        //input : superSamplingLvl := # rays shoot each pixel, i := row Index, j := col Index ;
        //Do    : divide each color component (r comp, g comp, and b comp) by superSampling ;

    }


    private void fillpixelColor(byte[] rgbData, Vector rgbValues, int i, int j) {
        // fix heights (if we pass the value 1)
        double red = Math.min(1, rgbValues.cartesian(0));
        double green = Math.min(1, rgbValues.cartesian(1));
        double blue = Math.min(1, rgbValues.cartesian(2));


        rgbData[(j * this.imageWidth + i) * 3] = (byte) ((int) (255 * red));
        rgbData[(j * this.imageWidth + i) * 3 + 1] = (byte) ((int) (255 * green));
        rgbData[(j * this.imageWidth + i) * 3 + 2] = (byte) ((int) (255 * blue));
    }

    public List<Intersection> findAllIntersectionOnRay(Ray ray) {
        List<Intersection> allIntersectionOnRay = new ArrayList<>();
        List<Surface> allSurfaces = new ArrayList<>();
        //Collections.copy(allSurfaces, scene.getSurfaces());
        allSurfaces.addAll(scene.getSurfaces());

        boolean stopCondition = true;
        Intersection hit;

        while (stopCondition) {
            hit = getIntersection(ray, allSurfaces);
            if (hit == null ) {
                break;
            }
            Material material = scene.getMaterials().get(hit.getSurface().getMaterialIndex() - 1);

            if (hit != null && allIntersectionOnRay.size() == 0 && material.getTrans() == 0) {
                allIntersectionOnRay.add(hit);
                stopCondition = false;

            }
            else if (hit == null || material.getTrans() == 0) {
                stopCondition = false;
            }
            else {
                allSurfaces.remove(hit.getSurface());
                allIntersectionOnRay.add(hit);
            }
        }
        Comparator<Intersection> cmp = new Comparator<Intersection>() {
            @Override
            public int compare(Intersection o1, Intersection o2) {
                if(camera.getPosition().distanceTo(o1.getPoint()) > camera.getPosition().distanceTo(o2.getPoint())){
                    return 1;
                }else if (camera.getPosition().distanceTo(o1.getPoint()) == camera.getPosition().distanceTo(o2.getPoint())){
                    return 0;
                }
                return -1;
            }
        };

        Collections.sort(allIntersectionOnRay, cmp);
        return allIntersectionOnRay;
    }

    public Intersection getIntersection(Ray ray, List<Surface> allSurface) {
        Vector tempHit, hitPoint = null;
        Vector p0 = ray.getStart();
        Surface hitSurface = null;
        double tempDist, minDist = Double.POSITIVE_INFINITY;

        // searching for minimum intersection point
        for (Surface surface : allSurface) {
            tempHit = surface.findIntersection(ray);
            if (tempHit != null) {
                tempDist = tempHit.distanceTo(p0);
                if (tempDist < minDist || hitSurface == null) {
                    minDist = tempDist;
                    hitSurface = surface;
                    hitPoint = new Vector(tempHit);
                }
            }
        }
        if (hitSurface == null) {       // no intersection...
            return null;
        }
        else {
            return new Intersection(hitSurface, hitPoint, ray);
        }
    }

    // GETTERS:
    public int getImageWidth() {
        return imageWidth;
    }

    public int getImageHeight() {
        return imageHeight;
    }

    public Vector getBg_color() {
        return bg_color;
    }

    public double getSh_rays() {
        return sh_rays;
    }

    public int getRec_max() {
        return rec_max;
    }

    public int getSs_level() {
        return ss_level;
    }

    public Camera getCamera() {
        return camera;
    }

}

