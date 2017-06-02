//package RayTracing;

import java.awt.Transparency;
import java.awt.color.*;
import java.awt.image.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

/**
 *  Main class for ray tracing exercise.
 */
public class RayTracer {

    public int imageWidth;
    public int imageHeight;
    private Vector bg_color;    // background color (r, g, b)
    private int sh_rays;     // root number of shadow rays (N^2 rays will be shot)
    private int rec_max;        // maximum number of recursions
    private int ss_level;       // super sampling level
    private Camera camera;
    private ArrayList<Surface> surfaces;
    private ArrayList<Material> materials;
    private ArrayList<Light> lights;


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

        this.materials = new ArrayList<>();
        this.lights = new ArrayList<>();
        this.surfaces = new ArrayList<>();

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

                    System.out.println(String.format("Parsed general settings (line %d)", lineNum));
                }
                else if (code.equals("mtl"))
                {
                    Material mtr = new Material(new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2])),
                            new Vector(Double.parseDouble(params[3]), Double.parseDouble(params[4]), Double.parseDouble(params[5])),
                            new Vector(Double.parseDouble(params[6]), Double.parseDouble(params[7]), Double.parseDouble(params[8])),
                            Float.parseFloat(params[9]), Float.parseFloat(params[10]));
                    this.materials.add(mtr);

                    System.out.println(String.format("Parsed material (line %d)", lineNum));
                }
                else if (code.equals("sph"))
                {
                    Sphere sph = new Sphere(new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2])),
                            Float.parseFloat(params[3]), Integer.parseInt(params[4]));
                    this.surfaces.add(sph);

                    System.out.println(String.format("Parsed sphere (line %d)", lineNum));
                }
                else if (code.equals("pln"))
                {
                    Plane pln = new Plane(new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2])),
                            Integer.parseInt(params[3]), Integer.parseInt(params[4]));
                    this.surfaces.add(pln);

                    System.out.println(String.format("Parsed plane (line %d)", lineNum));
                }
                else if (code.equals("trg"))
                {
                    Triangle trg = new Triangle(new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2])),
                            new Vector(Double.parseDouble(params[3]), Double.parseDouble(params[4]), Double.parseDouble(params[5])),
                            new Vector(Double.parseDouble(params[6]), Double.parseDouble(params[7]), Double.parseDouble(params[8])),
                            Integer.parseInt(params[9]));
                    this.surfaces.add(trg);

                    System.out.println(String.format("Parsed cylinder (line %d)", lineNum));
                }
                else if (code.equals("lgt"))
                {
                    Light lgt = new Light(new Vector(Double.parseDouble(params[0]), Double.parseDouble(params[1]), Double.parseDouble(params[2])),
                            new Vector(Double.parseDouble(params[3]), Double.parseDouble(params[4]), Double.parseDouble(params[5])),
                            Float.parseFloat(params[6]), Float.parseFloat(params[7]), Float.parseFloat(params[8]));
                    this.lights.add(lgt);

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
        Intersection hit;
        double imageRatio = imageHeight / imageWidth;
        double pixelWidth = camera.getScreenWidth() / imageWidth;
        double pixelHeight = imageRatio * pixelWidth;
        double shadowValue;
        for (int i = 0; i < this.imageWidth; i++) {
            for (int j = 0; j < this.imageHeight; j++) {
                Ray ray = Ray.constructRayThroughPixel(camera, i, j, imageWidth, imageHeight, pixelWidth, pixelHeight);
                // find intersection and find the closest intersection
                hit = getIntersection(ray);
                if (hit == null) {       // no intersection, setting background color
                    rgbData[(j * this.imageWidth + i) * 3] = (byte) ((int) (255 * bg_color.cartesian(0)));
                    rgbData[(j * this.imageWidth + i) * 3 + 1] = (byte) ((int) (255 * bg_color.cartesian(1)));
                    rgbData[(j * this.imageWidth + i) * 3 + 2] = (byte) ((int) (255 * bg_color.cartesian(2)));
                }
                else {
                    // get color of pixel (i,j) using rbgData
                    shadowValue = getShadowValue(hit);
                    int mat_idx = hit.getSurface().getMaterialIndex();
                    Material mat = materials.get(mat_idx - 1);
                    rgbData[(j * this.imageWidth + i) * 3] = (byte) ((int) (255 * (mat.getDiff().cartesian(0) * shadowValue)));
                    rgbData[(j * this.imageWidth + i) * 3 + 1] = (byte) ((int) (255 * (mat.getDiff().cartesian(1)  * shadowValue)));
                    rgbData[(j * this.imageWidth + i) * 3 + 2] = (byte) ((int) (255 * (mat.getDiff().cartesian(2))  * shadowValue));
                }
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
        for (Surface surface : surfaces) {
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
            return new Intersection(hitSurface, hitPoint);
        }
    }

    private double getShadowValue(Intersection hit) {
        double result = 0, temp;

        for (Light light: this.lights) {
            // getting the area light grid according number of shadow rays
            Light[] lightsGrid = light.getAreaLight(hit.getPoint(), this.sh_rays);
            // computing the shadow value at hit point
            temp = light.computeSoftShadow(lightsGrid, hit, this);
            // TODO: check how to calculate the shadow value on hit point, by all lights. I currently do multiplication
            if (result == 0) {
                result = temp;
            } else {
                result *= temp;
            }
        }

        return result;
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

    public ArrayList<Surface> getSurfaces() {
        return surfaces;
    }

    public ArrayList<Material> getMaterials() {
        return materials;
    }

    public ArrayList<Light> getLights() {
        return lights;
    }
}

