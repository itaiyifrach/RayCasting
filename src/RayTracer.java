//package RayTracing;

import java.awt.*;
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

    private Camera camera;
    private Scene scene;
    public int imageWidth;
    public int imageHeight;
    private Vector bg_color;    // background color (r, g, b)
    private int sh_rays;     // root number of shadow rays (N^2 rays will be shot)
    private int rec_max;        // maximum number of recursions
    private int ss_level;       // super sampling level

    byte[] rgbData = new byte[this.imageWidth * this.imageHeight * 3];

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
        Color pixelColor;
        for (int i = 0; i < this.imageWidth; i++) {
            for (int j = 0; j < this.imageHeight; j++) {

                List<Ray> mListRays = Ray.constructRayThroughPixel(camera, i, j, imageWidth, imageHeight, pixelWidth, pixelHeight, this.getSs_level());
                // find intersection and find the closest intersection
                this.superSamplingMehod(mListRays, i,j);
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

    private void AVGPixelVal(int superSamplingLvl, int i, int j) {
        //input : superSamplingLvl := # rays shoot each pixel, i := row Index, j := col Index ;
        //Do    : divide each color component (r comp, g comp, and b comp) by superSampling ;

    }


    private void fillpixelColor(Vector rgbValues, int i, int j){

        rgbData[(j * this.imageWidth + i) * 3] = (byte) ((int) (255 * (rgbValues.cartesian(0))));
        rgbData[(j * this.imageWidth + i) * 3 + 1] = (byte) ((int) (255 * (rgbValues.cartesian(1))));
        rgbData[(j * this.imageWidth + i) * 3 + 2] = (byte) ((int) (255 * (rgbValues.cartesian(2))));
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


    private Vector getRayColor(Vector ray, Vector source, List<Obj> intersectedObjects, int transIndex, int recIndex) {

//		if(DEBUG){
//			if(intersectedObjects.size()==0){
//				return this.backgroundColor;
//			}
//			else{
//				if(false) System.out.println("OOOO t="+intersectedObjects.get(0).t);
////				double level = 1-(intersectedObjects.get(0).t)/20;
////				double level = 5/(intersectedObjects.get(0).t+1);
////				return new Vector(level, level, level);
//				return list_materials.get(intersectedObjects.get(0).matIndex).diffuseColor;
//			}
//		}


        //if all elements in List was handled - return background
        if(transIndex >= intersectedObjects.size() || recIndex>=this.recurtionDepth){
            return this.backgroundColor;
        }

        //there's a hit
        else{
            Vector colorBehind = getRayColor(ray, source, intersectedObjects, transIndex+1, recIndex+1);

            Obj intersected = intersectedObjects.get(transIndex);
            Material material = list_materials.get(intersected.matIndex);
            Vector hitPoint = ray.cloneVector().multiplyBy(intersected.t-eps).add(source);
            Vector N = intersected.getNormal(hitPoint);
            //set normal towards source position
            if(Vector.dotProduct(ray,N)>0) N.multiplyBy(-1);
            Vector mirror = ray.projection(N).multiplyBy(-2).add(ray).normalize();
            Vector L;
            double p;
            double cos;
            Vector specularColor;
            Vector difColor;
            Vector IL;
            Screen lightScreen;
            Vector totalColor = new Vector(0, 0, 0);
            for (Light light: this.list_lights){
                L = Vector.fromPoints(hitPoint, light.position).normalize();
                specularColor = new Vector(0, 0, 0);
                if(Vector.dotProduct(L, N)<=0){
                    //if light is behind object
//					difColor = material.diffuseColor.multiplyBy(-Vector.dotProduct(N, L));
                    difColor = new Vector(0,0,0); //TODO ??
//					p=0;
                }
                else{
                    //if light in front of object
                    if(!material.diffuseColor.isZero() || !material.specularColor.isZero()){
                        difColor = material.diffuseColor.cloneVector();
                        difColor.multiplyBy(Vector.dotProduct(N, L));
                        cos = Vector.dotProduct(L, mirror);
                        if(cos>0){
                            specularColor = material.specularColor.cloneVector().multiplyBy(light.specularIntensity);
                            specularColor.multiplyBy(Math.pow(cos,material.phongCoef));
                        }
                        lightScreen = new Screen(L, light.position, light.lightWidth, light.lightWidth, this.shadowRays, this.shadowRays, 1);
                        p = lightScreen.howMuchIlluminated(hitPoint);
//						if(p==0) p=(1-light.shadowIntensity);
                        IL = light.color.cloneVector().multiplyBy((1+light.shadowIntensity*(p-1)));
                        totalColor = IL.multiplyColor(difColor.add(specularColor)).add(totalColor);
                    }
                }
//				IL = light.color.multiplyBy(Math.max(1-light.shadowIntensity, p));
                //TODO check IL correctness
//				IL = light.color.multiplyBy(p);


            }
//			totalColor.fixHighs(); //TODO
            //here

            //reflection
            Vector reflectionColor;
            if( !material.reflectionColor.isZero()){
                List<Obj> reflected_intersectedObjects = findIntersectedObjects(mirror,hitPoint,0);
                reflectionColor = getRayColor(mirror, hitPoint, reflected_intersectedObjects, 0, recIndex+1);
                reflectionColor = reflectionColor.cloneVector().multiplyColor(material.reflectionColor);
            }
            else reflectionColor = new Vector(0, 0, 0);

            double tr = material.transparencyValue;
            Vector res = colorBehind.cloneVector().multiplyBy(tr).add(totalColor.multiplyBy(1-tr)).add(reflectionColor);

            res.fixHighs();
            return res;
        }
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

    public Color getTheColor(Intersection intersection, int maxRecLvl){

        if(maxRecLvl == 0 || intersection.get){
            return new Color();
        }
        Color pixelColor = new Color();



        float transperencyVal = intersection.getSurface().getMaterial().getTrans();
        Color reflectionColor = new Color(intersection.getSurface().getMaterial().getRefl());

        Color colorByDirectedIllumination = new Color();

        for (Light lightSource : this.scene.getLights()) {
            int maxRecPerLight = this.rec_max;

            colorByDirectedIllumination = calculateColorByDirectedIllumination(intersection, lightSource);

            if(transperencyVal > 0){

            }

        }
        pixelColor.addColor(colorByDirectedIllumination);

        return null;
    }

    private Color calculateColorByDirectedIllumination(Intersection intersection,  Light lightSource){
        Ray lightRay = buildRayByPoints(intersection.getPoint(), lightSource.getPosition());

        Vector diffuseCoefficient = intersection.getSurface().getMaterial().getDiff();
        double diffuseVal = this.calculateDiffuseVal(intersection, lightRay);

        Vector specularCoefficient = intersection.getSurface().getMaterial().getSpec();
        double specularVal = this.calculateSpecularVal(intersection, lightRay);

        diffuseCoefficient = diffuseCoefficient.scale(diffuseVal);      // = K_d*(N.dot(L))
        specularCoefficient = specularCoefficient.scale(specularVal);   // = K_s*((V.dot(R))^n)

        int shadowIndicator = ShadowTermIndicator(lightSource.getPosition(), intersection.getPoint());

        Color lightSourceColor = lightSource.getColor();

        Color diffuseAndSpecularColorSum = new Color(diffuseCoefficient.plus(specularCoefficient));

        Color colorValueByDirectedIllumination = Color.multipleColor(diffuseAndSpecularColorSum, lightSourceColor.multipleByScalar(shadowIndicator)); // = I
        return colorValueByDirectedIllumination;
    }


    private int ShadowTermIndicator(Vector lightSourcePoint, Vector hitPoint){
        Ray inverseLightRay = buildRayByPoints(lightSourcePoint, hitPoint);
        Intersection lightSourceIntersection = this.getIntersection(inverseLightRay);
        if(hitPoint.isEqualTo(lightSourceIntersection.getPoint())){
            return 1;
        }
        return 0;
    }

    private Ray buildLightSourceRay(Light sourceLight, Intersection intersection){
        Vector dirVector;

        dirVector = (intersection.getPoint().minus(sourceLight.getPosition())).direction();
        Ray lightRay = new Ray(sourceLight.getPosition(), dirVector);

        return lightRay;
    }

    private Ray buildRayByPoints(Vector fromPoint, Vector toPoint){

        Vector dirVector = toPoint.minus(fromPoint).direction();
        Ray lightRay = new Ray(fromPoint, dirVector);
        return  lightRay;
    }

    public Color getColor(Intersection intersection) {
        int maxRecForLight = rec_max;
        Color pixelColor = new Color();
        Surface surface = intersection.getSurface();

        //System.out.println("hey from getColorMain method, intersection point point is :=" + intersection.getPoint().toString());
        for (Light sourceLight : this.scene.getLights()) {
            //System.out.println("hey from getColorMain method, source_light point is :=" + sourceLight.getPosition().toString());
            // constructing the light ray to hit point, and finding first intersection
            Ray lightRay = buildLightSourceRay(sourceLight, intersection);
            Intersection lightSourceIntersection = getIntersection(lightRay);

            // if the hit surface is the checked surface from the pixel ray, then the light hits the point directly
            // in that case, we calculate the color value it induced on the surface
            // otherwise, there is no effect from that light (but there is a shadow effect), we continue to the next one
            if (lightSourceIntersection.getSurface() == surface) {
                Color colorValPerLight = this.getColor(sourceLight, intersection, maxRecForLight);
                pixelColor.addColor(colorValPerLight);
            }
        }

        return pixelColor;
    }

    private Color getColor(Light light, Intersection intersection, int maxRecLevel) {
        if (maxRecLevel == 0) {
            return new Color();
        }

        // if the surface isn't null
        if (intersection != null) {
            Material material = this.scene.getMaterials().get(intersection.getSurface().getMaterialIndex() - 1);
            // calculate the diffuse and specular values of this surface...
            Color colorValBySurface = getColorBySurface(light, intersection, material);

            // if the surface is transference, the we continue to search background surfaces
            if (material.getTrans() > 0){
                // construct ray with same direction, but the start point is the hit point
                Ray transRay = new Ray(intersection.getPoint(), intersection.getHitRay().direction);
                Intersection transIntersection = getIntersection(transRay, intersection.getSurface());

                // if the surface isn't null
                if (transIntersection != null) {
                    maxRecLevel = maxRecLevel - 1;
                    Color transColor = this.getColor(light, transIntersection, maxRecLevel);
                    transColor.multypleByScalar(material.getTrans());
                    colorValBySurface.addColor(transColor);
                }
                // 1. build new ray from most far away hit point with same direction ;
                // 2. find inter section to this ray ;
                // 3. if we have inter section then
                // 4.

            }

            // if the surface is reflective, the we we continue to search by reflective ray
            if (material.isReflectence()) {
                Ray reflectionRay = intersection.getReflectionRay();
                Intersection reflectionIntersection = getIntersection(reflectionRay, intersection.getSurface());

                if (reflectionIntersection != null) {
                    maxRecLevel = maxRecLevel - 1;
                    Color reflectionColor = this.getColor(light, reflectionIntersection, maxRecLevel);
                    reflectionColor.multColor(new Color(material.getRefl()));
                    colorValBySurface.addColor(reflectionColor);
                }
            }
            return colorValBySurface;
        }
        else {
            return new Color(scene.getmBackGroundColor());
        }
    }

    /**
     *
     * (N.dot(L))
     *
     *
     * @param intersection
     * @param lightRay
     * @return
     */
    private double calculateDiffuseVal(Intersection intersection, Ray lightRay){
        Vector N = intersection.getNormal();
        Vector L = lightRay.getDirection();
        double cosTeta = N.dot(L);

        return cosTeta;
    }

    /**
     *
     * (V.dot(R))^n
     *
     * @param intersection
     * @param lightRay
     * @return
     */
    private double calculateSpecularVal(Intersection intersection, Ray lightRay){
        Vector V = intersection.getHitRay().getDirection().scale(-1);

        Intersection lightInter = new Intersection(intersection.getSurface(), intersection.getPoint(), lightRay);
        Vector R = lightInter.getReflectionRay().getDirection();

        double cosTeta = V.dot(R);
        double phong_specularity_coefficient = intersection.getSurface().getMaterial().getPhong();

        double specularVal = Math.pow(cosTeta, phong_specularity_coefficient);

        return specularVal;
    }

    // intersection parameter is from pixel ray!!
    private Color getColorBySurface(Light light, Intersection intersection, Material material) {
        Surface hitSurface = intersection.getSurface();
        Material hitMaterial = this.scene.getMaterials().get(hitSurface.getMaterialIndex() - 1);
        double cos, dotProduct, shadowValue;
        Color totalColor;

        //Vector specularColor = new Vector(0, 0, 0);
        Vector specularColor, diffColor;

        double diffuseVal, specularVal;

        // iterating over all scene lights
        for (Light light: this.scene.getLights()) {
            specularColor = new Vector(0, 0, 0);

            // constructing the light ray to hit point
            Ray lightRay = buildLightSourceRay(light, intersection);
            diffuseVal = this.calculateDiffuseVal(intersection, lightRay);
            specularVal = this.calculateSpecularVal(intersection,lightRay);

            if (diffuseVal <= 0) {    // if the light is behind the object
                // TODO: check what to return in that case
                diffColor = new Vector(0,0,0);
                return (new Color(diffColor));
            }
            else {                  // if light in front of the object
                if (!hitMaterial.getDiff().isZeroVector() || !hitMaterial.getSpec().isZeroVector()) {
                    diffColor = (hitMaterial.getDiff()).scale(diffuseVal);

                    if (specularVal > 0) {
                        specularColor = hitMaterial.getSpec().scale(light.getSpec());
                        specularColor = specularColor.scale(Math.pow(cos, hitMaterial.getPhong()));
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
                    Light[] lightsGrid = light.getAreaLight(intersection.getPoint(), this.sh_rays);
                    // computing the shadow value at hit point
                    shadowValue = light.computeSoftShadow(lightsGrid, intersection, this.scene);
                    Vector illuminateLight = light.getColor().scale(1 + (light.getShadow()*(shadowValue - 1)));
                    Color lightColor = new Color(illuminateLight);
                    /////

                    totalColor = new Color(lightColor);
                    totalColor.multColor(new Color(diffColor.plus(specularColor)));
                    return totalColor;
                }
            }
        }

        return null;
    }

    private double getShadowValue(Intersection hit) {
        double result = 0, temp;

        for (Light light: this.scene.getLights()) {
            // getting the area light grid according number of shadow rays
            Light[] lightsGrid = light.getAreaLight(hit.getPoint(), this.sh_rays);
            // computing the shadow value at hit point
            temp = light.computeSoftShadow(lightsGrid, hit, this.scene);
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

}

