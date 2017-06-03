import java.util.ArrayList;

/**
 * Created by oshriamir on 6/3/17.
 */
public class Scene {


    private Vector mBackGroundColor;
    private int shadowRaysNumber;
    private int maxRecursionLevel;
    private int superSamplingLevel;

    private ArrayList<Surface> surfaces;
    private ArrayList<Material> materials;
    private ArrayList<Light> lights;

    public Scene(Vector mBackGroundColor, int shadowRaysNumber, int maxRecursionLevel, int superSamplingLevel) {

        this.mBackGroundColor = mBackGroundColor;
        this.shadowRaysNumber = shadowRaysNumber;
        this.maxRecursionLevel = maxRecursionLevel;
        this.superSamplingLevel = superSamplingLevel;

        this.surfaces = new ArrayList<>();
        this.materials = new ArrayList<>();
        this.lights = new ArrayList<>();
    }

    public Scene() {
        this.surfaces = new ArrayList<>();
        this.materials = new ArrayList<>();
        this.lights = new ArrayList<>();
    }

    public void buildScene(Vector mBackGroundColor, int shadowRaysNumber, int maxRecursionLevel, int superSamplingLevel) {

        this.mBackGroundColor = mBackGroundColor;
        this.shadowRaysNumber = shadowRaysNumber;
        this.maxRecursionLevel = maxRecursionLevel;
        this.superSamplingLevel = superSamplingLevel;

    }

    public Vector getmBackGroundColor() {
        return mBackGroundColor;
    }
    public void setmBackGroundColor(Vector mBackGroundColor) {
        this.mBackGroundColor = mBackGroundColor;
    }

    public int getShadowRaysNumber() {
        return shadowRaysNumber;
    }
    public void setShadowRaysNumber(int shadowRaysNumber) {
        this.shadowRaysNumber = shadowRaysNumber;
    }

    public int getMaxRecursionLevel() {
        return maxRecursionLevel;
    }
    public void setMaxRecursionLevel(int maxRecursionLevel) {
        this.maxRecursionLevel = maxRecursionLevel;
    }

    public int getSuperSamplingLevel() {
        return superSamplingLevel;
    }
    public void setSuperSamplingLevel(int superSamplingLevel) {
        this.superSamplingLevel = superSamplingLevel;
    }

    public ArrayList<Surface> getSurfaces() {
        return surfaces;
    }
    public void setSurfaces(ArrayList<Surface> surfaces) {
        this.surfaces = surfaces;
    }

    public ArrayList<Material> getMaterials() {
        return materials;
    }
    public void setMaterials(ArrayList<Material> materials) {
        this.materials = materials;
    }

    public ArrayList<Light> getLights() {
        return lights;
    }
    public void setLights(ArrayList<Light> lights) {
        this.lights = lights;
    }
}

