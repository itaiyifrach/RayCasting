import java.util.ArrayList;

public class Scene {
    private Vector mBackGroundColor;

    private ArrayList<Surface> surfaces;
    private ArrayList<Material> materials;
    private ArrayList<Light> lights;

    public Scene() {
        this.surfaces = new ArrayList<>();
        this.materials = new ArrayList<>();
        this.lights = new ArrayList<>();
    }


    public Vector getmBackGroundColor() {
        return mBackGroundColor;
    }

    public void setmBackGroundColor(Vector mBackGroundColor) {
        this.mBackGroundColor = mBackGroundColor;
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

