

public class Material {
    private Vector diff;
    private Vector spec;
    private Vector refl;
    private float phong;
    private float trans;

    public Material(Vector diff, Vector spec, Vector refl, float phong, float trans) {
        this.diff = diff;
        this.spec = spec;
        this.refl = refl;
        this.phong = phong;
        this.trans = trans;
    }

    public Vector getDiff() {
        return diff;
    }

    public Vector getSpec() {
        return spec;
    }

    public Vector getRefl() {
        return refl;
    }

    public float getPhong() {
        return phong;
    }

    public float getTrans() {
        return trans;
    }

    public boolean isReflectence(){
        return !(refl.isZeroVector());
    }
}
