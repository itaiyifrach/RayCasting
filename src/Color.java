/**
 * Created by oshriamir on 6/3/17.
 */
public class Color {
    private Vector rgbValues;

    public Vector getRgbValues() {
        return rgbValues;
    }

    public Color(Vector rgbValues) {
        this.rgbValues = new Vector(rgbValues);
    }

    public Color() {
        this.rgbValues = new Vector(3);
    }

    public Color(Color otherColor) {
        this.rgbValues = new Vector(otherColor.getRgbValues());
    }

    public void setRgbValues(Vector rgbValues) {
        this.rgbValues = rgbValues;
    }

    protected void addColor(Color other){
        this.setRgbValues(getRgbValues().plus(other.getRgbValues()));
    }

    protected void multColor(Color other) {
        double[] multRes = new double[3];

        for (int i = 0; i < multRes.length; i++ ) {
            multRes[i] = this.rgbValues.cartesian(i) * other.rgbValues.cartesian(i);
        }

        this.setRgbValues(new Vector(multRes));
    }

    protected Color multipleByScalar(double scalar){
        this.setRgbValues(this.getRgbValues().scale(scalar));
        return this;
    }

    public static Color multipleColor(Color rgb_a, Color rgb_b){
        Vector rgb_1 = rgb_a.rgbValues;
        Vector rgb_2 = rgb_b.rgbValues;

        Vector result_rgb = Vector.multipleByCoordinates(rgb_1, rgb_2);
        return new Color(result_rgb);
    }

    protected Color diffAndSpecValues(Vector diffVector, Vector specVector){
        Color diffuseVal = new Color(diffVector);
        Color specularVal =  new Color(specVector);

        diffuseVal.addColor(specularVal);
        return new Color(diffuseVal);
    }

    public boolean isWhite(){
        return (this.getRgbValues().getData()[0] == 0 && this.getRgbValues().getData()[1] == 0 && this.getRgbValues().getData()[2] == 0);
    }
}
