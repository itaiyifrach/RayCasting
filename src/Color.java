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

    protected Color multypleByScalar(double scalar){
        this.setRgbValues(this.getRgbValues().scale(scalar));
        return this;
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
