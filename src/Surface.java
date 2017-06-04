/**
 * Created by Admin on 29/05/2017.
 */
public interface Surface {

    Vector findIntersection(Ray ray);

    Vector getNormalVector(Vector vec);

    int getMaterialIndex();

    Material getMaterial();

}
