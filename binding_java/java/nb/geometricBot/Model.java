/**
 * Model: Definition of the Geometry. This class is used as input for most
 * of the Geometric Bot routines.
 * 2011-2015 Victor Eduardo Cardoso Nungaray
 * Twitter: @victore_cardoso
 * email: victorc@cimat.mx
 *
 */

package nb.geometricBot;

public class Model {
    protected float[] vertices;   /* Concatenated coordinates of vertices */
    protected int[]  edges;       /* Concatenated IDs of vertices forming edges */
    protected float[] holes;      /* Concatenated coordinates of holes */

    protected Model() {
	vertices = null;
	edges = null;
	holes = null;
    }

    static {
	System.loadLibrary("nb_geometric_bot-jni");
    }

    public static native Model combine(Model modelA, Model modelB);
    public static native Model unify(Model modelA, Model modelB);
    public static native Model intersect(Model modelA, Model modelB);
    public static native Model substract(Model modelA, Model modelB);
    public static native Model difference(Model modelA, Model modelB);
    
    public native ModelStatus verify();

    public int getNVertices()
    {
	if (null == vertices)
	    return 0;
	else
	    return vertices.length / 2;
	
    }

    public int getNEdges()
    {
	if (null == edges)
	    return 0;
	else
	    return edges.length / 2;
    }

    public int getNHoles()
    {
	if (hasHoles())
	    return holes.length / 2;
	else
	    return 0;
    }

    public boolean hasHoles()
    {
	return (null != holes);
    }

    public float[] getVerticesRef()
    {
	return vertices;
    }

    public int[] getEdgesRef()
    {
	return edges;
    }

    public float[] getHolesRef()
    {
	return vertices;
    }
}
