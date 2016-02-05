/**
 * Mesh: Mesh produced by the Geometric Bot.
 * 2011-2015 Victor Eduardo Cardoso Nungaray
 * Twitter: @victore_cardoso
 * email: victorc@cimat.mx
 *
 */

package vcn.geometricBot;

public class Mesh {
    protected float[] vertices; /* Concatenated Vertices coordinates */
    protected int[] connEdges;  /* Concatenated IDs of vtx forming edges */
    protected int[] connMtx;    /* Concatemated IDs of vtx forming elems */
    protected int[] connAdj;    /* Concatenated IDs of elems adjacent to elems.
				 * If the ID is out of bounds, then the element
				 * is part of the boundary; there is not an 
				 * adjacency.
				 */
    protected int[] modelVtx;   /* Concatenated IDs of mesh vtx corresponding
				 * to model vertices.
				 */
    protected int[][] modelSgm; /* An array of concatenated mesh vtx IDs for
				 * each input segment (of the model).
				 */

    protected Mesh()
    {
	vertices = null;
	connEdges = null;
	connMtx = null;
	connAdj = null;
	modelVtx = null;
	modelSgm = null;
    }


    /* Public methods */
    public int getNVertices()
    {
	if (null == vertices)
	    return 0;
	else
	    return vertices.length / 2;
	
    }

    public int getNEdges()
    {
	if (null == vertices)
	    return 0;
	else
	    return connEdges.length / 2;
    }

    public int getNElements()
    {
	if (null == connMtx)
	    return 0;
	else
	    return connMtx.length / 3;
    }

    public float[] getVerticesRef()
    {
	return vertices;
    }

    public int[] getConnMtxRef()
    {
	return connMtx;
    }

    public int[] getConnEdgesRef()
    {
	return connEdges;
    }

    public float[] getEdges()
    {
	if (null == connEdges)
	    return null;
	float[] edges = new float[connEdges.length];
	for (int i = 0; i < connEdges.length; i++)
	    edges[i] = vertices[connEdges[i]];
	return edges;
    }

    public int[] getConnAdjRef()
    {
	return connAdj;
    }

    public int[] getModelVtxRef()
    {
	return modelVtx;
    }

    public int[][] getModelSgmRef()
    {
	return modelSgm;
    }
}
