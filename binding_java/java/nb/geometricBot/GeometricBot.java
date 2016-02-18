/**
 * Geometric Bot: Geometric tesselations for numerical analysis.
 * 2011-2015 Victor Eduardo Cardoso Nungaray
 * Twitter: @victore_cardoso
 * email: victorc@cimat.mx
 *
 */

package nb.geometricBot;

public class GeometricBot {
    private GeometricBot() {}

    static {
	System.loadLibrary("nb_geometric_bot-jni");
    }
    
    public static native Mesh generateMesh(Model model);
    private static native int JNICreateMesh
	(int NVtx, float vtx[], int NEdges, int edge[],
	 int NHoles, float hole[],
	 int maxNVtx, int maxNElem, float minAngle,
	 float maxEdgeLength, float maxSgmLength,
	 boolean includeEdges, boolean includeElements,
	 boolean includeElemAdj, boolean includeInputVtx,
	 boolean includeInputSgm,
	 float meshVtx[][], int meshConnEdge[][],
	 int meshConnMtx[][], int meshConnAdj[][],
	 int meshModelVtx[][], int meshModelSgm[][][]);

    private static native int JNICreateMesh
	(int NVtx, float vtx[], int NEdges, int edge[],
	 int NHoles, float hole[], String imageFile,
	 int maxNVtx, int maxNElem, float minAngle,
	 float scale, float xDisp, float yDisp,
	 boolean includeEdges, boolean includeElements,
	 boolean includeElemAdj, boolean includeInputVtx,
	 boolean includeInputSgm,
	 float meshVtx[][], int meshConnEdge[][],
	 int meshConnMtx[][], int meshConnAdj[][],
	 int meshModelVtx[][], int meshModelSgm[][][]);

    public static Box getEnvelopingBox(float vertices[])
    {
	float xMin = vertices[0];
	float yMin = vertices[1];
	float xMax = vertices[0];
	float yMax = vertices[1];

	for (int i = 1; i < vertices.length / 2; i++) {
	    if (vertices[i * 2] < xMin)
		xMin = vertices[i * 2];
	    else if (vertices[i * 2] > xMax)
		xMax = vertices[i * 2];

	    if (vertices[i*2+1] < yMin)
		yMin = vertices[i*2+1];
	    else if (vertices[i*2+1] > yMax)
		yMax = vertices[i*2+1];
	}
	
	return new Box(xMin, yMin, xMax, yMax);
    }

    public static float getCenterAndZoom(int width, int height, 
					 Box box, float center[])
    {
	if(null == box)
	    return 0.0f;

	center[0] = (box.vA[0] + box.vB[0]) / 2.0f;
	center[1] = (box.vA[1] + box.vB[1]) / 2.0f;
	    
	float zoom = width / (box.vB[0] - box.vA[0]);
	if (zoom > height / (box.vB[1] - box.vA[1]))
	    zoom = height / (box.vB[1] - box.vA[1]);
	return zoom * 0.9f;

    }

    /* Public functions */
    public static float ANGLE_MAX = 0.46163958715250017309f; /* 26.45 degrees */
    public static float ANGLE_UNC = 0.0f;

    public static Model getCircle(float x, float y, float r, float sideLength)
    {
	float p = (float)(2.0 * Math.PI * r);

	if (sideLength <= 0.0f)
	    /* Prevent division by zero */
	    sideLength = p / 10;

	float n = p / sideLength;
	if (10 > n)
	    n = 10;
	return getPolygon(x, y, r, Math.round((float)n));
    }

    public static Model getRect(float x1, float y1, float x2, float y2)
    {	
	Model model = new Model();

	model.vertices = new float[8];
	model.edges = new int[8];
	model.holes = null;

	/* Set vertices */
	model.vertices[0] = x1;
	model.vertices[1] = y1;
	model.vertices[2] = x2;
	model.vertices[3] = y1;
	model.vertices[4] = x2;
	model.vertices[5] = y2;
	model.vertices[6] = x1;
	model.vertices[7] = y2;

	/* Set edges */
	model.edges[0] = 0;
	model.edges[1] = 1;
	model.edges[2] = 1;
	model.edges[3] = 2;
	model.edges[4] = 2;
	model.edges[5] = 3;
	model.edges[6] = 3;
	model.edges[7] = 0;

	return model;
    }

    public static Model getPolygon(float x, float y, float r, int n)
    {
	Model model = new Model();
	
	float angleStep = (float)(Math.PI * 2.0 / n);
	
	model.vertices = new float[2 * n];
	model.edges = new int[2 * n];
	model.holes = null;

	for (int i = 0; i < n; i++) {
	    float angle = i * angleStep;
	    model.vertices[i * 2] = x + r * (float)Math.cos(angle);
	    model.vertices[i*2+1] = y + r * (float)Math.sin(angle);
	    model.edges[i * 2] = i;
	    model.edges[i*2+1] = (i+1) % n;
	}

	return model;
    }

    public static Mesh createMesh(Model model,
				  int maxNVtx, int maxNElem, float minAngle,
				  float maxEdgeLength, float maxSgmLength,
				  boolean includeEdges, boolean includeElements,
				  boolean includeElemAdj, boolean includeInputVtx,
				  boolean includeInputSgm)
    {
	float meshVertices[][] = new float[1][];
	int meshConnEdges[][] = new int[1][];
	int meshConnMtx[][] = new int[1][];
	int meshConnAdj[][] = new int[1][];	
	int meshInputVtx[][] = new int[1][];
	int meshInputSgm[][][] = new int[1][][];

	int status = 
	    JNICreateMesh(model.getNVertices(), model.vertices,
			  model.getNEdges(), model.edges,
			  model.getNHoles(), model.holes,
			  maxNVtx, maxNElem, minAngle, 
			  maxEdgeLength, maxSgmLength,
			  includeEdges, includeElements,
			  includeElemAdj, includeInputVtx,
			  includeInputSgm, meshVertices, 
			  meshConnEdges, meshConnMtx,
			  meshConnAdj, meshInputVtx,
			  meshInputSgm);
	if (0 != status) 
	    return null;

	Mesh mesh = new Mesh();
	mesh.vertices = meshVertices[0];
	
	if (includeEdges)
	    mesh.connEdges = meshConnEdges[0];

	if (includeElements)
	    mesh.connMtx = meshConnMtx[0];

	if (includeElemAdj)
	    mesh.connAdj = meshConnAdj[0];

	if (includeInputVtx)
	    mesh.modelVtx = meshInputVtx[0];

	if (includeInputSgm)
	    mesh.modelSgm = meshInputSgm[0];

	return mesh;
    }

    public static Mesh createMesh(Model model)
    {
	return createMesh(model, 0, 0, ANGLE_MAX, 0.0f, 0.0f,
			  true, true, true, true, true);
    }

    public static Mesh createMesh(Model model, int maxNVtx)
    {
	return createMesh(model, maxNVtx, 0, ANGLE_MAX, 0.0001f, 0.0f,
			  true, true, true, true, true);
    }

    public static Mesh createMesh(Model model, float minAngle)
    {
	return createMesh(model, 0, 0, minAngle, 0.0f, 0.0f,
			  true, true, true, true, true);
    }

    public static Mesh createMesh(Model model, float maxEdgeLength, 
				  float maxSgmLength)
    {
	return createMesh(model, 0, 0, ANGLE_MAX,
			  maxEdgeLength, maxSgmLength,
			  true, true, true, true, true);	
    }

    public static Mesh createMesh(Model model, boolean includeEdges,
				  boolean includeElements, boolean includeElemAdj,
				  boolean includeInputVtx, boolean includeInputSgm)
    {
	return createMesh(model, 0, 0, ANGLE_MAX, 0.0f, 0.0f,
			  includeEdges, includeElements, includeElemAdj,
			  includeInputVtx, includeInputSgm);
	
    }

    public static Mesh createMesh(Model model, float maxEdgeLength, 
				  float maxSgmLength, boolean includeEdges,
				  boolean includeElements, boolean includeElemAdj,
				  boolean includeInputVtx, boolean includeInputSgm)
    {
	return createMesh(model, 0, 0, ANGLE_MAX,
			  maxEdgeLength, maxSgmLength,
			  includeEdges, includeElements, includeElemAdj,
			  includeInputVtx, includeInputSgm);
	
    }

    public static Mesh createMesh(Model model, String imageFile,
				  int maxNVtx, int maxNElem, float minAngle,
				  float scale, float xDisp, float yDisp,
				  boolean includeEdges, boolean includeElements,
				  boolean includeElemAdj, boolean includeInputVtx,
				  boolean includeInputSgm)
    {
	float meshVertices[][] = new float[1][];
	int meshConnEdges[][] = new int[1][];
	int meshConnMtx[][] = new int[1][];
	int meshConnAdj[][] = new int[1][];	
	int meshInputVtx[][] = new int[1][];
	int meshInputSgm[][][] = new int[1][][];

	int status = 
	    JNICreateMesh(model.getNVertices(), model.vertices,
			  model.getNEdges(), model.edges,
			  model.getNHoles(), model.holes,
			  imageFile,
			  maxNVtx, maxNElem, minAngle, 
			  scale, xDisp, yDisp,
			  includeEdges, includeElements,
			  includeElemAdj, includeInputVtx,
			  includeInputSgm, meshVertices, 
			  meshConnEdges, meshConnMtx,
			  meshConnAdj, meshInputVtx,
			  meshInputSgm);
	if (0 != status) 
	    return null;

	Mesh mesh = new Mesh();
	mesh.vertices = meshVertices[0];
	
	if (includeEdges)
	    mesh.connEdges = meshConnEdges[0];

	if (includeElements)
	    mesh.connMtx = meshConnMtx[0];

	if (includeElemAdj)
	    mesh.connAdj = meshConnAdj[0];

	if (includeInputVtx)
	    mesh.modelVtx = meshInputVtx[0];

	if (includeInputSgm)
	    mesh.modelSgm = meshInputSgm[0];

	return mesh;
    }

    public static Mesh createMesh(String imageFile,
				  int maxNVtx, int maxNElem, float minAngle,
				  float scale, float xDisp, float yDisp,
				  boolean includeEdges, boolean includeElements,
				  boolean includeElemAdj, boolean includeInputVtx,
				  boolean includeInputSgm)
    {
	float meshVertices[][] = new float[1][];
	int meshConnEdges[][] = new int[1][];
	int meshConnMtx[][] = new int[1][];
	int meshConnAdj[][] = new int[1][];	
	int meshInputVtx[][] = new int[1][];
	int meshInputSgm[][][] = new int[1][][];

	int status = 
	    JNICreateMesh(0, null, 0, null, 0, null,
			  imageFile,
			  maxNVtx, maxNElem, minAngle, 
			  scale, xDisp, yDisp,
			  includeEdges, includeElements,
			  includeElemAdj, includeInputVtx,
			  includeInputSgm, meshVertices, 
			  meshConnEdges, meshConnMtx,
			  meshConnAdj, meshInputVtx,
			  meshInputSgm);
	/* The zero in the number of vertices of the model indicates to the
	 * JNI-binding that the model must be auto-adjusted to the picture.
	 */

	if (0 != status) 
	    return null;

	Mesh mesh = new Mesh();
	mesh.vertices = meshVertices[0];
	
	if (includeEdges)
	    mesh.connEdges = meshConnEdges[0];

	if (includeElements)
	    mesh.connMtx = meshConnMtx[0];

	if (includeElemAdj)
	    mesh.connAdj = meshConnAdj[0];

	if (includeInputVtx)
	    mesh.modelVtx = meshInputVtx[0];

	if (includeInputSgm)
	    mesh.modelSgm = meshInputSgm[0];

	return mesh;
    }

    public static Mesh createMesh(Model model, String imageFile, int maxNVtx,
				  float scale, float xDisp, float yDisp)
    {
	return createMesh(model, imageFile, maxNVtx, 0, ANGLE_MAX,
			  scale, xDisp, yDisp,
			  true, false, false, false, false);
    }

    public static Mesh createMesh(String imageFile, int maxNVtx, float minAngle)
    {
	return createMesh(imageFile, maxNVtx, 0, minAngle, 1.0f, 0.0f, 0.0f,
			  true, false, false, false, false);
    }

    public static Mesh createMesh(String imageFile)
    {
	return createMesh(imageFile, 15000, 0.2f);
    }
}
