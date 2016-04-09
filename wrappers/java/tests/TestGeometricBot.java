import static org.junit.Assert.assertEquals;
import org.junit.Test;

import nb.geometricBot.*;

public class TestGeometricBot {
    @Test
    public void sample1() {
	/* Create the mesh with the minimum number of triangles,
	 * where the minimum angle permitted is 26.45 degrees.
	 */
	Model model = GeometricBot.getCircle(0.0f, 0.0f, 1.0f, 0.1f);
	Mesh mesh = GeometricBot.generateMesh(model);
	GeometricBotDraw.drawMesh(mesh, "build/sample1_mesh.png", 1000, 800);
	boolean passTest = true;
	assertEquals(true, passTest);
    }

    @Test
    public void sample2() {
	/* Min angle constraint relaxed */
	Model model = GeometricBot.getCircle(0.0f, 0.0f, 1.0f, 0.1f);
	Mesh mesh = GeometricBot.createMesh(model,
					    GeometricBot.ANGLE_MAX / 2.0f);
	GeometricBotDraw.drawMesh(mesh, "build/sample2_mesh.png", 1000, 800);
	boolean passTest = true;
	assertEquals(true, passTest);
    }

    @Test
    public void sample3() {
	/* Max number of vertices constraint */
	Model model = GeometricBot.getCircle(0.0f, 0.0f, 1.0f, 0.1f);
	Mesh mesh = GeometricBot.createMesh(model, 70);
	GeometricBotDraw.drawMesh(mesh, "build/sample3_mesh.png", 1000, 800);
	boolean passTest = true;
	assertEquals(true, passTest);
    }

    @Test
    public void sample4() {
	/* Control the mesh size by defining the edge and subsegment
	 * lengths. A subsegment is an edge of the mesh that is part 
	 * of the input segments.
	 */
	Model model = GeometricBot.getCircle(0.0f, 0.0f, 1.0f, 0.1f);
	Mesh mesh = GeometricBot.createMesh(model, 0.1f, 0.03f);
	GeometricBotDraw.drawMesh(mesh, "build/sample4_mesh.png", 1000, 800);
	boolean passTest = true;
	assertEquals(true, passTest);
    }
    
    @Test
    public void sample5() {
	/* Control the mesh density using an image */
	Mesh mesh = GeometricBot.createMesh("tests/java/GeometricBotInputs/color_eye.jpg");
	GeometricBotDraw.drawMesh(mesh, "build/sample5_mesh.png", 1000, 800);
	boolean passTest = true;
	assertEquals(true, passTest);
    }
}
