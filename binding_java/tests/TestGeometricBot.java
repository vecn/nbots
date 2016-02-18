import nb.geometricBot.*;

public class TestGeometricBot {
    public static void main(String[] args) {
	sample1();
	sample2();
	sample3();
	sample4();
	sample5();
    }
    public static void sample1() {
	/* Create the mesh with the minimum number of triangles,
	 * where the minimum angle permitted is 26.45 degrees.
	 */
	Model model = GeometricBot.getCircle(0.0f, 0.0f, 1.0f, 0.1f);
	Mesh mesh = GeometricBot.generateMesh(model);
	GeometricBotDraw.drawMesh(mesh, "sample1_mesh.png", 1000, 800);
    }
    public static void sample2() {
	/* Min angle constraint relaxed */
	Model model = GeometricBot.getCircle(0.0f, 0.0f, 1.0f, 0.1f);
	Mesh mesh = GeometricBot.createMesh(model, GeometricBot.ANGLE_MAX / 2.0f);
	GeometricBotDraw.drawMesh(mesh, "sample2_mesh.png", 1000, 800);
    }
    public static void sample3() {
	/* Max number of vertices constraint */
	Model model = GeometricBot.getCircle(0.0f, 0.0f, 1.0f, 0.1f);
	Mesh mesh = GeometricBot.createMesh(model, 70);
	GeometricBotDraw.drawMesh(mesh, "sample3_mesh.png", 1000, 800);
    }
    public static void sample4() {
	/* Control the mesh size by defining the edge and subsegment
	 * lengths. A subsegment is an edge of the mesh that is part 
	 * of the input segments.
	 */
	Model model = GeometricBot.getCircle(0.0f, 0.0f, 1.0f, 0.1f);
	Mesh mesh = GeometricBot.createMesh(model, 0.1f, 0.03f);
	GeometricBotDraw.drawMesh(mesh, "sample4_mesh.png", 1000, 800);
    }
    public static void sample5() {
	/* Control the mesh density using an image */
	Mesh mesh = GeometricBot.createMesh("../tests/java/GeometricBotInputs/color_eye.jpg");
	GeometricBotDraw.drawMesh(mesh, "sample5_mesh.png", 1000, 800);
    }
}
