import static org.junit.Assert.assertEquals;
import org.junit.Test;

import nb.geometricBot.*;

public class TestModel {
    @Test
    public void checkCombine() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.combine(model1, model2);
	GeometricBotDraw.drawModel(model, "build/java_combine.png", 1000, 800);
	boolean passTest = true;
	assertEquals(true, passTest);
    }

    @Test
    public void checkUnify() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.unify(model1, model2);
	GeometricBotDraw.drawModel(model, "build/java_unify.png", 1000, 800);
	boolean passTest = true;
	assertEquals(true, passTest);
    }

    @Test
    public void checkIntersect() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.intersect(model1, model2);
	GeometricBotDraw.drawModel(model, "build/java_intersect.png",
				   1000, 800);
	boolean passTest = true;
	assertEquals(true, passTest);
    }

    @Test
    public void checkSubstractA() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.substract(model1, model2);
	GeometricBotDraw.drawModel(model, "build/java_substract_A.png",
				   1000, 800);
    }

    @Test
    public void checkSubstractB() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.substract(model2, model1);
	GeometricBotDraw.drawModel(model, "build/java_substract_B.png",
				   1000, 800);
	boolean passTest = true;
	assertEquals(true, passTest);
    }

    @Test
    public void checkDifference() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.difference(model2, model1);
	GeometricBotDraw.drawModel(model, "build/java_difference.png",
				   1000, 800);
	boolean passTest = true;
	assertEquals(true, passTest);
    }

    @Test
    public void checkVerify() {
	Model model = GeometricBot.getCircle(0, 0, 7, 1);
	ModelStatus status = model.verify();
	boolean passTest = (status == ModelStatus.OK);
	assertEquals(true, passTest);
    }

    @Test
    public void checkIsContinuum() {
	Model model = GeometricBot.getCircle(0, 0, 7, 1);
	boolean passTest = model.isContinuum();
	assertEquals(true, passTest);
    }

    @Test
    public void checkIsPointInside() {
	Model model = GeometricBot.getCircle(0, 0, 7, 1);
	boolean passTest = model.isPointInside(0, 0);
	assertEquals(true, passTest);
    }
}
