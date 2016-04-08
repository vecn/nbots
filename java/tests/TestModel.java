import nb.geometricBot.*;

public class TestModel {
    public static void main(String[] args) {
	int N_success = 0;
	int N_total = 9;
	checkCombine();
	checkUnify();
	checkIntersect();
	checkSubstractA();
	checkSubstractB();
	checkDifference();
	boolean ok1 = checkVerify();
	if (ok1)
	    N_success += 1;
	else
	    System.out.println(" > Test fails: 'Check verify()'");

	boolean ok2 = checkIsContinuum();
	if (ok2)
	    N_success += 1;
	else
	    System.out.println(" > Test fails: 'Check isContinuum()'");
	
	boolean ok3 = checkIsPointInside();
	if (ok3)
	    N_success += 1;
	else
	    System.out.println(" > Test fails: 'Check isContinuum()'");

	System.out.println("["+N_success+","+N_total+"] " + "Success");
    }

    public static void checkCombine() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.combine(model1, model2);
	GeometricBotDraw.drawModel(model, "java_combine.png", 1000, 800);
    }

    public static void checkUnify() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.unify(model1, model2);
	GeometricBotDraw.drawModel(model, "java_unify.png", 1000, 800);
    }

    public static void checkIntersect() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.intersect(model1, model2);
	GeometricBotDraw.drawModel(model, "java_intersect.png", 1000, 800);
    }

    public static void checkSubstractA() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.substract(model1, model2);
	GeometricBotDraw.drawModel(model, "java_substract_A.png", 1000, 800);
    }
    
    public static void checkSubstractB() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.substract(model2, model1);
	GeometricBotDraw.drawModel(model, "java_substract_B.png", 1000, 800);
    }
    
    public static void checkDifference() {
	Model model1 = GeometricBot.getCircle(-5, 0, 7, 1);
	Model model2 = GeometricBot.getCircle(5, 0, 7, 1);
	Model model = Model.difference(model2, model1);
	GeometricBotDraw.drawModel(model, "java_difference.png", 1000, 800);
    }

    public static boolean checkVerify() {
	Model model = GeometricBot.getCircle(0, 0, 7, 1);
	ModelStatus status = model.verify();
	return status == ModelStatus.OK;
    }

    public static boolean checkIsContinuum() {
	Model model = GeometricBot.getCircle(0, 0, 7, 1);
	return model.isContinuum();
    }

    public static boolean checkIsPointInside() {
	Model model = GeometricBot.getCircle(0, 0, 7, 1);
	return model.isPointInside(0, 0);
    }
}
