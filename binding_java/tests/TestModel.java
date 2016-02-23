import nb.geometricBot.*;

public class TestModel {
    public static void main(String[] args) {
	int N_success = 0;
	int N_total = 6;
	checkCombine();
	checkUnify();
	checkIntersect();
	checkSubstractA();
	checkSubstractB();
	boolean ok = checkVerify();
	if (ok)
	    N_success += 1;
	else
	    System.out.println(" > Test fails: 'Check verify()'");
	
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
    
    public static boolean checkVerify() {
	Model model = GeometricBot.getCircle(0, 0, 7, 1);
	ModelStatus status = model.verify();
	return status == ModelStatus.OK;
    }
}
