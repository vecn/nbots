import vcn.geometricBot.*;
import vcn.pdeBot.*;
import vcn.pdeBot.finiteElement.*;
import vcn.pdeBot.finiteElement.solidMechanics.*;

public class TestPdeBot {
    public static void main(String[] args) {
	test1();
    }
    public static void test1() {
	Model model = GeometricBot.getRect(-3, 0, 3, 1);
	Material material = new Material(2.11e11, 0.3);
	BoundaryConditions bCondVtx = new BoundaryConditions(BCType.DIRICHLET);
	BoundaryConditions bCondSgm = new BoundaryConditions(BCType.DIRICHLET);
	MeshResults mesh =
	    StaticElasticity2D.solve(model, material, bCondVtx,
				     bCondSgm, 1.0);
    }
}
