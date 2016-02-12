import nb.geometricBot.*;
import nb.pdeBot.*;
import nb.pdeBot.finiteElement.*;
import nb.pdeBot.finiteElement.solidMechanics.*;

public class TestPdeBot {
    public static void main(String[] args) {
	test1();
    }
    public static void test1() {
	Model model = GeometricBot.getRect(-3, 0, 3, 1);
	Material material = new Material(2.11e11, 0.3);
	BoundaryConditions DirichletVtx = getDirichletVtx();

	BoundaryConditions NewmannVtx = 
	    new BoundaryConditions(BCType.NEUMANN);

	BoundaryConditions DirichletSgm =
	    new BoundaryConditions(BCType.DIRICHLET);

	BoundaryConditions NewmannSgm =
	    new BoundaryConditions(BCType.NEUMANN);

	MeshResults mesh =
	    StaticElasticity2D.solve(model, material,
				     DirichletVtx, NewmannVtx,
				     DirichletSgm, NewmannSgm,
				     1.0, 200);
	MeshResultsDraw.draw(mesh, "TEMPORAL_Java.png", 1000, 800);
    }
    private static BoundaryConditions getDirichletVtx() {
	BoundaryConditions BC = new BoundaryConditions(BCType.DIRICHLET);
	int ids[] = new int[3];
	ids[0] = 0;
	ids[1] = 2;
	ids[2] = 3;
	
	char dof[] = new char[3];
	dof[0] = 'x';
	dof[1] = 'y';
	dof[2] = 'a';
       
	double val[] = new double[6];
	val[0] = 0;
	val[3] = -0.1;
	val[4] = 0;
	val[5] = 0;
	BC.set(ids, dof, val);
	return BC;
    }
}
