package nb.pdeBot.finiteElement.solidMechanics;

import nb.pdeBot.*;
import nb.pdeBot.finiteElement.*;
import nb.geometricBot.Model;

public class StaticElasticity2D {
    private StaticElasticity2D() {}

    static {
	System.loadLibrary("nb_bots");
	System.loadLibrary("nb_pde_bot-jni");
    }
    
    public static native MeshResults solve(Model model, Material material,
					   BoundaryConditions DirichletVtx,
					   BoundaryConditions NeumannVtx,
					   BoundaryConditions DirichletSgm,
					   BoundaryConditions NeumannSgm,
					   double thickness);
}
