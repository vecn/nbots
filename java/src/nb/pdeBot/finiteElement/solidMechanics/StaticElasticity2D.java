package nb.pdeBot.finiteElement.solidMechanics;

import nb.pdeBot.*;
import nb.pdeBot.finiteElement.*;
import nb.geometricBot.Model;

public class StaticElasticity2D {
    private StaticElasticity2D() {}

    static {
	System.loadLibrary("nb_pde_bot_jni");
    }
    
    public static native MeshResults solve(Model model, Material material,
					   Analysis2D analysis2D,
					   BoundaryConditions DirichletVtx,
					   BoundaryConditions NeumannVtx,
					   BoundaryConditions DirichletSgm,
					   BoundaryConditions NeumannSgm,
					   int N_nodes);
}
