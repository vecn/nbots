package vcn.pdeBot.finiteElement.solidMechanics;

import vcn.pdeBot.*;
import vcn.pdeBot.finiteElement.*;
import vcn.geometricBot.Model;

public class StaticElasticity2D {
    private StaticElasticity2D() {}

    static {
	System.loadLibrary("vcn_bots");
	System.loadLibrary("vcn_pde_bot-jni");
    }
    
    public static native MeshResults solve(Model model, Material material,
					   BoundaryConditions DirichletVtx,
					   BoundaryConditions NeumannVtx,
					   BoundaryConditions DirichletSgm,
					   BoundaryConditions NeumannSgm,
					   double thickness);
}