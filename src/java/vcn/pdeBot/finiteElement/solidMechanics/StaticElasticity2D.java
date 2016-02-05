package vcn.pdeBot.finiteElement.solidMechanics;

public class StaticElasticity2D {
    private StaticElasticity2D() {}

    static {
	System.loadLibrary("vcn_bots");
	System.loadLibrary("vcn_pde_bot-jni");
    }
    
    private static native int JNISolve(Model model, Material material,
				       BoundaryCond bCond,
				       bool enablePlaneStress,
				       double thickness,
				       MeshResults results,
				       Mesh mesh);

    public static MeshResults solve(Model model, Material material,
				    BoundaryCond bCond_vtx,
				    BoundaryCond bCond_sgm,
				    double thickness);
			     
}
