package vcn.pdeBot.finiteElement;

import vcn.geometricBot.Mesh;

public class MeshResults extends Mesh {
    protected float[] displacement;
    protected float[] strain;
    
    public float[] getDisplacementRef() {
	return displacement;
    }
    public float[] getStrainRef() {
	return strain;
    }
}
