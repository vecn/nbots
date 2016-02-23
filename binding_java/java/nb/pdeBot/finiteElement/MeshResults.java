package nb.pdeBot.finiteElement;

import nb.geometricBot.Mesh;

public class MeshResults extends Mesh {
    protected float[] displacement;   /* 2 values (u, v) */
    protected float[] strain;         /* 3 values (dux, dvy, [dvx + duy]/2) */
    protected float[] VonMisesStress; /* 1 value */

    public float[] getDisplacementRef() {
	return displacement;
    }
    public float[] getStrainRef() {
	return strain;
    }
    public float[] getVonMisesStressRef()
    {
	return VonMisesStress;
    }
}
