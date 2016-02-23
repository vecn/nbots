package nb.pdeBot.finiteElement;

import nb.geometricBot.Mesh;

public class MeshResults extends Mesh {
    protected float[] displacement; /* 2 values x node (u, v) */
    protected float[] strain; /* 3 values x elem (dux, dvy, 0.5*[dvx + duy]) */
    protected float[] VonMisesStress; /* 1 value x elem */

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
