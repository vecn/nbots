package nb.pdeBot.finiteElement;

import nb.geometricBot.Mesh;

public class MeshResults extends Mesh {
    private int status;
    protected float[] displacement;   /* 2 values (u, v) */
    protected float[] strain;         /* 3 values (dux, dvy, [dvx + duy]/2) */
    protected float[] VonMisesStress; /* 1 value */

    private MeshResults() {
	int status = -1;
	displacement = null;
	strain = null;
	VonMisesStress = null;
    }

    public float[] getDisplacementRef() {
	return displacement;
    }
    
    public float[] getStrainRef() {
	return strain;
    }

    public float[] getVonMisesStressRef() {
	return VonMisesStress;
    }

    public boolean FEMSuccess() {
	return status == 0;
    }

    private void setStatus(int status) {
	this.status = status;
    }
}
