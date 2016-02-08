package vcn.pdeBot;

public class BoundaryConditions {
    protected BCType type;
    protected int N;
    protected int id[];
    protected char dof[]; /* 'x', 'y' or 'a' */
    protected double val[];
    private BoundaryConditions(){}
    public BoundaryConditions(BCType type) {
	this.type = type;
    }
    public boolean isDirichlet() {
	return type == BCType.DIRICHLET;
    }
    public boolean isNewmann() {
	return type == BCType.NEWMANN;
    }
    public int getN() {
	return N;
    }
    public int[] getIdsRef() {
	return id;
    }
    public char[] getDofRef() {
	return dof;
    }
    public double[] getValuesRef() {
	return val;
    }
}
