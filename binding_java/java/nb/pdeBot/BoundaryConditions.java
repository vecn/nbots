package nb.pdeBot;

public class BoundaryConditions {
    protected BCType type;
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
	return type == BCType.NEUMANN;
    }
    public int getN() {
	if (null == id)
	    return 0;
	else
	    return id.length;
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
    public void set(int id[], char dof[], double val[]) {
	this.id = id;
	this.dof = dof;
	this.val = val;
    }
}
