package vcn.pdeBot;

public class BoundaryConditions {
    protected BCType type;
    protected int N;
    protected int id[];
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
}
