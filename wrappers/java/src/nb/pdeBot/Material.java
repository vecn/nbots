package nb.pdeBot;

public class Material {
    protected double E;
    protected double v;
    public Material(double youngModulus,
		    double poissonModulus) {
	E = youngModulus;
	v = poissonModulus;
    }
    public void setYoungModulus(double E) {
	this.E = E;
    }
    public double getYoungModulus() {
	return E;
    }
    public void setPoissonModulus(double v) {
	this.v = v;
    }
    public double getPoissonModulus() {
	return v;
    }			     
}
