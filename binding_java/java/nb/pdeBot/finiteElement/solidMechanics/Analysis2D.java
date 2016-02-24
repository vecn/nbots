package nb.pdeBot.finiteElement.solidMechanics;

public enum Analysis2D {
    PLANE_STRESS,
    PLANE_STRAIN,
    SOLID_OF_REVOLUTION;

    double thickness;

    private Analysis2D() {
	thickness = 1.0;
    }

    public void setThickness(double thickness) {
	this.thickness = thickness;
    }

    public double getThickness() {
	return thickness;
    }
}
