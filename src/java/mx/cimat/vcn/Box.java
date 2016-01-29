/**
 * Box used by the Geometric Bot.
 * 2011-2015 Victor Eduardo Cardoso Nungaray
 * Twitter: @victore_cardoso
 * email: victorc@cimat.mx
 *
 */

package mx.cimat.vcn;

public class Box {
    protected float[] vA;
    protected float[] vB;

    protected Box()
    {
	vA = new float[2];
	vB = new float[2];
    }

    protected Box(float xA, float yA, float xB, float yB)
    {
	this();
	vA[0] = xA;
	vA[1] = yA;
	vB[0] = xB;
	vB[1] = yB;
    }

    public float getWidth()
    {
	return vB[0] - vA[0];
    }

    public float getHeight()
    {
	return vB[1] - vA[1];
    }
}