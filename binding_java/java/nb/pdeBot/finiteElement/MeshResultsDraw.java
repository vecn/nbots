/**
 * Drawing functions of the Geometric Bot
 * 2011-2015 Victor Eduardo Cardoso Nungaray
 * Twitter: @victore_cardoso
 * email: victorc@cimat.mx
 *
 */
package nb.pdeBot.finiteElement;

import nb.geometricBot.*;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.Line2D;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class MeshResultsDraw {
    private MeshResultsDraw() {}

    public static void draw(MeshResults mesh, String imageFile,
			    int width, int height)
    {
	try {
	    /* TYPE_INT_ARGB specifies the image format: 
	     * 8-bit RGBA packed into integer pixels 
	     */
	    BufferedImage bi = 
		new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

	    Graphics2D ig2 = bi.createGraphics();

	    Box box = GeometricBot.getEnvelopingBox(mesh.getVerticesRef());
	    
	    float[] center = new float[2];
	    float zoom = GeometricBot.getCenterAndZoom(width, height,
						       box, center);

	    /* Draw background */
	    ig2.setColor(Color.WHITE);
	    ig2.fillRect(0, 0, width, height);
	    
	    ig2.setColor(Color.BLACK);
	    drawTrg(ig2, mesh, width, height, center, zoom);
	    
	    ImageIO.write(bi, "PNG", new File(imageFile));
	    /* You could use "JPEG", "gif" and "BMP" instead of "PNG" */      
	} catch (IOException ie) {
	    ie.printStackTrace();
	}
    }

    private static void drawTrg(Graphics2D g2, MeshResults mesh,
				int width, int height, float center[], 
				float zoom)
    {
	float d = (float)Math.sqrt(Math.pow(mesh.displacement[0], 2) +
				   Math.pow(mesh.displacement[1], 2));;
	float min = d;
	float max = d;
	for (int i = 1; i < mesh.getNVertices(); i++) {
	    d = (float)Math.sqrt(Math.pow(mesh.displacement[i * 2], 2) +
				 Math.pow(mesh.displacement[i*2+1], 2));
	    if (min > d)
		min = d;
	    if (max < d)
		max = d;
	}
	double scale = 1/max;

	for (int i = 0; i < mesh.getNElements(); i++) {
	    int id1 = mesh.getConnMtxRef()[i * 3];
	    int id2 = mesh.getConnMtxRef()[i*3+1];
	    int id3 = mesh.getConnMtxRef()[i*3+2];
	    double[] v1 = new double[2];
	    v1[0] = mesh.getVerticesRef()[id1 * 2] + scale * mesh.displacement[id1 * 2];
	    v1[1] = mesh.getVerticesRef()[id1*2+1] + scale * mesh.displacement[id1*2+1];
	    double[] v2 = new double[2];
	    v2[0] = mesh.getVerticesRef()[id2 * 2] + scale * mesh.displacement[id2 * 2];
	    v2[1] = mesh.getVerticesRef()[id2*2+1] + scale * mesh.displacement[id2*2+1];
	    double[] v3 = new double[2];
	    v3[0] = mesh.getVerticesRef()[id3 * 2] + scale * mesh.displacement[id3 * 2];
	    v3[1] = mesh.getVerticesRef()[id3*2+1] + scale * mesh.displacement[id3*2+1];
	    int[] x = new int[3];
	    int[] y = new int[3];
	    x[0] = (int)(zoom * v1[0] - zoom * center[0] + width / 2.0f);
	    x[1] = (int)(zoom * v2[0] - zoom * center[0] + width / 2.0f);
	    x[2] = (int)(zoom * v3[0] - zoom * center[0] + width / 2.0f);
	    y[0] = (int)(zoom * v1[1] - zoom * center[1] + height / 2.0f);
	    y[1] = (int)(zoom * v2[1] - zoom * center[1] + height / 2.0f);
	    y[2] = (int)(zoom * v3[1] - zoom * center[1] + height / 2.0f);
	    Color trgColor = getTrgColor(mesh, i, min, max);
	    g2.setColor(trgColor);
	    g2.fillPolygon(x, y, 3);
	    g2.setColor(Color.BLACK);
	    g2.drawPolygon(x, y, 3);
	}
    }

    private static Color getTrgColor(MeshResults mesh, int i,
				     float min, float max) {
	int id1 = mesh.getConnMtxRef()[i * 3];
	int id2 = mesh.getConnMtxRef()[i*3+1];
	int id3 = mesh.getConnMtxRef()[i*3+2];
	float d1 =  (float)Math.sqrt(Math.pow(mesh.displacement[id1 * 2], 2) +
				     Math.pow(mesh.displacement[id1*2+1], 2));
	float d2 =  (float)Math.sqrt(Math.pow(mesh.displacement[id2 * 2], 2) +
				     Math.pow(mesh.displacement[id2*2+1], 2));
	float d3 =  (float)Math.sqrt(Math.pow(mesh.displacement[id3 * 2], 2) +
				     Math.pow(mesh.displacement[id3*2+1], 2));
	float dist = (float)((d1 + d2 + d3)/3.0);
	float factor = (dist - min) / (max - min);
	return new Color(getRed(factor), getGreen(factor), getBlue(factor));
    }

    private static float getRed(float factor) {
	return factor;
    }
    private static float getGreen(float factor) {
	return factor * factor;
    }
    private static float getBlue(float factor) {
	return 1.0f - factor;
    }
}
